"""Independent path identities, with their mathematical authority stated per test.

These tests check the exploratory tools and explicitly derived laws. They do not
adopt a physical requirement for what incompatible fixed readings must report.
"""
from copy import deepcopy
import math

import numpy as np
import pytest

from met.inspection.case import Case, Observation
from met.inspection.geometry import profile
from met.inspection.paths import expand_paths, scale_noise
from met.noise_path import angle_reference, validate


def isotropic_observation(force, acceleration=(1.0, 0.0), sigma=1.0):
    return Observation("r0", force, acceleration, sigma ** 2 * np.eye(4))


def declared_case(name="correlated"):
    mixing = np.array([[1.0, 0.2, 0.1, 0.0],
                       [0.0, 0.8, 0.0, -0.1],
                       [0.1, 0.0, 0.6, 0.1],
                       [0.0, 0.1, 0.0, 0.9]])
    covariance = mixing @ mixing.T + 0.1 * np.eye(4)
    return Case.from_dict(dict(
        name=name, question="Vary the declared precision of fixed readings.",
        prior={"family": "half_cauchy", "center": 2.0},
        metadata={"units": {"force": "N"}}, future_note={"values": [1, 2]},
        truth={"mass": 2.0, "readings": [{"force_error": [0.1, -0.2]}]},
        readings=[dict(id="measured", force=[2.1, -0.2], acceleration=[1.0, 0.3],
                       covariance=covariance.tolist(), subspace=[[1.0], [0.0]],
                       window=[0.0, 1.0], metadata={"instrument": {"name": "sensor"}},
                       future_reading_note={"values": [3]})]))


def test_noise_scaling_preserves_declarations_truth_and_independent_copies():
    """A common SD factor changes every covariance entry by its square only."""
    source = declared_case()
    before = source.to_dict(include_truth=True)
    changed = scale_noise(source, 0.25)
    expected = deepcopy(before)
    expected["readings"][0]["covariance"] = (
        0.25 ** 2 * source.inputs.readings[0].covariance).tolist()
    assert changed.to_dict(include_truth=True) == expected
    assert source.to_dict(include_truth=True) == before
    assert "truth" not in changed.inputs.to_dict()

    changed.inputs.readings[0].force[0] = 999
    changed.inputs.readings[0].covariance[0, 0] = 999
    changed.inputs.readings[0].subspace[0, 0] = 0
    changed.inputs.readings[0].metadata["instrument"]["name"] = "changed"
    changed.inputs.readings[0].extras["future_reading_note"]["values"].append(4)
    changed.inputs.prior["center"] = 999
    changed.metadata["units"]["force"] = "changed"
    changed.extras["future_note"]["values"].append(3)
    changed.investigator_truth["readings"][0]["force_error"][0] = 999
    assert source.to_dict(include_truth=True) == before


def test_common_covariance_scaling_changes_cost_not_the_compatible_projection():
    """Multiplying the quadratic objective by a constant preserves its minimizer."""
    source = declared_case()
    scale = 0.2
    changed = scale_noise(source, scale)
    masses = np.array([0.2, 1.0, 4.0])
    original = profile(source.inputs.readings[0], masses)
    scaled = profile(changed.inputs.readings[0], masses)
    for key in ("candidate_force", "candidate_acceleration", "force_correction",
                "acceleration_correction", "conditional_mean"):
        np.testing.assert_allclose(scaled[key], original[key], rtol=3e-13, atol=3e-14)
    for key in ("residual_cost", "declared_cost", "within_subspace_cost",
                "outside_subspace_cost"):
        np.testing.assert_allclose(scaled[key], original[key] / scale ** 2, rtol=3e-13)
    np.testing.assert_allclose(scaled["conditional_covariance"],
                               scale ** 2 * original["conditional_covariance"], rtol=3e-13)


@pytest.mark.parametrize("scale", [0.0, -1.0, float("nan"), float("inf"), -float("inf")])
def test_noise_scale_requires_a_finite_positive_declaration(scale):
    with pytest.raises(ValueError):
        scale_noise(declared_case(), scale)


def test_path_expansion_retains_source_parameter_order_and_distinct_replays():
    sources = [declared_case("first"), declared_case("second")]
    before = [case.to_dict(include_truth=True) for case in sources]
    scales = [1.0, 0.3, 0.1]
    expanded = expand_paths(sources, scales)
    assert len(expanded) == 6
    assert len({case.name for case in expanded}) == 6
    for source in sources:
        selected = [case for case in expanded if case.metadata["path"]["source"] == source.name]
        assert [case.metadata["path"]["value"] for case in selected] == scales
        assert [case.metadata["path"]["index"] for case in selected] == [0, 1, 2]
        for case, scale in zip(selected, scales):
            assert case.metadata["path"]["parameter"] == "noise_scale"
            assert case.investigator_truth == source.investigator_truth
            assert case.metadata["units"] == source.metadata["units"]
            np.testing.assert_allclose(case.inputs.readings[0].covariance,
                                       scale ** 2 * source.inputs.readings[0].covariance)
    assert [case.to_dict(include_truth=True) for case in sources] == before


def test_path_transform_is_replaceable_without_changing_case_storage():
    """The exploration axis is a caller-supplied transform, not a fixed schema."""
    source = declared_case()

    def scale_force(case, scale):
        changed = Case.from_dict(case.to_dict(include_truth=True))
        changed.inputs.readings[0].force *= scale
        return changed

    expanded = expand_paths([source], [0.5, 2.0], transform=scale_force, parameter="force_scale")
    for case, scale in zip(expanded, [0.5, 2.0]):
        assert case.metadata["path"]["parameter"] == "force_scale"
        np.testing.assert_array_equal(case.inputs.readings[0].force,
                                      scale * source.inputs.readings[0].force)
        np.testing.assert_array_equal(case.inputs.readings[0].covariance,
                                      source.inputs.readings[0].covariance)


def test_nearby_path_coordinates_can_refine_a_threshold():
    expanded = expand_paths([declared_case()], [1.0000001, 1.0000002])
    assert len({case.name for case in expanded}) == 2
    assert [case.metadata["path"]["value"] for case in expanded] == [1.0000001, 1.0000002]


@pytest.mark.parametrize("scales", [[], [1.0, 1.0], [0.0], [-0.1], [np.nan], [np.inf]])
def test_path_values_must_be_distinct_finite_positive_scales(scales):
    with pytest.raises(ValueError):
        expand_paths([declared_case()], scales)


def test_noise_path_validation_expands_ordinary_saved_cases():
    """Supported saved cases retain annotations and investigator truth on expansion."""
    source = declared_case().to_dict(include_truth=True)
    source["readings"][0]["covariance"] = np.eye(4).tolist()
    source["readings"][0]["subspace"] = np.eye(2).tolist()
    document = dict(kind="noise_path", name="small-noise-path", seed=0,
                    question="How do fixed readings change with declared precision?",
                    noise_scales=[1.0, 0.5],
                    cases=[source])
    before = deepcopy(document)
    spec = validate(document)
    assert len(spec["cases"]) == 2 and len(spec["cells"]) == 2
    assert all(isinstance(case, Case) for case in spec["cases"])
    assert document == before
    for case in spec["cases"]:
        assert case.investigator_truth == document["cases"][0]["truth"]
        assert case.extras["future_note"] == {"values": [1, 2]}
    del document["noise_scales"]
    with pytest.raises(ValueError):
        validate(document)


@pytest.mark.parametrize("unsupported", ["multiple", "correlated", "restricted"])
def test_noise_recipe_rejects_inputs_outside_its_analytic_summary_scope(unsupported):
    """The generic path builder is broader than this one-reading isotropic recipe."""
    source = declared_case().to_dict(include_truth=True)
    original_covariance = source["readings"][0]["covariance"]
    source["readings"][0]["covariance"] = np.eye(4).tolist()
    source["readings"][0]["subspace"] = np.eye(2).tolist()
    if unsupported == "multiple":
        second = deepcopy(source["readings"][0])
        second["id"] = "another-reading"
        source["readings"].append(second)
    elif unsupported == "correlated":
        source["readings"][0]["covariance"] = original_covariance
    else:
        source["readings"][0]["subspace"] = [[1.0], [0.0]]
    document = dict(kind="noise_path", name="out-of-scope", seed=0,
                    noise_scales=[1.0, 0.5], cases=[source])
    # Case construction and expansion remain valid for these scientific inputs.
    assert len(expand_paths([Case.from_dict(source)], [1.0, 0.5])) == 2
    with pytest.raises(ValueError):
        validate(document)


@pytest.mark.parametrize("sigma,weight", [(0.8, 1.0), (0.1, 1.0), (0.1, 0.01)])
def test_perpendicular_equal_energy_reference_is_exactly_the_prior(sigma, weight):
    """H is constant in mass, so normalization leaves the half-Cauchy prior."""
    masses = np.array([0.03, 0.2, 1.0, 5.0, 30.0])
    got = angle_reference(isotropic_observation([0.0, 1.0], sigma=sigma),
                          prior_center=1.0, weight=weight, masses=masses)
    np.testing.assert_allclose(got["cdf"], 2 / math.pi * np.arctan(masses),
                               rtol=2e-11, atol=2e-12)
    assert math.isfinite(got["normalizer"]) and got["normalizer"] > 0
    assert 0 <= got["quadrature_error_bound"] < 1e-7


def test_aligned_reference_concentrates_inside_a_fixed_band_as_noise_shrinks():
    """The unique zero-cost positive mass is 1 for the aligned unit pair."""
    masses = np.array([0.5, 1.0, 2.0])
    broad = angle_reference(isotropic_observation([1.0, 0.0], sigma=0.8),
                            prior_center=1.0, weight=1.0, masses=masses)["cdf"]
    narrow = angle_reference(isotropic_observation([1.0, 0.0], sigma=0.08),
                             prior_center=1.0, weight=1.0, masses=masses)["cdf"]
    assert narrow[2] - narrow[0] > 0.999
    assert broad[2] - broad[0] < 0.9
    assert broad[1] == pytest.approx(0.5, abs=2e-11)
    assert narrow[1] == pytest.approx(0.5, abs=2e-11)


def test_acute_reference_freezes_when_inverse_misfit_weight_cancels_noise_scale():
    """For cos(theta)=0.6, Q=0.4/sigma² and g=sigma²/0.4 in the saturated regime."""
    masses = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
    references = []
    for sigma in (0.2, 0.1):
        references.append(angle_reference(
            isotropic_observation([0.6, 0.8], sigma=sigma), prior_center=1.0,
            weight=sigma ** 2 / 0.4, masses=masses)["cdf"])
    np.testing.assert_allclose(references[0], references[1], rtol=2e-11, atol=2e-12)


def test_opposite_reference_has_reciprocal_symmetry_while_leaving_a_fixed_band():
    """Reflection phi -> pi/2-phi fixes the median, but the endpoints dominate."""
    masses = np.array([0.1, 0.5, 1.0, 2.0, 10.0])
    references = []
    for sigma in (0.3, 0.15):
        got = angle_reference(isotropic_observation([-1.0, 0.0], sigma=sigma),
                              prior_center=1.0, weight=1.0, masses=masses)["cdf"]
        np.testing.assert_allclose(got + got[::-1], 1.0, rtol=0, atol=2e-11)
        assert got[2] == pytest.approx(0.5, abs=2e-11)
        references.append(got)
    wide_band = references[0][3] - references[0][1]
    narrow_band = references[1][3] - references[1][1]
    assert 0 <= narrow_band < wide_band


def test_complementary_angles_have_opposite_preference_despite_equal_unsigned_misfit():
    """Changing cos(theta)'s sign reverses the excess-cost tilt around mass 1."""
    masses = np.array([0.5, 2.0])
    sigma = 0.2
    weight = sigma ** 2 / 0.4
    acute = angle_reference(isotropic_observation([0.6, 0.8], sigma=sigma),
                            prior_center=1.0, weight=weight, masses=masses)["cdf"]
    obtuse = angle_reference(isotropic_observation([-0.6, 0.8], sigma=sigma),
                             prior_center=1.0, weight=weight, masses=masses)["cdf"]
    prior_band = 2 / math.pi * (math.atan(2.0) - math.atan(0.5))
    assert acute[1] - acute[0] > prior_band > obtuse[1] - obtuse[0]


def test_angle_reference_prior_jacobian_with_unequal_instruments():
    """With zero evidence weight the CDF is that of the declared prior at any s."""
    observation = Observation("r", [2, 0], [1, 0], np.diag([0.16, 0.16, 0.01, 0.01]))
    masses = np.array([0.01, 0.5, 3, 10, 100.0])
    got = angle_reference(observation, prior_center=3, weight=0, masses=masses)
    np.testing.assert_allclose(got["cdf"], 2 / math.pi * np.arctan(masses / 3), atol=2e-12)


@pytest.mark.parametrize("center", [1e-8, 1e8])
def test_null_reference_resolves_a_prior_far_from_the_instrument_ratio(center):
    observation = isotropic_observation([0, 0], acceleration=[0, 0])
    masses = center * np.array([0.1, 1, 10])
    got = angle_reference(observation, prior_center=center, weight=1, masses=masses)
    np.testing.assert_allclose(got["cdf"], 2 / math.pi * np.arctan(masses / center), atol=2e-7)
    assert got["normalizer"] == pytest.approx(math.pi / 2, rel=2e-7)


def test_noise_recipe_exports_finite_comparison_and_replay(tmp_path):
    """The complete route must retain native laws, unavailable scalars and saved inputs."""
    import json
    import yaml
    from met.noise_path import run
    from met.inspection.case import load_cases

    case = dict(name="null-force", prior=dict(family="half_cauchy", center=1),
                readings=[dict(force=[0, 0], acceleration=[1, 0], covariance=np.eye(4).tolist())])
    document = dict(kind="noise_path", name="replay", seed=0, noise_scales=[0.5, 0.2],
                    detail_scales=[], mass_band=[0.5, 4], checks=[], cases=[case],
                    numerics=dict(nodes=11, log_half_span=2))
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    out, bundle = run(path, tmp_path, progress=None)
    comparison = json.loads((out / "comparison.json").read_text(encoding="utf-8"))
    assert len(comparison["rows"]) == 2
    assert all(row["angle_degrees"] is None for row in comparison["rows"])
    assert all(row["laws"]["unweighted"]["curvature"] is None for row in comparison["rows"])
    assert (out / "recipe.yaml").read_bytes() == path.read_bytes()
    assert "Path identities" in (out / "comparison.html").read_text(encoding="utf-8")
    assert len(bundle["cases"]) == 2
    replay, cases = load_cases(out / "case_001" / "case.yaml")
    assert replay["checks"] == [] and len(replay["numerics"]["mass_grid"]) == 11
    np.testing.assert_allclose(cases[0].inputs.readings[0].covariance, 0.04 * np.eye(4))


def test_noise_recipe_refuses_parallel_workers_and_replicate_override():
    from met.runner import run
    path = "presets/inspection_noise_paths.yaml"
    for kwargs in (dict(workers=2), dict(workers="max"), dict(replicates=2)):
        with pytest.raises(ValueError, match="serial"):
            run(path, **kwargs)

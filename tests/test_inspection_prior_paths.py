"""Declared-prior changes preserve the experiment; final answers may depend on it."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest
import yaml

from met.inspection.case import Case, load_cases
from met.inspection.paths import expand_paths, set_prior_center
from met.inspection.references import equal_perpendicular_point
from met.inspection.run import inspect_case
from met.prior_path import compare_paths, run, summarize, validate


def case():
    return Case.from_dict(dict(name="perpendicular", question="Separate law and point.",
        prior=dict(family="half_cauchy", center=1), metadata={"units": {"mass": "kg"}},
        notes={"values": [1]}, truth={"private": [2]}, readings=[dict(
            force=[1, 0], acceleration=[0, 1], covariance=(0.01 * np.eye(4)).tolist(),
            metadata={"instrument": "fixed"}, subspace=np.eye(2).tolist(), window=[0, 1])]))


def recipe():
    return dict(kind="prior_path", name="prior-test", seed=0, cases=[case().to_dict(include_truth=True)],
                prior_centers=[0.1, 1, 10], detail_centers=[], mass_band=[0.5, 4], checks=[],
                numerics=dict(nodes=31, log_half_span=5))


def test_prior_transform_changes_only_the_explicit_center_and_copies_truth():
    original = case()
    before = original.to_dict(include_truth=True)
    changed = set_prior_center(original, 10)
    expected = deepcopy(before)
    expected["prior"]["center"] = 10
    assert changed.to_dict(include_truth=True) == expected
    changed.inputs.readings[0].force[0] = 7
    changed.inputs.readings[0].metadata["instrument"] = "changed"
    changed.inputs.prior["center"] = 3
    changed.investigator_truth["private"].append(3)
    changed.extras["notes"]["values"].append(2)
    assert original.to_dict(include_truth=True) == before
    assert "truth" not in changed.inputs.to_dict()


@pytest.mark.parametrize("center", [0, -1, True, float("nan"), float("inf")])
def test_prior_center_must_be_a_valid_declared_scale(center):
    with pytest.raises(ValueError):
        set_prior_center(case(), center)


def test_custom_prior_transform_must_define_its_own_parameter_semantics():
    original = case()
    original.inputs.prior = dict(family="future", center=1, extra=2)
    with pytest.raises(ValueError, match="other families"):
        set_prior_center(original, 2)


def test_validation_keeps_one_physical_mass_grid_across_prior_centres():
    source = recipe()
    before = deepcopy(source)
    spec = validate(source)
    assert source == before
    assert [c.inputs.prior["center"] for c in spec["cases"]] == [0.1, 1, 10]
    assert [c.metadata["path"]["parameter"] for c in spec["cases"]] == ["prior_center"] * 3
    grid = spec["numerics"]["mass_grid"]
    assert len(grid) == 31 and grid[15] == 1
    assert all(c.investigator_truth == source["cases"][0]["truth"] for c in spec["cases"])


@pytest.mark.parametrize("centers", [[], [1, 1], [0], [float("nan")]])
def test_validation_rejects_invalid_prior_paths(centers):
    document = recipe()
    document["prior_centers"] = centers
    with pytest.raises(ValueError):
        validate(document)


def test_perpendicular_law_is_the_supplied_prior_while_cloud_point_has_its_own_reference():
    results, rows = [], []
    grid = np.exp(np.linspace(-7, 7, 31))
    for changed in expand_paths([case()], [0.1, 1, 10], transform=set_prior_center, parameter="prior_center"):
        result = inspect_case(changed, masses=grid, checks=[])
        center = changed.inputs.prior["center"]
        quoted = result["evaluation"]["laws"]["quoted"]
        # Analytic half-Cauchy density per log mass, independent of adapter code.
        expected = 1 / (math.pi * np.cosh(quoted["log_mass"] - math.log(center)))
        np.testing.assert_allclose(quoted["density"], expected, rtol=3e-9, atol=1e-12)
        assert quoted["median"] == pytest.approx(center, rel=2e-7)
        assert quoted["point"] == pytest.approx(equal_perpendicular_point(center), rel=2e-7)
        if center != 1:
            assert not np.isclose(quoted["point"], center, rtol=0.1)
        rows.append(summarize(result, [0.5, 4]))
        results.append(result)
    comparison = compare_paths(rows, results)[0]
    assert max(comparison["maximum_absolute_errors"].values()) < 1e-8
    # A future adapter cannot quietly make the evidence depend on the chosen prior.
    results[-1]["evaluation"]["laws"]["quoted"]["log_kernel"][0] += 0.1
    assert compare_paths(rows, results)[0]["maximum_absolute_errors"]["data_only_log_evidence"] > 0.09


def test_prior_recipe_saves_comparisons_truth_and_exact_replay_grid(tmp_path):
    document = recipe()
    document["prior_centers"] = [0.1, 10]
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    out, bundle = run(path, tmp_path, progress=None)
    comparison = json.loads((out / "comparison.json").read_text(encoding="utf-8"))
    assert len(bundle["cases"]) == 2 and len(comparison["rows"]) == 2
    assert comparison["rows"][0]["laws"]["quoted"]["absolute_log_point_error"] < 1e-7
    assert (out / "recipe.yaml").read_bytes() == path.read_bytes()
    assert "What does the declared prior determine?" in (out / "comparison.html").read_text(encoding="utf-8")
    replay_a, a = load_cases(out / "case_000" / "case.yaml")
    replay_b, b = load_cases(out / "case_001" / "case.yaml")
    assert replay_a["numerics"]["mass_grid"] == replay_b["numerics"]["mass_grid"]
    assert replay_a["truth"] == replay_b["truth"] == document["cases"][0]["truth"]
    assert a[0].inputs.prior["center"] == 0.1 and b[0].inputs.prior["center"] == 10


def test_prior_recipe_refuses_parallelism_and_replicate_overrides():
    from met.runner import run as runner
    for kwargs in (dict(workers=2), dict(workers="max"), dict(replicates=2)):
        with pytest.raises(ValueError, match="serial"):
            runner("presets/inspection_prior_paths.yaml", **kwargs)


@pytest.mark.xfail(strict=True, reason="Exact reciprocal symmetry fixes median 1. With the saved covariance 0.01 I, "
    "met.full instead reports 0.4230518228822659 on the current numerical CDF plateau. "
    "This is the existing O13 median defect; equation and implementation are retained.")
def test_saved_opposite_prior_case_has_mathematical_median_one():
    from met.full import Reading, ours29
    from met.preset import load
    chosen = next(c for c in load("presets/inspection_prior_paths.yaml")["cases"]
                  if c.metadata["path"]["source"] == "equal-opposite" and c.inputs.prior["center"] == 1)
    r = chosen.inputs.readings[0]
    result = ours29([Reading(r.force, r.acceleration, r.covariance, r.subspace)], c=1, readout=False)
    assert result["median"] == pytest.approx(1, rel=2e-7)

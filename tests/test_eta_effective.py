"""B04: fitted-score moments, the isolated eta-prime probe, and paired replay."""
from copy import deepcopy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
import yaml

import met.full as full
from met.eta_effective import intervals, linear_moments, statistics
from met.law import ReadingSet


def _batch(force, acceleration, sf=1.0, sa=1.0):
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    shape = force.shape[:2]
    return ReadingSet(force, acceleration, np.broadcast_to(sf, shape), np.broadcast_to(sa, shape))


def _full(readings, index):
    return [full.isotropic(f, a, sf, sa) for f, a, sf, sa in zip(
        readings.force[index], readings.acceleration[index],
        readings.force_sd[index], readings.acceleration_sd[index])]


@pytest.mark.parametrize("c", [[1.0], [2.0, 2.0, 2.0], [0.25, 2, 7], [1, 0, 3, 0, 8]])
def test_linear_gaussian_bias_is_an_exact_matrix_moment(c):
    c = np.asarray(c, float)
    result = linear_moments(c)
    total = c.sum()
    # Scores are A epsilon, with independent epsilon_i of variance 1/c_i.
    # Zero precisions give identically zero scores and need no random variable.
    residualizer = np.eye(len(c)) - np.ones((len(c), 1)) @ (c / total)[None]
    a = np.diag(c) @ residualizer
    variances = np.divide(1.0, c, out=np.zeros_like(c), where=c > 0)
    covariance = a @ np.diag(variances) @ a.T
    np.testing.assert_allclose(result["score_covariance"], covariance, atol=2e-14)
    assert result["expected_J"] == pytest.approx(np.trace(covariance), abs=2e-14)
    assert result["variance_J"] == pytest.approx(2 * np.trace(covariance @ covariance), abs=2e-13)
    assert result["C"] == total
    assert result["n_eff"] == pytest.approx(total ** 2 / np.dot(c, c))
    assert result["expected_J"] / total == pytest.approx(1 - 1 / result["n_eff"], abs=2e-14)


def test_null_precision_leaves_exact_linear_moments_unchanged():
    before, after = linear_moments([1, 4, 2]), linear_moments([1, 4, 2, 0])
    for key in ("C", "n_eff", "expected_J"):
        assert before[key] == after[key]
    np.testing.assert_array_equal(after["score_covariance"][-1], np.zeros(4))


@pytest.mark.parametrize("c", [[], [-1, 2], [0, 0], [np.nan, 1]])
def test_linear_moments_rejects_inputs_outside_the_gaussian_identity(c):
    with pytest.raises(ValueError):
        linear_moments(c)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_vectorized_fitted_scores_curvatures_and_eta_match_full_type(k):
    rng = np.random.default_rng(4040 + k)
    f, a = rng.normal(size=(3, 7, k)), rng.normal(size=(3, 7, k))
    sigma = rng.uniform(0.5, 2.0, size=(3, 7))
    given = _batch(f, a, 1.7 * sigma, sigma)
    result = statistics(given)
    for b in range(3):
        readings = _full(given, b)
        states = [full.reading_state(r, 1.7) for r in readings]
        g = [full.misfit_weight(s["Q_plus"], s["kappa"]) for s in states]
        eta, curvature, j, phi = full.bartlett_eta(readings, 1.7, g)
        terms = [full._h_terms(r.pre(), 1.7, np.array([phi]), 2) for r in readings]
        scores = np.array([gi * h[1][0] / 2 for gi, h in zip(g, terms)])
        curvatures = np.array([-gi * h[2][0] / 2 for gi, h in zip(g, terms)])
        np.testing.assert_allclose(result["g"][b], g, atol=2e-12)
        np.testing.assert_allclose(result["scores"][b], scores, atol=3e-11)
        np.testing.assert_allclose(result["curvatures"][b], curvatures, atol=3e-11)
        assert result["C"][b] == pytest.approx(curvature, abs=3e-11)
        assert result["J"][b] == pytest.approx(j, abs=3e-11)
        assert result["eta"][b] == pytest.approx(eta, abs=3e-11)
        assert abs(np.sin(result["phi"][b] - phi)) < 2e-11
        assert result["scores"][b].sum() == pytest.approx(0, abs=2e-12)
        assert result["n_eff"][b] == pytest.approx(curvature ** 2 / np.dot(curvatures, curvatures))


def test_candidate_is_applied_to_the_same_evidence_without_changing_g():
    given = _batch([[[3], [1]]], [[[1], [3]]])
    result = statistics(given)
    assert result["C"][0] == pytest.approx(12)
    assert result["J"][0] == pytest.approx(32)
    assert result["n_eff"][0] == pytest.approx(2)
    assert result["eta"][0] == pytest.approx(3 / 8)
    assert result["eta_prime"][0] == pytest.approx(3 / 16)
    np.testing.assert_allclose(result["g"], 1)


def test_candidate_corrects_the_uncapped_ratio_before_applying_the_cap():
    result = statistics(_batch([[[2], [1]]], [[[1], [2]]]))
    assert result["C"][0] == pytest.approx(8)
    assert result["J"][0] == pytest.approx(4.5)
    assert result["n_eff"][0] == pytest.approx(2)
    assert result["eta"][0] == 1
    assert result["eta_prime"][0] == pytest.approx(8 / 9)


def test_signed_curvatures_can_put_neff_below_one_and_trigger_stated_fallback():
    given = _batch([[[2], [0], [0.5]]], [[[0], [1], [0.5]]])
    result = statistics(given)
    assert np.min(result["curvatures"]) < 0
    assert result["J"][0] > 0
    assert 0 < result["n_eff"][0] < 1
    assert result["eta_prime"][0] == result["eta"][0]


def test_cancelling_readings_exhibit_inverse_square_target_and_inverse_radius_ratio():
    """An admitted open g=1 region approaches a nonzero-curvature cancellation.

    O14 supplies the probability-density argument for divergent target mean and
    ratio variance. This checks the exact residual-ray algebra, not sampled tails.
    """
    t, theta = 0.25, math.pi / 4
    epsilon = np.array([1e-3, 1e-4, 1e-5])
    force, acceleration = np.zeros((3, 5, 1)), np.zeros((3, 5, 1))
    force[:, 0, 0] = acceleration[:, 0, 0] = t
    second = np.sqrt(-2j * t ** 2 + epsilon * np.exp(1j * theta))
    acceleration[:, 1, 0], force[:, 1, 0] = second.real, second.imag
    result = statistics(_batch(force, acceleration))
    np.testing.assert_array_equal(result["g"], np.ones((3, 5)))
    np.testing.assert_allclose(result["C"], epsilon, rtol=5e-12)
    first_curvature = 2 * t ** 2 * math.sin(theta)
    exact_square_sum = first_curvature ** 2 + (epsilon - first_curvature) ** 2
    np.testing.assert_allclose(epsilon ** 2 * (1 - result["linear_target"]), exact_square_sum, rtol=2e-11)
    exact_j = 2 * t ** 4 * math.cos(theta) ** 2
    np.testing.assert_allclose(result["J"], exact_j, rtol=2e-11)
    np.testing.assert_allclose(epsilon ** 2 * result["ratio"] ** 2, exact_j ** 2, rtol=2e-11)
    assert result["linear_target"][-1] < 80 * result["linear_target"][0]
    assert result["ratio"][-1] > 99 * result["ratio"][0]


def test_all_zero_series_exposes_undefined_ratios_and_retains_bartlett_fallback():
    result = statistics(_batch(np.zeros((2, 3, 1)), np.zeros((2, 3, 1))))
    np.testing.assert_array_equal(result["C"], [0, 0])
    np.testing.assert_array_equal(result["J"], [0, 0])
    for key in ("n_eff", "linear_target", "ratio"):
        assert np.isnan(result[key]).all()
    np.testing.assert_array_equal(result["eta"], [1, 1])
    np.testing.assert_array_equal(result["eta_prime"], [1, 1])


def test_appended_exact_zero_readings_leave_statistics_and_both_laws_unchanged():
    given = _batch([[[3], [1]]], [[[1], [3]]])
    padded = _batch(np.concatenate([given.force, np.zeros((1, 5, 1))], axis=1),
                    np.concatenate([given.acceleration, np.zeros((1, 5, 1))], axis=1))
    before, after = statistics(given), statistics(padded)
    for key in ("C", "J", "n_eff", "linear_target", "ratio", "eta", "eta_prime"):
        np.testing.assert_allclose(before[key], after[key], atol=1e-13)
    np.testing.assert_array_equal(after["scores"][:, 2:], np.zeros((1, 5)))
    np.testing.assert_array_equal(after["curvatures"][:, 2:], np.zeros((1, 5)))
    for key, value in intervals(given, before).items():
        np.testing.assert_allclose(value, intervals(padded, after)[key], rtol=1e-11)


def test_units_and_swap_preserve_powers_and_transform_intervals():
    given = _batch([[[3], [1], [2]]], [[[1], [3], [1.2]]])
    changed = _batch(given.force * 7, given.acceleration * 2, 7, 2)
    swapped = _batch(given.acceleration, given.force)
    result, scaled, swap = statistics(given), statistics(changed), statistics(swapped)
    for key in ("C", "J", "n_eff", "eta", "eta_prime"):
        np.testing.assert_allclose(result[key], scaled[key], atol=2e-12)
        np.testing.assert_allclose(result[key], swap[key], atol=2e-12)
    ordinary, in_units, inverse = intervals(given, result), intervals(changed, scaled), intervals(swapped, swap)
    for prefix in ("interval", "prime"):
        np.testing.assert_allclose(in_units[prefix + "_lo"], 3.5 * ordinary[prefix + "_lo"], rtol=2e-7)
        np.testing.assert_allclose(in_units[prefix + "_hi"], 3.5 * ordinary[prefix + "_hi"], rtol=2e-7)
        np.testing.assert_allclose(inverse[prefix + "_lo"], 1 / ordinary[prefix + "_hi"], rtol=2e-7)
        np.testing.assert_allclose(inverse[prefix + "_hi"], 1 / ordinary[prefix + "_lo"], rtol=2e-7)


def _independent_interval(given, stats, power):
    f = given.force[0] / given.force_sd[0, :, None]
    a = given.acceleration[0] / given.acceleration_sd[0, :, None]
    g = stats["g"][0]
    def h(phi):
        return float(np.dot(g, np.sum((f * math.sin(phi) + a * math.cos(phi)) ** 2, axis=1)))
    search = np.linspace(0, math.pi / 2, 4001)
    maximum = max(h(phi) for phi in search)
    density = lambda phi: math.exp(power * (h(phi) - maximum) / 2)
    normalizer = quad(density, 0, math.pi / 2, epsabs=1e-12, epsrel=1e-12)[0]
    c = float(given.s[0, 0])
    return [c * math.tan(brentq(lambda phi: quad(density, 0, phi, epsabs=1e-12)[0]
                               - probability * normalizer, 0, math.pi / 2, xtol=1e-13))
            for probability in (0.025, 0.975)]


@pytest.mark.parametrize("force,acceleration", [
    ([[3], [1]], [[1], [3]]),
    ([[-1], [0.5], [0.2]], [[1], [-0.8], [0.1]]),
    ([[20], [18], [23], [17]], [[19], [20], [21], [20]]),
    ([[0], [0]], [[0], [0]]),
])
def test_paired_intervals_match_full_estimator_and_independent_integrals(force, acceleration):
    given = _batch([force], [acceleration])
    result = statistics(given)
    measured = intervals(given, result)
    oracle = full.ours30(_full(given, 0), 1.0, readout=False)
    actual = [measured["interval_lo"][0], measured["interval_hi"][0]]
    np.testing.assert_allclose(actual, oracle["interval"], rtol=3e-5)
    for prefix, power in (("interval", result["eta"][0]), ("prime", result["eta_prime"][0])):
        expected = _independent_interval(given, result, power)
        np.testing.assert_allclose([measured[prefix + "_lo"][0], measured[prefix + "_hi"][0]],
                                   expected, rtol=3e-5)


def _recipe():
    return dict(kind="eta_effective", name="eta-effective-test", question="Check the finite-N identity.",
                seed=20261029, replicates=8, world=dict(dimension=2, design="new_excitation",
                mass=1.0, acceleration_snr=1.0, direction="random", readings=5,
                noise={"force_sd": 1.0, "acceleration_sd": 1.0}, supplied_noise="exact"))


def test_saved_study_is_paired_reproducible_and_preserves_replicate_prefix(tmp_path):
    from met.eta_effective_study import run
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(_recipe()), encoding="utf-8")
    first, record = run(path, tmp_path, progress=None)
    second, _ = run(path, tmp_path, progress=None)
    longer, _ = run(path, tmp_path, replicates=12, progress=None)
    assert first != second
    with np.load(first / "values.npz") as one, np.load(second / "values.npz") as two, np.load(longer / "values.npz") as more:
        assert sorted(one.files) == sorted(two.files)
        for key in one.files:
            np.testing.assert_array_equal(one[key], two[key])
            np.testing.assert_array_equal(one[key], more[key][:8])
    assert (first / "preset.yaml").read_bytes() == path.read_bytes()
    assert yaml.safe_load((longer / "preset.yaml").read_text(encoding="utf-8"))["replicates"] == 12
    assert (longer / "source_preset.yaml").read_bytes() == path.read_bytes()
    for name in ("summary.md", "summary.csv", "run.json"):
        assert (first / name).exists()
    restored = json.loads((first / "run.json").read_text(encoding="utf-8"))
    assert restored["workers"] == record["workers"] == 1
    assert restored["code"]["preset_sha256"]


@pytest.mark.parametrize("change", [
    lambda d: d.update(replicates=True), lambda d: d.update(seed=-1),
    lambda d: d.update(name="../escape"), lambda d: d["world"].update(dimension=4),
    lambda d: d["world"].update(acceleration_snr={"uniform": [0, 1]}),
    lambda d: d["world"].update(supplied_noise={"force_scale": 2, "acceleration_scale": 2}),
])
def test_study_validates_its_scope(change):
    from met.eta_effective_study import validate
    data = deepcopy(_recipe())
    change(data)
    with pytest.raises(ValueError):
        validate(data)


def test_summary_distinguishes_mean_ratio_and_ratio_of_means_and_uses_paired_errors():
    from met.eta_effective_study import summarize
    c, j = np.array([1., 3., 2., 6.]), np.array([2., 3., 6., 6.])
    target = np.array([0.2, 0.3, 0.5, 0.7])
    covered, prime = np.array([0., 1., 1., 0.]), np.array([1., 1., 1., 0.])
    ratio, gap = j / c, j / c - target
    result = summarize(dict(C=c, J=j, linear_target=target, ratio=ratio, paired_gap=gap,
                            covered=covered, prime_covered=prime))
    assert result["summary"]["ratio"]["mean"] == pytest.approx(1.75)
    assert result["expectation_ratio"]["value"] == pytest.approx(17 / 12)
    assert result["summary"]["ratio"]["mc_se"] == pytest.approx(np.std(ratio, ddof=1) / 2)
    assert result["summary"]["paired_gap"]["mc_se"] == pytest.approx(np.std(gap, ddof=1) / 2)
    influence = (j - (17 / 12) * c) / np.mean(c)
    assert result["expectation_ratio"]["mc_se"] == pytest.approx(np.std(influence, ddof=1) / 2)
    weighted_target = np.dot(c, target) / c.sum()
    assert result["expectation_ratio"]["curvature_weighted_target"] == pytest.approx(weighted_target)
    assert result["coverage"]["ours30"]["coverage"] == 0.5
    assert result["coverage"]["eta_prime"]["coverage"] == 0.75
    assert result["coverage"]["paired_change"]["mean"] == 0.25
    assert result["coverage"]["paired_change"]["mc_se"] == pytest.approx(np.std(prime - covered, ddof=1) / 2)


def test_summary_reports_undefined_replicates_without_changing_the_population():
    from met.eta_effective_study import summarize
    result = summarize(dict(C=np.array([1., 0.]), J=np.array([2., 0.]),
                            linear_target=np.array([0.5, np.nan]), ratio=np.array([2., np.nan]),
                            covered=np.array([1., 0.]), prime_covered=np.array([1., 1.])))
    assert result["summary"]["ratio"]["mean"] is None
    assert result["summary"]["ratio"]["finite_count"] == 1
    assert result["summary"]["ratio"]["undefined_count"] == 1
    assert result["expectation_ratio"]["value"] is None
    assert result["expectation_ratio"]["undefined_count"] == 1
    assert result["coverage"]["ours30"]["count"] == 2


def test_summary_undefined_count_includes_nonfinite_scores_with_finite_target():
    from met.eta_effective_study import summarize
    result = summarize(dict(C=np.array([1., 2.]), J=np.array([2., np.nan]),
                            linear_target=np.array([0.5, 0.5]), ratio=np.array([2., np.nan]),
                            covered=np.array([1., 0.]), prime_covered=np.array([1., 1.])))
    assert result["expectation_ratio"]["value"] is None
    assert result["expectation_ratio"]["undefined_count"] == 1


def test_shipped_recipe_covers_requested_dimensions_signals_and_sample_sizes():
    from met.preset import load
    path = Path(__file__).resolve().parents[1] / "presets" / "eta_effective_readings.yaml"
    spec = load(path)
    assert spec["raw"]["seed"] == 20261029
    assert len(spec["cells"]) == 72
    assert {cell["world"]["dimension"] for cell in spec["cells"]} == {1, 2, 3}
    assert {cell["world"]["readings"] for cell in spec["cells"]} == {5, 20, 100, 400}


def test_cli_validates_registered_study_kind(capsys):
    from met.__main__ import main
    path = Path(__file__).resolve().parents[1] / "presets" / "eta_effective_readings.yaml"
    assert main(["validate", str(path)]) == 0
    output = capsys.readouterr().out
    assert "72" in output and "serial" in output


def test_production_runner_rejects_parallel_study_recipe(tmp_path):
    from met.runner import run
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(_recipe()), encoding="utf-8")
    for workers in (2, "max", True, 1.0):
        with pytest.raises(ValueError, match="serial"):
            run(path, tmp_path, workers=workers, progress=None)

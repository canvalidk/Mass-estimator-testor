"""T56/T57 full-domain quadrature and independent Gaussian scalar geometry."""
import math

import numpy as np
import pytest
from scipy import integrate

from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


def _quad(count, **controls):
    case = witness_constants()
    return ScalarScoreQuadrature(count*case["population_scatter"], case["covariance"],
        prior_center=1., anchor_exponent=.5, epsabs=controls.get("epsabs", 1e-11),
        epsrel=controls.get("epsrel", 1e-10), limit=500)


def test_positive_orthogonal_witness_and_two_laplace_weights():
    case = witness_constants()
    true = np.array([case["true_mass"], 1.])
    false = np.array([case["false_mass"], 1.])
    precision = np.linalg.inv(case["covariance"])
    assert true@precision@false == pytest.approx(0., abs=2e-15)
    assert true@precision@true == pytest.approx(28/19, rel=1e-13)
    assert case["false_limit_weight"] == pytest.approx(.753554709424963, rel=1e-13)
    split = math.sqrt(case["true_mass"]*case["false_mass"])
    errors = []
    for count in (100, 1000, 10000):
        q = _quad(count)
        errors.append(abs(q.cdf(split, "score")-case["false_limit_weight"]))
    assert errors[-1] < 1e-4 and errors[-1] < errors[0]
    assert _quad(10000).quantile(.5, "score") < .3
    assert abs(_quad(10000).quantile(.5, "raw")-case["true_mass"]) < .005


def test_anchor_keeps_global_profile_ordering_in_this_witness():
    case = witness_constants()
    split = math.sqrt(case["true_mass"]*case["false_mass"])
    false_probabilities = []
    for count in (100, 1000, 10000):
        q = _quad(count)
        false_probabilities.append(q.cdf(split, "anchored_score"))
    assert false_probabilities[0] > false_probabilities[1] > false_probabilities[2]
    assert false_probabilities[-1] < 1e-20
    assert abs(_quad(10000).quantile(.5, "anchored_score")-case["true_mass"]) < .005


def test_full_positive_domain_quadrature_refinement():
    ordinary, fine = _quad(10000), _quad(10000, epsabs=1e-13, epsrel=1e-12)
    for method in ("raw", "score", "anchored_score"):
        assert ordinary.total(method)[0] == pytest.approx(fine.total(method)[0], rel=1e-10)
        assert ordinary.quantile(.5, method) == pytest.approx(fine.quantile(.5, method), abs=2e-10)


def test_complete_scalar_laws_respect_reciprocal_dual_and_noise_units():
    case = witness_constants()
    scatter, covariance = 1000*case["population_scatter"], case["covariance"]
    original = _quad(1000)
    swap = np.array([[0., 1.], [1., 0.]])
    dual = ScalarScoreQuadrature(swap@scatter@swap, swap@covariance@swap,
        prior_center=1., anchor_exponent=.5, epsabs=1e-11, epsrel=1e-10, limit=500)
    units = np.diag([3., .7])
    changed = ScalarScoreQuadrature(units@scatter@units, units@covariance@units,
        prior_center=3/.7, anchor_exponent=.5, epsabs=1e-11, epsrel=1e-10, limit=500)
    assert dual.oscillation == pytest.approx(original.oscillation, rel=1e-13)
    assert changed.anchor_weight == pytest.approx(original.anchor_weight, rel=1e-13)
    for method in ("raw", "score", "anchored_score"):
        for mass in (.25, .7, 1.2):
            assert original.cdf(mass, method)+dual.cdf(1/mass, method) == pytest.approx(1., rel=0., abs=2e-9)
            assert changed.cdf((3/.7)*mass, method) == pytest.approx(original.cdf(mass, method), rel=0., abs=2e-9)


def test_closed_scatter_statistic_matches_independent_joint_score_reference():
    from met.joint_score_law import scalar_common_score_from_scatter
    case, count = witness_constants(), 317
    scatter = count*case["population_scatter"]
    q = _quad(count)
    masses = np.array([.03, case["false_mass"], .7, case["true_mass"], 3., 30.])
    expected = scalar_common_score_from_scatter(scatter, masses,
        covariance=case["covariance"], rank_rtol=0., rank_atol=0.)["statistic"]
    actual = np.array([q.components(math.atan(m))[1] for m in masses])
    np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=2e-25)


def test_anchor_does_not_change_the_local_godambe_limit():
    from met.joint_calibration import gaussian_joint_score_moments, joint_covariance_benchmark
    from met.joint_mass import common_acceleration
    case = witness_constants()
    event = common_acceleration([[0.]], [0.], case["covariance"])
    moments = gaussian_joint_score_moments([event], [case["true_mass"]], [[1.]])
    benchmark = joint_covariance_benchmark(moments["sensitivity"], moments["variability"])
    local_precision = 1/benchmark["sampling_covariance"][0, 0]
    h = .5
    expected = -.5*h*h*local_precision
    anchor_changes = []
    for count in (10000, 1000000):
        q = _quad(count)
        truth_angle = math.atan(case["true_mass"])
        shifted_angle = math.atan(case["true_mass"]+h/math.sqrt(count))
        log_score = math.log(q.kernel(shifted_angle, "score")/q.kernel(truth_angle, "score"))
        log_anchored = math.log(q.kernel(shifted_angle, "anchored_score")/q.kernel(truth_angle, "anchored_score"))
        assert log_score == pytest.approx(expected, rel=.02)
        anchor_changes.append(abs(log_anchored-log_score))
        expected_anchor_change = q.anchor_weight*.5*h*h*moments["sensitivity"][0, 0]
        assert anchor_changes[-1] == pytest.approx(expected_anchor_change, rel=.02)
    assert .08 < anchor_changes[1]/anchor_changes[0] < .12
    assert anchor_changes[-1] < .001


def test_conditional_density_determinant_has_a_nonintegrable_pole():
    # Y=(1,-1), unit covariance, theta=pi/4+delta:
    # U=sqrt(2)sin(delta), V=sqrt(2)cos(delta).
    def factor(delta):
        return math.exp(-math.cos(delta)**2)/(math.sqrt(2)*abs(math.sin(delta)))
    increments = []
    for cutoff in (1e-3, 1e-5, 1e-7):
        value = integrate.quad(factor, cutoff, 10*cutoff, epsabs=1e-12)[0]
        increments.append(value)
    assert increments[-1] == pytest.approx(math.exp(-1)*math.log(10)/math.sqrt(2), rel=1e-10)
    assert increments[-1] == pytest.approx(increments[-2], rel=1e-8)


def test_seeded_full_law_study_uses_standard_runner_and_saves_all_worlds(tmp_path):
    import yaml
    from pathlib import Path
    from met import preset, runner
    config = yaml.safe_load(Path("presets/score_law_sanity.yaml").read_text(encoding="utf-8"))
    config.update(replicates=2, counts=[100], population_counts=[100])
    path = tmp_path/"recipe.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    assert preset.load(path)["kind"] == "score_law_sanity"
    first, record = runner.run(path, tmp_path/"first", workers=1, progress=None)
    second, again = runner.run(path, tmp_path/"second", workers=1, progress=None)
    assert record == again
    assert len(record["summary"]) == 3
    assert (first/"summary.md").is_file()
    with np.load(first/"values.npz") as a, np.load(second/"values.npz") as b:
        assert a["observations_100"].shape == (2, 100, 2)
        np.testing.assert_array_equal(a["observations_100"], b["observations_100"])
    with pytest.raises(ValueError, match="serial"):
        runner.run(path, tmp_path/"invalid", workers=2, progress=None)
    config["counts"] = [100, 100]
    with pytest.raises(ValueError, match="distinct"):
        preset.validate(config)

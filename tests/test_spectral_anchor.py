"""Named spectral-anchor comparison, not an estimator or coverage ranking."""
import math

import numpy as np
import pytest
from scipy import stats

from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


CONTROLS = dict(prior_center=1., anchor_exponent=.5,
                epsabs=1e-11, epsrel=1e-10, limit=500)


def _path(count, *, schedule="spectral"):
    # A deterministic expected-scatter path with exact whitened eigenvalues
    # N and N+Delta. It checks a local experiment, not sample coverage.
    case = witness_constants()
    energy = count**.6
    squared = case["mean"]@np.linalg.solve(case["covariance"], case["mean"])
    scatter = count*case["covariance"]+energy/squared*np.outer(case["mean"], case["mean"])
    return ScalarScoreQuadrature(scatter, case["covariance"],
                                  anchor_schedule=schedule, **CONTROLS)


def test_historical_default_and_both_research_states_are_explicit():
    case = witness_constants()
    default = ScalarScoreQuadrature(100*case["population_scatter"], case["covariance"], **CONTROLS)
    explicit = ScalarScoreQuadrature(100*case["population_scatter"], case["covariance"],
                                     anchor_schedule="profile_oscillation", **CONTROLS)
    assert default.anchor_schedule == "profile_oscillation"
    assert default.anchor_weight == explicit.anchor_weight
    assert default.anchor_weight == pytest.approx((1+default.oscillation)**(-.5))
    assert default.spectral_anchor_weight == pytest.approx(1/math.sqrt(1+default.eigenvalues[-1]))
    assert default.anchor_state["spectral_prefactor"] == 1.
    assert default.anchor_state["anchor_weight_override"] is None
    assert default.anchor_state["effective_anchor_weight"] == default.derived_anchor_weight
    np.testing.assert_array_equal(default.eigenvalues, default.scatter_eigenvalues)
    assert default.gap == default.largest_eigenvalue-default.smallest_eigenvalue


def test_eigenbasis_profile_retains_order_one_loss_on_an_order_trillion_baseline():
    count, energy, theta = 1e12, 1e12**.6, math.pi/4
    direction = np.array([math.sin(theta), math.cos(theta)])
    scatter = count*np.eye(2)+energy*np.outer(direction, direction)
    q = ScalarScoreQuadrature(scatter, np.eye(2), anchor_schedule="spectral", **CONTROLS)
    assert abs(q.profile_log_factor(theta)) < 1e-20
    for h in (.5, 1., -1.):
        displacement = h/math.sqrt(energy)
        expected = -.5*energy*math.sin(displacement)**2
        assert q.profile_log_factor(theta+displacement) == pytest.approx(expected, rel=2e-9, abs=1e-12)
    # The positive chart contains the principal line and has its minimum
    # projection at a channel boundary, not at the negative minor line.
    assert q.oscillation == pytest.approx(energy/4, rel=2e-9)


def test_weak_identifying_path_separates_global_suppression_from_local_distortion():
    profile_ratios, spectral_ratios, penalties = [], [], []
    for count in (1e4, 1e8, 1e12):
        spectral = _path(count)
        profile = _path(count, schedule="profile_oscillation")
        energy = count**.6
        assert spectral.eigenvalues == pytest.approx([count, count+energy], rel=2e-12)
        assert spectral.gap == pytest.approx(energy, rel=2e-9)
        profile_ratios.append(profile.anchor_weight*spectral.largest_eigenvalue/spectral.gap)
        spectral_ratios.append(spectral.anchor_weight*spectral.largest_eigenvalue/spectral.gap)
        penalties.append(spectral.anchor_weight*spectral.gap)
    assert profile_ratios[0] < profile_ratios[1] < profile_ratios[2]
    assert profile_ratios[-1] > 20
    assert spectral_ratios[0] > spectral_ratios[1] > spectral_ratios[2]
    assert spectral_ratios[-1] < .07
    assert penalties[0] < penalties[1] < penalties[2]
    assert penalties[-1] > 15


def test_full_positive_weak_path_width_matches_the_local_sandwich_comparison():
    count, energy = 1e12, 1e12**.6
    case = witness_constants()
    precision = np.linalg.inv(case["covariance"])
    squared = case["mean"]@precision@case["mean"]
    # Derivative of whitened line angle with respect to physical mass.
    speed = math.sqrt(np.linalg.det(precision))/squared
    reference = stats.norm.ppf(.975)*math.sqrt(count+energy)/energy/speed
    ratios = {}
    for schedule in ("profile_oscillation", "spectral"):
        q = _path(count, schedule=schedule)
        lower, upper = (q.quantile(p, "anchored_score") for p in (.025, .975))
        ratios[schedule] = (upper-lower)/(2*reference)
    assert .15 < ratios["profile_oscillation"] < .25
    assert .95 < ratios["spectral"] < 1.08


def test_explicit_weight_override_is_recorded_and_uses_an_independent_cache():
    case = witness_constants()
    q = ScalarScoreQuadrature(100*case["population_scatter"], case["covariance"],
                              anchor_schedule="spectral", anchor_weight_override=0., **CONTROLS)
    assert q.anchor_weight == q.anchor_state["effective_anchor_weight"] == 0.
    assert q.anchor_state["anchor_weight_override"] == 0.
    assert q.derived_anchor_weight == q.spectral_anchor_weight > 0
    for mass in (.2, .6, 1.5):
        assert q.cdf(mass, "anchored_score") == q.cdf(mass, "score")
    assert q.total("anchored_score") == q.total("score")


def test_spectral_state_and_full_law_respect_units_and_reciprocal_coordinates():
    case = witness_constants()
    scatter, covariance = 100*case["population_scatter"], case["covariance"]
    original = ScalarScoreQuadrature(scatter, covariance, anchor_schedule="spectral", **CONTROLS)
    scale = 7.5
    transform = np.diag([3., .4])
    converted = ScalarScoreQuadrature(transform@scatter@transform.T,
        transform@covariance@transform.T, prior_center=scale,
        anchor_exponent=.5, anchor_schedule="spectral", epsabs=1e-11, epsrel=1e-10, limit=500)
    swap = np.array([[0., 1.], [1., 0.]])
    reciprocal = ScalarScoreQuadrature(swap@scatter@swap, swap@covariance@swap,
                                       anchor_schedule="spectral", **CONTROLS)
    for q in (converted, reciprocal):
        np.testing.assert_allclose(q.eigenvalues, original.eigenvalues, rtol=3e-14)
        assert q.anchor_weight == pytest.approx(original.anchor_weight, rel=2e-14)
        assert q.oscillation == pytest.approx(original.oscillation, rel=2e-14)
    for mass in (.2, .6, 1.5):
        expected = original.cdf(mass, "anchored_score")
        assert converted.cdf(scale*mass, "anchored_score") == pytest.approx(expected, abs=2e-11)
        assert reciprocal.cdf(1/mass, "anchored_score") == pytest.approx(1-expected, abs=2e-11)


def test_no_count_input_and_exact_zero_scatter_addition_preserve_the_state():
    data = np.array([[2., .7], [.4, 1.1], [-.3, .5], [1.2, -.8]])
    scatter = data.T@data
    padded = np.vstack([data, [0., 0.]])
    first = ScalarScoreQuadrature(scatter, np.eye(2), anchor_schedule="spectral", **CONTROLS)
    second = ScalarScoreQuadrature(padded.T@padded, np.eye(2), anchor_schedule="spectral", **CONTROLS)
    assert first.anchor_state == second.anchor_state
    assert first.total("anchored_score") == second.total("anchored_score")


def test_flat_profile_remains_the_supplied_prior_under_either_schedule():
    for schedule in ("profile_oscillation", "spectral"):
        q = ScalarScoreQuadrature(1000*np.eye(2), np.eye(2), anchor_schedule=schedule, **CONTROLS)
        assert q.gap == q.oscillation == 0
        for method in ("raw", "score", "anchored_score"):
            assert q.cdf(2., method) == pytest.approx(2*math.atan(2.)/math.pi, abs=2e-15)


def test_explicit_schedule_override_and_numerical_control_validation():
    for override in (-1., math.nan, math.inf):
        with pytest.raises(ValueError):
            ScalarScoreQuadrature(np.eye(2), np.eye(2), anchor_weight_override=override, **CONTROLS)
    with pytest.raises(ValueError):
        ScalarScoreQuadrature(np.eye(2), np.eye(2), anchor_schedule="unselected", **CONTROLS)
    invalid = dict(CONTROLS, epsabs=math.nan)
    with pytest.raises(ValueError):
        ScalarScoreQuadrature(np.eye(2), np.eye(2), **invalid)

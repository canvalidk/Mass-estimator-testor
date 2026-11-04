"""T63: trace dilution under degenerating, correctly declared Gaussian geometry.

The continuum calculation below isolates the limiting informative field;
it is not a simulation of the full heterogeneous history. Finite actual
heterogeneous score fields and independent Gaussian moments are checked
separately. No new estimator or tuning parameter is adopted.
"""
import math

import numpy as np
import pytest

from met.joint_calibration import gaussian_joint_score_moments
from met.joint_mass import common_acceleration, gaussian_profile_cost
from met.joint_score_law import (gaussian_joint_score_statistic,
                                 scalar_common_score_from_scatter)
from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _bad_coefficients(index, mass):
    epsilon = 4.**(-index)
    x = epsilon*mass
    denominator = 1+x*x
    projector = np.array([[x*x, x], [x, 1.]])/denominator
    first = epsilon*np.array([[2*x, 1-x*x], [1-x*x, -2*x]])/denominator**2
    second = (2*epsilon**2*np.array([[1-3*x*x, -x*(3-x*x)],
                                    [-x*(3-x*x), -1+3*x*x]])/denominator**3)
    conditional = epsilon**2*np.outer([x, 1.], [x, 1.])/denominator**3
    return epsilon, projector, first, second, conditional


def test_bad_event_candidate_fields_vanish_uniformly_on_compact_mass_intervals():
    null_projector = np.diag([0., 1.])
    for index in (1, 3, 7, 15):
        for mass in np.geomspace(.05, 5., 31):
            epsilon, projector, first, second, conditional = _bad_coefficients(index, mass)
            denominator = 1+(epsilon*mass)**2
            assert np.linalg.norm(projector-null_projector, "fro") <= math.sqrt(2)*5*epsilon
            assert np.linalg.norm(first, "fro") == pytest.approx(
                math.sqrt(2)*epsilon/denominator, rel=5e-15)
            assert np.linalg.norm(second, "fro") == pytest.approx(
                2*math.sqrt(2)*epsilon**2/denominator**1.5, rel=5e-15)
            assert np.linalg.norm(conditional, 2) == pytest.approx(
                epsilon**2/denominator**2, rel=5e-15)
            assert np.linalg.norm(conditional, 2) <= epsilon**2*(1+5e-15)


def test_bad_event_random_fields_have_summable_variance_but_noise_trace_does_not():
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    for index in (1, 3, 7):
        for mass in (.2, 1.2, 5.):
            epsilon, projector, first, _, conditional = _bad_coefficients(index, mass)
            denominator = 1+(epsilon*mass)**2
            sums = np.zeros(7)
            for j, k in np.ndindex(3, 3):
                probability = weights[j]*weights[k]/math.pi
                z = math.sqrt(2)*np.array([nodes[j], nodes[k]])
                # Subtract the constant null-candidate loss analytically,
                # retaining tiny represented contrasts at degenerate geometry.
                contrast = (epsilon*mass*z[0]*z[1]
                    +.5*(epsilon*mass)**2*(z[0]**2-z[1]**2))/denominator
                score = .5*z@first@z
                covariance = z@conditional@z
                energy = z@z
                sums += probability*np.array([contrast, contrast**2, score**2,
                                               covariance, covariance**2, energy, energy**2])
            assert sums[0] == pytest.approx(0., abs=2e-17)
            assert sums[1] == pytest.approx((epsilon*mass)**2/denominator, rel=3e-14)
            assert sums[2] == pytest.approx(epsilon**2/denominator**2, rel=3e-14)
            assert sums[3] == pytest.approx(epsilon**2/denominator**2, rel=3e-14)
            assert sums[4]-sums[3]**2 == pytest.approx(
                2*epsilon**4/denominator**4, rel=3e-14)
            np.testing.assert_allclose(sums[5:], [2., 8.], rtol=3e-15)
    # Sum over every positive index bounds the bad subsequence as well.
    assert sum(16.**(-i) for i in range(1, 40)) == pytest.approx(1/15, rel=3e-15)
    assert sum(256.**(-i) for i in range(1, 40)) == pytest.approx(1/255, rel=3e-15)
    # The per-event mean trace is always two and variance four, regardless
    # of the rapidly vanishing candidate sensitivity just verified.
    assert 2*4096 > 1000*(2/15)


def _physical_scalar_fields(point, covariance, mass):
    channel = np.array([mass, 1.])
    constraint = np.array([1., -mass])
    precision_channel = np.linalg.solve(covariance, channel)
    fitted = point@precision_channel/(channel@precision_channel)
    variance = constraint@covariance@constraint
    residual = constraint@point
    return fitted*residual/variance, fitted*fitted/variance, residual*residual/variance


def test_actual_heterogeneous_history_aggregate_fields_match_physical_residuals():
    case = witness_constants()
    events = []
    for index in range(1, 13):
        noise = np.array([math.sin(index), math.cos(.7*index)])
        if index > 1 and index & (index-1) == 0:
            covariance = case["covariance"]
            point = case["mean"]+np.linalg.cholesky(covariance)@noise
        else:
            covariance = np.diag([16.**index, 1.])
            point = np.array([4.**index, 1.])*noise
        events.append(common_acceleration([[point[0]]], [point[1]], covariance))
    for mass in (.03, case["false_mass"], .6, case["true_mass"], 5.):
        reference = np.sum([_physical_scalar_fields(event.data, event.covariance, mass)
                            for event in events], axis=0)
        actual = gaussian_joint_score_statistic(events, [mass], **RANK)
        assert actual["score"][0] == pytest.approx(reference[0], rel=5e-13, abs=3e-14)
        assert actual["covariance"][0, 0] == pytest.approx(reference[1], rel=5e-13, abs=3e-14)
        assert actual["statistic"] == pytest.approx(reference[0]**2/reference[1],
                                                    rel=8e-13, abs=3e-14)
        assert sum(gaussian_profile_cost(event, [mass]) for event in events) == pytest.approx(
            reference[2], rel=5e-13, abs=3e-14)


def test_good_information_diverges_while_the_full_trace_anchor_loses_global_ordering():
    case = witness_constants()
    declaration = common_acceleration([[0.]], [0.], case["covariance"])
    moments = gaussian_joint_score_moments([declaration], [case["true_mass"]], [[1.]])
    a = moments["sensitivity"][0, 0]
    b = moments["variability"][0, 0]
    previous_information, previous_anchor_gap = 0., math.inf
    for good_count in (8, 32, 128, 512):
        # Good event locations are 2,4,8,...,2^M; all other events
        # have the degenerate null geometry. This is one fixed history.
        total_count = 2.**good_count
        expected_trace = 2*total_count+(28/19)*good_count
        weight = 1/math.sqrt(1+expected_trace)
        information = good_count*a*a/b
        anchor_gap = weight*good_count*(28/19)/2
        assert information > previous_information
        assert anchor_gap < previous_anchor_gap
        previous_information, previous_anchor_gap = information, anchor_gap
    assert information > 100
    assert anchor_gap < 1e-70


def test_limiting_informative_field_retains_false_laplace_mass_when_trace_dilutes_anchor():
    case = witness_constants()
    split = math.sqrt(case["true_mass"]*case["false_mass"])
    errors = []
    for good_count in (400, 1000):
        total_count = 2.**good_count
        expected_trace = 2*total_count+(28/19)*good_count
        q = ScalarScoreQuadrature(good_count*case["population_scatter"], case["covariance"],
            prior_center=1., anchor_exponent=.5, anchor_schedule="spectral",
            anchor_weight_override=1/math.sqrt(1+expected_trace),
            epsabs=1e-11, epsrel=1e-10, limit=500)
        # This explicitly decoupled limiting-field check drops the bounded
        # bad-event random perturbation. It verifies the two-root Laplace
        # mechanism and dilution rate, not finite-history Monte Carlo coverage.
        fields = scalar_common_score_from_scatter(good_count*case["population_scatter"],
            np.array([case["false_mass"], case["true_mass"]]),
            covariance=case["covariance"], rank_rtol=0., rank_atol=0.)
        np.testing.assert_allclose(fields["statistic"], 0., atol=2e-24)
        assert q.cdf(split, "anchored_score") == pytest.approx(q.cdf(split, "score"), abs=2e-10)
        # T63 fixes this reciprocal truncation before data. Conditioning
        # scales both root weights by one common prior normalizer.
        low, high = (q.cdf(mass, "anchored_score") for mass in (.1, 10.))
        truncated_false = (q.cdf(split, "anchored_score")-low)/(high-low)
        assert truncated_false == pytest.approx(case["false_limit_weight"], abs=.005)
        errors.append(abs(q.cdf(split, "anchored_score")-case["false_limit_weight"]))
    assert errors[0] > errors[1]
    assert errors[-1] < .002
    assert q.quantile(.5, "anchored_score") < .3
    # The declared good subsequence still identifies the mass on its own.
    assert abs(q.quantile(.5, "raw")-case["true_mass"]) < .05

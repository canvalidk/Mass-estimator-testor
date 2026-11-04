"""T65: full-prior disagreement versus compact-region mass information.

The final continuum check isolates the informative population field and
its diluted weight. It is not a full heterogeneous-history simulation,
coverage study, global-mixture calculation or median-limit assertion.
"""
import math

import numpy as np
import pytest
from scipy import integrate

from met.mass_contrast_anchor import scalar_prior_projector_contrast
from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


def _closed_mean(epsilon):
    a = epsilon/(1+epsilon)
    b = (1/math.pi if epsilon == 1 else
         2*epsilon*math.log(1/epsilon)/(math.pi*(1-epsilon*epsilon)))
    return np.array([[a, b], [b, 1-a]])


def _closed_variance(epsilon):
    mean = _closed_mean(epsilon)
    return float(mean[0, 0]*mean[1, 1]-mean[0, 1]**2)


@pytest.mark.parametrize("epsilon", [.0625, .2, 1., 5., 16.])
def test_full_half_cauchy_projector_average_has_the_exact_tail_formula(epsilon):
    # An independent uniform half-Cauchy angle integration computes the
    # physical graph projector directly, without the closed mean formula.
    def entry(theta, row, column):
        vector = np.array([epsilon*math.sin(theta), math.cos(theta)])
        return 2/math.pi*vector[row]*vector[column]/(vector@vector)

    numerical = np.array([[integrate.quad(entry, 0., math.pi/2, args=(row, column),
        epsabs=1e-12, epsrel=1e-12)[0] for column in range(2)] for row in range(2)])
    expected = _closed_mean(epsilon)
    np.testing.assert_allclose(numerical, expected, rtol=3e-13, atol=3e-14)
    variance = _closed_variance(epsilon)
    np.testing.assert_allclose(numerical-numerical@numerical,
                               variance*np.eye(2), rtol=4e-13, atol=3e-14)
    assert 0 < variance <= .25
    assert _closed_variance(1/epsilon) == pytest.approx(variance, rel=3e-14)
    reference = scalar_prior_projector_contrast(np.diag([epsilon**-2, 1.]),
        prior_center=1., mass_bounds=[0., math.inf],
        epsabs=1e-12, epsrel=1e-12, limit=500)
    assert reference["coefficient"] == pytest.approx(variance, rel=4e-12, abs=3e-14)


def test_prior_tail_disagreement_has_harmonic_mean_growth_with_finite_noise_variance():
    previous_ratio = 0.
    for count in (100, 1000, 10000, 100000):
        indices = np.arange(1, count+1, dtype=float)
        epsilon = 1/indices
        a = epsilon/(1+epsilon)
        b = np.empty(count)
        b[0] = 1/math.pi
        b[1:] = (2*epsilon[1:]*np.log(indices[1:])
                 /(math.pi*(1-epsilon[1:]**2)))
        weights = a*(1-a)-b*b
        # Remove exactly the good indices 4,16,256,... in the theorem.
        good = np.array([2**(2**j) for j in range(1, 5)])
        weights[good[good <= count]-1] = 0.
        expected_state = 2*weights.sum()
        variance_state = 4*(weights@weights)
        ratio = expected_state/(2*math.log(count))
        assert previous_ratio < ratio < 1.
        previous_ratio = ratio
        assert 0 < 2*math.log(count)-expected_state < 6.
        # This Gaussian quadratic-state variance stays bounded as its
        # mean diverges, giving relative concentration around 2 log N.
        assert variance_state < 1.
        # Local conditional-C information under any fixed candidate has
        # a summable i^-2 envelope despite that state divergence.
        assert np.sum(epsilon*epsilon) < math.pi**2/6


def test_local_null_fields_have_a_square_summable_leading_noise_and_remainder():
    exchange = np.array([[0., 1.], [1., 0.]])
    base = np.diag([0., 1.])
    for index in (2, 5, 20, 1000):
        epsilon = 1/index
        for mass in (.1, .3, 1.2, 3.):
            x = epsilon*mass
            denominator = 1+x*x
            projector = np.outer([x, 1.], [x, 1.])/denominator
            derivative = epsilon*np.array([[2*x, 1-x*x],
                                          [1-x*x, -2*x]])/denominator**2
            conditional = epsilon**2*np.outer([x, 1.], [x, 1.])/denominator**3
            assert np.linalg.norm(projector-base-epsilon*mass*exchange, 2) <= 12*epsilon**2
            assert np.linalg.norm(derivative-epsilon*exchange, 2) <= 12*epsilon**2
            assert np.linalg.norm(conditional, 2) <= epsilon**2*(1+3e-15)
    # Leading score terms are xi*zeta/i: their independent variances
    # sum to pi^2/6; all remaining quadratic coefficients are O(i^-2).
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    variance_product = sum(weights[j]*weights[k]/math.pi
        *(2*nodes[j]*nodes[k])**2 for j, k in np.ndindex(3, 3))
    assert variance_product == pytest.approx(1., abs=2e-15)


@pytest.mark.parametrize("gamma", [.125, .5, 2.])
def test_every_fixed_positive_power_has_the_same_sparse_history_rate_obstruction(gamma):
    previous = math.inf
    for good_count in (16, 64, 256, 1000):
        # At N=2^(2^M), log N=(log 2)*2^M. Work in logs to avoid
        # representing that unneeded doubly exponential reading count.
        log_state = math.log(2*math.log(2))+good_count*math.log(2)
        weight_times_gap = good_count*math.exp(-gamma*log_state)
        assert weight_times_gap < previous
        previous = weight_times_gap
    assert previous < 1e-30


def test_local_false_to_true_laplace_ratio_survives_the_diluted_anchor():
    case = witness_constants()
    expected_ratio = case["false_limit_weight"]/(1-case["false_limit_weight"])
    errors = []
    for good_count in (400, 1000):
        log_state = math.log(2*math.log(2))+good_count*math.log(2)
        q = ScalarScoreQuadrature(good_count*case["population_scatter"], case["covariance"],
            prior_center=1., anchor_exponent=.5, anchor_schedule="spectral",
            anchor_weight_override=math.exp(-.5*log_state),
            epsabs=1e-11, epsrel=1e-10, limit=500)
        # Ratios of unnormalized basin integrals cancel any other global
        # posterior mass. This deliberately decoupled field drops the
        # bounded random null perturbation; it checks the local mechanism.
        false = q.cdf(.5, "anchored_score")-q.cdf(.1, "anchored_score")
        true = q.cdf(1.8, "anchored_score")-q.cdf(.8, "anchored_score")
        errors.append(abs(false/true-expected_ratio))
        assert false/true == pytest.approx(
            (q.cdf(.5, "score")-q.cdf(.1, "score"))
            /(q.cdf(1.8, "score")-q.cdf(.8, "score")), abs=2e-10)
    assert errors[0] > errors[1]
    assert errors[-1] < .04
    assert 1/(1+expected_ratio) < .25


def test_tail_threshold_example_does_not_assume_a_monotone_schedule():
    def schedule(value):
        return (2+math.sin(value))/math.sqrt(1+value)

    # This positive vanishing schedule increases across these points;
    # threshold selection must control its whole tail, not its value at
    # one threshold. Its envelope is 3/sqrt(1+u).
    assert schedule(2*math.pi+.05) > schedule(2*math.pi)
    previous_log_index = 0.
    for good_count in (1, 2, 3, 4, 8, 16):
        threshold = 9*good_count**4-1
        log_index = max(threshold, good_count**2, previous_log_index+math.log(2))
        for offset in (0., .7, math.pi, 10.):
            state = log_index+offset
            assert schedule(state)*good_count <= (1/good_count)*(1+3e-15)
        previous_log_index = log_index

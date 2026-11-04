"""T67 independent continuous full-prior integrals and normalized-law checks.

The n-scaled population field isolates the limiting mechanism. It is not
a Gaussian sample-history simulation or a coverage study.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.pairwise_anchor import ScalarPairwiseQuadrature
from met.score_law_obstruction import witness_constants


def _reference(count):
    case = witness_constants()
    return ScalarPairwiseQuadrature(count*case["population_scatter"], case["covariance"],
        prior_center=1., anchor_prefactor=1., epsabs=1e-10, epsrel=1e-9, limit=400)


@pytest.mark.parametrize("mass", [.26666666666666666, .7, 1.2, 3.])
def test_pairwise_continuous_penalty_matches_physical_projector_quadrature(mass):
    q, case = _reference(30), witness_constants()
    covariance = case["covariance"]
    w = np.linalg.inv(np.linalg.cholesky(covariance))
    scatter = w@(30*case["population_scatter"])@w.T
    def projector(m):
        v = w@np.array([m, 1.])
        return np.outer(v, v)/(v@v)
    candidate = projector(mass)
    def integrand(rank):
        rival = projector(math.tan(math.pi*rank/2))
        delta = rival-candidate
        gain = np.trace(delta@scatter)/2
        variance = np.trace(delta@delta@scatter)
        return max(gain, 0.)/math.sqrt(1+variance)
    expected = integrate.quad(integrand, 0., 1., epsabs=2e-9, epsrel=2e-9,
                              points=[2/math.pi*math.atan(mass)], limit=500)[0]
    assert q.penalty(math.atan(mass))[0] == pytest.approx(expected, abs=1e-8)


def test_false_root_ordering_and_local_penalty_rates_with_full_half_cauchy():
    false, local = [], []
    for count in (100, 400, 1600):
        q = _reference(count)
        assert q.penalty(math.atan(1.2))[0] < 1e-12
        false.append(q.penalty(math.atan(4/15))[0])
        local.append(q.penalty(math.atan(1.2+1/math.sqrt(count)))[0])
    np.testing.assert_allclose(np.array(false)/np.sqrt([100,400,1600]), .211,
                               rtol=.01)
    assert local[2] < local[1] < local[0] < .02


def test_full_normalization_suppresses_false_root_and_keeps_godambe_window():
    probabilities, errors = [], []
    case = witness_constants()
    m0, covariance = 1.2, case["covariance"]
    # Independent physical sensitivity A and graph contribution G for a=1.
    c = np.array([m0, 1.])
    omega = np.linalg.inv(covariance)
    residual_precision = omega-np.outer(omega@c, omega@c)/(c@omega@c)
    sensitivity = residual_precision[0,0]
    variability = sensitivity+sensitivity/(c@omega@c)
    for count in (100, 400, 1600):
        q = _reference(count)
        probabilities.append(q.cdf(.6, "pairwise"))
        width = math.sqrt(variability/count)/sensitivity
        central = q.cdf(m0+width, "pairwise")-q.cdf(m0-width, "pairwise")
        errors.append(abs(central-(stats.norm.cdf(1)-stats.norm.cdf(-1))))
    assert probabilities[2] < probabilities[1] < probabilities[0]
    assert probabilities[-1] < .001
    assert errors[-1] < .015
    assert errors[-1] < errors[0]

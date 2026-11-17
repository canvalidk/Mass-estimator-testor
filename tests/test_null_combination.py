"""Finite ingredients for E2's null-combination obstruction and scoped control.

No simulations: Gaussian quadrature, exact mixture probabilities and positive
scatter fixtures. Infinite-array and history quantifiers are proved in
explore/e2-null-combination.md, not inferred from finite tests.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats


@pytest.mark.parametrize("radius", [.15, .7, 1.2])
@pytest.mark.parametrize("external_mean", [-1.3, 2.1])
def test_joint_likelihood_ratio_second_moment_with_unchanged_external_experiment(
    radius, external_mean
):
    nodes, weights = np.polynomial.hermite_e.hermegauss(48)
    weights /= math.sqrt(2 * math.pi)
    x = external_mean + .8 * nodes[:, None]
    y = nodes[None, :]
    external = stats.norm.pdf(x, external_mean, .8)
    reference = external * stats.norm.pdf(y)
    alternative = external * (
        stats.norm.pdf(y, radius) + stats.norm.pdf(y, -radius)
    ) / 2
    ratio = alternative / reference
    second = np.einsum("i,j,ij->", weights, weights, ratio**2)
    first = np.einsum("i,j,ij->", weights, weights, ratio)
    assert first == pytest.approx(1, abs=2e-13)
    assert second == pytest.approx(math.cosh(radius**2), rel=3e-13)


@pytest.mark.parametrize("radius,threshold,coefficient", [
    (.1, 3., .4), (.7, 4., 1.2), (1.4, 2., -.6)
])
def test_event_bound_for_event_that_depends_on_both_experiments(
    radius, threshold, coefficient
):
    # A={X+kY>threshold}; X has a fixed informative Gaussian law in both
    # worlds. Integrate its conditional event probability independently.
    mean, sd = .3, .65
    conditional = lambda y: stats.norm.sf((threshold-coefficient*y-mean)/sd)
    alpha = integrate.quad(lambda y: conditional(y)*stats.norm.pdf(y),
                           -12, 12, epsabs=1e-13)[0]
    mixed = integrate.quad(
        lambda y: conditional(y) * (
            stats.norm.pdf(y-radius)+stats.norm.pdf(y+radius)
        )/2, -12, 12, epsabs=1e-13
    )[0]
    total_sd = math.hypot(sd, coefficient)
    assert alpha == pytest.approx(stats.norm.sf((threshold-mean)/total_sd), abs=1e-12)
    assert mixed <= math.sqrt(alpha*math.cosh(radius**2)) + 1e-13
    assert math.cosh(radius**2) <= math.exp(radius**4/2)


@pytest.mark.parametrize("n,log_inverse_alpha", [(100, 3.), (100, 20.), (10000, 15.), (10000, 200.)])
def test_slow_alternative_has_exact_energy_and_the_claimed_finite_bound(n, log_inverse_alpha):
    d2 = min(log_inverse_alpha, math.sqrt(n))
    delta = math.sqrt(n*d2)
    radius_squared = delta/n
    information = delta**2/(n+delta)
    assert information == pytest.approx(d2/(1+math.sqrt(d2/n)))
    exact_log_bound = -.5*log_inverse_alpha + .5*n*math.log(math.cosh(radius_squared))
    upper_log_bound = -.5*log_inverse_alpha + delta**2/(4*n)
    assert exact_log_bound <= upper_log_bound + 1e-12
    assert upper_log_bound <= -.25*log_inverse_alpha + 1e-12


@pytest.mark.parametrize("weight,block_point", [(0., 3.), (.4, 1.), (4., 30.)])
def test_inverse_weight_concentration_bound_with_positive_external_point(weight, block_point):
    # B lognormal; its density is bounded. Thus Q_B(h)<=h*sup f_B.
    # Check exact conditional success probability against this bound.
    w0, mass, tolerance = 2., 2., .01
    mu, sigma = math.log(mass), .35
    baseline = stats.lognorm(s=sigma, scale=math.exp(mu))
    density_sup = baseline.pdf(math.exp(mu-sigma*sigma))
    a = w0/(w0+weight)
    centre = (1-a)*block_point
    lower, upper = (mass-tolerance-centre)/a, (mass+tolerance-centre)/a
    success = baseline.cdf(upper)-baseline.cdf(lower)
    assert success <= 2*tolerance*density_sup/a + 1e-14
    assert 1-success >= 1-2*tolerance*density_sup/a - 1e-14


def test_vanishing_weight_does_not_by_itself_make_the_point_neutral():
    baseline, w0 = 2., 4.
    n = np.array([10., 100., 1000., 10000.])
    weight, point = 1/n, n
    combined = (w0*baseline+weight*point)/(w0+weight)
    assert combined[-1]-baseline == pytest.approx(1/w0, abs=6e-5)
    assert weight[-1] < weight[0]/100


@pytest.mark.parametrize("signals", [[0., 0.], [.2, -1.4], [0., 3.]])
def test_gram_control_exact_gaussian_moments_by_quadrature(signals):
    nodes, weights = np.polynomial.hermite_e.hermegauss(10)
    weights /= math.sqrt(2*math.pi)
    u = np.array([.6, .8])
    noise = np.stack(np.meshgrid(nodes, nodes, indexing="ij"), -1).reshape(-1, 2)
    quadrature = np.outer(weights, weights).ravel()
    matrix_variance, trace_variance = 0., 0.
    for r in signals:
        data = r*u + noise
        outer = np.einsum("ni,nj->nij", data, data)
        expected = np.eye(2)+r*r*np.outer(u, u)
        assert np.einsum("n,nij->ij", quadrature, outer) == pytest.approx(expected, abs=1e-12)
        matrix_variance += np.dot(quadrature, ((outer-expected)**2).sum(axis=(1, 2)))
        trace_variance += np.dot(quadrature, (np.trace(outer, axis1=1, axis2=2)-(2+r*r))**2)
    n, delta = len(signals), np.dot(signals, signals)
    assert matrix_variance == pytest.approx(6*(n+delta), abs=1e-11)
    assert trace_variance == pytest.approx(4*(n+delta), abs=1e-11)


def test_gram_control_zero_padding_bounds_and_deterministic_null_scale():
    # A deterministic scatter sequence verifies the CLT substitution, not
    # a sampling limit. Eigenvalues stay positive at every checked N.
    fluctuation = np.array([[1.1, -.7], [-.7, -.3]])
    target = ((fluctuation[0, 0]-fluctuation[1, 1])**2+4*fluctuation[0, 1]**2)/2
    errors = []
    for n in (100, 10000, 1000000):
        gram = n*np.eye(2)+math.sqrt(n)*fluctuation
        data = np.linalg.cholesky(gram).T
        padded = np.vstack([data[::-1], np.zeros((2, 2))])
        assert padded.T@padded == pytest.approx(gram, rel=1e-14, abs=1e-10)
        eigenvalues = np.linalg.eigvalsh(padded.T@padded)
        energy = eigenvalues.sum()
        weight = (eigenvalues[1]-eigenvalues[0])**2/energy
        assert 0 <= weight <= energy
        errors.append(abs(weight-target))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < .001

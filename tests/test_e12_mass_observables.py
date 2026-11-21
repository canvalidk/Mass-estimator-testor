"""Two-body Gaussian observables, exact innovation cost and finite ray odds."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from explore.e12_tangent_pairings import tangent_pairings
from met.full import Reading


def profiled_cost(y, covariance, graph):
    root = np.linalg.cholesky(covariance)
    whitened_y = np.linalg.solve(root, y)
    whitened_line = np.linalg.solve(root, graph)
    drive = (whitened_y @ whitened_line) / (whitened_line @ whitened_line)
    residual = whitened_y - drive * whitened_line
    return residual @ residual, drive


def test_collapse_of_same_total_depends_on_force_sensor_rows():
    masses, alternative = np.array([1., 3.]), np.array([2., 2.])
    for sensor, expected in [(np.eye(2), 2.), (np.ones((1, 2)), 0.)]:
        mean = np.r_[sensor @ masses, 1.]
        graph = np.r_[sensor @ alternative, 1.]
        twice_kl, nearest_drive = profiled_cost(mean, np.eye(mean.size), graph)
        assert nearest_drive == pytest.approx(1.)
        assert twice_kl == pytest.approx(expected)
        assert norm.cdf(-math.sqrt(max(0., twice_kl)) / 2) == pytest.approx(
            0.2397500611 if expected else 0.5
        )


@pytest.mark.parametrize("r,T", [(.17, .8), (.5, 1.9), (.73, 2.8)])
def test_innovation_sum_matches_the_full_profiled_cost_and_force_ratio_residual(r, T):
    root = np.array([[1.3, 0., 0.], [.2, .7, 0.], [-.1, .4, 1.1]])
    covariance = root @ root.T
    y = np.array([1.7, -.35, 1.1])
    contrasts = np.array([[1., 1., -T], [1-r, -r, 0.]])
    total, ray = contrasts @ y
    V = contrasts @ covariance @ contrasts.T
    p, c, q = V[0, 0], V[0, 1], V[1, 1]
    forward = total**2/p + (ray-c*total/p)**2/(q-c**2/p)
    reverse = ray**2/q + (total-c*ray/q)**2/(p-c**2/q)
    true_profile, _ = profiled_cost(y, covariance, np.array([T*r, T*(1-r), 1.]))
    assert true_profile == pytest.approx(forward)
    assert true_profile == pytest.approx(reverse)
    assert abs(forward - (total**2/p + ray**2/q)) > .01
    # Force ratio m1/m2 uses the two physical force channels with their own covariance.
    reading = Reading([y[0]], [y[1]], covariance[:2, :2])
    parts = tangent_pairings([reading], r/(1-r))
    assert parts.within + parts.between + parts.common == pytest.approx(ray**2/q)


def test_two_individual_candidate_statistics_are_both_chi_square_but_correlated():
    T, r, a = 2., .25, 1.3
    graph = np.array([T*r, T*(1-r), 1.])
    contrasts = np.array([[1., 1., -T], [1-r, -r, 0.]])
    np.testing.assert_allclose(contrasts @ (a*graph), 0, atol=1e-14)
    noise = contrasts @ contrasts.T
    np.testing.assert_allclose(noise, [[6, .5], [.5, .625]])
    innovations = np.array([[1, 0], [-noise[0, 1]/noise[0, 0], 1]])
    assert (innovations @ noise @ innovations.T)[0, 1] == pytest.approx(0.)


def test_smooth_common_positive_drive_priors_reverse_the_finite_bayes_order():
    y = np.array([1., 2., 1.])
    first, second = np.array([.5, 1.5, 1.]), np.array([1., 1., 1.])

    def likelihood(amp, graph):
        delta = y - amp*graph
        return math.exp(-0.5*(delta @ delta))

    for amp, expected in [(1., .25), (3., -.75)]:
        assert math.log(likelihood(amp, first)/likelihood(amp, second)) == pytest.approx(expected)
        # A common lognormal nuisance prior. Integrate in standardized log-drive coordinates.
        def integrated(graph):
            return quad(lambda u: norm.pdf(u)*likelihood(amp*math.exp(.01*u), graph), -12, 12)[0]
        log_odds = math.log(integrated(first)/integrated(second))
        assert log_odds * expected > 0
        assert log_odds == pytest.approx(expected, abs=.02)


@pytest.mark.parametrize("corr", [0., .35])
def test_known_total_force_exact_quadratic_likelihood_and_factorization(corr):
    y = np.array([[.6, 1.9, 1.15], [-.3, -1.1, -.6]])
    settings = [2., -1.3]
    covariance = np.array([[1.2, .3, corr], [.3, .9, 0.], [corr, 0., .7]])
    P = np.linalg.inv(covariance)
    design = np.array([[1., 0.], [-1., 0.], [0., 1.]])
    J = np.zeros((2, 2))
    z = np.zeros(2)
    baseline = 0.
    for row, s in zip(y, settings):
        D = s*design
        y0 = row - [0., s, 0.]
        J += D.T @ P @ D
        z += D.T @ P @ y0
        baseline -= .5 * (y0 @ P @ y0)
    assert np.linalg.eigvalsh(J).min() > 0
    assert bool(abs(J[0, 1]) < 1e-10) == (corr == 0.)

    def loglik(r, inv_total):
        value = 0.
        for row, s in zip(y, settings):
            residual = row - s*np.array([r, 1-r, inv_total])
            value -= .5*(residual @ P @ residual)
        beta = np.array([r, inv_total])
        assert value == pytest.approx(baseline + beta @ z - .5*(beta @ J @ beta))
        return value

    difference = loglik(.8, .7)-loglik(.8, .4)-loglik(.2, .7)+loglik(.2, .4)
    assert difference == pytest.approx(-(.8-.2)*(.7-.4)*J[0, 1])


def test_sharp_drive_bounds_and_an_interior_bayes_factor():
    y = np.array([1., 2., 1.])
    rays = np.array([[.5, 1.5, 1.], [1., 1., 1.]])
    drives = np.array([1., 3.])
    residuals = y - drives[:, None, None] * rays
    kernels = np.exp(-.5 * np.sum(residuals**2, axis=2))
    ratios = kernels[:, 0] / kernels[:, 1]
    # On the full positive drive range log R(a) = 1/4 - (a-1)^2/4.
    np.testing.assert_allclose(ratios, np.exp(.25 - .25*(drives-1)**2))
    target = .95
    weight = (target - ratios[1]) / (ratios[0] - ratios[1])
    prior = np.array([weight, 1-weight]) / kernels[:, 1]
    prior /= prior.sum()
    evidence = prior @ kernels
    assert evidence[0] / evidence[1] == pytest.approx(target)

"""Covariance-dependent positive-mass geometry, independent of the candidate law."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from met.full import Reading, reading_state


def cost_on_arc(matrix, beta):
    values, vectors = np.linalg.eigh(matrix @ matrix.T)
    line = math.atan2(vectors[1, -1], vectors[0, -1]) % math.pi
    distance = max(0.0, min(line - beta, math.pi - line))
    return values[0] + (values[1] - values[0]) * math.sin(distance) ** 2


@pytest.mark.parametrize("rho", [-math.sqrt(3) / 2, -.5, 0., .7])
def test_aperture_agrees_with_repo_full_covariance_cost(rho):
    chan = np.array([[4., 2 * rho], [2 * rho, 1.]])
    spatial = np.array([[1.5, .3], [.3, .8]])
    beta = math.acos(-rho)
    chol = np.linalg.cholesky(chan)
    endpoints = np.linalg.solve(chol, np.eye(2))
    frame, _ = np.linalg.qr(endpoints)
    frame *= np.sign(np.diag(frame.T @ endpoints))
    directions = frame.T @ endpoints
    assert math.atan2(*directions[::-1, 1]) == pytest.approx(beta)
    for raw in [np.array([[4., 1.], [1., 4.]]), np.array([[1., 2.], [-3., -2.]])]:
        whitened = np.linalg.solve(np.linalg.cholesky(spatial), raw.T).T
        matrix = frame.T @ np.linalg.solve(chol, whitened)
        reading = Reading(*raw, np.kron(chan, spatial))
        assert cost_on_arc(matrix, beta) == pytest.approx(reading_state(reading, 2.)["Q_plus"], abs=2e-9)


@pytest.mark.parametrize("beta", [math.pi / 6, math.pi / 2, .87 * math.pi])
def test_null_mean_from_angular_integral(beta):
    # Dimension one: radial squared length has mean 2; eigenline is uniform.
    integrand = lambda x: 2 * math.sin(max(0., min(x - beta, math.pi - x))) ** 2 / math.pi
    observed, _ = quad(integrand, 0, math.pi, points=[beta, (math.pi + beta) / 2])
    assert observed == pytest.approx(1 - (beta + math.sin(beta)) / math.pi, abs=1e-12)
    if beta == math.pi / 6:
        assert observed > .67


@pytest.mark.parametrize("rho,strength,threshold", [(-.65, 4., .4), (-.5, 6., .2), (-.32, 8., .08)])
def test_negative_correlation_exceeds_right_angle_mixture_tail(rho, strength, threshold):
    conditional_sd = math.sqrt(1 - rho * rho)
    gained, gain_error = quad(lambda z: norm.pdf(z) * norm.sf((threshold - rho*z) / conditional_sd),
                              -np.inf, -threshold - strength, epsabs=1e-22)
    lost, loss_error = quad(lambda z: norm.pdf(z) * norm.cdf((-threshold - rho*z) / conditional_sd),
                            -np.inf, threshold - strength, epsabs=1e-22)
    assert gained - lost > 100 * (gain_error + loss_error)


@pytest.mark.parametrize("k", [1, 2, 5])
def test_wide_arc_ray_bound_and_narrow_arc_limit(k):
    rng = np.random.default_rng(38 + k)
    for _ in range(12):
        x, y = rng.normal(size=(2, k))
        x[0] += 4
        matrix = np.stack((x, y))
        e = x / np.linalg.norm(x)
        ray_cost = np.sum((y - e * max(0., e @ y)) ** 2)
        assert cost_on_arc(matrix, .72 * math.pi) <= cost_on_arc(matrix, math.pi/2) + 1e-12
        assert cost_on_arc(matrix, math.pi/2) <= ray_cost + 1e-12
        assert cost_on_arc(matrix, 1e-9) == pytest.approx(y @ y, abs=1e-7)

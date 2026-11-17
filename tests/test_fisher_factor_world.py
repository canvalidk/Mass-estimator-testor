"""What the consistency of T72's counterexamples rests on (problems/O18.md, 17 Nov; Claude's check of B25).

T72's sign-dependent factors keep the true mass as the population optimum when the world's errors are exactly
Gaussian with the declared variance. That is arithmetic of the Gaussian's fourth moment, not symmetry: in a world
whose errors are a scale mixture with the same declared covariance the optimum moves, while a factor that is a
function of the reading's own residual alone stays at the truth. Deterministic quadrature; independent of
tests/test_fisher_factor_counterexample.py.
"""
import math

import numpy as np
import pytest

NODES, WEIGHTS = np.polynomial.hermite_e.hermegauss(60)
WEIGHTS = WEIGHTS / WEIGHTS.sum()

GAUSSIAN = ((1.0, 1.0),)                 # (weight, variance) in noise units
MIXTURE = ((0.5, 0.5), (0.5, 1.5))       # isotropic, covariance still the declared identity, not Gaussian


def _residual(x, y, theta):
    return (math.cos(theta) * x - math.sin(theta) * y) ** 2


def polynomial_loss(x, y, theta):
    """T72 section 1."""
    energy, sign_term, k = x * x + y * y, x * y, _residual(x, y, theta)
    return (energy ** 2 + 5 * energy) * k + sign_term * (k * k - 6 * k)


def damped_loss(x, y, theta, a=0.1):
    """T72 section 4's sign-dependent summand, the one that also meets bounded confidence."""
    v = 1 / (1 + 2 * a)
    energy, sign_term, k = x * x + y * y, x * y, _residual(x, y, theta)
    return np.exp(-a * energy) * ((energy ** 2 + 5 * v * energy) * k + sign_term * (k * k - 6 * v * k))


def residual_only_loss(x, y, theta, nu=3.0):
    """A function of the reading's own residual alone (the cut factor's form, T39)."""
    return np.log1p(_residual(x, y, theta) / nu)


def population_loss(loss, theta, true_theta, signal, world):
    mu_x, mu_y = signal * math.sin(true_theta), signal * math.cos(true_theta)
    total = 0.0
    for weight, variance in world:
        s = math.sqrt(variance)
        values = loss(mu_x + s * NODES[:, None], mu_y + s * NODES[None, :], theta)
        total += weight * float(np.sum(WEIGHTS[:, None] * WEIGHTS[None, :] * values))
    return total


def score_at_truth(loss, true_theta, signal, world, step=1e-5):
    up = population_loss(loss, true_theta + step, true_theta, signal, world)
    down = population_loss(loss, true_theta - step, true_theta, signal, world)
    return (up - down) / (2 * step)


def optimum(loss, true_theta, signal, world):
    grid = np.linspace(0.02, math.pi / 2 - 0.02, 801)
    values = [population_loss(loss, t, true_theta, signal, world) for t in grid]
    return math.tan(grid[int(np.argmin(values))])


TRUE = math.atan(2.0)


@pytest.mark.parametrize("signal", [1.0, 3.0])
def test_the_polynomial_is_stationary_at_the_truth_only_in_the_gaussian_world(signal):
    assert abs(score_at_truth(polynomial_loss, TRUE, signal, GAUSSIAN)) < 1e-5
    # with E tau^2 = 1 and E tau^4 = 5/4, E tau^6 = 7/4: the score is -3 (r^2 + 2) cos(2 theta0)
    expected = -3 * (signal ** 2 + 2) * math.cos(2 * TRUE)
    assert score_at_truth(polynomial_loss, TRUE, signal, MIXTURE) == pytest.approx(expected, rel=1e-5)


def test_the_bounded_confidence_counterexample_also_rests_on_the_gaussian_world():
    assert abs(score_at_truth(damped_loss, TRUE, 1.0, GAUSSIAN)) < 1e-6
    assert score_at_truth(damped_loss, TRUE, 1.0, MIXTURE) > 0.5


def test_the_optimum_moves_off_the_truth_under_the_mixture():
    assert optimum(polynomial_loss, TRUE, 1.0, GAUSSIAN) == pytest.approx(2.0, abs=0.01)
    assert 1.85 < optimum(polynomial_loss, TRUE, 1.0, MIXTURE) < 1.92      # 1.887
    assert 1.90 < optimum(damped_loss, TRUE, 1.0, MIXTURE) < 1.95          # 1.927


@pytest.mark.parametrize("signal", [1.0, 3.0])
def test_a_function_of_the_residual_alone_stays_at_the_truth_in_both_worlds(signal):
    for world in (GAUSSIAN, MIXTURE):
        assert abs(score_at_truth(residual_only_loss, TRUE, signal, world)) < 1e-8
        assert optimum(residual_only_loss, TRUE, signal, world) == pytest.approx(2.0, abs=0.01)

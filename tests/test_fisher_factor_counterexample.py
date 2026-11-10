"""T72: deterministic checks of the sign-dependent Fisher counterexample.

Tensor Gauss--Hermite quadrature integrates the degree-six losses and scores
exactly in exact arithmetic. These tests check a scalar Gaussian theorem,
not coverage, a production factor, or a full-type consistency claim.
"""
import math

import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import brentq


NODES, WEIGHTS = hermgauss(12)
NOISE_X, NOISE_Y = np.meshgrid(np.sqrt(2) * NODES, np.sqrt(2) * NODES)
NORMAL_WEIGHTS = np.outer(WEIGHTS, WEIGHTS) / math.pi


def residual(x, y, theta):
    return (math.cos(theta) * x - math.sin(theta) * y) ** 2


def loss_from_invariants(k, energy, cross):
    return (energy**2 + 5 * energy) * k + cross * (k**2 - 6 * k)


def loss(x, y, theta):
    return loss_from_invariants(residual(x, y, theta), x*x + y*y, x*y)


def angular_score(x, y, theta):
    t = math.sin(theta) * x + math.cos(theta) * y
    w = math.cos(theta) * x - math.sin(theta) * y
    energy, cross = x*x + y*y, x*y
    return -2 * t * w * (energy**2 + 5 * energy + cross * (2*w*w - 6))


def normal_expectation(function, mean):
    return float(np.sum(NORMAL_WEIGHTS * function(NOISE_X + mean[0],
                                                NOISE_Y + mean[1])))


def theoretical_gap(mean, theta):
    radius_squared = float(mean @ mean)
    normal = np.array([math.cos(theta), -math.sin(theta)])
    b = float(normal @ mean)
    c = mean[0] * normal[1] + mean[1] * normal[0]
    q = normal[0] * normal[1]
    return b*b * (radius_squared**2 + 21*radius_squared + 78
                  + mean[0]*mean[1]*b*b + 4*b*c + 12*q)


@pytest.mark.parametrize("radius", [-3., -.4, 0., .1, 1.3, 4.])
def test_population_score_vanishes_at_every_checked_truth_by_quadrature(radius):
    for theta in (.12, .5, math.pi/4, 1.2, 1.48):
        mean = radius * np.array([math.sin(theta), math.cos(theta)])
        score = normal_expectation(lambda x, y: angular_score(x, y, theta), mean)
        assert score == pytest.approx(0., abs=2e-10)


@pytest.mark.parametrize("radius", [-3., -.4, .1, 1.3, 4.])
def test_exact_population_gap_and_strict_global_bound_by_quadrature(radius):
    for truth in (.12, .5, math.pi/4, 1.2, 1.48):
        mean = radius * np.array([math.sin(truth), math.cos(truth)])
        for theta in (.03, .31, .84, 1.17, 1.54):
            # Integrate the difference, not two nearly equal expectations.
            measured = normal_expectation(
                lambda x, y: loss(x, y, theta) - loss(x, y, truth), mean)
            expected = theoretical_gap(mean, theta)
            b = radius * math.sin(truth - theta)
            lower = b*b * (radius**4 + 17*radius**2 + 72)
            assert measured == pytest.approx(expected, rel=2e-10, abs=2e-10)
            assert expected >= lower - 2e-10
            assert lower > 0


def test_true_null_population_is_flat_but_a_measured_zero_is_inert():
    values = [normal_expectation(lambda x, y: loss(x, y, theta), [0., 0.])
              for theta in (.01, .3, .8, 1.4, 1.56)]
    # E[E^2 K]=24 and E[5 E K]=20; the signed perturbation cancels.
    np.testing.assert_allclose(values, 44., rtol=0., atol=2e-12)
    for theta in (.01, .3, .8, 1.4, 1.56):
        assert loss(0., 0., theta) == 0.


def test_loss_derivative_agrees_with_candidate_angle_finite_difference():
    for x, y, theta in [(1., -2., .4), (-.7, 1.1, 1.3), (2., 3., .9)]:
        step = 1e-5
        numerical = (loss(x, y, theta + step) - loss(x, y, theta - step)) / (2*step)
        assert numerical == pytest.approx(angular_score(x, y, theta), rel=2e-9)


def test_loss_is_increasing_and_nonnegative_on_each_readings_residual_range():
    for energy in (0., .001, .2, 1., 7., 100.):
        for cross in np.linspace(-energy/2, energy/2, 11):
            for k in np.linspace(0., energy, 17):
                derivative = energy**2 + 5*energy + cross*(2*k - 6)
                assert derivative >= 2*energy - 1e-10
                assert loss_from_invariants(k, energy, cross) >= 0.


def test_sign_dependence_changes_mass_contrasts_and_cannot_be_a_constant():
    # Same E and absolute data; both K=1.5 and K=3.5 are attained at
    # positive candidate masses for each of (1,2) and (1,-2).
    targets = (1.5, 3.5)
    contrast = []
    for y in (2., -2.):
        bracket = (math.atan(.5), math.pi/2) if y > 0 else (0., math.atan(2.))
        angles = [brentq(lambda theta: residual(1., y, theta) - k, *bracket)
                  for k in targets]
        for theta, k in zip(angles, targets):
            assert residual(1., y, theta) == pytest.approx(k)
        contrast.append(loss(1., y, angles[1]) - loss(1., y, angles[0]))
    np.testing.assert_allclose(contrast, [96., 104.], atol=2e-10, rtol=0.)


def test_frame_swap_units_and_zero_padding_preserve_relative_factors():
    data = np.array([[1., -2.], [-.3, .6], [2., 3.]])
    theta, sigmas, mass = .43, np.array([.7, 1.9]), 1.3
    value = loss(data[:, 0], data[:, 1], theta)
    np.testing.assert_allclose(loss(-data[:, 0], -data[:, 1], theta), value)
    np.testing.assert_allclose(loss(data[:, 1], data[:, 0], math.pi/2-theta), value)
    with_zero = np.vstack([data, [0., 0.]])
    assert sum(loss(with_zero[:, 0], with_zero[:, 1], theta)) == pytest.approx(sum(value))
    physical = data * sigmas
    scales = np.array([2.3, .4])
    original_angle = math.atan(mass / (sigmas[0] / sigmas[1]))
    new_sigmas = sigmas * scales
    new_mass = mass * scales[0] / scales[1]
    converted = physical * scales / new_sigmas
    new_angle = math.atan(new_mass / (new_sigmas[0] / new_sigmas[1]))
    np.testing.assert_allclose(loss(converted[:, 0], converted[:, 1], new_angle),
                               loss(data[:, 0], data[:, 1], original_angle))


def test_antiparallel_noise_shrink_span_diverges_so_t30_is_not_claimed():
    # On positive angles (1,-2) has residual infimum 1 at theta=0 and
    # maximum 5 at atan(2). Monotonicity of rho makes these span endpoints.
    # Open-domain endpoint infima give the same evidence span.
    energy, cross, minimum, maximum = 5., -2., 1., 5.
    leading = energy**2*(maximum-minimum) + cross*(maximum**2-minimum**2)
    next_order = (5*energy-6*cross)*(maximum-minimum)
    assert leading == 52.
    assert next_order == 148.
    for noise_sd in (1., .5, .1, .02, .005):
        low = loss(1/noise_sd, -2/noise_sd, 0.)
        high = loss(1/noise_sd, -2/noise_sd, math.atan(2.))
        assert noise_sd**6*(high-low) == pytest.approx(
            leading + next_order*noise_sd**2, rel=3e-14)

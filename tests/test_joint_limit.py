"""T54 continuous-prior small-noise check on an explicit finite window.

This deterministic quadrature normalizes the continuous law *inside* an
h-box. It does not certify global posterior tails or total variation, and
does not apply an atomic-prior Gaussian-limit claim. T54 proves the global
concentration step. The quantities below are shape and quantile widths,
not moments of a finite-noise posterior or sampling distribution.
"""
from functools import lru_cache
import math

import numpy as np
import pytest
from scipy import optimize, stats

from met.joint_mass import common_acceleration, gaussian_profile_cost


_TRUE_MASSES = np.array([1., 2.])
_NOISE_DIRECTION = np.array([.3, -.4, .2])
_LOCAL_COVARIANCE = np.eye(2)+np.outer(_TRUE_MASSES, _TRUE_MASSES)
_PROBABILITY = .95
_HALF_WINDOW = 20.


def _physical_cost(data, first, second, epsilon):
    # The physical constraints have covariance I_2+m*m', rather than a
    # determinant factor in the profiled Gaussian likelihood.
    residual_first = data[0]-first*data[2]
    residual_second = data[1]-second*data[2]
    inner = first*residual_first+second*residual_second
    return ((residual_first**2+residual_second**2)
            -inner**2/(1+first**2+second**2))/epsilon**2


def _row_integral_to(axis, density, cumulative, upper):
    """Integrate the piecewise-linear row density up to each row's bound."""
    step = axis[1]-axis[0]
    position = np.clip((upper-axis[0])/step, 0., len(axis)-1.)
    index = np.minimum(np.floor(position).astype(int), len(axis)-2)
    fraction = position-index
    rows = np.arange(len(axis))
    return (cumulative[rows, index]+step*(fraction*density[rows, index]
        +fraction**2*(density[rows, index+1]-density[rows, index])/2))


@lru_cache(maxsize=4)
def _window_metrics(epsilon, nodes):
    data = np.r_[_TRUE_MASSES, 1.]+epsilon*_NOISE_DIRECTION
    fitted = data[:2]/data[2]  # Positive, saturated one-event Gaussian MLE.
    axis = np.linspace(-_HALF_WINDOW, _HALF_WINDOW, nodes)
    step = axis[1]-axis[0]
    first, second = fitted[0]+epsilon*axis[:, None], fitted[1]+epsilon*axis[None, :]
    assert np.min(first) > 0 and np.min(second) > 0
    cost = _physical_cost(data, first, second, epsilon)
    # Independent Exp(1) priors form a proper continuous joint mass prior.
    log_density = -.5*cost-first-second
    density = np.exp(log_density-np.max(log_density))
    edge = np.ones(nodes)
    edge[[0, -1]] = .5
    weights = np.outer(edge, edge)*step**2
    density /= np.sum(weights*density)

    precision = np.linalg.inv(_LOCAL_COVARIANCE)
    h1, h2 = axis[:, None], axis[None, :]
    gaussian = np.exp(-.5*(precision[0, 0]*h1**2
        +2*precision[0, 1]*h1*h2+precision[1, 1]*h2**2))
    gaussian /= np.sum(weights*gaussian)
    shape_error = .5*np.sum(weights*np.abs(density-gaussian))

    cumulative = np.column_stack([np.zeros(nodes), np.cumsum(
        .5*(density[:, :-1]+density[:, 1:])*step, axis=1)])
    marginal_weights = edge*step
    probabilities = [(1-_PROBABILITY)/2, (1+_PROBABILITY)/2]
    z = stats.norm.ppf(probabilities[1])
    results = dict(shape_error=float(shape_error))
    for operation, gradient in (("sum", np.ones(2)), ("parallel", np.array([4/9, 1/9]))):
        fitted_composite = (sum(fitted) if operation == "sum"
                            else fitted[0]*fitted[1]/sum(fitted))

        def cdf(scaled_value):
            value = fitted_composite+epsilon*scaled_value
            if operation == "sum":
                upper_mass = value-first[:, 0]
            elif value <= 0:
                return 0.
            else:
                denominator = first[:, 0]-value
                # Parallel sum is below value for all positive m2 when m1<=value.
                upper_mass = np.full(nodes, np.inf)
                positive = denominator > 0
                upper_mass[positive] = value*first[positive, 0]/denominator[positive]
            upper = (upper_mass-fitted[1])/epsilon
            return float(marginal_weights@_row_integral_to(axis, density, cumulative, upper))

        quantiles = [optimize.brentq(lambda value: cdf(value)-probability,
                                     -_HALF_WINDOW, _HALF_WINDOW, xtol=2e-10)
                     for probability in probabilities]
        width = (quantiles[1]-quantiles[0])/2
        gaussian_width = z*math.sqrt(gradient@_LOCAL_COVARIANCE@gradient)
        results[operation+"_width"] = float(width)
        results[operation+"_relative_width_error"] = float(abs(width/gaussian_width-1))
    return results


@pytest.mark.parametrize("epsilon", [.04, .01])
def test_vectorized_physical_residual_matches_the_independent_joint_qr_profile(epsilon):
    data = np.r_[_TRUE_MASSES, 1.]+epsilon*_NOISE_DIRECTION
    fitted = data[:2]/data[2]
    h = np.array([[0., 0.], [.7, -.5], [-3., 2.], [5., -4.]])
    masses = fitted+epsilon*h
    expected = _physical_cost(data, masses[:, 0], masses[:, 1], epsilon)
    reading = common_acceleration(data[:2, None], data[2:], epsilon**2*np.eye(3))
    actual = np.array([gaussian_profile_cost(reading, point) for point in masses])
    np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=2e-12)


def test_continuous_joint_window_shape_and_composite_widths_approach_the_gaussian_limit():
    metrics = [_window_metrics(epsilon, 401) for epsilon in (.04, .02, .01)]
    # Halving epsilon removes roughly half the shape/width error. A wrong
    # limiting covariance would instead leave a nonzero error plateau.
    for quantity in ("shape_error", "sum_relative_width_error", "parallel_relative_width_error"):
        assert metrics[0][quantity] > metrics[1][quantity] > metrics[2][quantity]
        for earlier, later in zip(metrics, metrics[1:]):
            assert .45 < later[quantity]/earlier[quantity] < .6
    assert metrics[-1]["shape_error"] < .009
    assert metrics[-1]["sum_relative_width_error"] < .0042
    assert metrics[-1]["parallel_relative_width_error"] < .0036


def test_window_grid_refinement_is_smaller_than_the_observed_small_noise_change():
    coarse, refined = _window_metrics(.01, 401), _window_metrics(.01, 601)
    previous = _window_metrics(.02, 401)
    for quantity in ("shape_error", "sum_width", "parallel_width"):
        numerical_change = abs(coarse[quantity]-refined[quantity])
        small_noise_change = abs(previous[quantity]-coarse[quantity])
        assert numerical_change < .02*small_noise_change

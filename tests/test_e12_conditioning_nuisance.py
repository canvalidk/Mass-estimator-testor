"""Ordinary Gaussian outcome kernels and the four-corner likelihood test."""
from fractions import Fraction

import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import logsumexp

from explore.e12_residual_join import residual_log_evidence


ANGLES = np.arctan([1., 2.])
LINES = np.column_stack((np.sin(ANGLES), np.cos(ANGLES)))
CORNERS = np.array([[2., 2.], [3., 1.], [1., 4.], [2., 3.]])


def mixed(values):
    first, second, third, fourth = values
    return fourth - third - second + first


@pytest.mark.parametrize("tau", [0., 64.])
def test_e12_singletons_violate_gaussian_mixture_odds_separability(tau):
    odds = []
    for y in CORNERS:
        factors = residual_log_evidence([np.outer(y, y)], ANGLES, tau=tau)
        odds.append(factors[0] - factors[1])
    expected = [Fraction(2, 5), Fraction(-9, 50), Fraction(53, 170), Fraction(27, 20)]
    assert mixed(expected) == Fraction(2751, 1700)
    np.testing.assert_allclose(odds, np.array(expected, float), atol=3e-14)


def test_separate_candidate_specific_drive_mixtures_obey_the_rectangle_identity():
    # These nuisance populations differ in location, shape and scale across the two candidates.
    priors = [([-2., .4, 2.1], [.05, .25, .70]), ([-4., .1, 1., 6.], [.4, .01, .49, .1])]
    log_densities = []
    for u, (r, weights) in zip(LINES, priors):
        errors = CORNERS[:, None, :] - np.outer(r, u)[None, :, :]
        log_densities.append(logsumexp(np.log(weights) - .5 * (errors**2).sum(axis=-1), axis=-1))
    odds = log_densities[0] - log_densities[1]
    assert np.ptp(odds) > 1
    assert mixed(odds) == pytest.approx(0., abs=4e-15)


@pytest.mark.parametrize("dimension", [1, 3])
@pytest.mark.parametrize("variance", [0., .4, 10.])
def test_gaussian_drive_yields_exactly_one_power_of_the_residual(dimension, variance):
    y = np.array([[1.4, -.2], [-1., 3.], [.34, .53]])[:dimension]
    exponent = variance / (1 + variance)
    direct, residual = [], []
    for u in LINES:
        # A pre-data spherical push law lives on the candidate's own k-plane.
        design = np.kron(np.eye(dimension), u[:, None])
        cov = np.eye(2 * dimension) + variance * design @ design.T
        direct.append(-.5 * (np.linalg.slogdet(cov)[1] + y.ravel() @ np.linalg.solve(cov, y.ravel())))
        transverse = y @ [u[1], -u[0]]
        residual.append(-.5 * exponent * (transverse @ transverse))
    assert np.subtract(*direct) == pytest.approx(np.subtract(*residual), abs=1e-14)


def test_bounded_pointwise_residual_loses_normalized_odds_structure():
    costs = np.array([np.square(CORNERS @ [u[1], -u[0]]) for u in LINES])
    affine_logs = -.5 * costs
    assert mixed(affine_logs[0] - affine_logs[1]) == pytest.approx(0., abs=6e-15)
    bounded_logs = np.expm1(-.5 * costs)
    assert abs(mixed(bounded_logs[0] - bounded_logs[1])) > .1


def test_nontrivial_e12_update_differs_from_null_conditionally_reliable_probability():
    # At a true zero drive the mass is independent of every report. For
    # a uniform angle prior the conditionally reliable lower-half probability is 1/2.
    angles = np.linspace(0, np.pi / 2, 2001)
    for y in [np.array([2., 1.]), np.array([2.01, .99])]:
        logp = residual_log_evidence([np.outer(y, y)], angles)
        density = np.exp(logp - logp.max())
        middle = len(angles) // 2
        left = trapezoid(density[:middle + 1], angles[:middle + 1])
        total = trapezoid(density, angles)
        assert left / total < .49

"""Calibrated bounded ray scores for the unknown-drive force experiment."""
import math

import numpy as np
import pytest
from scipy.special import roots_hermitenorm
from scipy.stats import norm

from explore.e12_mass_ray_diagnostic import mass_ray_diagnostic


@pytest.mark.parametrize("masses,drive", [
    ([3., 7.], 5.), ([5., 2.], 0.), ([5., 2.], -.7), ([5., 2.], 1.8),
])
def test_integrating_joint_force_noise_gives_the_declared_discrepancy(masses, drive):
    r = .3
    cov = np.array([[1.7, .43], [.43, .8]])
    c = np.array([1-r, -r])
    nodes, weights = roots_hermitenorm(48)
    independent = np.array(np.meshgrid(nodes, nodes)).reshape(2, -1).T
    mean = drive * np.array(masses)
    forces = mean + independent @ np.linalg.cholesky(cov).T
    scores = mass_ray_diagnostic(forces, cov, r).scores
    integral = scores @ np.outer(weights, weights).ravel() / (2*np.pi)
    target = -np.expm1(-(c @ mean)**2 / (4*(c @ cov @ c)))
    assert integral == pytest.approx(target, abs=1e-12)


def test_two_pushes_heterogeneous_covariance_exact_chi_square_tail_and_bound():
    forces = np.array([[1.2, -1.3], [-.65, 1.]])
    cov = np.array([[[.7, -.2], [-.2, 1.]], [[.3, .08], [.08, .4]]])
    r, confidence = .25, .1
    c = np.array([1-r, -r])
    z = np.array([c @ row / np.sqrt(c @ matrix @ c) for row, matrix in zip(forces, cov)])
    result = mass_ray_diagnostic(forces, cov, r, confidence=confidence)
    np.testing.assert_allclose(result.standardized_contrasts, z)
    expected = np.mean(1-np.sqrt(2)*np.exp(-z*z/2))
    assert result.mean_score == pytest.approx(expected)
    assert result.upper_discrepancy == pytest.approx(expected + np.sqrt(-np.log(.9)/2))
    assert result.chi_square == pytest.approx(z @ z)
    # Central chi-square(2) is exponential of mean two.
    assert result.null_p_value == pytest.approx(np.exp(-(z @ z)/2))
    reflected = mass_ray_diagnostic(-forces[:, ::-1], cov[:, ::-1, ::-1], 1-r)
    np.testing.assert_allclose(reflected.standardized_contrasts, z)


def test_single_push_tail_matches_the_exact_two_sided_gaussian_probability():
    result = mass_ray_diagnostic([2., -.4], np.eye(2), .5)
    assert result.null_p_value == pytest.approx(2*norm.sf(2.4/np.sqrt(2)))
    assert result.upper_discrepancy == 1.


def test_smallest_score_and_lower_clipping_of_bound():
    result = mass_ray_diagnostic(np.tile([3., 9.], (100, 1)), np.eye(2), .25)
    assert result.mean_score == pytest.approx(1-math.sqrt(2))
    assert result.upper_discrepancy == 0
    assert result.null_p_value == 1


def test_rank_one_instrument_with_positive_contrast_variance_and_huge_signal():
    result = mass_ray_diagnostic([1e300, 1.], [[1., -1.], [-1., 1.]], .5)
    assert result.standardized_contrasts[0] == pytest.approx(.5e300)
    assert result.scores[0] == 1.
    assert np.isinf(result.chi_square) and result.null_p_value == 0.
    noiseless_second = mass_ray_diagnostic([2., 10.], [[4., 0.], [0., 0.]], .25)
    assert noiseless_second.standardized_contrasts[0] == pytest.approx(-1/1.5)


@pytest.mark.parametrize("field,value", [
    ("forces", []), ("forces", [0., np.nan]), ("forces", [[1., 2., 3.]]),
    ("allocation", 0.), ("allocation", 1.), ("allocation", np.nan),
    ("confidence", 1.), ("confidence", 0.),
    ("covariances", [[1., 0.], [.1, 1.]]), ("covariances", [[1., 2.], [2., 1.]]),
    ("covariances", np.ones((2, 2))), ("covariances", [[1., .1], [.1, 0.]]),
    ("covariances", np.zeros((3, 2, 2))),
])
def test_invalid_inputs_and_unstandardizable_contrasts(field, value):
    arguments = dict(forces=[1., 1.], allocation=.5, confidence=.95, covariances=np.eye(2))
    arguments[field] = value
    with pytest.raises(ValueError):
        mass_ray_diagnostic(**arguments)

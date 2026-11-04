"""T55: independent Gaussian integration, profile derivatives and composition."""
import math

import numpy as np
import pytest

from met.composition_calibration import gaussian_linear_coverage
from met.joint_calibration import (gaussian_joint_observed_moments,
    gaussian_joint_profile_score, gaussian_joint_score_moments,
    joint_covariance_benchmark)
from met.joint_mass import (common_acceleration, efficient_mass_information,
                            separate_accelerations)


def _witness():
    return [common_acceleration([[1.], [1.]], [1.], covariance)
            for covariance in (np.eye(3), np.diag([1., 9., 1.]))]


def _constraint_log_factor(reading, masses):
    # Direct physical residual F_j-m_j*A, not QR/projector differentiation.
    constraint = np.column_stack([np.eye(2), -np.asarray(masses)])
    residual = constraint@reading.data
    covariance = constraint@reading.covariance@constraint.T
    return -.5*residual@np.linalg.solve(covariance, residual)


def test_observed_score_and_hessian_against_independent_constraint_derivatives():
    factor = np.array([[1.2, .1, -.3], [.2, 1.4, .1], [.1, -.2, .9]])
    reading = common_acceleration([[2.1], [.4]], [.7], factor@factor.T)
    masses, step = np.array([1.3, .8]), 1e-4
    score = gaussian_joint_profile_score([reading], masses)
    observed = gaussian_joint_observed_moments([reading], masses)
    reference = []
    hessian = np.empty((2, 2))
    for axis in range(2):
        direction = np.eye(2)[axis]
        values = [_constraint_log_factor(reading, masses+t*step*direction)
                  for t in (-2, -1, 1, 2)]
        reference.append((values[0]-8*values[1]+8*values[2]-values[3])/(12*step))
        # Independently differentiate the envelope score to check the
        # analytic least-squares Hessian.
        hessian[:, axis] = (gaussian_joint_profile_score([reading], masses+step*direction)
                           -gaussian_joint_profile_score([reading], masses-step*direction))/(2*step)
    assert score == pytest.approx(reference, rel=1e-11, abs=2e-11)
    assert observed["score"] == pytest.approx(score, rel=1e-13, abs=2e-14)
    assert observed["sensitivity"] == pytest.approx(-hessian, abs=2e-8)
    assert observed["variability"] == pytest.approx(np.outer(score, score), rel=1e-13, abs=2e-14)


@pytest.mark.parametrize("covariance", [np.eye(3), np.diag([1., 9., 1.]),
    np.array([[1.3, .2, -.1], [.2, 1.5, .3], [-.1, .3, 1.1]])])
def test_exact_score_moments_and_expected_hessian_by_gaussian_quadrature(covariance):
    masses, signal = np.array([1.2, .8]), np.array([.7])
    declaration = common_acceleration([[0.], [0.]], [0.], covariance)
    theoretical = gaussian_joint_score_moments([declaration], masses, [signal])
    assert theoretical["sensitivity"] == pytest.approx(
        efficient_mass_information([declaration], masses, [signal]), rel=1e-13, abs=2e-14)
    nodes, weights = np.polynomial.hermite.hermgauss(4)
    indices = np.array(list(np.ndindex(4, 4, 4)))
    noise = math.sqrt(2)*nodes[indices]
    probabilities = np.prod(weights[indices], axis=1)/math.pi**1.5
    points = declaration.design(masses)@signal+noise@np.linalg.cholesky(covariance).T
    mean_score, second_score, mean_hessian = np.zeros(2), np.zeros((2, 2)), np.zeros((2, 2))
    for point, probability in zip(points, probabilities):
        event = common_acceleration(point[:2, None], point[2:], covariance)
        observed = gaussian_joint_observed_moments([event], masses)
        mean_score += probability*observed["score"]
        second_score += probability*observed["variability"]
        mean_hessian += probability*observed["sensitivity"]
    assert mean_score == pytest.approx(np.zeros(2), abs=2e-15)
    assert second_score == pytest.approx(theoretical["variability"], rel=1e-13, abs=3e-14)
    assert mean_hessian == pytest.approx(theoretical["sensitivity"], rel=1e-13, abs=3e-14)


def test_t54_conditional_joint_score_by_residual_gaussian_integration():
    masses, fitted_signal = np.array([1.2, .8]), .7
    covariance = np.array([[1.3, .2, -.1], [.2, 1.5, .3], [-.1, .3, 1.1]])
    chol = np.linalg.cholesky(covariance)
    design = np.r_[masses, 1.][:, None]
    basis = np.linalg.solve(chol, design)
    rotation, _ = np.linalg.qr(basis, mode="complete")
    constraint = np.column_stack([np.eye(2), -masses])
    # Physical residual precision gives the conditional covariance separately
    # from the projector-score implementation under test.
    expected_covariance = fitted_signal**2*np.linalg.inv(constraint@covariance@constraint.T)
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    mean, second, statistic_moments = np.zeros(2), np.zeros((2, 2)), np.zeros(2)
    for i, j in np.ndindex(3, 3):
        probability = weights[i]*weights[j]/math.pi
        residual = math.sqrt(2)*np.array([nodes[i], nodes[j]])
        point = chol@(basis[:, 0]*fitted_signal+rotation[:, 1:]@residual)
        event = common_acceleration(point[:2, None], point[2:], covariance)
        score = gaussian_joint_profile_score([event], masses)
        statistic = score@np.linalg.solve(expected_covariance, score)
        mean += probability*score
        second += probability*np.outer(score, score)
        statistic_moments += probability*np.array([statistic, statistic**2])
    assert mean == pytest.approx(np.zeros(2), abs=2e-15)
    assert second == pytest.approx(expected_covariance, rel=1e-13, abs=2e-15)
    # First two moments of chi-square_2; linear Gaussianity proves its law.
    assert statistic_moments == pytest.approx([2., 8.], rel=1e-13, abs=2e-14)


def test_separate_latent_vectors_with_correlated_errors_have_the_same_moment_identity():
    masses, signal = np.array([1.2, .8]), np.array([.7, -.4])
    covariance = np.eye(4)+.1*np.ones((4, 4))
    declaration = separate_accelerations(np.zeros((2, 1)), np.zeros((2, 1)), covariance)
    theoretical = gaussian_joint_score_moments([declaration], masses, [signal])
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    mean_hessian, second_score = np.zeros((2, 2)), np.zeros((2, 2))
    for index in np.ndindex(3, 3, 3, 3):
        probability = np.prod(weights[list(index)])/math.pi**2
        noise = math.sqrt(2)*nodes[list(index)]
        point = declaration.design(masses)@signal+np.linalg.cholesky(covariance)@noise
        event = separate_accelerations(point[:2, None], point[2:, None], covariance)
        observed = gaussian_joint_observed_moments([event], masses)
        mean_hessian += probability*observed["sensitivity"]
        second_score += probability*observed["variability"]
    assert mean_hessian == pytest.approx(theoretical["sensitivity"], rel=1e-13, abs=3e-14)
    assert second_score == pytest.approx(theoretical["variability"], rel=1e-13, abs=3e-14)


@pytest.mark.parametrize("mass", [1e4, 1e8, 1e12])
def test_compatible_large_mass_retains_small_transverse_curvature(mass):
    event = common_acceleration([[mass]], [1.], np.eye(2))
    observed = gaussian_joint_observed_moments([event], [mass])
    expected = 1/(1+mass**2)
    assert abs(gaussian_joint_profile_score([event], [mass])[0]) < 5e-16/mass
    assert abs(observed["score"][0]) < 5e-16/mass
    assert observed["sensitivity"][0, 0] == pytest.approx(expected, rel=1e-13, abs=0.)


def test_explicit_heterogeneous_protocol_requires_distinct_scalar_inflations():
    first = np.array([[2/3, -1/3], [-1/3, 2/3]])
    second = np.array([[10/19, -1/19], [-1/19, 2/19]])
    moments = gaussian_joint_score_moments(_witness(), [1., 1.], [[1.], [1.]])
    assert moments["sensitivity"] == pytest.approx(first+second, rel=1e-13, abs=2e-15)
    assert moments["noise_geometry"] == pytest.approx(first/3+9*second/19, rel=1e-13, abs=2e-15)
    benchmark = joint_covariance_benchmark(moments["sensitivity"], moments["variability"])
    raw, actual = benchmark["raw_covariance"], benchmark["sampling_covariance"]
    assert raw == pytest.approx(np.array([[1., .5], [.5, 17/11]]), rel=1e-13, abs=2e-15)
    assert actual[0, 0] == pytest.approx(80/57, rel=1e-13, abs=2e-15)
    spectrum = benchmark["inflation_eigenvalues"]
    assert spectrum == pytest.approx([848/627, 80/57], rel=1e-13, abs=2e-15)
    required = []
    for gradient in (np.array([1., 0.]), np.array([0., 1.]), np.ones(2)):
        eta = (gradient@raw@gradient)/(gradient@actual@gradient)
        required.append(eta)
        assert gaussian_linear_coverage(actual, raw/eta, gradient, probability=.95)["coverage"] == pytest.approx(.95, abs=1e-15)
    assert required == pytest.approx([.7125, .7349007170435741, .7220942593905031], abs=2e-15)


def test_conservative_scalar_is_sharp_for_all_contrasts_but_does_not_select_exact_covariance():
    moments = gaussian_joint_score_moments(_witness(), [1., 1.], [[1.], [1.]])
    benchmark = joint_covariance_benchmark(moments["sensitivity"], moments["variability"])
    raw, actual = benchmark["raw_covariance"], benchmark["sampling_covariance"]
    inflation = benchmark["minimal_conservative_inflation"]
    assert np.min(np.linalg.eigvalsh(inflation*raw-actual)) >= -2e-15
    assert np.min(np.linalg.eigvalsh((inflation-1e-4)*raw-actual)) < -1e-5
    for gradient in ([1., 0.], [0., 1.], [1., 1.], [1., -1.]):
        coverage = gaussian_linear_coverage(actual, inflation*raw, gradient, probability=.95)["coverage"]
        assert coverage >= .95-1e-15
    assert gaussian_linear_coverage(actual, inflation*raw, [0., 1.], probability=.95)["coverage"] > .953


def test_scalar_reduction_and_fixed_protocol_homogeneity_allow_a_single_inflation():
    event = common_acceleration([[0., 0.]], [0., 0.], np.eye(4))
    moments = gaussian_joint_score_moments([event], [2.], [[.7, -.4]])
    assert moments["sensitivity"][0, 0] == pytest.approx(.65/5)
    assert moments["noise_geometry"][0, 0] == pytest.approx(2/25)
    # In a two-body common isotropic protocol every tangent has the same
    # inflation when latent energy and noise geometry are proportional.
    common = common_acceleration(np.zeros((2, 2)), np.zeros(2), np.eye(6))
    moments = gaussian_joint_score_moments([common], [1., 2.], [[.7, -.4]])
    benchmark = joint_covariance_benchmark(moments["sensitivity"], moments["variability"])
    assert benchmark["inflation_eigenvalues"] == pytest.approx(np.full(2, 1+2/(6*.65)))


def test_common_frame_and_unit_covariance_of_joint_benchmark():
    masses, signal = np.array([1.3, .8]), np.array([.7, -.4])
    factor = np.eye(6)+.03*np.arange(36).reshape(6, 6)
    covariance = factor@factor.T
    event = common_acceleration(np.zeros((2, 2)), np.zeros(2), covariance)
    original = gaussian_joint_score_moments([event], masses, [signal])
    rotation = np.array([[0., -1.], [1., 0.]])
    transform = np.kron(np.eye(3), rotation)
    rotated = common_acceleration(np.zeros((2, 2)), np.zeros(2), transform@covariance@transform.T)
    changed = gaussian_joint_score_moments([rotated], masses, [rotation@signal])
    for key in ("sensitivity", "variability", "noise_geometry"):
        assert changed[key] == pytest.approx(original[key], abs=2e-14)
    scale = 3.
    transform = np.diag([scale]*4+[1.]*2)
    scaled = common_acceleration(np.zeros((2, 2)), np.zeros(2), transform@covariance@transform.T)
    changed = gaussian_joint_score_moments([scaled], scale*masses, [signal])
    for key in ("sensitivity", "variability", "noise_geometry"):
        assert changed[key] == pytest.approx(original[key]/scale**2, abs=2e-14)


def test_observable_sums_respect_zero_padding_and_replication_without_signal_inputs():
    events = [common_acceleration([[8.+x], [12.+y]], [4.], np.eye(3))
              for x, y in ((.2, -.3), (-.4, .1), (.3, .4))]
    masses = [2., 3.]
    observed = gaussian_joint_observed_moments(events, masses)
    zero = common_acceleration([[0.], [0.]], [0.], np.diag([2., 3., 4.]))
    padded = gaussian_joint_observed_moments(events+[zero], masses)
    repeated = gaussian_joint_observed_moments(events*5, masses)
    for key in ("score", "sensitivity", "variability"):
        np.testing.assert_array_equal(padded[key], observed[key])
        assert repeated[key] == pytest.approx(5*observed[key], abs=2e-13)
    original_cov = joint_covariance_benchmark(observed["sensitivity"], observed["variability"])
    repeated_cov = joint_covariance_benchmark(repeated["sensitivity"], repeated["variability"])
    assert repeated_cov["sampling_covariance"] == pytest.approx(original_cov["sampling_covariance"]/5, abs=2e-15)


def test_empty_and_singular_protocols_are_not_silently_regularized():
    zero = gaussian_joint_observed_moments([], [1., 2.])
    assert np.all(zero["sensitivity"] == 0) and np.all(zero["variability"] == 0)
    with pytest.raises(ValueError):
        joint_covariance_benchmark(zero["sensitivity"], np.eye(2))
    event = _witness()[0]
    for signals in ([], [[1., 2.]], [[math.nan]]):
        with pytest.raises(ValueError):
            gaussian_joint_score_moments([event], [1., 1.], signals)
    with pytest.raises(ValueError):
        joint_covariance_benchmark(np.eye(2), np.ones((3, 3)))

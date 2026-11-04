"""Independent physical-protocol geometry and finite joint-law checks for T54."""
from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest

from met.composition_calibration import gaussian_linear_coverage
from met.full import Reading
from met.joint_mass import (FiniteJointMassLaw, common_acceleration,
                            efficient_mass_information, gaussian_joint_geometry,
                            gaussian_joint_posterior, gaussian_profile_cost,
                            group_joint_law, scalar_composite_law,
                            separate_accelerations)
from met.operator_score import gaussian_operator_profile


def _covariance(size):
    matrix = np.eye(size)+np.arange(size*size).reshape(size, size)/(7*size*size)
    return matrix@matrix.T


@pytest.mark.parametrize("protocol", ["common", "separate"])
def test_full_covariance_profile_matches_independent_physical_constraints(protocol):
    masses = np.array([1.7, .8])
    forces = np.array([[2., -.3], [.4, 1.2]])
    if protocol == "common":
        reading = common_acceleration(forces, [.7, -.2], _covariance(6))
        constraint = np.c_[np.eye(4), -np.vstack([m*np.eye(2) for m in masses])]
    else:
        reading = separate_accelerations(forces, [[.7, -.2], [.1, .4]], _covariance(8))
        constraint = np.c_[np.eye(4), -np.diag(np.repeat(masses, 2))]
    residual = constraint@reading.data
    expected = residual@np.linalg.solve(constraint@reading.covariance@constraint.T, residual)
    assert gaussian_profile_cost(reading, masses) == pytest.approx(expected, rel=2e-14)
    np.testing.assert_allclose(constraint@reading.design(masses), 0., atol=1e-15)


def test_high_signal_cost_uses_residual_coordinates_instead_of_energy_subtraction():
    radius = 1e8
    masses, error = np.array([2., 3.]), np.array([.5, -.7, .3])
    data = radius*np.r_[masses, 1.]+error
    reading = common_acceleration(data[:2, None], data[2:], np.eye(3))
    constraint = np.c_[np.eye(2), -masses]
    physical = constraint@error
    expected = physical@np.linalg.solve(constraint@constraint.T, physical)
    assert gaussian_profile_cost(reading, masses) == pytest.approx(expected, rel=2e-6)


def test_one_body_both_protocols_reduce_to_the_existing_scalar_graph_profile():
    force, acceleration, covariance = [2., -.4], [.7, .2], _covariance(4)
    legacy = Reading(force, acceleration, covariance)
    for mass in (.2, 2., 7.):
        expected = -2*gaussian_operator_profile([legacy], mass*np.eye(2))
        for reading in (common_acceleration([force], acceleration, covariance),
                        separate_accelerations([force], [acceleration], covariance)):
            assert gaussian_profile_cost(reading, [mass]) == pytest.approx(expected, abs=2e-14)


@pytest.mark.parametrize("protocol", ["common", "separate"])
def test_public_design_derivatives_and_explicit_body_selectors(protocol):
    masses = np.array([.8, 2., 1.4])
    if protocol == "common":
        reading = common_acceleration([[1., .3], [.4, -.2]], [.7, .1], np.eye(6), body_indices=[2, 0])
    else:
        reading = separate_accelerations([[1., .3], [.4, -.2]], [[.7, .1], [-.3, .2]], np.eye(8), body_indices=[2, 0])
    derivatives = reading.design_derivatives(masses)
    assert derivatives.shape == (3, len(reading.data), reading.latent_dimension)
    assert not np.any(derivatives[1])
    for index in range(3):
        delta = np.eye(3)[index]*1e-5
        expected = (reading.design(masses+delta)-reading.design(masses-delta))/(2e-5)
        np.testing.assert_allclose(derivatives[index], expected, atol=1e-10)


def test_efficient_information_matches_independent_full_fisher_schur_complement():
    masses = np.array([1.7, .8])
    reading = common_acceleration([[2., -.3], [.4, 1.2]], [.7, -.2], _covariance(6))
    signal = np.array([1.2, -.5])
    design = reading.design(masses)
    derivatives = np.zeros((6, 2))
    derivatives[:2, 0], derivatives[2:4, 1] = signal, signal
    jacobian = np.c_[derivatives, design]
    fisher = jacobian.T@np.linalg.solve(reading.covariance, jacobian)
    expected = fisher[:2, :2]-fisher[:2, 2:]@np.linalg.solve(fisher[2:, 2:], fisher[2:, :2])
    result = efficient_mass_information([reading], masses, [signal])
    np.testing.assert_allclose(result, expected, atol=2e-15)
    changed_data = common_acceleration([[100., -30.], [40., 120.]], [7., -20.], reading.covariance)
    np.testing.assert_allclose(efficient_mass_information([changed_data], masses, [signal]), result, atol=0.)
    np.testing.assert_allclose(efficient_mass_information([reading], masses, [np.zeros(2)]), 0., atol=0.)


def test_shared_acceleration_yields_the_positive_cross_covariance_and_sum_coverage_gap():
    masses, signal = np.array([2., 3.]), np.array([1.2, -.5])
    force_variances, acceleration_variance = np.array([.7, 1.3]), .4
    covariance = np.diag(np.repeat(np.r_[force_variances, acceleration_variance], 2))
    reading = common_acceleration([[2., -.3], [.4, 1.2]], [.7, -.2], covariance)
    information = efficient_mass_information([reading], masses, [signal])
    expected = (np.diag(force_variances)+acceleration_variance*np.outer(masses, masses))/(signal@signal)
    actual = np.linalg.inv(information)
    np.testing.assert_allclose(actual, expected, rtol=3e-14)
    assert actual[0, 1] > 0
    correct = gaussian_linear_coverage(actual, actual, [1., 1.], probability=.95)
    wrong = gaussian_linear_coverage(actual, np.diag(np.diag(actual)), [1., 1.], probability=.95)
    assert correct["coverage"] == pytest.approx(.95, abs=1e-15)
    assert wrong["coverage"] < .9


def test_local_ratio_error_covariance_by_independent_physical_gaussian_quadrature():
    # This integrates the linearized sampling error, not the finite-noise
    # ratio's moments or the moments of a heavy-tailed continuous mass law.
    masses, signal = np.array([2., 3.]), 1.7
    covariance = np.array([[.7, .1, -.03], [.1, 1.3, .04], [-.03, .04, .4]])
    reading = common_acceleration([[2.], [3.]], [1.], covariance)
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    indices = np.array(list(np.ndindex(3, 3, 3)))
    errors = math.sqrt(2)*nodes[indices]@np.linalg.cholesky(covariance).T
    probabilities = np.prod(weights[indices], axis=1)/math.pi**1.5
    local_mass_errors = (errors[:, :2]-errors[:, 2, None]*masses)/signal
    integrated = np.einsum("n,ni,nj->ij", probabilities, local_mass_errors, local_mass_errors)
    expected = np.linalg.inv(efficient_mass_information([reading], masses, [[signal]]))
    np.testing.assert_allclose(integrated, expected, atol=2e-15)
    for gradient in (np.ones(2), masses[::-1]**2/masses.sum()**2):
        projection = local_mass_errors@gradient
        assert probabilities@(projection**2) == pytest.approx(gradient@expected@gradient, abs=2e-15)


def test_separately_observed_accelerations_have_their_declared_free_latent_vectors():
    masses, signals = np.array([2., 3.]), np.array([1.2, -.5])
    force_variances, acceleration_variances = np.array([.7, 1.3]), np.array([.4, .6])
    covariance = np.diag(np.r_[force_variances, acceleration_variances])
    reading = separate_accelerations([[2.], [4.]], [[1.], [.8]], covariance)
    expected = np.diag(signals**2/(force_variances+masses**2*acceleration_variances))
    np.testing.assert_allclose(efficient_mass_information([reading], masses, [signals]), expected, atol=2e-15)
    # A shared observation cannot be duplicated through this factory.
    with pytest.raises(ValueError):
        separate_accelerations([[2.], [4.]], [1.], covariance)


def test_declared_units_frame_and_body_permutation_preserve_the_model():
    masses, signal = np.array([1.7, .8]), np.array([1.2, -.5])
    reading = common_acceleration([[2., -.3], [.4, 1.2]], [.7, -.2], _covariance(6))
    angle = .4
    rotation = np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]])
    factors = np.array([3., 7., .4])
    transform = np.zeros((6, 6))
    for j, factor in enumerate(factors):
        transform[2*j:2*j+2, 2*j:2*j+2] = factor*rotation
    changed_data = transform@reading.data
    changed = common_acceleration(changed_data[:4].reshape(2, 2), changed_data[4:],
                                  transform@reading.covariance@transform.T)
    changed_masses = masses*factors[:2]/factors[2]
    assert gaussian_profile_cost(changed, changed_masses) == pytest.approx(gaussian_profile_cost(reading, masses), abs=2e-14)
    information = efficient_mass_information([reading], masses, [signal])
    inverse_jacobian = np.diag(factors[2]/factors[:2])
    expected = inverse_jacobian@information@inverse_jacobian
    np.testing.assert_allclose(efficient_mass_information([changed], changed_masses, [factors[2]*rotation@signal]), expected, rtol=3e-14)
    permutation = [2, 3, 0, 1, 4, 5]
    permuted = common_acceleration(reading.data[permutation[:4]].reshape(2, 2), reading.data[4:],
                                   reading.covariance[np.ix_(permutation, permutation)], body_indices=[1, 0])
    assert gaussian_profile_cost(permuted, masses) == pytest.approx(gaussian_profile_cost(reading, masses), abs=2e-14)


def _product_prior():
    first, second = np.array([1., 2., 4.]), np.array([2., 3.])
    atoms = np.array([[x, y] for x in first for y in second])
    weights = np.outer([.2, .5, .3], [.4, .6]).ravel()
    return FiniteJointMassLaw(atoms, weights)


def test_finite_profile_update_is_proper_retains_log_evidence_and_is_zero_invariant():
    prior = _product_prior()
    reading = common_acceleration([[1000.], [-1000.]], [.2], np.eye(3))
    posterior = gaussian_joint_posterior([reading], prior)
    assert np.all(posterior.probabilities > 0)
    assert math.isfinite(posterior.log_evidence) and posterior.log_evidence < -1000
    assert posterior.probabilities.sum() == pytest.approx(1., abs=1e-15)
    # Ordinary factors exp(-cost/2) underflow, while the normalized law exists.
    assert all(math.exp(-gaussian_profile_cost(reading, atom)/2) == 0 for atom in prior.atoms)
    zero = common_acceleration([[0.], [0.]], [0.], np.array([[2., .1, .2], [.1, 1., -.1], [.2, -.1, 1.5]]))
    padded = gaussian_joint_posterior([reading, zero], prior)
    np.testing.assert_array_equal(padded.atoms, posterior.atoms)
    np.testing.assert_array_equal(padded.probabilities, posterior.probabilities)
    assert padded.log_evidence == posterior.log_evidence
    no_data = gaussian_joint_posterior([zero], prior)
    np.testing.assert_allclose(no_data.probabilities, prior.probabilities, atol=1e-16)
    assert no_data.log_evidence == pytest.approx(0., abs=3e-16)


def test_huge_constant_cost_keeps_the_supplied_unequal_prior_weights():
    reading = common_acceleration([[1e100], [-1e100]], [0.], np.eye(3))
    # For every candidate (s,s), the physical datum is exactly orthogonal
    # to (s,s,1). The true cost is the same 2e200, not evidence about s.
    for atoms, probabilities in (([[1., 1.], [2., 2.]], [.2, .8]),
                                 ([[1., 1.], [2., 2.], [3., 3.]], [.1, .2, .7])):
        prior = FiniteJointMassLaw(atoms, probabilities)
        costs = np.array([gaussian_profile_cost(reading, mass) for mass in prior.atoms])
        np.testing.assert_array_equal(costs, np.full(len(atoms), 2e200))
        posterior = gaussian_joint_posterior([reading], prior)
        np.testing.assert_allclose(posterior.probabilities, prior.probabilities, atol=1e-16)
        assert posterior.log_evidence == -1e200


def test_product_separation_is_derived_only_for_an_explicit_product_experiment_and_prior():
    prior = _product_prior()
    first = common_acceleration([[2.5]], [1.2], [[.8, .1], [.1, .4]], body_indices=[0])
    second = common_acceleration([[3.7]], [.9], [[1.2, -.2], [-.2, .6]], body_indices=[1])
    actual = gaussian_joint_posterior([first, second], prior)
    first_weights = np.array([.2, .5, .3])*np.exp(-np.array([gaussian_profile_cost(first, [x, 1.]) for x in [1., 2., 4.]])/2)
    second_weights = np.array([.4, .6])*np.exp(-np.array([gaussian_profile_cost(second, [1., x]) for x in [2., 3.]])/2)
    expected = np.outer(first_weights/first_weights.sum(), second_weights/second_weights.sum()).ravel()
    np.testing.assert_allclose(actual.probabilities, expected, atol=2e-16)


def test_supplied_joint_prior_and_shared_noise_coupling_survive_the_update_and_pushforward():
    coupled = FiniteJointMassLaw([[1., 1.], [2., 2.]], [.5, .5])
    independent = FiniteJointMassLaw([[1., 1.], [1., 2.], [2., 1.], [2., 2.]], [.25]*4)
    first = scalar_composite_law(gaussian_joint_posterior([], coupled), operation="addition")
    second = scalar_composite_law(gaussian_joint_posterior([], independent), operation="addition")
    np.testing.assert_allclose(first.support, [2., 4.], atol=0.)
    np.testing.assert_allclose(first.probabilities, [.5, .5], atol=0.)
    np.testing.assert_allclose(second.probabilities, [.25, .5, .25], atol=0.)
    reading = common_acceleration([[2.], [3.]], [1.], np.diag([.2, .3, 2.]))
    posterior = gaussian_joint_posterior([reading], independent)
    table = posterior.probabilities.reshape(2, 2)
    assert abs(table[0, 0]*table[1, 1]-table[0, 1]*table[1, 0]) > 1e-4


@pytest.mark.parametrize("operation", ["addition", "parallel"])
def test_grouped_joint_pushforwards_are_associative_and_preserve_weights(operation):
    law = FiniteJointMassLaw([[1., 2., 4.], [2., 3., 5.], [3., 4., 6.]], [.2, .3, .5])
    grouped = group_joint_law(law, [[0, 1], [2]], operation=operation)
    twice = scalar_composite_law(grouped, operation=operation)
    direct = scalar_composite_law(law, operation=operation)
    np.testing.assert_allclose(twice.support, direct.support, rtol=2e-15)
    np.testing.assert_allclose(twice.probabilities, direct.probabilities, atol=0.)
    if operation == "addition":
        reference = law.atoms.sum(axis=1)
    else:
        reference = 1/np.sum(1/law.atoms, axis=1)
    np.testing.assert_allclose(direct.support, reference, rtol=2e-15)
    marginal = group_joint_law(law, [[1]], operation=operation)
    np.testing.assert_allclose(marginal.atoms[:, 0], law.atoms[:, 1], atol=0.)


def test_readings_and_joint_laws_copy_inputs_and_reject_invalid_protocols():
    forces, covariance = np.array([[2.], [3.]]), np.eye(3)
    reading = common_acceleration(forces, [1.], covariance)
    forces[:] = 99.
    covariance[:] = 99.
    np.testing.assert_allclose(reading.data, [2., 3., 1.], atol=0.)
    with pytest.raises(ValueError):
        reading.data[0] = 99.
    with pytest.raises(FrozenInstanceError):
        reading.protocol = "separate_accelerations"
    law = FiniteJointMassLaw([[2., 3.], [1., 2.], [2., 3.]], [1., 2., 3.])
    np.testing.assert_allclose(law.probabilities, [1/3, 2/3], atol=1e-16)
    with pytest.raises(ValueError):
        law.atoms[0, 0] = 99.
    with pytest.raises(ValueError):
        common_acceleration([[2.], [3.]], [1.], np.eye(4))
    with pytest.raises(ValueError):
        common_acceleration([[2.], [3.]], [1.], np.eye(3), body_indices=[0, 0])
    with pytest.raises(ValueError):
        gaussian_profile_cost(reading, [1.])
    with pytest.raises(ValueError):
        efficient_mass_information([reading], [1., 2.], [[1., 2.]])
    with pytest.raises(ValueError):
        group_joint_law(law, [[0], [0, 1]], operation="addition")
    with pytest.raises(ValueError):
        group_joint_law(law, [[0, 1]], operation="unspecified")

"""Independent scalar/tensor compatibility checks for O17's declared world."""
import math

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from met.full import Reading, isotropic
from met.tensor_geometry import (common_scalar_minimum, scalar_profile, scalar_representability,
                                 scalar_gaussian_compatibility_bound, tensor_known_cost, validate_tensor)


@pytest.mark.parametrize("tensor", [np.diag([1., 2., 4.]), 2.0 * np.eye(3)])
def test_principal_axis_is_scalar_but_mixed_eigenvalues_are_not(tensor):
    one = scalar_representability(tensor, [[2, 0, 0], [-1, 0, 0]])
    assert one["status"] == "one" and one["unique_scalar"] == pytest.approx(tensor[0, 0])
    mixed = scalar_representability(tensor, np.eye(3))
    assert mixed["scalar_exact"] == bool(np.all(np.diag(tensor) == tensor[0, 0]))
    if not mixed["scalar_exact"]:
        assert mixed["status"] == "none" and mixed["unique_scalar"] is None


def test_repeated_eigenspace_can_hide_anisotropy_in_a_plane():
    hidden = scalar_representability(np.diag([2., 2., 5.]), [[1, 0, 0], [1, 2, 0]])
    assert hidden["excitation_rank"] == 2
    assert hidden["scalar_exact"] and hidden["unique_scalar"] == pytest.approx(2)
    oblique = scalar_representability(np.diag([1., 2., 4.]), [[1, 1, 0], [2, 2, 0]])
    assert oblique["excitation_rank"] == 1
    assert not oblique["scalar_exact"]
    assert oblique["least_squares_scalar"] == pytest.approx(1.5)


def test_zero_excitation_admits_every_scalar_and_tolerance_is_relative():
    zero = scalar_representability(np.diag([1., 2., 3.]), np.zeros((2, 3)))
    assert zero["status"] == "every" and zero["unique_scalar"] is None
    for amplitude in (1.0, 1e-9, 1e9):
        result = scalar_representability(np.diag([1., 2., 3.]), amplitude * np.array([[1., 1., 0.]]))
        assert not result["scalar_exact"]
        assert result["relative_residual"] == pytest.approx(1 / math.sqrt(10))


def test_common_scalar_profile_matches_residual_covariance_cost_and_constraints():
    factor = np.eye(6) + .08 * np.tril(np.ones((6, 6)), -1)
    covariance = factor @ factor.T
    readings = [Reading([2., -.3, 1.], [1., 1.4, -.2], covariance)]
    masses = np.array([.2, 1.3, 5.])
    profile = scalar_profile(readings, masses)
    row = profile["readings"][0]
    for index, mass in enumerate(masses):
        residual = readings[0].force - mass * readings[0].acceleration
        transform = np.c_[np.eye(3), -mass * np.eye(3)]
        v = transform @ covariance @ transform.T
        expected = residual @ np.linalg.solve(v, residual)
        assert row["correction_cost"][index] == pytest.approx(expected, rel=1e-12)
    np.testing.assert_allclose(row["candidate_force"], masses[:, None] * row["candidate_acceleration"], atol=1e-13)


@pytest.mark.parametrize("forces,accelerations", [
    ([[1., 0, 0], [0, 2., 0]], [[1., 0, 0], [0, 1., 0]]),
    ([[2., 0, 0], [0, 0, 0]], [[-1., 0, 0], [0, 1., 0]]),
    ([[0., 0, 0]], [[1., 0, 0]]),
    ([[1., 0, 0]], [[0., 0, 0]]),
    ([[0., 0, 0]], [[0., 0, 0]]),
])
def test_common_scalar_infimum_matches_independent_angle_qr_fits(forces, accelerations):
    sf, sa = .3, .7
    readings = [isotropic(f, a, sf, sa) for f, a in zip(forces, accelerations)]
    stacked = np.concatenate([np.r_[np.asarray(f) / sf, np.asarray(a) / sa]
                              for f, a in zip(forces, accelerations)])

    def cost(phi):
        # Independent latent candidate matrix at the compact angle, including endpoints.
        blocks = []
        for index in range(len(readings)):
            block = np.zeros((6, 3 * len(readings)))
            block[:, 3 * index:3 * (index + 1)] = np.vstack((math.sin(phi) * np.eye(3),
                                                          math.cos(phi) * np.eye(3)))
            blocks.append(block)
        design = np.vstack(blocks)
        fitted = np.linalg.lstsq(design, stacked, rcond=None)[0]
        residual = stacked - design @ fitted
        return float(residual @ residual)

    interior = minimize_scalar(cost, bounds=(0, math.pi / 2), method="bounded", options={"xatol": 1e-13})
    minimum = min(cost(0), cost(math.pi / 2), interior.fun)
    found = common_scalar_minimum(readings)
    assert found["positive_min_cost"] == pytest.approx(minimum, abs=1e-9)
    if found["positive_min_location"] != "interior":
        assert found["positive_min_mass"] is None


def test_axes_have_zero_individual_q_but_positive_common_scalar_cost():
    readings = [isotropic([1, 0, 0], [1, 0, 0], .1, .1),
                isotropic([0, 2, 0], [0, 1, 0], .1, .1)]
    found = common_scalar_minimum(readings)
    assert found["sum_individual_projective_min_cost"] == pytest.approx(0, abs=1e-25)
    assert found["positive_min_cost"] > 0
    assert found["extra_common_mass_cost"] == pytest.approx(found["positive_min_cost"])


def test_common_minimum_is_bounded_by_every_fixed_mass_cost_for_the_conservative_reference():
    readings = [isotropic([1.1, .1, -.1], [1.0, -.1, .1], .1, .1),
                isotropic([0, 2.1, .2], [.1, 1., -.1], .1, .1)]
    minimum = common_scalar_minimum(readings)["positive_min_cost"]
    fixed = scalar_profile(readings, [.2, 1., 1.5, 2., 5.])["total_cost"]
    assert np.all(minimum <= fixed + 1e-10)
    bounded = scalar_gaussian_compatibility_bound(readings)
    assert bounded["reference_degrees_of_freedom"] == 6
    assert 0 <= bounded["conservative_p"] <= 1


def test_independent_compatibility_diagnostic_distinguishes_exact_tensor_and_scalar_fits():
    tensor = np.diag([1., 2., 4.])
    acceleration = np.array([1., 1., 0.]) / np.sqrt(2)
    readings = [isotropic(tensor @ acceleration, acceleration, .03, .03)] * 8
    assert tensor_known_cost(readings, tensor)["total_cost"] < 1e-20
    assert scalar_gaussian_compatibility_bound(readings)["exceeds_cutoff"]
    control = [isotropic(2 * acceleration, acceleration, .03, .03)] * 8
    assert scalar_gaussian_compatibility_bound(control)["conservative_p"] == pytest.approx(1)


def test_known_tensor_exact_cost_is_zero_while_oblique_scalar_cost_is_positive():
    tensor = np.diag([1., 2., 4.])
    accelerations = np.tile([1., 1., 0.], (3, 1))
    readings = [isotropic(tensor @ a, a, .1, .1) for a in accelerations]
    known = tensor_known_cost(readings, tensor)
    assert known["total_cost"] == pytest.approx(0, abs=1e-25)
    assert known["residual_degrees_of_freedom"] == 9
    assert common_scalar_minimum(readings)["positive_min_cost"] > 0


def test_known_tensor_cost_corrections_commute_with_a_common_frame_rotation():
    tensor = np.array([[2., .2, 0], [.2, 1., .1], [0, .1, 3.]])
    factor = np.eye(6) + .07 * np.tril(np.ones((6, 6)), -1)
    covariance = factor @ factor.T
    readings = [Reading([2.3, -.2, 1.5], [1., .7, -.1], covariance)]
    phi = .51
    change = np.array([[math.cos(phi), -math.sin(phi), 0], [math.sin(phi), math.cos(phi), 0], [0, 0, 1]])
    block = np.block([[change, np.zeros((3, 3))], [np.zeros((3, 3)), change]])
    rotated = [Reading(change @ r.force, change @ r.acceleration, block @ r.cov @ block.T) for r in readings]
    first = tensor_known_cost(readings, tensor)
    second = tensor_known_cost(rotated, change @ tensor @ change.T)
    assert second["total_cost"] == pytest.approx(first["total_cost"], rel=1e-13)
    np.testing.assert_allclose(second["force_correction"], first["force_correction"] @ change.T, atol=1e-13)
    np.testing.assert_allclose(second["candidate_force"], second["candidate_acceleration"] @ (change @ tensor @ change.T).T, atol=1e-13)


def test_known_tensor_gaussian_quadratic_has_three_unfitted_degrees_per_reading():
    """The independent Gaussian quadratic is a rank-three orthogonal projection.

    Its mean trace(P)=3 and variance 2 trace(P²)=6 establish the unfitted
    chi-square reference without sampling or estimating the tensor from data.
    """
    tensor = np.array([[2., .2, 0], [.2, 1., .1], [0, .1, 3.]])
    factor = np.eye(6) + .07 * np.tril(np.ones((6, 6)), -1)
    covariance = factor @ factor.T
    constraint = np.c_[np.eye(3), -tensor]
    propagated = constraint @ covariance @ constraint.T
    projection = factor.T @ constraint.T @ np.linalg.solve(propagated, constraint @ factor)
    np.testing.assert_allclose(projection @ projection, projection, atol=1e-13)
    assert np.trace(projection) == pytest.approx(3)
    assert 2 * np.trace(projection @ projection) == pytest.approx(6)
    acceleration = np.array([1., -.3, .8])
    center = np.r_[tensor @ acceleration, acceleration]
    sum_cost = 0.0
    for noise in factor.T:
        observed = center + noise
        sum_cost += tensor_known_cost([Reading(observed[:3], observed[3:], covariance)], tensor)["total_cost"]
    assert sum_cost == pytest.approx(3, rel=1e-13)


@pytest.mark.parametrize("tensor", [np.eye(2), np.diag([1., 0., 2.]), np.diag([1., -1., 2.]),
                                     [[1., .3, 0], [0, 2., 0], [0, 0, 3.]], np.full((3, 3), np.nan)])
def test_non_spd_tensors_are_rejected(tensor):
    with pytest.raises(ValueError):
        validate_tensor(tensor)


def test_common_minimum_rejects_noise_pairs_it_cannot_reduce_to_one_gram():
    readings = [isotropic([1, 0, 0], [1, 0, 0], .1, .1),
                isotropic([0, 1, 0], [0, 1, 0], .2, .1)]
    with pytest.raises(ValueError, match="one shared"):
        common_scalar_minimum(readings)


def test_restricted_candidate_subspace_requires_different_geometry():
    reading = Reading([1., 0, 0], [1., 0, 0], np.eye(6), subspace=np.array([[1.], [0.], [0.]]))
    with pytest.raises(ValueError, match="full 3D"):
        scalar_profile([reading], [1.0])

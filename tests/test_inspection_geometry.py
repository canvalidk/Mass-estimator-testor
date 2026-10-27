"""Independent constrained-fit identities for the case inspector's geometry.

Expectations come from linear constraints, covariance propagation and orthogonal
projection, not from a second implementation of the mass equation.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from met.inspection.geometry import profile


def observation(force, acceleration, covariance, subspace=None):
    dimension = len(force)
    return SimpleNamespace(force=np.asarray(force, dtype=float),
                           acceleration=np.asarray(acceleration, dtype=float),
                           covariance=np.asarray(covariance, dtype=float),
                           subspace=np.eye(dimension) if subspace is None else np.asarray(subspace, dtype=float))


def correlated_observation():
    mixing = np.array([[1.1, 0.2, 0.0, 0.3, 0.0, 0.1],
                       [0.1, 0.9, 0.2, 0.0, -0.2, 0.0],
                       [0.0, 0.1, 1.3, 0.1, 0.0, 0.2],
                       [0.2, 0.0, 0.1, 0.7, 0.1, 0.0],
                       [0.0, -0.1, 0.0, 0.2, 0.8, 0.1],
                       [0.1, 0.0, 0.2, 0.0, 0.1, 0.6]])
    covariance = mixing @ mixing.T + 0.2 * np.eye(6)
    plane = np.array([[1.0, 0.0], [0.0, 0.6], [0.0, 0.8]])
    return observation([3.1, -0.7, 1.2], [0.9, 0.4, -0.3], covariance, plane)


def test_correlated_restricted_fit_matches_independent_linear_constraints():
    """KKT projection onto F*=m a*, B.T a*=0 is an independent reference."""
    observed = correlated_observation()
    masses = np.array([0.2, 0.7, 2.0, 8.0])
    got = profile(observed, masses)
    complement = np.array([[0.0], [-0.8], [0.6]])
    y = np.r_[observed.force, observed.acceleration]
    for index, mass in enumerate(masses):
        constraints = np.block([[np.eye(3), -mass * np.eye(3)],
                                [np.zeros((1, 3)), complement.T]])
        discrepancy = constraints @ y
        multiplier = np.linalg.solve(constraints @ observed.covariance @ constraints.T, discrepancy)
        correction = observed.covariance @ constraints.T @ multiplier
        fitted = y - correction
        np.testing.assert_allclose(got["candidate_force"][index], fitted[:3], atol=2e-14)
        np.testing.assert_allclose(got["candidate_acceleration"][index], fitted[3:], atol=2e-14)
        assert got["declared_cost"][index] == pytest.approx(discrepancy @ multiplier, rel=2e-13)
        candidate = np.r_[mass * observed.subspace, observed.subspace]
        np.testing.assert_allclose(candidate.T @ np.linalg.solve(observed.covariance, correction), 0, atol=1e-13)
        expected_covariance = np.linalg.inv(candidate.T @ np.linalg.solve(observed.covariance, candidate))
        np.testing.assert_allclose(got["conditional_covariance"][index], expected_covariance, rtol=2e-13, atol=2e-14)


def test_zero_newton_residual_does_not_establish_subspace_compatibility():
    """The latent pair can satisfy Newton II while violating declared directions."""
    observed = observation([0, 2], [0, 1], np.eye(4), [[1], [0]])
    got = profile(observed, [2])
    np.testing.assert_array_equal(got["residual"], [[0, 0]])
    assert got["residual_cost"][0] == 0
    assert got["declared_cost"][0] == pytest.approx(5)
    assert got["outside_subspace_cost"] == pytest.approx(5)
    assert got["within_subspace_cost"][0] == pytest.approx(0)


def test_full_subspace_fit_equals_residual_covariance_projection():
    """For N=I the residual constraint contains every restriction on a pair."""
    source = correlated_observation()
    observed = observation(source.force, source.acceleration, source.covariance)
    masses = np.array([0.1, 1.0, 6.0])
    got = profile(observed, masses)
    for index, mass in enumerate(masses):
        constraint = np.c_[np.eye(3), -mass * np.eye(3)]
        residual = observed.force - mass * observed.acceleration
        covariance = constraint @ observed.covariance @ constraint.T
        correction = observed.covariance @ constraint.T @ np.linalg.solve(covariance, residual)
        np.testing.assert_allclose(got["force_correction"][index], correction[:3], atol=2e-14)
        np.testing.assert_allclose(got["acceleration_correction"][index], correction[3:], atol=2e-14)
    np.testing.assert_allclose(got["declared_cost"], got["residual_cost"], rtol=2e-13, atol=2e-14)
    assert got["outside_subspace_cost"] < 1e-25


def test_projection_energies_and_constraints_are_consistent():
    """Pythagoras and constraint inclusion are algebraic authorities."""
    observed = correlated_observation()
    masses = np.array([0.15, 0.6, 2.5, 9.0])
    got = profile(observed, masses)
    np.testing.assert_allclose(got["declared_cost"],
                               got["outside_subspace_cost"] + got["within_subspace_cost"], rtol=2e-13)
    np.testing.assert_allclose(got["declared_cost"] + got["explained_energy"], got["data_energy"], rtol=2e-13)
    assert got["span_energy"] + got["outside_subspace_cost"] == pytest.approx(got["data_energy"], rel=2e-13)
    assert np.all(got["declared_cost"] >= got["residual_cost"] - 1e-12)
    np.testing.assert_allclose(got["candidate_force"], masses[:, None] * got["candidate_acceleration"], atol=2e-14)
    np.testing.assert_allclose(got["candidate_acceleration"], got["conditional_mean"] @ observed.subspace.T, atol=2e-14)


def test_frame_rotation_preserves_costs_and_rotates_vectors():
    observed = correlated_observation()
    rotation = np.array([[0.0, -1.0, 0.0], [0.6, 0.0, -0.8], [0.8, 0.0, 0.6]])
    transform = np.kron(np.eye(2), rotation)
    rotated = observation(rotation @ observed.force, rotation @ observed.acceleration,
                          transform @ observed.covariance @ transform.T, rotation @ observed.subspace)
    a, b = profile(observed, [0.3, 1.5, 4.0]), profile(rotated, [0.3, 1.5, 4.0])
    for name in ("residual_cost", "declared_cost", "within_subspace_cost", "explained_energy",
                 "data_energy", "span_energy", "outside_subspace_cost", "conditional_mean", "conditional_covariance"):
        np.testing.assert_allclose(b[name], a[name], rtol=2e-13, atol=2e-14)
    for name in ("residual", "candidate_force", "candidate_acceleration", "force_correction", "acceleration_correction"):
        np.testing.assert_allclose(b[name], a[name] @ rotation.T, atol=2e-14)
    np.testing.assert_allclose(b["residual_covariance"], rotation @ a["residual_covariance"] @ rotation.T, atol=2e-14)


def test_changing_channel_units_transforms_mass_and_covariance_together():
    observed = correlated_observation()
    force_scale, acceleration_scale = 7.0, 0.4
    mass_scale = force_scale / acceleration_scale
    transform = np.diag([force_scale] * 3 + [acceleration_scale] * 3)
    changed = observation(force_scale * observed.force, acceleration_scale * observed.acceleration,
                          transform @ observed.covariance @ transform.T, observed.subspace)
    masses = np.array([0.2, 1.0, 3.0])
    a, b = profile(observed, masses), profile(changed, mass_scale * masses)
    for name in ("residual_cost", "declared_cost", "within_subspace_cost", "explained_energy",
                 "data_energy", "span_energy", "outside_subspace_cost"):
        np.testing.assert_allclose(b[name], a[name], rtol=2e-13, atol=2e-14)
    for name in ("residual", "candidate_force", "force_correction"):
        np.testing.assert_allclose(b[name], force_scale * a[name], atol=2e-13)
    for name in ("candidate_acceleration", "acceleration_correction", "conditional_mean"):
        np.testing.assert_allclose(b[name], acceleration_scale * a[name], atol=2e-14)
    np.testing.assert_allclose(b["residual_covariance"], force_scale ** 2 * a["residual_covariance"], rtol=2e-13)
    np.testing.assert_allclose(b["conditional_covariance"], acceleration_scale ** 2 * a["conditional_covariance"], rtol=2e-13)


def test_geometry_is_not_restricted_to_estimator_dimensions():
    observed = observation([2, 1, 0, -1], [1, 0, 0.5, -0.5], np.eye(8))
    got = profile(observed, [0.5, 2.0])
    assert got["candidate_acceleration"].shape == (2, 4)
    assert got["conditional_covariance"].shape == (2, 4, 4)
    np.testing.assert_allclose(got["declared_cost"], got["residual_cost"], rtol=2e-13)


@pytest.mark.parametrize("masses", [[], 1.0, [[1.0]], [0.0], [-1.0], [np.inf], [np.nan], [1.0, 0.0]])
def test_invalid_candidate_grid_is_rejected(masses):
    observed = observation([2.0], [1.0], np.eye(2))
    with pytest.raises(ValueError, match="masses"):
        profile(observed, masses)

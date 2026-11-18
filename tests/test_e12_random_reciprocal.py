"""Geometry checks for the shared reciprocal direction law in the E12 exploration."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtri


def direction_density(gram, mass):
    force_sq, dot, accel_sq = gram[0, 0], gram[0, 1], gram[1, 1]
    return math.sqrt(np.linalg.det(gram)) / (math.pi * (
        accel_sq * mass**2 - 2 * dot * mass + force_sq))


def directions_in_plane(count):
    angles = np.pi * (np.arange(count) + 0.5) / count
    return np.column_stack((np.cos(angles), np.sin(angles)))


@pytest.mark.parametrize("force", [[4., 1.], [-4., -1.], [4., -1.]])
def test_direction_sweep_matches_positive_probability_shape_and_median(force):
    force = np.array(force)
    accel = np.array([1., 4.])
    pair = np.column_stack((force, accel))
    gram = pair.T @ pair
    directions = directions_in_plane(100003)
    measured = (directions @ force) / (directions @ accel)

    rho = (force @ accel) / (np.linalg.norm(accel) * np.linalg.norm(force))
    positive_prob = .5 + math.asin(rho) / math.pi
    assert np.mean(measured > 0) == pytest.approx(positive_prob, abs=2 / len(measured))
    assert quad(lambda m: direction_density(gram, m), 0, np.inf)[0] == pytest.approx(positive_prob)
    median = np.linalg.norm(force) / np.linalg.norm(accel)
    half_prob = quad(lambda m: direction_density(gram, m), 0, median)[0]
    assert half_prob == pytest.approx(positive_prob / 2)
    for cut in [-.7, .2, 2.5]:
        analytic = quad(lambda m: direction_density(gram, m), -np.inf, cut)[0]
        assert np.mean(measured <= cut) == pytest.approx(analytic, abs=2 / len(measured))


@pytest.mark.parametrize("change", [
    np.array([[2., 0.], [0., .25]]),
    np.array([[0., 1.], [1., 0.]]),
    np.array([[2., -1.], [1., 2.]]),
])
def test_full_line_channel_transform_uses_projective_jacobian(change):
    gram = np.array([[17., 8.], [8., 13.]])
    mass = np.array([-.3, .2, 1.7])
    transformed_gram = change @ gram @ change.T
    denominator = change[1, 0] * mass + change[1, 1]
    new_mass = (change[0, 0] * mass + change[0, 1]) / denominator
    jacobian = np.abs(np.linalg.det(change)) / denominator**2
    np.testing.assert_allclose(direction_density(transformed_gram, new_mass) * jacobian,
                               direction_density(gram, mass), rtol=1e-13)


def test_independent_direction_law_keeps_the_conflict_lost_in_relative_evidence():
    directions = directions_in_plane(4099)
    aligned = np.column_stack(([2., 2.], [2., 2.]))
    conflict = np.column_stack(([4., 1.], [1., 4.]))
    np.testing.assert_array_equal(conflict.T @ conflict - aligned.T @ aligned, np.eye(2) * 9)
    aligned_scores = directions @ aligned
    np.testing.assert_array_equal(aligned_scores[:, 0] / aligned_scores[:, 1], 1)
    conflict_scores = directions @ conflict
    positive = conflict_scores[conflict_scores[:, 0] * conflict_scores[:, 1] > 0]
    positive_masses = positive[:, 0] / positive[:, 1]
    assert np.quantile(positive_masses, .25) < .8 < 1.2 < np.quantile(positive_masses, .75)


def test_reciprocal_lift_respects_frame_units_swap_and_minimum_norm():
    pair = np.array([[4., 1.], [1., 4.], [0., 0.]])
    normal = np.array([0., 0., 1.])
    rotation = np.eye(3) - 2 * np.outer([1., 2., 2.], [1., 2., 2.]) / 9

    def lift(data, mass):
        return data @ np.linalg.solve(data.T @ data, [mass, 1.])

    for mass in [.03, 1., 3.]:
        covector = lift(pair, mass)
        np.testing.assert_allclose(covector @ pair, [mass, 1.])
        np.testing.assert_allclose((covector + 2 * normal) @ pair, [mass, 1.])
        assert np.linalg.norm(covector + 2 * normal) > np.linalg.norm(covector)
        np.testing.assert_allclose(lift(rotation @ pair, mass), rotation @ covector)
        np.testing.assert_allclose(lift(pair[:, ::-1], 1 / mass), covector / mass)
        np.testing.assert_allclose(lift(pair * [2., 3.], 2 * mass / 3), covector / 3)


def test_antiparallel_positive_limit_and_axes():
    rho = -1 + 1e-9
    normalizer = math.sqrt((1 - rho) * (1 + rho)) / math.acos(-rho)
    for x in [.01, .5, 1., 3., 100.]:
        density = normalizer / (1 - 2 * rho * x + x**2)
        assert density == pytest.approx(1 / (1 + x)**2, rel=3e-9)
    directions = directions_in_plane(100)
    accel = np.array([1., 0.])
    anti = -3 * accel
    np.testing.assert_allclose((directions @ anti) / (directions @ accel), -3)
    np.testing.assert_array_equal(directions @ np.zeros(2), 0)


def test_identifying_gaussian_array_has_inconsistent_geometric_limit():
    # A,B,C per reading in the strong-law limit with truth (2,1) and identity noise.
    gram = np.outer([2., 1.], [2., 1.]) + np.eye(2)
    mass_limit = math.sqrt(gram[0, 0] / gram[1, 1])
    alive = quad(lambda m: direction_density(gram, m), 0, np.inf)[0]
    left = quad(lambda m: direction_density(gram, m), 0, mass_limit)[0]
    assert left == pytest.approx(alive / 2)
    assert mass_limit == pytest.approx(math.sqrt(5 / 2))
    assert mass_limit < 2
    assert math.sqrt(np.linalg.det(gram)) / gram[1, 1] == pytest.approx(math.sqrt(1.5))
    for n in [10., 10000.]:
        assert direction_density(n * gram, 1.) == pytest.approx(direction_density(gram, 1.))


def test_concentrating_direction_law_has_cauchy_width_ratios_despite_positive_cut():
    mu = 2.
    target = math.tan(.475 * math.pi) / math.tan(.25 * math.pi)
    normal_ratio = ndtri(.975) / ndtri(.75)
    ranks = np.array([.025, .25, .75, .975])
    errors = []
    for gamma in [.05, .005, .0005, .00005]:
        negative_mass = .5 + math.atan(-mu / gamma) / math.pi
        probabilities = negative_mass + (1 - negative_mass) * ranks
        quantiles = mu + gamma * np.tan(np.pi * (probabilities - .5))
        observed = (quantiles[3] - quantiles[0]) / (quantiles[2] - quantiles[1])
        errors.append(abs(observed - target))
    assert np.all(np.diff(errors) < 0)
    assert errors[-1] < .003
    assert target > 4 * normal_ratio


@pytest.mark.parametrize("angle", [.13, .74])
def test_transport_of_two_angular_atoms_with_the_instrument_tensor(angle):
    def angular_atoms(gram, covariance):
        root = np.linalg.cholesky(covariance)
        white = np.linalg.inv(root)
        eigenvalues, axes = np.linalg.eigh(white @ gram @ white.T)
        principal, transverse = axes[:, 1], axes[:, 0]
        directions = np.array([np.cos(angle) * principal + sign * np.sin(angle) * transverse
                               for sign in [-1, 1]])
        return eigenvalues, directions @ root.T

    def line_projectors(rows):
        unit = rows / np.linalg.norm(rows, axis=1, keepdims=True)
        return np.einsum("ij,ik->ijk", unit, unit)

    gram = np.array([[17., 8.], [8., 13.]])
    covariance = np.array([[3., 1.], [1., 2.]])
    change = np.array([[-2., .3], [1., 1.7]])
    original_values, original_lines = angular_atoms(gram, covariance)
    transformed_values, transformed_lines = angular_atoms(change @ gram @ change.T,
                                                          change @ covariance @ change.T)
    np.testing.assert_allclose(original_values, transformed_values)
    expected = line_projectors(original_lines @ change.T)
    observed = line_projectors(transformed_lines)
    assert np.allclose(expected, observed) or np.allclose(expected, observed[::-1])
    for strength in [.1, 100.]:
        np.testing.assert_allclose(direction_density(strength * covariance, .8),
                                   direction_density(covariance, .8))

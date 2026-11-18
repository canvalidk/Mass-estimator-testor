"""Finite checks of the covector reconstruction in explore/e12-covector-information.md."""

import itertools

import numpy as np
import pytest
from scipy.linalg import eig


def correlated_instrument():
    return (np.array([[2., .3], [.3, 3.7]]),
            np.array([[.4, -.15], [.8, .1]]),
            np.array([[1., .1], [.1, .8]]))


def residual_cost(blocks, reading, mass):
    a, b, d = blocks
    force, accel = np.split(reading, 2)
    residual = force - mass * accel
    covariance = a - mass * (b + b.T) + mass ** 2 * d
    return residual @ np.linalg.solve(covariance, residual)


def contraction_coordinates(blocks):
    """A linearization gives w in the first k entries of each eigenvector."""
    a, b, d = blocks
    k = len(d)
    eye = np.eye(k)
    zero = np.zeros_like(d)
    left = np.block([[zero, eye], [-a, b + b.T]])
    right = np.block([[eye, zero], [zero, d]])
    roots, vectors = eig(left, right)
    upper = np.flatnonzero(roots.imag > 0)
    upper = upper[np.argsort(roots[upper].imag)]
    assert len(upper) == k
    rows = []
    for root, w in zip(roots[upper], vectors[:k, upper].T):
        np.testing.assert_allclose((a - root * (b + b.T) + root**2 * d) @ w, 0, atol=1e-12)
        # Polynomial continuation uses bilinear transpose even for complex roots.
        covector = np.concatenate([w, -root * w])
        rows.extend([covector.real, covector.imag])
    return roots[upper], np.stack(rows)


def sign_alternatives(coordinates, reading):
    values = coordinates @ reading
    for signs in itertools.product((1, -1), repeat=len(reading) // 2):
        signed = values * np.repeat(signs, 2)
        yield signs, np.linalg.solve(coordinates, signed)


def physical_pairings(reading):
    force, accel = np.split(reading, 2)
    product = accel @ force
    return np.array([product / (accel @ accel), product / (force @ force)])


@pytest.mark.parametrize("mass", [.03, 1., 3., 29.])
def test_directional_covector_attains_the_full_residual(mass):
    blocks = correlated_instrument()
    a, b, d = blocks
    reading = np.array([1., -2., .7, 1.2])
    force, accel = np.split(reading, 2)
    r = force - mass * accel
    v = a - mass * (b + b.T) + mass**2 * d
    witness = np.linalg.solve(v, r)
    k = residual_cost(blocks, reading, mass)
    assert (witness @ r)**2 / (witness @ v @ witness) == pytest.approx(k)
    angles = np.arange(100) * np.pi / 100
    directions = np.column_stack([np.cos(angles), np.sin(angles)])
    projected = (directions @ r)**2 / np.einsum("ij,jk,ik->i", directions, v, directions)
    assert np.all(projected <= k + 1e-13)


def test_simple_correlated_spectrum_and_all_four_relative_alternatives():
    blocks = correlated_instrument()
    a, b, d = blocks
    assert np.linalg.eigvalsh(np.block([[a, b], [b.T, d]])).min() > .62
    roots, coordinates = contraction_coordinates(blocks)
    np.testing.assert_allclose(roots, [.474586 + 1.428980j, -.025219 + 2.020063j], atol=1e-6)
    reading = np.array([1., -2., .7, 1.2])
    offsets = []
    candidates = []
    for signs, other in sign_alternatives(coordinates, reading):
        delta = residual_cost(blocks, other, 0) - residual_cost(blocks, reading, 0)
        offsets.append(delta)
        candidates.append(other)
        for t in [0.02, .3, 1., 4., 11., 1000.]:
            assert residual_cost(blocks, other, t) - residual_cost(blocks, reading, t) == pytest.approx(delta, abs=1e-12)
        if signs == (-1, 1):
            np.testing.assert_allclose(other, [-1.049856, -5.130173, -2.654676, .429211], atol=1e-6)
    assert len(np.unique(np.round(candidates, 7), axis=0)) == 4
    np.testing.assert_allclose(offsets, [0, 5.551836, 5.551836, 0], atol=1e-6)
    x_map = np.linalg.inv(coordinates)[2:]
    gram = x_map.T @ np.linalg.solve(d, x_map)
    assert np.linalg.norm(gram[:2, 2:]) > .1


def test_sheared_instrument_absolute_curve_fails_to_determine_spatial_pairings():
    shear = np.array([[1., 1.], [0., 1.]])
    blocks = (shear @ np.diag([1., 4.]) @ shear.T, np.zeros((2, 2)), shear @ shear.T)
    first = np.array([3., 2., 2., 1.])
    second = np.array([-1., -2., 0., -1.])
    np.testing.assert_allclose(physical_pairings(first), [8 / 5, 8 / 13])
    np.testing.assert_allclose(physical_pairings(second), [2., 2 / 5])
    for mass in [0., .01, .5, 1., 2., 5., 1e3]:
        expected = (1 - mass)**2 / (1 + mass**2) + (2 - mass)**2 / (4 + mass**2)
        for reading in [first, second]:
            assert residual_cost(blocks, reading, mass) == pytest.approx(expected)
    _, coordinates = contraction_coordinates(blocks)
    x_map = np.linalg.inv(coordinates)[2:]
    gram = x_map.T @ np.linalg.solve(blocks[2], x_map)
    np.testing.assert_allclose(gram[:2, 2:], 0, atol=1e-12)
    alternatives = [candidate for _, candidate in sign_alternatives(coordinates, first)]
    assert any(np.allclose(second, candidate) for candidate in alternatives)


def test_repeated_noise_profile_hides_one_free_gram_trace():
    blocks = (np.diag([1., 1., 4.]), np.zeros((3, 3)), np.eye(3))
    force = np.array([2., 2., 3.])
    accel = np.array([2., 2., -1.])
    first = np.r_[force, accel]
    for increment in [.125, 1., 9.]:
        # Realize the entire family of positive trace increments in the two equal modes.
        gram = np.array([[8., 8.], [8., 8.]]) + increment * np.eye(2)
        new_force, new_accel = np.linalg.cholesky(gram)
        other = np.r_[new_force, force[2], new_accel, accel[2]]
        for mass in [.05, .5, 1., 2., 6.]:
            assert residual_cost(blocks, other, mass) - residual_cost(blocks, first, mass) == pytest.approx(increment)

"""Finite identities supporting the geometry and the two-experiment obstruction."""
import math

import numpy as np
import pytest


def shape(b):
    return np.array([b[0, 0] - b[1, 1], 2 * b[0, 1]]) / np.trace(b)


@pytest.mark.parametrize("p", [[-0.6, 0.2], [0.0, 0.1], [0.28, 0.9]])
def test_every_positive_cross_shape_is_both_a_fixed_conflict_and_a_population_scatter(p):
    radius = np.linalg.norm(p)
    angle = math.atan2(p[1], p[0]) / 2
    z = math.sqrt(2 * radius / (1 - radius)) * np.array([math.cos(angle), math.sin(angle)])
    b = np.eye(2) + np.outer(z, z)
    rows = np.linalg.cholesky(b).T
    assert np.linalg.det(rows) != 0
    np.testing.assert_allclose(shape(b), p, atol=2e-15)
    np.testing.assert_allclose(rows.T @ rows, b)


def test_regular_noisy_design_and_shrinking_two_reading_conflict_have_same_exact_scatter():
    rows = np.array([[3, 2], [3, 0], [1, 2], [1, 0]])
    pair = np.linalg.cholesky(rows.T @ rows).T
    for n in [1, 16, 81]:
        repeated = np.tile(rows, (n, 1))
        shrinking_conflict = math.sqrt(n) * pair
        padded = np.pad(shrinking_conflict, ((0, 4 * n - 2), (0, 0)))
        np.testing.assert_allclose(repeated.T @ repeated, padded.T @ padded)


def test_planar_rotations_of_scalar_pairs_concentrate_any_full_rank_block_into_two_rows():
    # Thin QR using only local two-row rotations, which are spatial frame changes when packaged.
    raw = np.array([[3, 2], [-1, 1], [4, 0], [2, -3]], dtype=float)
    reduced = raw.copy()
    for col in range(2):
        for row in range(len(reduced) - 1, col, -1):
            x, y = reduced[[col, row], col]
            if x == y == 0:
                continue
            rotation = np.array([[x, y], [-y, x]]) / math.hypot(x, y)
            reduced[[col, row]] = rotation @ reduced[[col, row]]
    np.testing.assert_allclose(reduced[2:], 0, atol=1e-14)
    np.testing.assert_allclose(reduced[:2].T @ reduced[:2], raw.T @ raw)
    np.testing.assert_allclose(reduced[:2], np.linalg.cholesky(raw.T @ raw).T, atol=1e-14)

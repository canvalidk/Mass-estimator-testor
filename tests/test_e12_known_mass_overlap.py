"""Known-standard paths and agreement of evidence on their common mass interval."""
from fractions import Fraction

import numpy as np
import pytest

from explore.e12_compact_join import JoinState
from explore.e12_residual_join import residual_log_evidence


def upper(b):
    return np.array([[1, b], [0, 1]])


def lower(c):
    return upper(c).T


def mass_map(matrix, mass):
    f, a = matrix @ [mass, 1]
    return f / a


P = upper(1) @ lower(1)
Q = lower(Fraction(1, 5)) @ upper(2)
R = Fraction(1, 5) * np.array([[4, 3], [-3, 4]])
MASSES = (Fraction(1, 4), Fraction(1, 2))


def test_identical_observation_and_covariance_with_overlapping_output_masses():
    np.testing.assert_array_equal(P @ R, Q)
    np.testing.assert_array_equal(R @ R.T, np.eye(2))
    np.testing.assert_array_equal(Q @ Q.T, P @ P.T)
    point = np.array([1, 1])
    rotated = R @ point
    np.testing.assert_array_equal(rotated, [Fraction(7, 5), Fraction(1, 5)])
    np.testing.assert_array_equal(P @ rotated, Q @ point)
    assert tuple(mass_map(R, m) for m in MASSES) == (Fraction(16, 13), 2)
    for m in MASSES:
        target_mass = mass_map(Q, m)
        assert target_mass == mass_map(P, mass_map(R, m))
        assert Fraction(10, 7) < target_mass < 2


def singleton_factors(channels, masses, family, tau):
    f, a = np.asarray(channels, dtype=float)
    angles = np.arctan(np.array(masses, dtype=float))
    if family == "compact":
        return JoinState.from_reading(f, a).log_evidence(angles, tau=tau)
    matrix = np.stack([np.atleast_1d(f), np.atleast_1d(a)])
    return residual_log_evidence([matrix @ matrix.T], angles,
                                 dimension=matrix.shape[1], tau=tau)


@pytest.mark.parametrize("family", ["compact", "all-candidate"])
@pytest.mark.parametrize("tau", [0, 64])
@pytest.mark.parametrize("balanced", [False, True])
def test_incompatible_contrasts_even_with_identical_final_records(family, tau, balanced):
    channels_prime = np.c_[np.eye(2), np.zeros(2)] if balanced else np.ones(2)
    channels = np.asarray(R, dtype=float) @ channels_prime
    masses_rotated = [mass_map(R, m) for m in MASSES]
    shifted = singleton_factors(channels, masses_rotated, family, tau)
    original = singleton_factors(channels_prime, MASSES, family, tau)
    exact = Fraction(36, 425) if balanced else Fraction(4368, 53125)
    assert np.subtract(*(shifted - original)) == pytest.approx(float(exact), abs=2e-14)


def test_uniform_angle_priors_agree_on_overlap_as_measures():
    # Density of g_A(m) at the image of a unit angular prior point.
    def density_at(matrix, target_mass):
        x, y = np.linalg.solve(matrix, [target_mass, 1])
        return (2 / np.pi) / (np.linalg.det(matrix) * (x*x + y*y))

    for target in np.linspace(1.44, 1.99, 15):
        assert density_at(P.astype(float), target) == pytest.approx(
            density_at(Q.astype(float), target))


@pytest.mark.parametrize("phi", [.005, .09, .20])
def test_small_channel_rotations_realizable_on_overlap_by_positive_standards(phi):
    cosine, sine = np.cos(phi), np.sin(phi)
    rotation = np.array([[cosine, sine], [-sine, cosine]])
    base = upper(1) @ lower(1) @ upper(1)
    q = base @ rotation
    left = (q[0, 0] - 1) / q[1, 0]
    right = (q[1, 1] - 1) / q[1, 0]
    assert min(left, q[1, 0], right) > 0
    np.testing.assert_allclose(upper(left) @ lower(q[1, 0]) @ upper(right), q)
    # Rowwise spatial vectors are F=e1 and a=e2. A spatial rotation realizes the channel action.
    channels = np.c_[np.eye(2), np.zeros(2)]
    spatial = np.eye(3)
    spatial[:2, :2] = rotation
    np.testing.assert_allclose(channels @ spatial.T, rotation.T @ channels)
    half_turn = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    assert np.linalg.det(half_turn) == 1
    np.testing.assert_array_equal(channels @ half_turn.T, channels[::-1])


def angular_kernel(theta, orientation, strength, secondary):
    """An even projective kernel for which overlaps can be checked without fitting."""
    offset = theta - orientation
    return np.exp(-strength * np.sin(offset)**2) + secondary * np.exp(-strength * np.cos(offset)**2)


@pytest.mark.parametrize("secondary", [0, 1])
def test_scalar_interval_lemma_and_alternatives_at_increasing_amplitude(secondary):
    endpoint = np.pi / 2
    grid = np.linspace(0, endpoint, 20001)
    middle = len(grid) // 2
    angle = 0.31
    for strength in (30, 400):
        axis = angular_kernel(grid, 0, strength, secondary)
        anti = angular_kernel(grid, -endpoint / 2, strength, secondary)
        compatible = angular_kernel(grid, angle, strength, secondary)
        # The half-arc conditional laws in the proof coincide after shifting.
        np.testing.assert_allclose(anti[:middle + 1] / anti[:middle + 1].sum(),
                                   axis[middle:] / axis[middle:].sum(), atol=1e-15)
        np.testing.assert_allclose(anti, anti[::-1], atol=1e-14)
        if strength == 400:
            assert compatible[np.abs(grid - angle) < 0.2].sum() / compatible.sum() > 0.999
            assert axis[(grid > 0.2) & (grid < endpoint - 0.2)].sum() / axis.sum() < 1e-6
            if secondary == 1:
                assert axis[middle:].sum() / axis.sum() == pytest.approx(0.5)
                assert anti[np.abs(grid - endpoint / 2) < 0.2].sum() / anti.sum() > 0.999
            else:
                assert axis[grid < 0.2].sum() / axis.sum() > 0.999
                assert anti[grid < 0.2].sum() / anti.sum() == pytest.approx(0.5)

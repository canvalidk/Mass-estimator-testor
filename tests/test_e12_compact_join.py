"""Finite checks for the four-coordinate theorem and its comparison of laws."""
import math
from functools import reduce

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid
from scipy.stats import ncx2

from explore.e12_compact_join import JoinState, _residual_count


ANGLES = np.linspace(0, math.pi / 2, 151)


def state(rows):
    return reduce(JoinState.merge, (JoinState.from_reading(f, a) for f, a in rows))


def test_associativity_and_zero_identity():
    readings = [([3, 1], [2, 0]), ([0, -2], [5, 2]), ([-1, 4], [2, 3])]
    a, b, c = (state([row]) for row in readings)
    left, right = a.merge(b).merge(c), c.merge(b.merge(a))
    assert (left.force_energy, left.cross, left.acceleration_energy) == (31, 12, 46)
    assert left.count == pytest.approx(right.count)
    assert left == left.merge(state([([0, 0], [0, 0])]))
    np.testing.assert_allclose(left.log_evidence(ANGLES), right.log_evidence(ANGLES))


@pytest.mark.parametrize("row", [(1, 2), (-1, 3), ([1, 2], [2, -3]), ([0, 0, 0], [0, 0, 0])])
def test_same_singleton_law_for_every_member(row):
    obs = state([row])
    np.testing.assert_array_equal(obs.log_evidence(ANGLES, tau=1e6), obs.log_evidence(ANGLES))


def test_distinct_same_prior_log_odds_on_two_positive_pushes():
    obs = state([(2, 1), (1, 2)])
    angles = np.arctan([1.0, 0.5])
    unit = 1 / (1 - math.sqrt(0.5))
    c = 2 * unit * (1 - math.exp(-0.25))
    penalty = 36 / 11**3 * (c - unit * (1 - math.exp(-0.5))) / (1 + c)
    odds = [np.subtract(*obs.log_evidence(angles, tau=t)) for t in (0, 64)]
    assert odds[0] - odds[1] == pytest.approx(64 / 180 * (1 - 1 / (1 + 64 * penalty)))
    assert odds[0] > odds[1]


def test_fourth_coordinate_distinguishes_data_with_exactly_the_same_gram():
    rows = np.array([[3, 2], [-1, 1]])
    rotated = (np.array([[1, 1], [-1, 1]]) / math.sqrt(2)) @ rows
    obs, other = state(rows), state(rotated)
    np.testing.assert_allclose([obs.force_energy, obs.cross, obs.acceleration_energy],
                               [other.force_energy, other.cross, other.acceleration_energy])
    assert abs(obs.count - other.count) > 0.1
    assert np.ptp(obs.log_evidence(ANGLES, tau=64) - other.log_evidence(ANGLES, tau=64)) > 0.001


def test_covariance_under_channel_interchange_and_independent_spatial_frame_changes():
    rows = [([0, 2], [1, 3]), ([4, 2], [-1, 3])]
    obs = state(rows)
    rotations = [np.array([[0, -1], [1, 0]]), np.diag([-1, 1])]
    transformed = state([(r @ f, r @ a) for (f, a), r in zip(rows, rotations)])
    assert transformed == obs
    swapped = state([(a, f) for f, a in rows])
    np.testing.assert_allclose(swapped.log_evidence(math.pi / 2 - ANGLES, tau=64),
                               obs.log_evidence(ANGLES, tau=64), atol=1e-14)


@pytest.mark.parametrize("k,noncentrality", [(1, 0), (2, 5), (3, 1)])
def test_count_matches_declared_gaussian_experiment(k, noncentrality):
    actual = ncx2.expect(lambda cost: _residual_count(k, cost), args=(k, noncentrality))
    expected = k * (1 - 0.5**(0.5 * k) * math.exp(-noncentrality / 4)) / (1 - 0.5**(0.5 * k))
    assert actual == pytest.approx(expected, abs=2e-10)


@pytest.mark.parametrize("tau", [0, 64])
def test_incompatible_and_compatible_shape_limits(tau):
    scale = 1e4
    pair = state([(2 * scale, scale), (scale, 2 * scale)])
    limiting_precision = 2 / (1 - math.sqrt(0.5)) * (64 / 9) * (9 ** 0.125)
    sine2 = np.sin(ANGLES - math.pi / 4) ** 2
    np.testing.assert_allclose(pair.log_evidence(ANGLES, tau=tau), -0.5 * limiting_precision * sine2,
                               atol=1e-6)
    balanced = state([(scale, 0), (0, scale)])
    np.testing.assert_allclose(balanced.log_evidence(ANGLES, tau=tau), -sine2, atol=1e-8)
    axis = state([(scale, 0)])
    np.testing.assert_allclose(axis.log_evidence(ANGLES, tau=tau), -0.5 * np.cos(ANGLES) ** 2,
                               atol=1e-12)
    aligned = state([(scale, scale)])
    np.testing.assert_allclose(aligned.log_evidence(ANGLES, tau=tau), -scale**2 * sine2, atol=1e-7)


def test_regular_fixed_push_has_growing_precision_with_same_population_shape_as_conflict():
    # Balanced signs of a standard-sized transverse/parallel disturbance give B/N = I + zz^T exactly.
    z = np.array([2, 1])
    rows = z + np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]])
    repeats = 256
    repeated = state(np.tile(rows, (repeats, 1)))
    fixed = state(np.linalg.cholesky(repeats * (rows.T @ rows)).T)
    np.testing.assert_allclose([repeated.force_energy, repeated.cross, repeated.acceleration_energy],
                               [fixed.force_energy, fixed.cross, fixed.acceleration_energy])
    assert repeated.count > 100 * fixed.count
    assert np.ptp(repeated.log_evidence(ANGLES)) > 50 * np.ptp(fixed.log_evidence(ANGLES))


@pytest.mark.parametrize("forces,accelerations,s,c,tau,expected", [
    ([2, 1], [1, 2], 1, 1, 0, [0.113638, 1, 8.799860]),
    ([2, 1], [1, 2], 1, 1, 64, [0.102600, 1, 9.746570]),
    ([55, 56], [500, 500], 50, 5, 0, [5.409620, 5.549021, 5.688507]),
    ([50, 60], [500, 500], 50, 5, 0, [5.287793, 5.498295, 5.708992]),
])
def test_reported_central_mass_intervals(forces, accelerations, s, c, tau, expected):
    obs = state(zip(forces, accelerations))
    grid = np.linspace(0, math.pi / 2, 200001)
    logp = obs.log_evidence(grid, tau=tau)
    prior = 1 / ((c * np.cos(grid)) ** 2 + (s * np.sin(grid)) ** 2)
    cdf = cumulative_trapezoid(np.exp(logp - logp.max()) * prior, grid, initial=0)
    actual = s * np.tan(np.interp([0.025, 0.5, 0.975], cdf / cdf[-1], grid))
    np.testing.assert_allclose(actual, expected, atol=2e-6)


def test_fourth_coordinate_uses_only_defects_at_the_noise_scale_ratio():
    # Adding equal standardized values to the two channels changes the drive, keeping this residual.
    first = state([([2, 4], [1, 2]), ([1, 1], [0, 0])])
    second = state([([102, -46], [101, -48]), ([-2, 5], [-3, 4])])
    assert first.force_energy != second.force_energy
    assert first.count == pytest.approx(second.count)
    aligned = state([(1, 1), (2, 2)])
    assert aligned.count == 0
    almost = state([(1, 1 + 1e-5), (2, 2 - 1e-5)])
    np.testing.assert_allclose(aligned.log_evidence(ANGLES), almost.log_evidence(ANGLES), atol=3e-5)


def test_invalid_declarations_and_angles():
    with pytest.raises(ValueError):
        state([([2, 1], [1])])
    with pytest.raises(ValueError):
        state([(1, 2)]).merge(state([([1, 0], [1, 1])]))
    with pytest.raises(ValueError):
        state([(1, 2)]).log_evidence([-0.1])
    with pytest.raises(ValueError):
        state([(1, 2)]).log_evidence([0], tau=-1)


@pytest.mark.parametrize("f,a,relative_precision", [(2, 2, 1.0), (2, 1, 2 / 5)])
def test_singletons_against_liseo_reference_equation_14(f, a, relative_precision):
    # With uniform angle prior the reciprocal Jacobian in mass is already accounted for.
    liseo_log = -0.5 * (f * np.cos(ANGLES) - a * np.sin(ANGLES)) ** 2
    candidate = state([(f, a)]).log_evidence(ANGLES)
    np.testing.assert_allclose(candidate - candidate[0], relative_precision * (liseo_log - liseo_log[0]),
                               atol=1e-14)

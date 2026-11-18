"""Finite checks for reciprocal tangent pooling and the candidate-wise Gaussian theorem."""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from explore.e12_tangent_pairings import tangent_pairings
from met.conditional_score import gaussian_score_test
from met.full import Reading, _h_terms, isotropic


def energies(state):
    return np.array([state.within, state.between, state.common])


def scalar(force, accel):
    return isotropic([force], [accel], 1., 1.)


def test_reading_declaration_locates_conflict_before_normalization():
    blocks = [[scalar(2, 2), scalar(2, 2)],
              [scalar(4, 1), scalar(1, 4)],
              [isotropic([4, 1], [1, 4], 1., 1.)]]
    for rows, expected in zip(blocks, [[0, 0, 0], [0, 9, 0], [9, 0, 0]]):
        np.testing.assert_allclose(energies(tangent_pairings(rows, 1.)), expected, atol=1e-12)


def test_full_correlated_reciprocal_agrees_with_profile_and_T45():
    covariance_root = np.array([[1.2, .1, .3, 0.], [.2, .8, -.1, .3],
                                [.1, .3, 1.5, .1], [-.2, .1, .3, 1.]])
    cov = covariance_root @ covariance_root.T
    reading = Reading([2.2, -.7], [.5, 1.1], cov)
    for mass in [.2, 1.3, 4.]:
        d = np.hstack([np.eye(2), -mass * np.eye(2)])
        residual = reading.force - mass * reading.acceleration
        variance = d @ cov @ d.T
        precision_residual = np.linalg.solve(variance, residual)
        fitted = reading.acceleration - (cov[2:, :2] - mass * cov[2:, 2:]) @ precision_residual
        weight_m = fitted @ np.linalg.solve(variance, fitted)
        score_m = fitted @ precision_residual
        state = tangent_pairings([reading], mass)
        chart_scale = 1. + mass**2
        assert state.score == pytest.approx(chart_scale * score_m, abs=1e-13)
        assert state.weight == pytest.approx(chart_scale**2 * weight_m, abs=1e-13)
        assert state.correction == pytest.approx(score_m / weight_m / chart_scale)
        assert sum(energies(state)) == pytest.approx(residual @ precision_residual)
        assert state.degrees == (1, 0, 1)
        assert state.common == pytest.approx(gaussian_score_test([reading], mass)["statistic"])


def _restricted_example():
    matrix = np.array([[1., .2, .1, -.1], [-.3, 1.4, .1, .3],
                       [.1, .2, 1.5, 0.], [.2, .3, .1, .8]])
    return Reading([1.5, -.8], [.9, .5], matrix @ matrix.T, [[.8], [.6]])


def test_declared_subspace_drops_unrelated_residual_dimensions():
    reading = _restricted_example()
    mass = 1.2
    state = tangent_pairings([reading], mass)
    pre = reading.pre()
    expected = pre["energy"] - _h_terms(pre, 1., np.array([math.atan(mass)]))[0]
    assert state.degrees == (0, 0, 1)
    assert sum(energies(state)) == pytest.approx(expected, abs=1e-12)
    assert state.common == pytest.approx(gaussian_score_test([reading], mass)["statistic"])
    # A physical displacement orthogonal in the noise metric to both declared
    # channels is absent from K. A large one would dominate the ambient cost.
    extra = reading.cov @ np.r_[[-.6, .8], [1.2, -1.6]]
    moved = Reading(reading.force + 1e3 * extra[:2],
                    reading.acceleration + 1e3 * extra[2:], reading.cov, reading.subspace)
    np.testing.assert_allclose(energies(tangent_pairings([moved], mass)), energies(state), atol=1e-10)


@pytest.mark.parametrize("perturbation,active", [(0., 0), (1e-12, 1)])
def test_restricted_span_preserves_full_coordinate_tangent_rank(perturbation, active):
    cov = np.diag([2., 4., 1., 8.])
    y = cov @ np.array([1., 2., -2., -4.])
    y[0] += perturbation
    reading = Reading(y[:2], y[2:], cov, [[.6], [.8]])
    result = tangent_pairings([reading], 2.)
    reference = gaussian_score_test([reading], 2.)
    pre = reading.pre()
    cost = pre["energy"] - _h_terms(pre, 1., np.array([math.atan(2.)]))[0]
    assert result.active == active
    assert result.score == reference["score"]
    assert result.weight == reference["variance"]
    assert result.common == reference["statistic"]
    assert sum(energies(result)) == pytest.approx(cost)
    if active:
        assert 0. < result.weight < 1e-21
        assert result.common > 58.
    else:
        assert result.within > 58.


def test_candidate_merger_agrees_in_multiple_orders_and_parenthesizations():
    mass = 1.7
    rows = [scalar(4., 1.), _restricted_example(), isotropic([3., 1.], [1., 5.], .4, .8),
            scalar(-4., 2.)]
    states = [tangent_pairings([row], mass) for row in rows]
    one_pass = tangent_pairings(rows, mass)
    for combined in [states[0].merge(states[1].merge(states[2])).merge(states[3]),
                     states[3].merge(states[0]).merge(states[1].merge(states[2])),
                     tangent_pairings(reversed(rows), mass)]:
        np.testing.assert_allclose(energies(combined), energies(one_pass), rtol=5e-14)
        assert combined.correction == pytest.approx(one_pass.correction)
        assert combined.degrees == (1, 3, 1)
    left, right = states[0].merge(states[1]), states[2].merge(states[3])
    difference = (left.correction - right.correction)**2
    increment = left.weight * right.weight / (left.weight + right.weight) * difference
    assert one_pass.between == pytest.approx(left.between + right.between + increment)
    assert one_pass.common == pytest.approx(gaussian_score_test(rows, mass)["statistic"])
    with pytest.raises(ValueError, match="same candidate"):
        left.merge(tangent_pairings([], 2.))


def test_zero_padding_and_nonzero_zero_tangent_conventions():
    mass = 2.
    empty = tangent_pairings([], mass)
    assert empty.degrees == (0, 0, 0) and empty.correction is None
    assert tangent_pairings([scalar(0., 0.)], mass) == empty
    inactive = scalar(1., -2.)
    state = tangent_pairings([inactive], mass)
    assert state.degrees == (1, 0, 0) and state.correction is None
    np.testing.assert_allclose(energies(state), [5., 0., 0.])
    active = scalar(2., 3.)
    combined = tangent_pairings([active, inactive], mass)
    assert combined == tangent_pairings([scalar(0., 0.), active, inactive], mass)
    assert combined.degrees == (1, 0, 1)
    np.testing.assert_allclose(energies(combined), [5., 0., gaussian_score_test([active], mass)["statistic"]])


def test_spatial_rotation_unit_and_reciprocal_chart_changes():
    reading, mass = _restricted_example(), 2.3
    rows = [reading, scalar(-1., 2.)]
    original = tangent_pairings(rows, mass)
    rotation = np.array([[.6, -.8], [.8, .6]])
    transform = np.zeros((4, 4))
    transform[:2, :2] = 3. * rotation
    transform[2:, 2:] = rotation
    data = transform @ np.r_[reading.force, reading.acceleration]
    altered = Reading(data[:2], data[2:], transform @ reading.cov @ transform.T,
                      rotation @ reading.subspace)
    scaled = tangent_pairings([altered, isotropic([-3.], [2.], 3., 1.)], 3. * mass)
    derivative = 3 * (1 + mass**2) / (1 + (3 * mass)**2)
    np.testing.assert_allclose(energies(scaled), energies(original), atol=1e-12)
    assert scaled.correction == pytest.approx(derivative * original.correction)
    swapped = [Reading(rd.acceleration, rd.force, np.roll(np.roll(rd.cov, rd.dimension, 0),
                       rd.dimension, 1), rd.subspace) for rd in rows]
    other = tangent_pairings(swapped, 1. / mass)
    np.testing.assert_allclose(energies(other), energies(original), atol=1e-12)
    assert other.correction == pytest.approx(-original.correction)


@pytest.mark.parametrize("within,between,common", [(-.8, .2, .7), (.2, -1.9, -.5), (3., .4, -1.1)])
def test_explicit_orthogonal_normal_coordinates_at_fixed_graph_projections(within, between, common):
    # Holding the graph projections at (3,0) in 2D and 4 in 1D leaves
    # three independent Gaussian residual coordinates at the true candidate.
    residuals = np.array([3 * common + 4 * between, 5 * within, 4 * common - 3 * between]) / 5
    along = np.array([3., 0., 4.])
    force, accel = (along + residuals) / math.sqrt(2), (along - residuals) / math.sqrt(2)
    state = tangent_pairings([isotropic(force[:2], accel[:2], 1., 1.), scalar(force[2], accel[2])], 1.)
    assert state.degrees == (1, 1, 1)
    np.testing.assert_allclose(energies(state), np.square([within, between, common]), atol=2e-14)
    assert state.score / math.sqrt(state.weight) == pytest.approx(common)


@pytest.mark.parametrize("mass", [.3, 1., 4.])
def test_exterior_tensor_cost_agrees_with_candidate_tangent_partition(mass):
    # Fixed channel ratio s=2, different channel scales and spatial dimensions.
    pairs = [(np.array([4., 1.]), np.array([1., 4.])),
             (np.array([-.7, 2.3, .4]), np.array([1., -.2, .3])),
             (np.array([1.8]), np.array([.4]))]
    within = sum((f @ f) * (a @ a) - (f @ a)**2 for f, a in pairs)
    cross = sum(np.linalg.norm(np.outer(f, b) - np.outer(a, g))**2
                for index, (f, a) in enumerate(pairs) for g, b in pairs[index + 1:])
    data = np.concatenate([np.column_stack((f, a)) for f, a in pairs], axis=0)
    gram = data.T @ data
    assert within + cross == pytest.approx(np.linalg.det(gram))
    values = np.linalg.eigvalsh(gram)
    assert (within + cross) / values[-1] == pytest.approx(values[0])
    rows = [isotropic(2 * sd * f, sd * a, 2 * sd, sd) for sd, (f, a) in zip([.3, 1., .7], pairs)]
    state = tangent_pairings(rows, mass)
    # Convert from the helper's derivative at atan(m) to the noise angle.
    derivative = 2 * (1 + (mass / 2)**2) / (1 + mass**2)
    angle_weight = derivative**2 * state.weight
    assert angle_weight * (state.within + state.between) == pytest.approx(within + cross)


def test_widening_on_an_independent_between_reading_contrast_strictly_increases_coverage():
    p = .9
    critical = stats.norm.ppf((1 + p) / 2)
    # Integrate over the unit circle directions of two independent N(0,1)s
    # (between coordinate, common score); the squared radius is chi-square(2).
    def angle_probability(angle):
        c, s = abs(math.cos(angle)), abs(math.sin(angle))
        switch = 1 / c
        before = min(switch, critical / s)
        after = 1.5 * critical / s
        return -math.expm1(-before**2 / 2) + max(0., math.exp(-switch**2 / 2) - math.exp(-after**2 / 2))
    breaks = [math.atan(critical), math.atan(1.5 * critical)]
    actual = 2 / math.pi * integrate.quad(angle_probability, 0., math.pi / 2, points=breaks)[0]
    expected = p + 2 * stats.norm.sf(1.) * (2 * stats.norm.cdf(1.5 * critical) - 1 - p)
    assert actual == pytest.approx(expected, abs=1e-11)
    assert actual > p

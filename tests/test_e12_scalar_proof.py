"""Finite consequences of the composition obstruction and completed scalar witness."""
import math

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

from explore.e12_scalar_completion import completed_log_evidence
from explore.e12_scalar_witness import pooled_witness
from met.full import _h_terms, isotropic, reading_state
from tools.leftover_at_series_best_mass import law as historical_law


ALPHA = 1 / 8
QUADRANT = math.pi / 2


def gaussian_log_evidence(rows, angles):
    normals = np.stack((np.cos(angles), -np.sin(angles)))
    return -0.5 * np.sum((np.asarray(rows) @ normals) ** 2, axis=0)


def cumulative_mass(angles, log_factor, prior_density):
    density = prior_density * np.exp(log_factor - np.max(log_factor))
    result = cumulative_trapezoid(density, angles, initial=0)
    return result / result[-1]


def completed_quantiles(forces, sigma_f, sigma_a, c, alpha):
    rows = np.column_stack((np.asarray(forces) / sigma_f, np.full(len(forces), 1 / sigma_a)))
    angles = np.linspace(0, QUADRANT, 150001)
    ratio = sigma_f / sigma_a
    # The half-Cauchy prior on m, pulled back along m = ratio * tan(angle).
    prior = 1 / ((c * np.cos(angles)) ** 2 + (ratio * np.sin(angles)) ** 2)
    cdf = cumulative_mass(angles, completed_log_evidence(rows, angles, alpha=alpha), prior)
    return ratio * np.tan(np.interp([.025, .5, .975], cdf, angles))


@pytest.mark.parametrize("alpha,far_median,far_half_width", [
    (1 / 8, 5.492869, .389155), (1 / 16, 5.486878, .519363),
])
def test_completed_law_on_both_o14_pairs_with_its_mass_prior(alpha, far_median, far_half_width):
    for forces, median, half_width in [
        ([5.5, 5.6], 5.549021, .139444), ([5., 6.], far_median, far_half_width),
    ]:
        lo, med, hi = completed_quantiles(forces, .1, .002, c=5, alpha=alpha)
        assert (med, (hi - lo) / 2) == pytest.approx((median, half_width), abs=3e-6)


def test_o14_three_readings_distinguish_the_original_rules_and_the_completion():
    forces = [5., 5.05, 6.]
    for kind, median, half_width in [
        ("gauss", 5.350, .114), ("pooled", 5.337, .521), ("per reading", 5.126, .422),
    ]:
        med, lo, hi = historical_law(forces, kind=kind)
        assert (med, (hi - lo) / 2) == pytest.approx((median, half_width), abs=.00051)
    lo, med, hi = completed_quantiles(forces, .1, .002, c=5, alpha=ALPHA)
    assert (med, (hi - lo) / 2) == pytest.approx((5.346192, .289352), abs=3e-6)


@pytest.mark.parametrize("sigma,per_reading_report", [
    (.1, (1.038789, .802607, 1.347937)), (.01, (1.000419, .973090, 1.028517)),
])
def test_a_reading_exactly_on_the_series_fit(sigma, per_reading_report):
    forces = [.2, 1., 1.4]
    # The raw least-residual line is the middle row for both declared errors.
    for kind, expected in [
        ("per reading", per_reading_report), ("pooled", (.999970, .497171, 2.011259)),
    ]:
        assert historical_law(forces, sigma, sigma, c=1, kind=kind) == pytest.approx(expected, abs=2e-6)
    assert completed_quantiles(forces, sigma, sigma, c=1, alpha=ALPHA) == pytest.approx(
        (.549513, 1., 1.819794), abs=2e-6)


@pytest.mark.parametrize("scale", [1., 10., 100.])
def test_exact_product_duplicate_cancellation_and_oscillation_bound(scale):
    x, y = scale * np.array([[1., 1.], [4., 1.]])
    angles = np.linspace(0., QUADRANT, 501)
    one = gaussian_log_evidence([x], angles)
    mix = gaussian_log_evidence([x, y], angles)
    duplicate = gaussian_log_evidence([x, x, y], angles)
    np.testing.assert_allclose(duplicate - mix, one, atol=2e-11)
    assert np.ptp(one) <= np.ptp(duplicate) + np.ptp(mix)
    assert np.ptp(one) == pytest.approx(0.5 * scale ** 2)


def test_completed_joint_law_escapes_product_cancellation_on_the_same_lists():
    x, y = np.array([[1., 1.], [4., 1.]])
    angles = np.linspace(0., QUADRANT, 501)
    for rows in [[x, y], [x, x, y]]:
        first = completed_log_evidence(20 * np.array(rows), angles, alpha=ALPHA)
        second = completed_log_evidence(200 * np.array(rows), angles, alpha=ALPHA)
        np.testing.assert_allclose(first, second, atol=5e-11)
        assert np.ptp(first) > 0.1
    mixed = completed_log_evidence(200 * np.array([x, y]), angles, alpha=ALPHA)
    duplicated = completed_log_evidence(200 * np.array([x, x, y]), angles, alpha=ALPHA)
    singleton = completed_log_evidence(200 * np.array([x]), angles, alpha=ALPHA)
    assert np.ptp(singleton) > np.ptp(mixed) + np.ptp(duplicated)


def test_bounded_contrast_controls_posterior_interval_mass():
    angles = np.linspace(0., QUADRANT, 4097)
    factor = completed_log_evidence([[4., 0.]], angles, alpha=ALPHA)
    cdf = cumulative_mass(angles, factor, 1.)
    # A short interval on the favoured side gains mass; the bound is nontrivial.
    start, end = 3700, 3980
    prior_mass = (angles[end] - angles[start]) / QUADRANT
    bound = math.exp(np.ptp(factor)) * prior_mass
    assert prior_mass < cdf[end] - cdf[start] <= bound < 1


@pytest.mark.parametrize("row,endpoint", [([1., 0.], QUADRANT), ([0., 1.], 0.)])
def test_axis_gaussian_continuity_forces_concentration_and_the_cap_limits_it(row, endpoint):
    angles = np.linspace(0., QUADRANT, 4097)
    row = np.array(row)
    for sigma in [.1, .01]:
        gaussian = gaussian_log_evidence([row / sigma], angles)
        nearby = gaussian_log_evidence([(row + (1 - row) * 1e-11) / sigma], angles)
        np.testing.assert_allclose(nearby, gaussian, atol=1e-7)
        completed = completed_log_evidence([row / sigma], angles, alpha=ALPHA)
        np.testing.assert_allclose(completed, -.5 * np.sin(angles - endpoint) ** 2, atol=1e-14)
        assert np.ptp(gaussian) == pytest.approx(.5 / sigma ** 2)
        assert np.ptp(completed) == pytest.approx(.5)
    # Probability of the endpoint's fixed 0.05-radian neighbourhood.
    def end_mass(factor):
        cdf = cumulative_mass(angles, factor, 1.)
        return np.interp(.05, angles, cdf) if endpoint == 0 else 1 - np.interp(endpoint - .05, angles, cdf)

    assert end_mass(gaussian) > .99999
    assert end_mass(completed) < .1


def test_equal_eigenvalues_have_zero_gaussian_contrast_and_a_nonconstant_completed_limit():
    angles = np.array([0., QUADRANT / 2, QUADRANT])
    for scale in [16., 1024.]:
        rows = scale * np.eye(2)
        assert np.ptp(gaussian_log_evidence(rows, angles)) < 5e-10
        completed = completed_log_evidence(rows, angles, alpha=ALPHA)
        loading = scale ** 2 / (1 + scale ** 2)
        expected = loading / (1 + 2 * loading) * (2 * math.log1p(loading / 2) - math.log1p(loading))
        assert completed[0] - completed[1] == pytest.approx(expected, abs=2e-14)
        assert completed[2] == pytest.approx(completed[0])
    limit = (2 * math.log(1.5) - math.log(2.)) / 3
    assert completed[0] - completed[1] == pytest.approx(limit, abs=2e-7)
    assert limit > .03


def test_rotating_a_balanced_pair_changes_its_reading_information():
    angles = np.linspace(0., QUADRANT, 501)
    axis = np.array([[5., 0.], [0., 5.]])
    oblique = np.array([[3., 4.], [-4., 3.]])
    np.testing.assert_array_equal(oblique.T @ oblique, axis.T @ axis)
    left, right = [completed_log_evidence(rows, angles, alpha=ALPHA) for rows in (axis, oblique)]
    assert np.ptp(left - right) > .04


def test_negative_slope_shrinking_ray_has_a_spread_informative_limit():
    angles = np.linspace(0., QUADRANT, 1001)
    row = np.array([[-3., 2.]])
    first, second = [completed_log_evidence(scale * row, angles, alpha=ALPHA) for scale in [30., 300.]]
    np.testing.assert_allclose(first, second, atol=1e-13)
    assert 0.1 < np.ptp(first) < 10.
    assert first[-1] > first[0] > np.min(first)


@pytest.mark.parametrize("rows", [
    [[3., 4.], [-4., 3.]], [[4., 1.], [1., -3.], [0., 1.]], [[-3., 2.]], [[0., 3.]],
])
def test_zero_padding_row_symmetries_and_channel_swap(rows):
    rows = np.asarray(rows)
    angles = np.linspace(0., QUADRANT, 79)
    original = completed_log_evidence(rows, angles, alpha=ALPHA)
    changed = rows[::-1].copy()
    changed[::2] *= -1
    padded = np.vstack(([0., 0.], changed, [0., 0.]))
    np.testing.assert_allclose(completed_log_evidence(padded, angles, alpha=ALPHA), original, atol=1e-13)
    swapped = completed_log_evidence(padded[:, ::-1], QUADRANT - angles, alpha=ALPHA)
    np.testing.assert_allclose(swapped, original, atol=1e-13)


@pytest.mark.parametrize("rows", [[], [[0., 0.], [0., 0.]]])
def test_no_information_preserves_every_prior(rows):
    angles = np.linspace(0., QUADRANT, 101)
    np.testing.assert_array_equal(completed_log_evidence(rows, angles, alpha=ALPHA), 0.)


def test_identifying_interior_data_retain_the_local_pooled_kernel():
    true_angle = .48
    u = np.array([math.sin(true_angle), math.cos(true_angle)])
    normal = np.array([math.cos(true_angle), -math.sin(true_angle)])
    # The two equally displaced measurements retain a nonzero residual.
    rows = 5000 * u + np.array([3., -3.])[:, None] * normal
    fit = pooled_witness(rows, alpha=ALPHA)
    assert fit.angle == pytest.approx(true_angle) and fit.positive_residual == pytest.approx(18.)
    x = np.array([-3., -1., 0., .5, 2.])
    angles = fit.angle + x / math.sqrt(fit.precision)
    completed = completed_log_evidence(rows, angles, alpha=ALPHA)
    pooled = -.5 * fit.precision * np.sin(angles - fit.angle) ** 2
    np.testing.assert_array_equal(completed, pooled)
    np.testing.assert_allclose(completed, -.5 * x ** 2, atol=3e-7)


@pytest.fixture
def colliding_lists():
    compatible = math.sqrt(13 / 15) * np.ones((3, 2))
    conflicting = np.column_stack(([.2, 1., 1.4], np.ones(3)))
    return compatible, conflicting


def test_two_and_three_row_collisions_add_only_isotropic_scatter(colliding_lists):
    two_row_pair = (np.full((2, 2), 2), np.array([[4, 1], [1, 4]]))
    for (compatible, conflicting), increment in [(colliding_lists, .4), (two_row_pair, 9)]:
        aligned_gram = compatible.T @ compatible
        conflict_gram = conflicting.T @ conflicting
        np.testing.assert_allclose(conflict_gram - aligned_gram, increment * np.eye(2), atol=1e-15)
        assert np.linalg.matrix_rank(aligned_gram) == 1
        assert np.linalg.matrix_rank(conflict_gram) == 2


def test_one_declared_vector_reading_has_a_t36_curve_collision_and_distinct_t32_minima():
    angles = np.linspace(0., QUADRANT, 257)
    for sigma in [1., .1, .01]:
        compatible = isotropic([2., 2.], [2., 2.], sigma, sigma)
        conflicting = isotropic([4., 1.], [1., 4.], sigma, sigma)
        np.testing.assert_array_equal(compatible.cov, conflicting.cov)
        residuals = []
        for reading, expected_min in [(compatible, 0.), (conflicting, 9 / sigma ** 2)]:
            assert reading.k == reading.dimension == 2
            pre = reading.pre()
            residuals.append(pre["energy"] - _h_terms(pre, 1., angles))
            state = reading_state(reading, 1.)
            assert state["Q"] == pytest.approx(expected_min, abs=2e-9)
            assert state["Q_plus"] == pytest.approx(expected_min, abs=2e-9)
        left, right = residuals
        np.testing.assert_allclose(right - left, 9 / sigma ** 2, atol=2e-9)
        # Gaussian log-evidence relative to the common m=1 candidate.
        centre = len(angles) // 2
        np.testing.assert_allclose(-.5 * (left - left[centre]), -.5 * (right - right[centre]), atol=2e-9)


def test_reciprocal_pairings_split_the_vector_reading_cost_in_each_spatial_frame():
    rotation = np.array([[.6, -.8], [.8, .6]])
    angles = np.array([.15, .5, 1.1])
    masses = np.tan(angles)
    sigma_f, sigma_a = .7, .4
    for transform in (np.eye(2), rotation, rotation @ np.diag([-1., 1.])):
        for observed_force, observed_accel, expected in (
            ([2., 2.], [2., 2.], 1.),
            ([4., 1.], [1., 4.], 8 / 17),
            ([-4., -1.], [1., 4.], -8 / 17),
        ):
            force, accel = transform @ observed_force, transform @ observed_accel
            forward = (accel @ force) / (accel @ accel)
            reverse = (force @ accel) / (force @ force)
            assert forward == pytest.approx(expected)
            assert reverse == pytest.approx(expected)
            transverse = force - forward * accel
            assert transverse @ accel == pytest.approx(0, abs=1e-13)
            cost = ((accel @ accel) * (masses - forward) ** 2 + transverse @ transverse)
            cost /= sigma_f ** 2 + masses ** 2 * sigma_a ** 2
            pre = isotropic(force, accel, sigma_f, sigma_a).pre()
            np.testing.assert_allclose(cost, pre["energy"] - _h_terms(pre, 1., angles), atol=1e-12)


def test_normalized_aggregate_curves_coincide_at_each_noise_scale(colliding_lists):
    compatible, conflicting = colliding_lists
    angles = np.linspace(0., QUADRANT, 2001)
    for sigma in [1., .1, .01]:
        left, right = [gaussian_log_evidence(rows / sigma, angles) for rows in (compatible, conflicting)]
        np.testing.assert_allclose(right - left, -.2 / sigma ** 2, atol=2e-11)
        np.testing.assert_allclose(left - left.max(), right - right.max(), atol=2e-11)
        # Even after adding a nonuniform angular prior, the posterior CDFs agree.
        prior = 1 + np.sin(angles)
        np.testing.assert_allclose(cumulative_mass(angles, left, prior),
                                   cumulative_mass(angles, right, prior), atol=2e-13)


def test_completion_separates_the_colliding_lists_as_the_errors_shrink(colliding_lists):
    compatible, conflicting = colliding_lists
    angles = np.linspace(0., QUADRANT, 20001)
    near = np.array([math.pi / 4 - .015, math.pi / 4 + .015])
    masses = []
    for sigma in [.02, .01, .005]:
        pair = []
        for rows in (compatible, conflicting):
            factor = completed_log_evidence(rows / sigma, angles, alpha=ALPHA)
            cdf = cumulative_mass(angles, factor, 1.)
            lo, hi = np.interp(near, angles, cdf)
            pair.append(hi - lo)
        masses.append(pair)
    compatible_mass, conflicting_mass = np.array(masses).T
    assert compatible_mass[0] < .95 and compatible_mass[-1] > .999999
    assert np.all(np.diff(compatible_mass) > 0.)
    assert .07 < conflicting_mass[0] < .10
    np.testing.assert_allclose(conflicting_mass, conflicting_mass[0], atol=1e-12)


def test_a_single_normalized_curve_recovers_its_rank_one_gram_including_sign_and_zero():
    angles = np.array([0., math.pi / 4, QUADRANT])
    for row in [[2., 3.], [-2., 3.], [2., -3.], [-2., -3.], [0., 3.], [2., 0.], [0., 0.]]:
        for offset in [0., 19.]:
            q0, qmid, qend = -2 * gaussian_log_evidence([row], angles) + offset
            difference = q0 - qend
            twice_cross = q0 + qend - 2 * qmid
            trace = math.hypot(difference, twice_cross)
            gram = .5 * np.array([[trace + difference, twice_cross], [twice_cross, trace - difference]])
            np.testing.assert_allclose(gram, np.outer(row, row), atol=1e-14)


def test_orthogonal_row_mixing_preserves_complex_square_sum_in_each_realization():
    v = np.array([1., 2., -2.]) / 3
    transform = np.eye(3) - 2 * np.outer(v, v)
    np.testing.assert_allclose(transform.T @ transform, np.eye(3), atol=1e-15)
    means = np.array([1. + 2.j, 1. - .4j, 3. + .8j])
    rng = np.random.default_rng(1204)
    noise = rng.normal(size=(3, 20)) + 1j * rng.normal(size=(3, 20))
    for sigma in [1., .01]:
        observations = means[:, None] + sigma * noise
        mixed = transform @ observations
        statistic = np.sum(observations ** 2, axis=0)
        np.testing.assert_allclose(np.sum(mixed ** 2, axis=0), statistic, atol=1e-13)
        assert np.max(np.abs(statistic.imag)) > 1
        assert np.max(np.abs(mixed ** 2 - observations ** 2)) > 1

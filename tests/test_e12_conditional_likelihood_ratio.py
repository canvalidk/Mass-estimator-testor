"""Exact and independent numerical checks for the conditional graph likelihood-ratio test."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erfc, ndtri
from scipy.stats import chi2

from explore.e12_conditional_likelihood_ratio import conditional_lr_pvalue, gaussian_rank_weight


def test_far_tail_preserves_the_known_normal_lower_bound():
    # A beta-quantile integral with absolute tolerances can lose this tail.
    angle = np.arcsin(np.sqrt(500./1e6))
    value = conditional_lr_pvalue(np.diag([0., 1e6]), angle, 5)
    floor = erfc(np.sqrt(500./2))
    assert floor <= value <= 2*floor


def direction(angle):
    return np.array([np.sin(angle), np.cos(angle)])


def exact_one_dimensional_tail(largest, loss, n):
    # Integrate over the single normal coordinate instead of the chi-square
    # variable used in the reported formula.
    if loss == 0:
        return 1.0
    normal_tail = erfc(np.sqrt(loss / 2))
    if n == 1:
        return normal_tail
    limit = np.sqrt(loss)
    def integrand(z):
        threshold = largest * (1 - z*z/loss)
        return np.sqrt(2/np.pi) * np.exp(-z*z/2) * chi2.sf(threshold, n-1)
    return normal_tail + quad(integrand, 0., limit, epsabs=1e-12)[0]


@pytest.mark.parametrize("n", [2, 3, 16])
def test_agrees_with_independent_normal_coordinate_quadrature(n):
    scatter = np.array([[4.5, .8], [.8, 1.1]])
    largest = np.linalg.eigvalsh(scatter)[-1]
    angles = np.array([-.65, .1, .72, 1.35])
    expected = [exact_one_dimensional_tail(largest, largest-direction(t)@scatter@direction(t), n)
                for t in angles]
    np.testing.assert_allclose(conditional_lr_pvalue(scatter, angles, n), expected,
                               atol=2e-9, rtol=2e-8)


@pytest.mark.parametrize("profile", [
    np.array([2., 0., 0., 0., 0.]),
    np.array([1., -1., 1., -1., 0.]),
])
def test_conditional_tail_matches_full_gaussian_scatter_recomputation(profile):
    # Fix observed longitudinal coordinates; redraw only the independent
    # transverse coordinates, then refit a whole line at every draw.
    n = len(profile)
    longitudinal = profile @ profile
    observed_noise = np.array([1.3, -.7, 1.1, .6, -.4])
    observed = np.column_stack([observed_noise, profile])
    scatter = observed.T @ observed
    observed_largest = np.linalg.eigvalsh(scatter)[-1]
    rng = np.random.default_rng(203)
    errors = rng.standard_normal((90000, n))
    baseline = np.broadcast_to(profile, errors.shape)
    full = np.stack([errors, baseline], axis=2)
    draws = np.linalg.eigvalsh(full.transpose(0, 2, 1) @ full)[:, -1]
    estimate = np.mean(draws - longitudinal >= observed_largest - longitudinal)
    assert abs(conditional_lr_pvalue(scatter, 0., n)-estimate) < .006


def test_one_row_is_exact_squared_transverse_normal_test():
    row = np.array([2., 3.])
    angles = np.array([-.7, .1, .5, 1.3])
    expected = erfc(np.abs(row[0]*np.cos(angles)-row[1]*np.sin(angles)) / np.sqrt(2))
    np.testing.assert_allclose(conditional_lr_pvalue(np.outer(row, row), angles, 1), expected)
    with pytest.raises(ValueError):
        conditional_lr_pvalue(np.eye(2), 0., 1)


def test_level_sets_have_one_best_line_and_respect_metric_frame_changes():
    scatter = np.diag([1., 12.])
    angles = np.linspace(0., np.pi/2, 20)
    values = conditional_lr_pvalue(scatter, angles, 5)
    assert values[0] == pytest.approx(1)
    assert np.all(np.diff(values) < 0)
    shift = .63
    rotation = np.array([[np.cos(shift), np.sin(shift)], [-np.sin(shift), np.cos(shift)]])
    np.testing.assert_allclose(conditional_lr_pvalue(rotation@scatter@rotation.T, angles+shift, 5),
                               values, atol=2e-9)
    np.testing.assert_allclose(conditional_lr_pvalue(np.zeros((2, 2)), angles, 5), 1.)
    np.testing.assert_allclose(conditional_lr_pvalue(np.eye(2), angles, 5), 1., atol=5e-8)


@pytest.mark.parametrize("n", [1, 3, 6])
def test_orthogonal_positive_mass_has_central_chi_square_tail_instead_of_perfect_fit(n):
    precision = np.array([[2., -3.], [-3., 5.]])
    whitening = np.linalg.cholesky(precision).T
    u1, u2 = [whitening @ np.array([mass, 1.]) for mass in [1., 2.]]
    u1, u2 = u1 / np.linalg.norm(u1), u2 / np.linalg.norm(u2)
    assert abs(u1 @ u2) < 1e-14
    angles = [np.arctan2(u[0], u[1]) for u in [u1, u2]]
    length_squared = 25.
    scatter = length_squared*np.outer(u1, u1)
    values = conditional_lr_pvalue(scatter, angles, n)
    assert values[0] == pytest.approx(1., abs=1e-7)
    assert values[1] == pytest.approx(chi2.sf(length_squared, n), rel=1e-7, abs=1e-13)


def test_weak_floor_includes_every_angle_and_is_sharp_at_rank_one():
    n, largest = 8, 5.
    baseline = chi2.sf(largest, n)
    angles = np.linspace(0., np.pi, 31)
    values = conditional_lr_pvalue(np.diag([0., largest]), angles, n)
    assert values.min() == pytest.approx(baseline, abs=1e-9)
    assert np.all(values > .05)
    assert np.all(conditional_lr_pvalue(np.diag([1., largest]), angles, n) >= baseline-1e-10)


@pytest.mark.parametrize("r", [1., 2.])
def test_local_width_on_regular_scatter_tracks_actual_sampling_variance(r):
    n, level = 1600, .15
    # Population scatter at this design gives L=n(1+r^2), gap=n*r^2;
    # the estimator's angular standard deviation is sqrt(1+r^2)/(sqrt(n)*r^2).
    angle = ndtri(1-level/2)*np.sqrt(1+r*r)/(np.sqrt(n)*r*r)
    scatter = n*np.diag([1., 1+r*r])
    assert conditional_lr_pvalue(scatter, angle, n) == pytest.approx(level, abs=.005)


def test_two_sided_gaussian_rank_matching_has_strictly_decreasing_discount():
    scatter = np.diag([1., 6.])
    angles = np.array([.09, .3, .6, .95, 1.45])
    losses = 5*np.sin(angles)**2
    ranks = conditional_lr_pvalue(scatter, angles, 3)
    discount = ndtri(1-ranks/2)**2 / losses
    assert np.all((0 < discount) & (discount < 1))
    assert np.all(np.diff(discount) < 0)


def test_an_appended_zero_changes_nonpeak_weight_and_preserves_the_peak():
    block = np.array([[4., 1.], [1., 4.]])
    with_zero = np.vstack([block, [0., 0.]])
    angles = np.array([np.pi/4, np.arctan(2.)])
    first = conditional_lr_pvalue(block.T@block, angles, len(block))
    second = conditional_lr_pvalue(with_zero.T@with_zero, angles, len(with_zero))
    assert first[0] == pytest.approx(1., abs=1e-7)
    assert second[0] == pytest.approx(1., abs=1e-7)
    assert 0 < first[1] < second[1] < 1


def test_gaussian_weight_matches_single_row_likelihood_and_softens_on_multiple_rows():
    row = np.array([2., 3.])
    angles = np.array([-.25, .2, .6, 1.2])
    transverse = row[0]*np.cos(angles)-row[1]*np.sin(angles)
    np.testing.assert_allclose(gaussian_rank_weight(np.outer(row, row), angles, 1),
                               np.exp(-.5*transverse**2), rtol=1e-10)
    scatter = np.diag([2., 7.])
    values = gaussian_rank_weight(scatter, angles, 5)
    pvalues = conditional_lr_pvalue(scatter, angles, 5)
    np.testing.assert_allclose(values, np.exp(-.5*ndtri(pvalues/2)**2))
    raw = np.exp(-.5*(7-2)*np.sin(angles)**2)
    assert np.all((values > raw) & (values < 1))
    assert gaussian_rank_weight(scatter, 0., 5) == pytest.approx(1)

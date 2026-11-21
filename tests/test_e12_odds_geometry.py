"""Finite checks for pairwise minimax preferences and the scope of odds multiplication."""

from fractions import Fraction
from itertools import combinations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import hyp0f1, i0, ndtr


def unit(mass):
    return np.array([mass, 1.0]) / np.hypot(mass, 1.0)


def projection_energy(rows, mass):
    return np.linalg.norm(rows @ unit(mass))**2


def sphere_factor(dimension, t):
    return hyp0f1(dimension / 2, t**2 / 4)


@pytest.mark.parametrize("length", [.2, 1.7, 3.])
def test_sphere_mixtures_match_direct_gaussian_density_integration(length):
    y = np.array([[.3, -.7], [1.4, 2.1]])
    u = unit(1.7)
    def integrand(angle):
        b = length * np.array([np.sin(angle), np.cos(angle)])
        return np.exp(-.5 * np.sum((y - np.outer(b, u))**2))
    integrated = quad(integrand, 0., 2*np.pi)[0] / (2*np.pi)
    radius = length*np.linalg.norm(y @ u)
    relative_factor = np.exp(-.5 * length**2) * sphere_factor(2, radius)
    assert integrated == pytest.approx(np.exp(-.5*np.sum(y**2)) * relative_factor)
    assert sphere_factor(1, radius) == pytest.approx(np.cosh(radius))
    assert sphere_factor(2, radius) == pytest.approx(i0(radius))
    assert sphere_factor(3, radius) == pytest.approx(np.sinh(radius) / radius)


def test_sphere_ratio_has_strength_independent_order_but_crosses_cost_threshold():
    y = np.array([[1., 3.], [2., .1], [.1, 1.]])
    t1, t2 = (np.linalg.norm(y @ unit(m)) for m in [.5, 2.])
    assert t1 > t2
    ratios = [sphere_factor(3, delta*t1)/sphere_factor(3, delta*t2)
              for delta in [.01, .1, 1., 10., 40.]]
    assert all(r > 1 for r in ratios)
    assert ratios[0] < 1.0001 and ratios[-1] > 2
    # Every positive likelihood decision at threshold one picks the lower covector cost.
    assert np.sum(y**2) - t1*t1 < np.sum(y**2) - t2*t2


def test_class_and_profile_symmetries_preserve_conditional_errors():
    rng = np.random.default_rng(18)
    profile = np.array([.2, .9, -.3, 1.5])
    rotation, _ = np.linalg.qr(rng.normal(size=(4, 4)))
    u, v = unit(.4), unit(2.)
    bisector = (u + v)/np.linalg.norm(u + v)
    reflection = 2*np.outer(bisector, bisector) - np.eye(2)
    np.testing.assert_allclose(reflection @ u, v)
    y = rng.normal(size=(30, 4, 2)) + np.outer(profile, u)
    def choose_first(rows):
        return np.linalg.norm(rows @ u, axis=-1) > np.linalg.norm(rows @ v, axis=-1)
    np.testing.assert_array_equal(choose_first(rotation @ y), choose_first(y))
    np.testing.assert_array_equal(choose_first(y @ reflection), ~choose_first(y))


def test_single_push_error_matches_exact_two_independent_bisector_signs():
    u, v = unit(.5), unit(2.)
    half_angle = np.arccos(u @ v)/2
    b0, b1 = (u+v)/np.linalg.norm(u+v), (u-v)/np.linalg.norm(u-v)
    for r in [-1.3, 1.3]:
        p0, p1 = ndtr(r * np.cos(half_angle)), ndtr(r * np.sin(half_angle))
        exact_error = (1-p0)*p1 + p0*(1-p1)
        # Distribution projected onto the two independent bisector coordinates.
        p, q = ndtr(r*u @ b0), ndtr(r*u @ b1)
        assert exact_error == pytest.approx(p*(1-q) + q*(1-p))
        assert 0 < exact_error < .5


def test_positive_readings_reverse_existing_singleton_product_preference():
    def scores(f, a):
        diff = Fraction((f+a)**2, 2) - Fraction((2*f+a)**2, 5)
        cap = min(1, Fraction(1+min(f*f, a*a)**2, f*f+a*a))
        return diff, cap*diff/2
    first, second = scores(1, 2), scores(3, 1)
    assert first[0] + second[0] == Fraction(-1, 2)
    assert first[1] + second[1] == Fraction(2, 25)


def test_two_row_ordinal_score_can_fail_a_three_row_comparison():
    u, v = unit(1), unit(2)
    values, axes = np.linalg.eigh(np.outer(u, u)-np.outer(v, v))
    def row_at(d):
        index = int(d > 0)
        return np.sqrt(d/values[index])*axes[:, index]
    rows = np.array([row_at(d) for d in [2., 2., -6.]])
    differences = (rows @ u)**2 - (rows @ v)**2
    np.testing.assert_allclose(differences, [2, 2, -6])
    for i, j in combinations(range(3), 2):
        pair = differences[[i, j]]
        assert np.sign(sum(pair)) == np.sign(np.arctan(pair).sum())
    assert sum(differences) < 0 < np.arctan(differences).sum()


def test_telescope_on_three_distinct_graphs_fixes_the_whole_slope_vector():
    projectors = [np.outer(unit(m), unit(m)) for m in [.35, 1.2, 3.3]]
    changes = np.column_stack([(projectors[i]-projectors[j]).ravel()
                               for i, j in [(0, 1), (1, 2), (2, 0)]])
    assert np.linalg.matrix_rank(changes) == 2
    np.testing.assert_allclose(changes @ np.ones(3), 0., atol=1e-14)


def test_repeated_valid_pushes_require_different_scales_at_the_same_instrument():
    cdf_width = lambda r, c: 1 / np.sqrt(2*c*r*r)
    sampling_width = lambda r: np.sqrt(1+r*r)/(r*r)
    for r, c in [(1., .25), (2., .4)]:
        assert cdf_width(r, c) == pytest.approx(sampling_width(r))
    assert cdf_width(1., .4) < sampling_width(1.)
    assert cdf_width(2., .25) > sampling_width(2.)


def test_conflict_witnesses_keep_all_decisions_with_different_finite_odds():
    rows = np.array([[4., 1.], [1., 4.]])
    gram = rows.T @ rows
    small, large = np.linalg.eigvalsh(gram)
    h = np.exp(-np.sum(rows**2, axis=1)/2)
    count = np.sum((1-h)*(1+3*h))
    def log_factor_difference(alpha):
        discount = min(1., (count/small)*(large/small)**alpha)
        span = min(discount*(large-small)**2/large, 1+min(np.diag(gram))**2)
        coeff = span/(2*(large-small))
        return coeff*(projection_energy(rows, 1.)-projection_energy(rows, 2.))
    assert 0 < log_factor_difference(1.) < log_factor_difference(2.)

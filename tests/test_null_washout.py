"""Independent finite identities and controls for E10's washout continuation."""
import math
from itertools import permutations, product

import numpy as np
import pytest
from scipy.special import roots_hermitenorm

from met.null_washout import (fixed_signal_log_transfer_bound,
                             gated_precision, homogeneous_log_transfer_bound,
                             mixture_log_second_moment, ordered_scatter,
                             symmetric_mixture_log_product_bound,
                             washout_fit)


@pytest.mark.parametrize('radius', [.01, .3, .8, 1.2])
def test_sign_mixture_second_moment_by_independent_gaussian_quadrature(radius):
    nodes, weights = roots_hermitenorm(48)
    likelihood = np.exp(-radius**2/2)*np.cosh(radius*nodes)
    moment = np.dot(weights, likelihood**2)/math.sqrt(2*math.pi)
    assert mixture_log_second_moment([radius]) == pytest.approx(math.log(moment), abs=8e-15)
    assert np.dot(weights, likelihood)/math.sqrt(2*math.pi) == pytest.approx(1., abs=3e-15)


def test_exact_product_moment_and_homogeneous_upper_bound():
    radii = np.array([.4, .7, 1.1])
    exact = np.prod(np.cosh(radii*radii))
    assert mixture_log_second_moment(radii) == pytest.approx(math.log(exact))
    equal = np.full(12, .3)
    log_alpha = -8.
    actual = fixed_signal_log_transfer_bound(log_alpha, equal)
    bound = homogeneous_log_transfer_bound(log_alpha, count=12, energy=np.dot(equal, equal))
    assert actual <= bound < 0
    assert math.exp(actual) == pytest.approx(math.sqrt(math.exp(log_alpha)*np.prod(np.cosh(equal**2))))


def test_log_bound_retains_tiny_probabilities_and_caps_at_one():
    assert homogeneous_log_transfer_bound(-1e6, count=100, energy=10) == pytest.approx(-499999.75)
    assert homogeneous_log_transfer_bound(-1, count=1, energy=10) == 0
    assert fixed_signal_log_transfer_bound(-math.inf, [.5]) == -math.inf


def test_diagonal_schedule_has_growing_information_and_vanishing_power_bound():
    # Synthetic null rejection rates check the algebra, not a procedure's rate.
    ds, bounds = [], []
    for exponent in [4, 8, 16, 32, 64]:
        n = 10**exponent
        log_alpha = -math.sqrt(math.log(n))
        d_squared = min(-log_alpha, math.sqrt(n))
        energy = math.sqrt(n*d_squared)
        bound = homogeneous_log_transfer_bound(log_alpha, count=n, energy=energy)
        assert bound <= .25*log_alpha+1e-14
        ds.append(d_squared)
        bounds.append(bound)
    assert np.all(np.diff(ds) > 0)
    assert np.all(np.diff(bounds) < 0)


def test_gate_exact_zero_padding_and_permutation_invariance():
    y = np.array([[1.3, .4], [-.2, .8], [.9, 1.1], [.1, -.3]])
    reference = washout_fit(y)
    for transformed in [y[::-1], np.vstack((y, np.zeros((5, 2))))]:
        fit = washout_fit(transformed)
        assert fit['washout_precision'] == pytest.approx(reference['washout_precision'], abs=2e-15)
        assert fit['energy'] == pytest.approx(reference['energy'], abs=2e-15)
    assert 0 < reference['washout_precision'] < reference['precision']
    assert washout_fit(np.zeros((4, 2)))['washout_precision'] == 0


def test_gate_suppresses_bounded_null_strength_and_preserves_strong_margin():
    energies = [1e3, 1e6, 1e12, 1e24]
    null = [gated_precision(2., e) for e in energies]
    assert np.all(np.diff(null) < 0)
    ratios = []
    for e in energies:
        b = math.log(math.e+e)
        k = b*b
        result = gated_precision(k, e)
        ratios.append(result/k)
        assert result/k == pytest.approx(b/(b+1))
    assert np.all(np.diff(ratios) > 0)


def test_gate_continuity_at_eigenvalue_tie():
    base = np.eye(2)
    strengths = []
    for epsilon in [1e-2, 1e-4, 1e-6, 0.]:
        y = base.copy()
        y[0, 0] += epsilon
        fit = washout_fit(y)
        strengths.append(fit['washout_precision'])
        assert 0 <= fit['washout_precision'] <= fit['curv']
    assert np.all(np.diff(strengths) < 0)


def test_gate_loses_consistency_on_an_explicit_slow_information_sequence():
    # Delta=sqrt(N)*(log N)^1/4 gives K~sqrt(log N), log(e+E)~log N.
    # Use log N as the independent variable to avoid constructing huge arrays.
    values = []
    for log_n in [10., 30., 100., 300., 700.]:
        energy = 2*math.exp(log_n)
        k = math.sqrt(log_n)
        values.append(gated_precision(k, energy))
    assert np.all(np.diff(values) > 0)
    assert .96 < values[-1] < 1
    # At arbitrarily large log N, the exact asymptotic expression tends to 1.
    log_n = 1e8
    k = math.sqrt(log_n)
    assert k*k/(k+log_n+math.log(2)) == pytest.approx(1., abs=1.1e-4)


def test_permuted_sign_mixture_second_moment_and_negative_association_bound():
    a = np.array([.1, .7, 1.5])
    perms = list(permutations(range(len(a))))
    # Direct Gaussian MGF for every pair of signed/permuted mean vectors.
    means = np.array([np.sqrt(a[list(p)])*signs
                      for p in perms for signs in product([-1., 1.], repeat=len(a))])
    direct = np.exp(means@means.T).mean()
    cosh_expression = np.mean([np.prod(np.cosh(np.sqrt(a*a[list(p)]))) for p in perms])
    assert direct == pytest.approx(cosh_expression, rel=3e-15)
    row_bound = symmetric_mixture_log_product_bound(a)
    assert math.log(direct) <= row_bound
    expanded_bound = a.sum()**2/(2*len(a)) + math.exp(a.max()**2/2)*(a@a)**2/(8*len(a))
    assert row_bound <= expanded_bound


def _angular_contrast(result):
    # Identity multiples cancel in the normalized mass law.
    matrix = result['matrix']
    return np.array([matrix[0, 0]-matrix[1, 1], 2*matrix[0, 1]])


def test_ordered_control_changes_under_permutation_and_zero_padding():
    y = np.array([[3., 0.], [0., 1.], [.3, .1]])
    reference = _angular_contrast(ordered_scatter(y))
    permuted = _angular_contrast(ordered_scatter(y[::-1]))
    padded = _angular_contrast(ordered_scatter(np.vstack((y, np.zeros((3, 2))))))
    assert np.linalg.norm(reference-permuted) > 1
    assert np.linalg.norm(reference-padded) > 1


def test_ordered_scatter_fourth_moment_by_gaussian_quadrature():
    # E||YY'-E YY'||_F²=6(1+r²) for a two-dimensional unit-noise row.
    nodes, weights = roots_hermitenorm(16)
    x, z = np.meshgrid(nodes+1.4, nodes, indexing='ij')
    w = np.outer(weights, weights)/(2*math.pi)
    squared = (x*x-(1+1.4**2))**2 + (z*z-1)**2 + 2*(x*z)**2
    assert np.sum(w*squared) == pytest.approx(6*(1+1.4**2), rel=4e-15)


@pytest.mark.parametrize('precision,energy', [(-1., 2.), (1., -1.), (math.inf, 2.)])
def test_gate_rejects_invalid_arguments(precision, energy):
    with pytest.raises(ValueError):
        gated_precision(precision, energy)

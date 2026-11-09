"""The identities and bounds behind the fixed-history half of exploration E10 (problems/O15.md, 9 Nov).

Claude's check of Codex's argument on `explore/e10-null-limit-prior`: a law that returns to the starting law on pure
noise, and does not depend on the order of the readings, is not consistent on every series with Delta_N/sqrt(N) -> oo.
The proof is in O15; these are its finite ingredients, computed independently of the branch's code. Deterministic.
"""
import itertools
import math

import numpy as np
import pytest
from scipy.special import logsumexp


def _log_second_moment_bound(a):
    """log of prod_i mean_k exp(a_i a_k / 2): the negative-association bound on E_0 L_sym^2."""
    a = np.asarray(a, float)
    return float((logsumexp(0.5 * np.outer(a, a), axis=1) - math.log(len(a))).sum())


def _log_permutation_mean(a):
    """log E_sigma exp(sum_i a_i a_sigma(i) / 2), by enumerating every permutation."""
    a = np.asarray(a, float)
    perms = np.array(list(itertools.permutations(range(len(a)))))
    return float(logsumexp(0.5 * (a[None, :] * a[perms]).sum(1)) - math.log(len(perms)))


def test_second_moment_of_the_sign_and_permutation_mixture():
    # E_0 L_sym^2 = E_sigma prod_i cosh(sqrt(a_i a_sigma(i))): left side by Gauss-Hermite quadrature on the data,
    # right side by enumeration. L_sym averages fair signs (analytically, a cosh) and all placements of the amplitudes.
    a = np.array([0.5, 0.3, 0.8])
    x, w = np.polynomial.hermite_e.hermegauss(60)
    w = w / w.sum()
    grid = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3)
    weight = np.einsum("i,j,k->ijk", w, w, w).ravel()
    perms = list(itertools.permutations(range(3)))
    ratio = np.zeros(len(grid))
    for placement in perms:
        placed = a[list(placement)]
        ratio += np.prod(np.exp(-placed / 2) * np.cosh(np.sqrt(placed) * grid), axis=1)
    ratio /= len(perms)
    by_quadrature = float(np.sum(weight * ratio * ratio))
    by_enumeration = float(np.mean([np.prod(np.cosh(np.sqrt(a * a[list(s)]))) for s in perms]))
    assert by_quadrature == pytest.approx(by_enumeration, rel=1e-11)
    # cosh(x) <= exp(x^2/2), so the enumeration is below E_sigma exp(sum a_i a_sigma(i)/2)
    assert math.log(by_enumeration) <= _log_permutation_mean(a)


@pytest.mark.parametrize("population", [
    [0.2, 0.5, 0.9, 1.4, 2.0, 2.7, 3.5],          # spread
    [1.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0],          # repeats and zeros
    [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 6.0],    # one dominant amplitude
    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],          # constant: equality
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7],          # one nonzero: equality
])
def test_negative_association_bound_for_a_permuted_population(population):
    # Joag-Dev and Proschan (1983): a uniform permutation of a fixed population is negatively associated, so
    # E prod_i f_i(a_sigma(i)) <= prod_i E f_i(a_sigma(i)) for nonnegative increasing f_i; here f_i(x) = exp(a_i x/2).
    assert _log_permutation_mean(population) <= _log_second_moment_bound(population) + 1e-12


def test_the_bound_is_strict_for_a_spread_population():
    spread = [0.2, 0.5, 0.9, 1.4, 2.0, 2.7, 3.5]
    assert _log_permutation_mean(spread) < _log_second_moment_bound(spread) - 0.1


@pytest.mark.parametrize("j", [1, 2, 3])
@pytest.mark.parametrize("n", [10, 100, 1000, 4000])
def test_the_stage_bound(j, n):
    # With a_i <= b_i = j/sqrt(i): log prod_i mean_k exp(b_i b_k/2) <= 2 j^2 + e^{j^2/2} j^4 (1 + log n)^2 / (8 n).
    b = j / np.sqrt(np.arange(1, n + 1))
    exact = _log_second_moment_bound(b)
    assert exact <= 2 * j * j + math.exp(j * j / 2) * j ** 4 * (1 + math.log(n)) ** 2 / (8 * n)


def test_the_history_identifies_while_every_reading_weakens():
    # a_i = j/sqrt(i) on (N_{j-1}, N_j] with N_j >= (j+1)^4: Delta_N/sqrt(N) >= (3/4) kappa_{ceil(N/4)}, and a_i <= 1/j.
    ends = [(j + 1) ** 4 * 3 for j in range(1, 10)]
    i = np.arange(1, ends[-1] + 1)
    stage = 1 + np.searchsorted(ends, i, side="left")
    a = stage / np.sqrt(i)
    delta = np.cumsum(a)
    for n in range(64, ends[-1], 211):
        assert delta[n - 1] / math.sqrt(n) >= 0.75 * stage[math.ceil(n / 4) - 1]
    for j in range(2, 10):
        assert a[stage == j].max() <= 1.0 / j
    assert delta[-1] / math.sqrt(len(i)) > 10

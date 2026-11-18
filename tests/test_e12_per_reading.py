"""Finite checks for the scalar per-reading proof in explore/e12-per-reading.md."""

import math

import numpy as np
import pytest
from scipy import integrate, special, stats

from tools.leftover_at_series_best_mass import law as o14_law


def _scatter(rows, weights):
    return rows.T @ (rows * weights[:, None])


def _principal(scatter):
    eigenvalues, vectors = np.linalg.eigh(scatter)
    return vectors[:, -1], eigenvalues[-1] - eigenvalues[0]


def _pilot(rows, *, positive=False):
    scatter = rows.T @ rows
    line, _ = _principal(scatter)
    if positive and line[0] * line[1] < 0:
        return np.eye(2)[np.argmax(np.diag(scatter))]
    return line


def _weights(rows, pilot, kappa):
    normal = np.array([pilot[1], -pilot[0]])
    squared = (rows @ normal) ** 2
    return kappa / np.maximum(kappa, squared)


def _weighted_fit(rows, *, positive=False, kappa=1.0):
    pilot = _pilot(rows, positive=positive)
    weighted = _scatter(rows, _weights(rows, pilot, kappa))
    line, gap = _principal(weighted)
    return pilot, weighted, line, gap


def _distance(first, second):
    return math.acos(np.clip(abs(first @ second), 0, 1))


def _closed_gap(kappa):
    a = math.sqrt(kappa)
    return 4 * (a * stats.norm.pdf(a) - kappa * special.ndtr(-a))


@pytest.mark.parametrize("kappa", [0.1, 0.5, 1, 4, 9])
def test_population_gap_matches_independent_gaussian_quadrature(kappa):
    a = math.sqrt(kappa)

    def gap_integrand(z):
        weight = kappa / max(kappa, z * z)
        return 2 * weight * (1 - z * z) * stats.norm.pdf(z)

    gap = integrate.quad(gap_integrand, 0, a)[0] + integrate.quad(gap_integrand, a, np.inf)[0]
    assert gap == pytest.approx(_closed_gap(kappa), rel=2e-10, abs=1e-14)
    assert gap > 0


def test_out_of_quadrant_whole_line_pilot_and_positive_endpoint():
    rows = np.array([[-3.0, 1.0], [-4.0, 1.0], [-2.0, 1.0]])
    whole_line = _pilot(rows)
    positive = _pilot(rows, positive=True)
    assert whole_line[0] * whole_line[1] < 0
    np.testing.assert_array_equal(positive, [1, 0])
    np.testing.assert_allclose(_weights(rows, positive, 0.5), [0.5, 0.5, 0.5])
    assert _distance(whole_line, positive) > 0.2


def test_continuous_weighted_gram_matches_o14_existing_interval():
    forces = [0.2, 1.0, 1.4]
    sigma = 0.1
    rows = np.column_stack((forces, np.ones(len(forces)))) / sigma
    pilot, weighted, _, _ = _weighted_fit(rows, positive=True)
    np.testing.assert_allclose(abs(pilot), [2**-0.5, 2**-0.5])
    # c = s = 1 gives a uniform angular prior. Compute the continuous quantiles
    # independently of O14's finite log-mass grid and rectangle summation.
    theta = np.linspace(0, math.pi / 2, 40_001)
    u = np.column_stack((np.sin(theta), np.cos(theta)))
    log_density = 0.5 * np.einsum("ij,jk,ik->i", u, weighted, u)
    cdf = integrate.cumulative_trapezoid(np.exp(log_density - log_density.max()), theta, initial=0)
    quantiles = np.tan(np.interp([0.5, 0.025, 0.975], cdf / cdf[-1], theta))
    scratch = o14_law(forces, sigma_f=sigma, sigma_a=sigma, c=1.0, kind="per reading")
    np.testing.assert_allclose(quantiles, scratch, rtol=7e-5)


@pytest.mark.parametrize("allocation", ["dense signed", "single high", "mixed"])
def test_identifying_examples_with_different_signal_allocations(allocation):
    n = 8192
    truth = np.array([2, 1]) / math.sqrt(5)
    rng = np.random.default_rng(121826)
    amplitude = np.zeros(n)
    if allocation == "dense signed":
        amplitude[::2], amplitude[1::2] = 0.7, -0.7
    elif allocation == "single high":
        amplitude[0] = -math.sqrt(n)
    else:
        amplitude[:4] = [80, -60, -90, 40]
        amplitude[4:] = 0.2
    rows = rng.standard_normal((n, 2)) + amplitude[:, None] * truth
    _, _, weighted_line, gap = _weighted_fit(rows)
    assert _distance(weighted_line, truth) < 0.085
    assert gap / n > 0.27


def test_bounded_reading_count_with_a_growing_signal():
    truth = np.array([2, 3]) / math.sqrt(13)
    noise = np.array([[0.3, 0.8], [1.4, 0.0], [-0.7, 1.2]])
    for radius in (1e2, 1e4):
        rows = noise.copy()
        rows[0] -= radius * truth
        _, _, line, gap = _weighted_fit(rows)
        assert _distance(line, truth) < 2 / radius
        assert gap > 0.7 * radius**2


def test_null_weighted_gap_tracks_the_random_pilot():
    rng = np.random.default_rng(92312)
    rows = rng.standard_normal((60_000, 2))
    for positive in (False, True):
        pilot, _, line, gap = _weighted_fit(rows, positive=positive)
        assert _distance(pilot, line) < 0.035
        assert gap / len(rows) == pytest.approx(_closed_gap(1), abs=0.025)
    unweighted = np.linalg.eigvalsh(rows.T @ rows)
    share = 1 - unweighted[0] / unweighted[1]
    assert share * gap > 10
    assert (unweighted[1] - unweighted[0]) / len(rows) < 0.03


def test_weak_array_scatter_comparison_uniform_in_pilot_on_finite_fixture():
    rng = np.random.default_rng(91916)
    null = rng.standard_normal((4000, 2))
    truth = np.array([2, 1]) / math.sqrt(5)
    amplitude = np.zeros(len(null))
    amplitude[::53] = (1 + rng.random(len(amplitude[::53]))) * 0.2
    rows = null + amplitude[:, None] * truth
    norm = np.linalg.norm(null, axis=1)
    for kappa in (0.5, 2):
        upper = np.sum(2 / math.sqrt(kappa) * abs(amplitude) * norm**2
                       + 2 * abs(amplitude) * norm + amplitude**2)
        for angle in np.linspace(0, math.pi, 11, endpoint=False):
            pilot = np.array([math.cos(angle), math.sin(angle)])
            difference = (_scatter(rows, _weights(rows, pilot, kappa))
                          - _scatter(null, _weights(null, pilot, kappa)))
            assert np.linalg.norm(difference, ord=2) <= upper


def test_gram_share_warp_cannot_reach_an_exponentially_small_raw_tail():
    # Deterministic scale comparison from the proof, evaluated in log space.
    # At these gaps an ordinary probability computation underflows.
    for n in (10_000, 1_000_000):
        h = 1 / math.sqrt(n)
        gap = _closed_gap(1) * n
        log_warp_rank = special.log_ndtr(special.ndtri(0.025) / math.sqrt(h))
        assert -0.1 * gap < log_warp_rank < 0
        assert h * gap > 25


@pytest.mark.parametrize("a", [0.0, 0.5, 3.0])
def test_smooth_gaussian_residual_weights_and_constant_control(a):
    rows = np.random.default_rng(71219).normal(size=(45_000, 2))
    pilot = _pilot(rows)
    normal = np.array([pilot[1], -pilot[0]])
    weights = np.exp(-0.5 * a * (rows @ normal) ** 2)
    _, gap = _principal(_scatter(rows, weights))
    # Gaussian integration gives A=(1+a)^(-1/2), B=(1+a)^(-3/2).
    assert gap / len(rows) == pytest.approx(a / (1 + a)**1.5, abs=0.021)


@pytest.mark.parametrize("positive", [False, True])
def test_pilot_residual_budget_controls_low_weight_count(positive):
    rng = np.random.default_rng(12033)
    truth = np.array([1.0, 2.0]) / math.sqrt(5)
    rows = rng.normal(size=(5000, 2)) + rng.normal(0, 80, (5000, 1)) * truth
    pilot = _pilot(rows, positive=positive)
    residual = rows @ np.array([pilot[1], -pilot[0]])
    truth_residual = rows @ np.array([truth[1], -truth[0]])
    assert residual @ residual <= truth_residual @ truth_residual + 1e-8
    assert 4 * np.count_nonzero(np.abs(residual) > 2) <= truth_residual @ truth_residual

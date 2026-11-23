"""Two readings that disagree, under two ways of combining separate readings' errors
(problems/O14.md, 23 Nov; tools/how_errors_combine.py). Measurements, not a rule."""
import math

import numpy as np
import pytest
from scipy import stats

from met.evidence import Curve
from met.full import isotropic


def _interval(xs, sigma, lengths, grid):
    xs = np.asarray(xs, float)[:, None]
    log = -np.sum(np.abs(xs - grid), 0) / (sigma / math.sqrt(2)) if lengths else -0.5 * np.sum((xs - grid) ** 2, 0) / sigma ** 2
    w = np.exp(log - log.max())
    c = np.cumsum(w / w.sum())
    return float(grid[np.searchsorted(c, 0.025)]), float(grid[np.searchsorted(c, 0.975)])


def test_direct_readings_squares_ignore_the_gap_and_lengths_spread_over_it():
    grid = np.linspace(4.0, 7.0, 300001)
    for pair in ((5.5, 5.6), (5.0, 6.0)):
        lo, hi = _interval(pair, 0.1, False, grid)
        assert 0.5 * (hi - lo) == pytest.approx(0.139, abs=2e-3)               # the same at any gap (Kac-Bernstein)
    lo, hi = _interval((5.5, 5.6), 0.1, True, grid)
    assert (lo, hi) == pytest.approx((5.425, 5.675), abs=2e-3)
    for sigma, half in ((0.1, 0.510), (0.01, 0.478)):                            # no longer follows sigma
        lo, hi = _interval((5.0, 6.0), sigma, True, grid)
        assert 0.5 * (hi - lo) == pytest.approx(half, abs=3e-3)


def test_the_lengths_statistic_under_the_declared_gaussian_is_exact_overall_and_pays_on_close_pairs():
    q = stats.norm.isf(math.sqrt(0.025))
    assert q == pytest.approx(1.002, abs=1e-3)
    e = np.random.default_rng(1).normal(size=(400000, 2))
    covered = (e.min(1) - q <= 0) & (0 <= e.max(1) + q)
    gap = np.abs(e[:, 0] - e[:, 1])
    assert covered.mean() == pytest.approx(0.95, abs=2e-3)
    assert covered[gap < 0.5].mean() == pytest.approx(0.887, abs=6e-3) and covered[gap > 2].mean() > 0.995


def test_lengths_tell_apart_two_series_with_one_gaussian_curve_without_a_count():
    # E14-6's pair: two readings a gap apart at small sigma; fifty honest readings at sigma 0.5
    n, gap = 50, 1.0
    s_a, s_b = gap / math.sqrt(2 * (n - 1)), gap / 2 * math.sqrt(n / (n - 1))
    z = np.random.default_rng(2).normal(size=n)
    honest = 5.5 + s_b * (z - z.mean()) / z.std(ddof=1)
    grid = np.linspace(4.5, 6.5, 200001)
    widths = {}
    for name, xs, s in (("far", [5.0, 6.0], s_a), ("honest", honest, s_b)):
        for lengths in (False, True):
            lo, hi = _interval(xs, s, lengths, grid)
            widths[name, lengths] = 0.5 * (hi - lo)
    assert widths["far", False] == pytest.approx(widths["honest", False], abs=2e-3)
    assert widths["far", True] > 0.45 and widths["honest", True] < 0.2


def test_on_force_and_acceleration_readings_lengths_tip_onto_the_stronger_reading_as_sigma_shrinks():
    def tilt(sigma):
        pair = [isotropic([5.0], [1.0], sigma, sigma), isotropic([6.0], [1.0], sigma, sigma)]
        length = lambda m: sum(math.sqrt(max(Curve([r], scale=1.0)(m), 0.0)) for r in pair)
        return length(5.0) - length(6.0), length(5.5) - length(6.0)
    assert tilt(0.1) == pytest.approx((0.3172, 0.1449), abs=2e-3)
    assert tilt(0.001)[0] == pytest.approx(31.72, abs=0.05)
    zeros = [isotropic([0.0], [0.0], 3.0, 0.02)] * 5
    pair = [isotropic([5.0], [1.0], 0.1, 0.1), isotropic([6.0], [1.0], 0.1, 0.1)]
    total = lambda rs: sum(math.sqrt(max(Curve([r], scale=1.0)(5.3), 0.0)) for r in rs)
    assert total(pair + zeros) == total(pair)

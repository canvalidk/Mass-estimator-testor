"""The curve's least under errors that are not Gaussian, and on the whole arc of masses (T84; L11,
G25, O19). Harvested 24 Nov 2026 from the exploration branch `explore/e17-l25-and-g25` (its
E17-G1, E17-G2 and E17-G3) and recomputed by Claude on main's own code, at the sizes the branch's
summary quotes and beyond. The simulations are in `tools/point_under_other_errors.py`, which
prints T84's tables.

What is measured: with errors that have mean zero and the declared covariance, the curve's
expected rise from the true mass is exactly the noiseless rise, whatever the errors' shape; the
least then sits within a fixed multiple of sqrt(D + Delta) / Delta of the true angle on the full
type; on the whole arc a far-off instrument's noise-only readings can carry the least away while
B27's weighted condition holds; and without a fourth moment the least fails at a signal growth
where Gaussian errors would let it converge."""
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import point_under_other_errors as study                 # noqa: E402

from met.evidence import Curve, leftover                 # noqa: E402
from met.full import isotropic                           # noqa: E402


@pytest.mark.parametrize("shape", ["skewed", "two-point"])
def test_the_curves_expected_rise_from_the_true_mass_needs_only_the_declared_covariance(shape):
    # thirty readings of the full type (mixed dimensions, declared subspaces, correlated covariances)
    noiseless, mean, se = study.expected_rise(shape, [0.3, 1.0, 2.0, 4.0, 12.0], reps=800)
    assert noiseless[2] == pytest.approx(0.0, abs=1e-9)
    for j in (0, 1, 3, 4):
        assert abs(mean[j] - noiseless[j]) < 4.0 * se[j]
        assert se[j] < 0.02 * noiseless[j] + 0.5


def test_the_least_keeps_the_rate_under_skewed_errors_on_the_full_type():
    # twenty series a size: 0.57 and 1.39 at the median; the longer run below is the one T84 quotes
    for n in (60, 180):
        median, high = study.normalised_errors("skewed", n, reps=20, seed=1000 + n)
        assert 0.3 < median < 2.0 and high < 4.0


@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="forty series of up to 720 full-type readings per shape")
@pytest.mark.parametrize("shape", ["gaussian", "skewed", "two-point"])
def test_the_least_keeps_the_rate_as_the_series_grows(shape):
    # the error times Delta / sqrt(D + Delta) does not grow with the series, and is no worse than Gaussian's
    for n in (60, 240, 720):
        median, high = study.normalised_errors(shape, n, reps=40, seed=1000 + n)
        assert 0.3 < median < 1.5 and high < 4.0


def test_on_the_whole_arc_a_far_instruments_noise_can_carry_the_least_away():
    # the formula used for speed, against met.evidence on a small series of the same kind
    rng = np.random.default_rng(0)
    f1, a1, x, y = 3.0 + rng.normal(), 3.0 + rng.normal(), rng.normal(size=7), rng.normal(size=7)
    series = [isotropic([f1], [a1], 1.0, 1.0)] + [isotropic([7 * xi], [yi], 7.0, 1.0) for xi, yi in zip(x, y)]
    for m in (0.0, 0.3, 1.0, 5.0, 40.0, math.inf):
        assert float(Curve(series, scale=1.0)(m)) == pytest.approx(float(study.far_instrument_curve(f1, a1, x, y, 7, m)), rel=1e-10)

    away, near_on_the_interval = {}, {}
    for n in (100, 2500):
        away[n], near_on_the_interval[n] = study.far_instrument_shares(n, 400, 1720 + n)
        assert study.far_instrument_conditions(n)[0] > 8.0                   # B27's weighted condition grows
    assert away[100] == pytest.approx(0.25, abs=0.06) and away[2500] == pytest.approx(0.48, abs=0.06)
    assert near_on_the_interval[2500] > near_on_the_interval[100] > 0.8      # on a declared interval it converges
    # L11's own, unweighted condition falls: by the premise list this is not an identifying sequence
    assert study.far_instrument_conditions(2500)[1] < study.far_instrument_conditions(100)[1] < 4.0


def test_without_a_fourth_moment_the_least_fails_where_gaussian_errors_let_it_converge():
    # the grid point used for speed, against met.evidence
    rng = np.random.default_rng(3)
    f, a = 2.0 + rng.normal(size=6), 1.0 + rng.normal(size=6)
    series = [isotropic([x], [y], 1.0, 1.0) for x, y in zip(f, a)]
    grid_point = study.MASSES[study.one_instrument_curve(f, a).argmin()]
    assert grid_point == pytest.approx(leftover(series, scale=1.0)["mass"], rel=0.02)

    # Delta = 8 N^0.65: grows faster than sqrt(N), slower than N^0.8 (the 2.5-degree threshold)
    gaussian = [study.misses(None, n, 0.65, 8.0, 200, 1730 + n) for n in (100, 1600)]
    five = [study.misses(5.0, n, 0.65, 8.0, 200, 1731 + n) for n in (100, 1600)]
    heavy = [study.misses(2.5, n, 0.65, 8.0, 200, 1732 + n) for n in (100, 1600)]
    assert gaussian[1] < 0.12 < gaussian[0]
    assert five[1] < five[0] and five[1] < 0.25
    assert heavy[1] > heavy[0] - 0.02 and heavy[1] > 0.4
    # above that threshold it comes back
    assert study.misses(2.5, 1600, 0.95, 6.0, 200, 1733) < 0.12
    # a point read from sum log(1 + K_i / 4) converges under the same heavy errors
    log_heavy = [study.misses(2.5, n, 0.65, 8.0, 150, 1734 + n, log_curve=True) for n in (100, 1600)]
    assert log_heavy[1] < 0.05 < log_heavy[0]

"""Two plain cuts of the curve, and the cut read in the unit of the fluke that happened
(R-2026-11-24-4; problems/O14.md, 24 Nov evening). Measurements behind that section, not a rule.

A cut at a fixed height is right at its level and reacts to disagreement the wrong way: the set
narrows, then is empty. A cut at a fixed rise above the curve's least is never empty and ignores
disagreement, and with many weak readings it is not right at its level. Reading the rise in the
larger of sigma's unit and the leftover per direction that can be seen gives Can's two scenes of
17 Nov on direct readings (the 18 Nov pair weight's numbers, T82 (7))."""
import math

import numpy as np
import pytest
from scipy import stats

from met.evidence import Curve, leftover
from met.full import isotropic

MASSES = np.linspace(3.0, 8.0, 500001)
RISE_CUT = float(stats.chi2.ppf(0.95, 1))            # 3.84: one direction, the mass
HEIGHT_CUT = float(stats.chi2.ppf(0.95, 2))          # 5.99: two readings, two noise components


def _pair(first, second):
    """Two direct readings of a mass, each +- 0.1 (the acceleration read as 1 with next to no error)."""
    return [isotropic([first], [1.0], 0.1, 1e-6), isotropic([second], [1.0], 0.1, 1e-6)]


def _half_width(keep):
    kept = MASSES[keep]
    return 0.0 if kept.size == 0 else 0.5 * float(kept.max() - kept.min())


def test_a_cut_at_a_fixed_height_narrows_under_disagreement_and_then_keeps_nothing():
    widths = {}
    for name, (first, second) in {"agree": (5.5, 5.6), "mild": (5.5, 5.75), "far": (5.0, 6.0)}.items():
        widths[name] = _half_width(Curve(_pair(first, second), scale=1.0)(MASSES) <= HEIGHT_CUT)
    assert widths["agree"] == pytest.approx(0.166, abs=0.002)
    assert widths["mild"] == pytest.approx(0.120, abs=0.002)      # agreeing less, and narrower
    assert widths["far"] == 0.0                                    # nothing is kept
    assert leftover(_pair(5.0, 6.0), scale=1.0)["least"] == pytest.approx(50.0, rel=1e-4)


def test_a_cut_at_a_fixed_rise_is_never_empty_and_does_not_see_disagreement():
    widths = []
    for first, second in ((5.5, 5.6), (5.5, 5.75), (5.0, 6.0)):
        curve = Curve(_pair(first, second), scale=1.0)(MASSES)
        widths.append(_half_width(curve - curve.min() <= RISE_CUT))
    assert widths == pytest.approx([0.139, 0.139, 0.139], abs=0.002)


def test_the_rise_read_in_the_unit_of_the_fluke_gives_the_two_scenes():
    # the unit: the larger of 1 (sigma's own) and the leftover per direction that can be seen
    # (two noise components, one of them the mass: one seen direction)
    widths = []
    for first, second in ((5.5, 5.6), (5.5, 5.75), (5.0, 6.0)):
        curve = Curve(_pair(first, second), scale=1.0)(MASSES)
        unit = max(1.0, float(curve.min()) / 1.0)
        widths.append(_half_width(curve - curve.min() <= RISE_CUT * unit))
    assert widths[0] == pytest.approx(0.139, abs=0.002)            # an agreeing pair is left alone
    assert widths[1] == pytest.approx(0.245, abs=0.003)
    assert widths[2] == pytest.approx(0.980, abs=0.003)            # 5 and 6: about +- 1


ANGLES = np.linspace(0.0, 0.5 * math.pi, 2001)
COS, SIN = np.cos(ANGLES), np.sin(ANGLES)


def _held(n, energy, reps, seed):
    """Share of series in which each plain cut keeps the true mass 2, and in which the fixed-height
    cut keeps nothing; n scalar readings on one instrument, total signal energy `energy`."""
    rng = np.random.default_rng(seed)
    true_angle = math.atan(2.0)
    size = math.sqrt(energy / n)
    height = float(stats.chi2.ppf(0.95, n))
    rise = fixed = empty = 0
    for _ in range(reps):
        f = size * math.sin(true_angle) + rng.normal(size=n)
        a = size * math.cos(true_angle) + rng.normal(size=n)
        curve = COS * COS * (f @ f) - 2.0 * SIN * COS * (f @ a) + SIN * SIN * (a @ a)
        at_truth = (math.cos(true_angle) ** 2 * (f @ f) - 2.0 * math.sin(true_angle) * math.cos(true_angle) * (f @ a)
                    + math.sin(true_angle) ** 2 * (a @ a))
        rise += at_truth - curve.min() <= RISE_CUT
        fixed += at_truth <= height
        empty += curve.min() > height
    return rise / reps, fixed / reps, empty / reps


def test_the_fixed_rise_is_not_right_at_its_level_on_many_weak_readings_and_the_fixed_height_is():
    rng = np.random.default_rng(1)
    f, a = 2.0 + rng.normal(size=5), 1.0 + rng.normal(size=5)
    series = [isotropic([x], [y], 1.0, 1.0) for x, y in zip(f, a)]
    mine = COS * COS * (f @ f) - 2.0 * SIN * COS * (f @ a) + SIN * SIN * (a @ a)
    assert mine[1:-1] == pytest.approx(Curve(series, scale=1.0)(np.tan(ANGLES[1:-1])), rel=1e-9)

    no_signal = {n: _held(n, 0.0, 3000, 2411 + n) for n in (1, 5, 20, 100)}
    assert [no_signal[n][0] for n in (1, 5, 20, 100)] == pytest.approx([0.975, 0.82, 0.54, 0.33], abs=0.03)
    assert all(no_signal[n][1] == pytest.approx(0.95, abs=0.015) for n in (1, 5, 20, 100))
    strong = _held(20, 100.0, 3000, 2499)
    assert strong[0] == pytest.approx(0.93, abs=0.02) and strong[1] == pytest.approx(0.95, abs=0.015)
    assert 0.01 < strong[2] < 0.07                                 # and it keeps nothing in a few series in a hundred

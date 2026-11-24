""""Same curve, same result" (called L17 or reading splitting on the exploration branches E14 to
E17), against the pages (`problems/O20.md`, 24 Nov). Two things are measured here. The idealised
split, one reading against n readings 1/sqrt(n) the size (equally: the same values at n times the
variance), leaves the curve, the score, the point and the Gaussian score set alone. Cutting a real
window into sub-windows does not: each sub-window carries its own noise, and the Gaussian rule's own
point and set move, with no discount in play. That discount rules move under the idealised split,
and that two series can share a curve and differ in scatter, is in
`tests/test_common_scale_discount.py`."""
import math

import numpy as np
import pytest

from met.conditional_score import gaussian_score_test
from met.evidence import Curve, leftover
from met.full import Reading, isotropic

Z95 = 1.959964


def test_the_idealised_split_leaves_the_curve_the_score_the_point_and_the_set_alone():
    other = [isotropic([3.0], [2.0], 1.0, 1.0)]
    whole = [isotropic([5.0], [1.0], 1.0, 1.0)] + other
    masses = np.geomspace(0.1, 40.0, 25)
    for n in (10, 50):
        smaller = [isotropic([5.0 / math.sqrt(n)], [1.0 / math.sqrt(n)], 1.0, 1.0)] * n + other
        noisier = [Reading([5.0], [1.0], n * np.eye(2))] * n + other
        for parts in (smaller, noisier):
            assert Curve(parts, scale=1.0)(masses) == pytest.approx(Curve(whole, scale=1.0)(masses), rel=1e-10)
            for m in (0.7, 2.0, 6.0):
                assert (gaussian_score_test(parts, m)["normal_score"]
                        == pytest.approx(gaussian_score_test(whole, m)["normal_score"], rel=1e-9))
            assert leftover(parts, scale=1.0)["mass"] == pytest.approx(leftover(whole, scale=1.0)["mass"], rel=1e-6)
            assert leftover(parts, scale=1.0)["least"] == pytest.approx(leftover(whole, scale=1.0)["least"], abs=1e-8)


# One instrument at noise ratio 1, scalar readings in noise units (x = F / sigma_F, y = a / sigma_a).
# With m = tan(theta): the curve is sum (x cos - y sin)^2, the score sum U V over sqrt(sum U^2) with
# U = x sin + y cos along the candidate line and V = x cos - y sin across it. Both are checked against
# met.evidence and met.conditional_score in the test below before they are used.

ANGLES = np.linspace(0.0, 0.5 * math.pi, 20001)
GRID = np.geomspace(0.2, 20.0, 2001)


def _point_angle(x, y):
    curve = (np.cos(ANGLES) ** 2 * (x @ x) - 2.0 * np.sin(ANGLES) * np.cos(ANGLES) * (x @ y)
             + np.sin(ANGLES) ** 2 * (y @ y))
    return float(ANGLES[curve.argmin()])


def _kept(x, y, masses):
    angle = np.arctan(np.asarray(masses, dtype=float))
    s, c = np.sin(angle)[:, None], np.cos(angle)[:, None]
    along, across = x[None, :] * s + y[None, :] * c, x[None, :] * c - y[None, :] * s
    return np.abs((along * across).sum(1) / np.sqrt((along * along).sum(1))) <= Z95


def test_cutting_a_real_window_moves_the_gaussian_rules_own_point_and_set():
    rng = np.random.default_rng(2411)
    x, y = 3.0 * rng.normal(size=4), 3.0 * rng.normal(size=4)
    series = [isotropic([a], [b], 1.0, 1.0) for a, b in zip(x, y)]
    for m in (0.5, 2.0, 7.0):
        angle = math.atan(m)
        assert Curve(series, scale=1.0)(m) == pytest.approx(float(np.sum((x * math.cos(angle) - y * math.sin(angle)) ** 2)), rel=1e-9)
        assert bool(_kept(x, y, [m])[0]) == (gaussian_score_test(series, m)["p_value"] >= 0.05)
    assert math.tan(_point_angle(x, y)) == pytest.approx(leftover(series, scale=1.0)["mass"], rel=2e-3)

    # A body of mass 2 under a steady acceleration of 5 noise units over one window. Read once over
    # the window (unit sigma), or read in 50 sub-windows (each with sqrt(50) times the noise, and
    # so the same information in total). The whole-window reading is the mean of the sub-window ones.
    mass, acceleration, n, reps = 2.0, 5.0, 50, 600
    moved, error_whole, error_cut, width_whole, width_cut, held_whole, held_cut = [], [], [], [], [], 0, 0
    for _ in range(reps):
        force = mass * acceleration + math.sqrt(n) * rng.normal(size=n)
        acc = acceleration + math.sqrt(n) * rng.normal(size=n)
        xw, yw = np.array([force.mean()]), np.array([acc.mean()])
        xs, ys = force / math.sqrt(n), acc / math.sqrt(n)
        whole, cut = _point_angle(xw, yw), _point_angle(xs, ys)
        moved.append(abs(cut - whole))
        error_whole.append(abs(whole - math.atan(mass)))
        error_cut.append(abs(cut - math.atan(mass)))
        kept_whole, kept_cut = _kept(xw, yw, GRID), _kept(xs, ys, GRID)
        width_whole.append(np.ptp(np.log(GRID[kept_whole])))
        width_cut.append(np.ptp(np.log(GRID[kept_cut])))
        held_whole += bool(_kept(xw, yw, [mass])[0])
        held_cut += bool(_kept(xs, ys, [mass])[0])
    # the point moves by more than half of its own typical error, and the set is a fifth wider
    assert np.median(moved) == pytest.approx(0.037, abs=0.008)
    assert np.median(moved) > 0.5 * np.median(error_whole)
    assert np.median(error_cut) > 1.1 * np.median(error_whole)
    assert np.median(width_cut) / np.median(width_whole) == pytest.approx(1.22, abs=0.06)
    # and both are right at their level: two valid answers that differ
    assert held_whole / reps == pytest.approx(0.95, abs=0.025)
    assert held_cut / reps == pytest.approx(0.95, abs=0.025)

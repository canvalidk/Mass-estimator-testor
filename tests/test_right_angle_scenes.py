"""Readings at right angles in the force-acceleration plane (Claude, 18 Nov 2026; problems/O14.md, 18 Nov).

Scene C: one reading with force 2 and no acceleration, one with acceleration 1 and no force, each +- sigma. No mass
explains both (the leftover at the best line is 1/sigma^2), so by T30 as restated on 17 Nov the law should stay spread
as sigma shrinks. These tests record what the equations in met/ do with it.
"""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import right_angle_scenes as ras  # noqa: E402

EQUATIONS = ("ours29", "ours30", "ours33, nu = inf", "ours33, nu = 10")


@pytest.mark.parametrize("name", EQUATIONS)
def test_every_equation_sends_the_point_off_as_sigma_shrinks_when_the_pair_is_two_scalar_readings(name):
    # the point grows like 1/sigma and the interval's lower end rises with it: growing confidence in a huge mass
    # from readings no mass explains. eta is 1 throughout: both readings have zero score at the best line (J = 0).
    wide, narrow = ras.run(name, 1.0), ras.run(name, 0.05)
    assert wide["median"] == pytest.approx(1.9, abs=0.15)
    assert narrow["median"] > 50 and narrow["lo"] > 14
    assert wide["eta"] == pytest.approx(1.0) and narrow["eta"] == pytest.approx(1.0)
    if narrow["g"] is not None:
        assert narrow["g"] == pytest.approx([1.0, 1.0])


def test_the_point_scales_as_one_over_sigma():
    a, b = ras.run("ours30", 0.1), ras.run("ours30", 0.05)
    assert b["median"] / a["median"] == pytest.approx(2.0, rel=0.02)


@pytest.mark.parametrize("name,median", (("ours29", 1.97), ("ours30", 2.52)))
def test_the_same_numbers_declared_as_one_2d_reading_are_held_by_g(name, median):
    # the contradiction packaged inside one reading is discounted (T30); spread over two readings it is not
    for sigma in (0.5, 0.2, 0.05):
        r = ras.run(name, sigma, split=False)
        assert r["median"] == pytest.approx(median, abs=0.01)
        assert r["g"][0] < 1


def test_ours33_without_a_handed_over_number_does_not_hold_the_2d_reading_either():
    assert ras.run("ours33, nu = inf", 0.05, split=False)["median"] > 50


def test_one_weight_per_reading_grows_arbitrarily_sure_as_a_stray_nears_the_right_angle():
    # two readings on m = 1.73 and one eps degrees off the right angle to them, sigma -> 0. Read by count the interval
    # stays about [0.5, 13] whatever eps; read per reading its limit narrows with eps, to a point at eps = 0.
    widths = []
    for eps in (20.0, 5.0, 1.0):
        lo, hi = ras.near_right_angle(eps, 0.002, "one weight per reading", grid=40001)
        widths.append(math.log(hi / lo))
        lo, hi = ras.near_right_angle(eps, 0.002, "one weight, by count", grid=40001)
        assert lo < 0.6 and hi > 8
    assert widths[0] > 2 * widths[1] > 4 * widths[2]

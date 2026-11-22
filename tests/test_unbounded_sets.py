"""Why a true statement about mass has to be able to keep every large mass (T80). The curve's
set for one reading has four shapes, as a ratio's does (Fieller 1954); its ends move continuously
in the noise angle and reach infinite mass at finite data; and an always-finite interval is right
about a large true mass with probability going to zero when the signal is absent, where the
exact conditional test is right at its level."""
import math

import numpy as np
import pytest

from met.conditional_lr import conditional_test
from met.evidence import Curve
from met.full import isotropic
from met.scalar_witness import witness

CUT = 3.8415                                                    # the chi-square threshold at 95%


def _kept(force, acceleration, masses):
    return Curve([isotropic([force], [acceleration], 1.0, 1.0)])(np.asarray(masses, dtype=float)) <= CUT


@pytest.mark.parametrize("force,acceleration,low_end,high_end,shape", [
    (10.0, 5.0, False, False, "an interval"),
    (10.0, 1.0, False, True, "at least so much"),
    (1.0, 10.0, True, False, "at most so much"),
    (1.0, 1.0, True, True, "every positive mass"),
])
def test_one_reading_has_four_shapes_decided_by_which_channel_can_be_told_from_zero(force, acceleration, low_end, high_end, shape):
    # unit sigmas. The curve at mass 0 is F^2 / sigma_F^2 (the force is left unexplained) and at
    # infinite mass a^2 / sigma_a^2 (the acceleration is): an end is kept exactly when that channel
    # is within 1.96 sigma of zero.
    curve = Curve([isotropic([force], [acceleration], 1.0, 1.0)])
    assert curve(0.0) == pytest.approx(force ** 2) and curve(math.inf) == pytest.approx(acceleration ** 2)
    assert bool(_kept(force, acceleration, [0.0])[0]) is low_end, shape
    assert bool(_kept(force, acceleration, [math.inf])[0]) is high_end, shape
    assert _kept(force, acceleration, [force / acceleration])[0]          # the reading's own ratio is always kept
    if shape == "every positive mass":
        assert _kept(force, acceleration, np.geomspace(1e-6, 1e6, 200)).all()


def test_the_upper_end_reaches_infinite_mass_at_finite_data_and_continuously_in_the_angle():
    # force 10, acceleration falling towards 1.96 sigma: the largest kept mass grows without bound
    # and is infinite from 1.96 down, while its noise angle climbs smoothly to a quarter turn
    angles = []
    for acceleration in (3.0, 2.5, 2.1, 2.0, 1.97, 1.961):
        masses = np.geomspace(1.0, 1e9, 400001)
        top = masses[_kept(10.0, acceleration, masses)].max()
        assert top < 1e9
        angles.append(math.atan(top))
    assert all(b > a for a, b in zip(angles, angles[1:])) and angles[-1] > 0.5 * math.pi - 0.01
    steps = np.diff(angles)
    assert steps.max() < 0.25                                             # no jump on the way
    assert _kept(10.0, 1.959, [math.inf])[0] and _kept(10.0, 1.0, [1e12])[0]


def test_an_always_finite_interval_is_wrong_about_large_masses_when_there_is_no_signal():
    # five noise-only scalar readings. Their law is the same whatever the true mass is, so an
    # interval that is always finite cannot be right about every mass: the witness's 95% interval
    # (a proper law, centre 1) holds the mass 1 in nearly every series and the mass 10^4 in none.
    # No widening helps while the interval stays finite; the exact test keeps 10^4 at its level.
    rng = np.random.default_rng(80)
    reps, small, large, kept_large = 300, 0, 0, 0
    for rep in range(reps):
        f, a = rng.normal(size=5), rng.normal(size=5)
        lo, hi = witness(f, a, 1.0, 1.0, alpha=1 / 8, centre=1.0)["interval"]
        small += lo <= 1.0 <= hi
        large += lo <= 1e4 <= hi
        series = [isotropic([x], [y], 1.0, 1.0) for x, y in zip(f, a)]
        kept_large += conditional_test(series, 1e4, draws=199, seed=rep)["p_value"] > 0.05
    assert small / reps > 0.9 and large == 0
    assert kept_large / reps == pytest.approx(0.95, abs=0.04)

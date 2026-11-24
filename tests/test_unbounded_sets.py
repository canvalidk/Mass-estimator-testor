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


# T80 (6), 24 Nov (R-2026-11-24-3): a rule that is valid only where there is a signal, behind a
# check on the readings. One scalar reading, unit sigmas, true mass 10 (ten times the instruments'
# ratio), true acceleration `signal` noise units. The textbook rule is F/a with the propagated
# error, sigma_m/m = sqrt((sigma_F/F)^2 + (sigma_a/a)^2), at 1.96. The switched rule returns that
# interval when the measured acceleration is more than `threshold` sigma from zero and every
# positive mass otherwise.

SIGNALS = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0)


def _one_reading_worlds(mass, signal, seed, n=200_000):
    rng = np.random.default_rng(seed)
    return mass * signal + rng.normal(size=n), signal + rng.normal(size=n)


def _textbook_holds(force, acceleration, mass):
    estimate = force / acceleration
    half = 1.96 * np.abs(estimate) * np.sqrt(1.0 / force ** 2 + 1.0 / acceleration ** 2)
    return np.abs(estimate - mass) <= half


def test_the_textbook_interval_is_right_only_where_there_is_a_signal():
    held = {}
    for signal in SIGNALS:
        f, a = _one_reading_worlds(10.0, signal, seed=8060 + int(10 * signal))
        held[signal] = _textbook_holds(f, a, 10.0).mean()
    assert held[0.0] == pytest.approx(0.31, abs=0.01)             # no signal: right 31% of the time
    assert held[2.0] == pytest.approx(0.888, abs=0.005)
    assert held[12.0] == pytest.approx(0.950, abs=0.003)          # strong signal: it works
    values = [held[s] for s in SIGNALS]
    assert all(b > a - 0.004 for a, b in zip(values, values[1:]))  # and it gets there gradually


@pytest.mark.parametrize("threshold,worst,at_signal", [(1.96, 0.845, 1.0), (3.0, 0.889, 2.0), (5.0, 0.926, 4.0)])
def test_a_check_on_the_readings_lets_through_the_series_the_rule_is_wrong_on(threshold, worst, at_signal):
    # The textbook interval misses a large mass when the acceleration was read too HIGH (the mass
    # comes out too small, with too small an error). Those readings pass any "is the acceleration
    # far enough from zero" check. So the switched rule's worst case is close to the textbook
    # rule's own value at that signal (equal to it for the two higher thresholds); a higher
    # threshold moves where the worst case is and does not remove it.
    lowest, where = 1.0, None
    for signal in SIGNALS:
        f, a = _one_reading_worlds(10.0, signal, seed=8061 + int(10 * signal))
        held = np.where(np.abs(a) > threshold, _textbook_holds(f, a, 10.0), True).mean()
        if held < lowest:
            lowest, where = held, signal
    assert where == at_signal
    assert lowest == pytest.approx(worst, abs=0.006) and lowest < 0.93


def test_among_the_series_that_pass_the_check_the_rule_can_be_wrong_every_time():
    # true acceleration 3 sigma, check at 5 sigma: 2.3% of series pass, each with the acceleration
    # read at least 2 sigma too high, and the textbook interval holds the true mass in none of them
    f, a = _one_reading_worlds(10.0, 3.0, seed=8062, n=400_000)
    passed = np.abs(a) > 5.0
    assert passed.mean() == pytest.approx(0.0228, abs=0.002)
    assert _textbook_holds(f, a, 10.0)[passed].mean() < 0.002


def test_the_curves_set_is_right_at_its_level_at_every_signal_with_no_check():
    # the curve of one scalar reading at unit sigmas is (F - m a)^2 / (1 + m^2) (checked against
    # met.evidence on a few readings); cut at the 95% threshold it holds the true mass in 95% of
    # series at every signal, and it is the same rule at every signal
    rng = np.random.default_rng(8063)
    for _ in range(5):
        f, a, m = rng.normal(scale=4.0), rng.normal(scale=4.0), float(rng.uniform(0.2, 20.0))
        assert Curve([isotropic([f], [a], 1.0, 1.0)])(m) == pytest.approx((f - m * a) ** 2 / (1.0 + m * m), rel=1e-9, abs=1e-9)
    for signal in SIGNALS:
        f, a = _one_reading_worlds(10.0, signal, seed=8064 + int(10 * signal))
        assert ((f - 10.0 * a) ** 2 / 101.0 <= CUT).mean() == pytest.approx(0.95, abs=0.003)

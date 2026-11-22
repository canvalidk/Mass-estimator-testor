"""The witness family on one scalar instrument (T78; exploration E12, Codex, 18 Nov). Written again
from the formulas in met/scalar_witness.py; these tests reproduce E12's reported numbers on that
independent code, check the four properties its theorem states, and record two things the check
found: the discount goes away as the acceleration alone is made more precise, and the cap for a
reading on an axis misfires at a large noise ratio. Coverage runs are seeded; the long ones need
MET_SLOW=1."""
import math
import os
import sys

import numpy as np
import pytest
from scipy import integrate, stats

from met.evidence import Curve
from met.full import isotropic
from met.scalar_witness import summary, witness

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import leftover_at_series_best_mass as nov17                    # noqa: E402

SLOW = pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="hundreds of simulated series")
NEAR, FAR, ONES = (5.5, 5.6), (5.0, 6.0), (1.0, 1.0)


def _half(interval):
    return 0.5 * (interval[1] - interval[0])


# ---- E12's reported numbers, on independent code --------------------------------------------

@pytest.mark.parametrize("alpha,far_discount,far_half", [(1.0, 1.0, 0.13944), (1 / 8, 0.12841, 0.38914),
                                                         (1 / 16, 0.07210, 0.51934)])
def test_the_two_scenes_under_three_members_of_the_family(alpha, far_discount, far_half):
    # prior uniform in the noise angle, b = 1 (e12-pooled-witness.md, its table). The near pair is
    # never discounted; the far pair is, by an amount that depends on the member.
    near, far = (witness(scene, ONES, 0.1, 0.002, alpha=alpha) for scene in (NEAR, FAR))
    assert near["discount"] == 1.0 and _half(near["interval"]) == pytest.approx(0.13944, abs=2e-5)
    assert far["discount"] == pytest.approx(far_discount, abs=1e-5)
    assert _half(far["interval"]) == pytest.approx(far_half, abs=2e-5)


def test_two_members_meet_every_property_and_disagree_on_the_far_pair():
    # e12-scalar-proof.md section 3, half-Cauchy centre 5: the nonuniqueness, on the project's scenes
    got = {alpha: witness(FAR, ONES, 0.1, 0.002, alpha=alpha, centre=5.0)["interval"] for alpha in (1 / 8, 1 / 16)}
    assert got[1 / 8] == pytest.approx((5.104, 5.882), abs=1e-3)
    assert got[1 / 16] == pytest.approx((4.968, 6.007), abs=1e-3)
    for alpha in (1 / 8, 1 / 16):
        assert witness(NEAR, ONES, 0.1, 0.002, alpha=alpha, centre=5.0)["interval"] == pytest.approx((5.410, 5.689), abs=1e-3)


# ---- property 1: exact symmetries -----------------------------------------------------------

def test_zero_rows_order_signs_units_and_the_swap():
    f, a = np.array([0.2, 1.0, 1.4, 2.2]), np.array([1.0, 1.0, 1.0, 1.9])
    base = witness(f, a, 0.3, 0.2, alpha=1 / 8, centre=1.3)
    padded = witness(np.r_[f, 0.0, 0.0], np.r_[a, 0.0, 0.0], 0.3, 0.2, alpha=1 / 8, centre=1.3)
    assert padded["point"] == base["point"] and padded["interval"] == base["interval"]
    assert witness(f[::-1], a[::-1], 0.3, 0.2, alpha=1 / 8, centre=1.3)["point"] == pytest.approx(base["point"], rel=1e-12)
    flip = np.array([1, -1, 1, -1])
    assert witness(f * flip, a * flip, 0.3, 0.2, alpha=1 / 8, centre=1.3)["point"] == pytest.approx(base["point"], rel=1e-12)
    assert witness(1000 * f, a, 300.0, 0.2, alpha=1 / 8, centre=1300.0)["point"] == pytest.approx(1000 * base["point"], rel=1e-10)
    swapped = witness(a, f, 0.2, 0.3, alpha=1 / 8, centre=1 / 1.3)
    assert 1 / swapped["point"] == pytest.approx(base["point"], rel=1e-9)
    assert 1 / swapped["interval"][1] == pytest.approx(base["interval"][0], rel=1e-9)
    # the all-zero experiment returns the starting law: its median is the declared centre
    assert witness([0.0, 0.0], [0.0, 0.0], 0.3, 0.2, alpha=1 / 8, centre=1.3)["point"] == pytest.approx(1.3, rel=1e-4)


# ---- property 2: fixed readings as both sigmas shrink ----------------------------------------

def test_three_readings_no_mass_fits_stay_spread_where_a_weight_per_reading_collapses():
    # forces 0.2, 1, 1.4 at acceleration 1 (e12-scalar-proof.md section 1). One reading sits exactly
    # on the series' best line, so a weight per reading leaves it at full weight and the law
    # collapses on it; the pooled weights do not. This is O14's "one weight on the result".
    forces, sigmas = [0.2, 1.0, 1.4], (0.1, 0.01)
    per_reading = [nov17.law(forces, s, s, c=1.0, kind="per reading") for s in sigmas]
    assert per_reading[0][1:] == pytest.approx((0.8026, 1.3479), abs=1e-3)
    assert per_reading[1][1:] == pytest.approx((0.9731, 1.0285), abs=1e-3)
    for s in sigmas:
        assert nov17.law(forces, s, s, c=1.0, kind="pooled")[1:] == pytest.approx((0.4972, 2.0113), abs=1e-3)
    for s in (0.1, 0.01, 0.001):
        got = witness(forces, [1.0] * 3, s, s, alpha=1 / 8, centre=1.0)
        assert (got["interval"][0], got["point"], got["interval"][1]) == pytest.approx((0.5495, 1.0, 1.8198), abs=1e-3)


def test_fixed_shapes_as_the_sigmas_shrink():
    limit = lambda f, a: [witness(f, a, s, s, alpha=1 / 8, centre=1.0)["interval"] for s in (0.1, 0.01, 0.001)]
    compatible = limit((2.0, 4.0), (1.0, 2.0))                  # both readings say 2: concentrates
    assert _half(compatible[2]) < 0.01 * _half(compatible[0]) + 0.001 and compatible[2][0] < 2.0 < compatible[2][1]
    for shape in (((2.0, 1.0), (-1.0, -2.0)),                   # opposed readings
                  ((1.0,), (0.0,)),                             # one reading on the force axis
                  ((1.0, 0.0), (0.0, 1.0))):                    # right angles: equal eigenvalues (scene C)
        widths = [math.log(hi / lo) for lo, hi in limit(*shape)]
        assert widths[2] == pytest.approx(widths[0], rel=1e-3) and widths[2] > 5.0
    # the limiting precision of the three-reading shape is 3 H^2 / (Q L) * (L / Q)^alpha: 14^alpha apart
    rows = np.array([[0.2, 1.0], [1.0, 1.0], [1.4, 1.0]]) / 1e-4
    ratio = summary(rows, 1 / 8)["precision"] / summary(rows, 1 / 16)["precision"]
    assert ratio == pytest.approx(14 ** (1 / 16), rel=1e-6)


# ---- what the check found --------------------------------------------------------------------

def test_the_discount_goes_away_as_the_acceleration_alone_is_made_more_precise():
    # the far pair, forces 5 and 6 each +- 0.1. E12 notes it; measured 22 Nov: 0.128, 0.225, 0.400,
    # 0.711, 1 at sigma_a = 2e-3 ... 2e-7. So for two direct readings of a mass (an acceleration
    # known exactly) this family discounts nothing: its discount is not a function of gap / sigma.
    discounts = [witness(FAR, ONES, 0.1, sa, alpha=1 / 8, centre=5.0)["discount"] for sa in (2e-3, 2e-4, 2e-5, 2e-6, 2e-7)]
    assert discounts == pytest.approx([0.1284, 0.2250, 0.4000, 0.7113, 1.0], abs=2e-4)


def test_the_cap_for_a_reading_on_an_axis_misfires_at_a_large_noise_ratio():
    # two agreeing readings, forces 5.5 and 5.6 each +- 0.1: at sigma_a = 0.002 the law keeps
    # 5.41 to 5.69, as the Gaussian does. With the acceleration known to 2e-6 the readings lie
    # close to an axis in noise units, the cap K' = min(K, 1 + min(B11, B22)^2) engages, and the
    # same law keeps 0.2 to 15.6. The theorem's coverage statement is a limit at a fixed noise
    # ratio; this is what it does not say. E12's own code gives the same log factor there.
    usual = witness(NEAR, ONES, 0.1, 2e-3, alpha=1 / 8, centre=5.0)
    assert usual["interval"] == pytest.approx((5.410, 5.689), abs=1e-3) and usual["precision"] < 1 + 6161.0 ** 2
    precise = witness(NEAR, ONES, 0.1, 2e-6, alpha=1 / 8, centre=5.0)
    assert precise["precision"] == pytest.approx(1 + 6161.0 ** 2) and precise["discount"] == 1.0
    assert precise["interval"][0] < 0.25 and precise["interval"][1] > 15.0


def test_what_normalising_erases_is_a_constant_the_curve_keeps():
    # two readings agreeing on mass 1, against readings of 4 and 1/4 (unit noise). Their curves
    # differ by exactly 9 at every mass, so their normalised Gaussian evidence is identical and no
    # function of it can tell them apart; T73's curve does, through its level (T73 (5), T77 (4)).
    agree = Curve([isotropic([2.0], [2.0], 1.0, 1.0)] * 2)
    conflict = Curve([isotropic([4.0], [1.0], 1.0, 1.0), isotropic([1.0], [4.0], 1.0, 1.0)])
    masses = np.array([0.0, 0.25, 0.7, 1.0, 4.0, math.inf])
    assert conflict(masses) - agree(masses) == pytest.approx(np.full(6, 9.0), abs=1e-12)
    assert agree.least()[0] == pytest.approx(0.0, abs=1e-12) and conflict.least()[0] == pytest.approx(9.0)


def test_an_interval_that_widens_with_the_gap_covers_more_than_its_level():
    # two readings at strong signal: the centre's error Z and the gap G are independent standard
    # normals. Any multiplier g(G) >= 1 on the half-width gives coverage E[2 Phi(z g(G)) - 1] > p:
    # widening by half when |G| > 2 gives 0.9521 at p = 0.95. So "exactly the stated level on
    # identifying series" and "wider when the readings disagree" cannot both be asked.
    z = stats.norm.ppf(0.975)
    widened = lambda g: 2 * stats.norm.cdf(z * (1.5 if abs(g) > 2 else 1.0)) - 1
    pieces = [integrate.quad(lambda g: widened(g) * stats.norm.pdf(g), a, b)[0] for a, b in ((-12, -2), (-2, 2), (2, 12))]
    expected = 0.95 + 2 * stats.norm.cdf(-2) * (2 * stats.norm.cdf(1.5 * z) - 1 - 0.95)
    assert sum(pieces) == pytest.approx(expected, abs=1e-9) and expected == pytest.approx(0.95212, abs=1e-5)


# ---- properties 3 and 4: coverage ------------------------------------------------------------

def _covers(rng, signals, reps, alpha=1 / 8, angle=1.0):
    signals = np.asarray(signals, dtype=float)
    truth, law, report = math.tan(angle), 0, 0
    for _ in range(reps):
        y = signals[:, None] * np.array([math.sin(angle), math.cos(angle)])[None, :] + rng.normal(size=(signals.size, 2))
        got = witness(y[:, 0], y[:, 1], 1.0, 1.0, alpha=alpha)
        law += got["interval"][0] <= truth <= got["interval"][1]
        report += got["report"][0] <= truth <= got["report"][1]
    return law / reps, report / reps


def test_coverage_on_two_identifying_series_short_run():
    rng = np.random.default_rng(78)
    for signals in ([30.0] + [0.0] * 100, [2.0] * 60):
        law, report = _covers(rng, signals, 150)
        assert 0.90 < law < 0.995 and report >= law


@SLOW
def test_coverage_on_identifying_series_and_below():
    # measured 22 Nov, 300 series each, uniform-angle prior, true angle 1.0 rad; law's interval
    # (alpha 1/8; 1/16) and the hull report: 200 at signal 2: 0.967; 0.977; 0.970. 2000 at 0.5:
    # 0.927; 0.963; 0.957. One at 30 among 100 noise-only: 0.947; 0.940; 0.953. Two at 50: 0.937;
    # 0.943; 0.950. Not identifying: 20 at 1: 0.950; 0.963; 0.990. 20 at 0.5: 0.957; 0.973; 0.983.
    # 20 noise-only: 0.940; 0.973; 0.997.
    rng = np.random.default_rng(78)
    for signals in ([2.0] * 200, [0.5] * 2000, [30.0] + [0.0] * 100, [50.0, 50.0], [1.0] * 20, [0.5] * 20, [0.0] * 20):
        law, _ = _covers(rng, signals, 300)
        law16, _ = _covers(rng, signals, 300, alpha=1 / 16)
        _, report = _covers(rng, signals, 300)
        assert 0.90 < law < 0.995 and 0.90 < law16 < 0.995 and report >= 0.93


def test_final_factors_that_stay_bounded_on_conflicts_cannot_also_multiply():
    # if final log factors added over disjoint lists, then for a reading x and a reading y that
    # conflicts with it, osc(l_x) <= osc(l_xxy) + osc(l_xy): bounded confidence on the two
    # conflicting lists would bound the single good reading's factor too, and it could never
    # concentrate (E12, product obstruction). The witness keeps both ends by not multiplying: it
    # recomputes its pooled discount for the joined list.
    from met.scalar_witness import log_factor
    angles = np.linspace(0.0, 0.5 * math.pi, 2001)
    osc = lambda rows, t: float(np.ptp(log_factor(np.asarray(rows) * t, angles, 1 / 8)))
    x, y = [2.0, 1.0], [1.0, 3.0]
    for t in (10.0, 100.0, 1000.0):
        bounded = osc([x, x, y], t) + osc([x, y], t)
        assert bounded < 60.0                                   # the two conflicting lists stay bounded
        assert osc([x], t) > 0.3 * t * t                        # the good reading alone concentrates
    assert osc([x], 1000.0) > 1000 * (osc([x, x, y], 1000.0) + osc([x, y], 1000.0))

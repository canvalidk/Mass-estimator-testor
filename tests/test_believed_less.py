"""How much less a result is believed: Can's two scenes (R-2026-11-17-3; problems/O14.md, 17 Nov).

Scene A: two readings of one mass at 5.5 and 5.6. Scene B: at 5 and 6. Each force is known to
sigma = 0.1 and each acceleration is 1 (to 0.002), so a reading's mass is known to about 0.1.
The sigma is believed as given. B must come out less sure than A; the amount is the open item.
These tests record what each mechanism on the pages does with the two scenes.
"""
import math

import numpy as np
import pytest

from met.full import Reading, ours30
from met.ours33 import ours33

A, B = (5.5, 5.6), (5.0, 6.0)


def _pair(scene, sigma=0.1):
    cov = np.diag([sigma ** 2, 0.002 ** 2])
    return [Reading(np.array([f]), np.array([1.0]), cov) for f in scene]


def _half(interval):
    return 0.5 * (interval[1] - interval[0])


def test_the_exact_gaussian_evidence_is_as_wide_in_the_far_scene_as_in_the_near_one():
    # the relative probability of the pair under each candidate mass falls off identically in A and
    # in B: the pair's improbability multiplies every candidate alike. So the Gaussian alone gives
    # no discount, whatever the gap.
    a, b = (ours33(_pair(s), 5.0, math.inf, calibrate=False) for s in (A, B))
    assert _half(a["raw_interval"]) == pytest.approx(1.96 * 0.1 / math.sqrt(2), rel=0.01)
    assert _half(b["raw_interval"]) == pytest.approx(_half(a["raw_interval"]), rel=2e-3)


def test_far_pairs_are_as_reliable_as_near_ones_under_the_exact_gaussian_and_only_there():
    # two readings of one mass, kept only when they land more than 4 sigma apart. How often is
    # "the average, with the combined uncertainty" right? Under the exact Gaussian as often as ever
    # (the average is independent of the gap); with a noise size that varies from reading to reading
    # (nu = 10) much less often. So the Gaussian itself calls for no discount: what forces one is
    # another premise.
    rng = np.random.default_rng(1117)
    n, sigma = 4_000_000, 0.1
    half = 1.96 * sigma / math.sqrt(2)

    def right_among_far(errors):
        x1, x2 = errors(), errors()                       # the true mass is 0 without loss
        far = np.abs(x1 - x2) > 4 * sigma
        return float(np.mean(np.abs((x1 + x2)[far] / 2) < half))

    assert right_among_far(lambda: rng.normal(0, sigma, n)) == pytest.approx(0.95, abs=0.006)
    assert right_among_far(lambda: rng.normal(0, sigma, n) / np.sqrt(rng.chisquare(10, n) / 10)) < 0.74


def test_g_is_blind_to_readings_that_disagree_with_each_other():
    # each reading alone fits a positive mass exactly, so its own misfit is zero and g = 1
    for scene in (A, B):
        assert ours30(_pair(scene), 5.0)["g"] == pytest.approx([1.0, 1.0])


def test_eta_is_what_widens_the_far_scene():
    for run in (lambda rd: ours30(rd, 5.0), lambda rd: ours33(rd, 5.0, math.inf)):
        a, b = run(_pair(A)), run(_pair(B))
        assert a["eta"] == pytest.approx(1.0)                      # 0.5 sigma off each: no widening
        assert b["eta"] == pytest.approx(1 / 25, rel=0.02)         # 5 sigma off each: 1 / 5**2
        assert _half(a["interval"]) == pytest.approx(0.139, abs=0.003)
        assert _half(b["interval"]) == pytest.approx(0.70, abs=0.02)


def test_the_gram_eta_does_not_widen_the_far_scene():
    # eta_G = 1 - Q/H (T70, B22; T43's spectral form on this slice) is the weak-signal correction
    # alone: at strong signal it is 1 whatever the gap, so B comes out as narrow as A
    a, b = (ours30(_pair(s), 5.0, eta_form="gram") for s in (A, B))
    assert a["eta"] > 0.9998 and b["eta"] > 0.9998
    assert _half(b["interval"]) == pytest.approx(_half(a["interval"]), rel=2e-3)


@pytest.mark.parametrize("sigma", [0.25, 0.1, 0.05])
def test_past_two_sigma_the_declared_sigma_no_longer_sets_the_quoted_width(sigma):
    # eta = min(1, (2 sigma / gap)**2), so past a gap of 2 sigma the quoted half-width is about
    # 0.7 * gap whatever sigma was declared
    b = ours33(_pair(B, sigma), 5.0, math.inf)
    assert _half(b["interval"]) == pytest.approx(0.707, abs=0.012)
    per_reading = math.hypot(sigma, 0.002 * 5.5)                 # the acceleration's share of a reading's error
    assert _half(b["raw_interval"]) == pytest.approx(1.96 * per_reading / math.sqrt(2), rel=0.02)


def test_a_handed_over_nu_widens_the_far_scene_with_no_eta():
    # the case study: the discount's parameter is given (nu_eff = 10) and eta is switched off
    a, b = (ours33(_pair(s), 5.0, 10.0, calibrate=False) for s in (A, B))
    assert _half(a["interval"]) == pytest.approx(0.15, abs=0.01)
    assert _half(b["interval"]) == pytest.approx(0.55, abs=0.02)


def test_nu_and_eta_together_discount_the_same_gap_twice():
    # ours33 as it stands: the handed-over nu widens the raw law for the gap, and eta then widens the
    # report for the same gap; the warp's ranks leave the data and land in pi's tails
    b = ours33(_pair(B), 5.0, 30.0)
    assert 5.1 < b["raw_interval"][0] < b["raw_interval"][1] < 5.9
    assert math.log10(b["interval"][1] / b["interval"][0]) > 70


def test_bounded_confidence_between_readings_as_the_declared_sigma_shrinks():
    # T30's requirement, read for the disagreement between readings: the two readings stay at 5 and 6
    # while the declared sigma shrinks. The Gaussian and the Gram eta become certain of 5.5 (width
    # proportional to sigma); the squared-score eta and a handed-over nu keep a width set by the gap,
    # and stay informative (the interval does not run off toward pi's).
    widths = {name: [] for name in ("gaussian", "gram", "eta30", "eta33", "nu")}
    for sigma in (0.2, 0.05):
        rd = _pair(B, sigma)
        runs = {"gaussian": ours33(rd, 5.0, math.inf, calibrate=False), "gram": ours30(rd, 5.0, eta_form="gram"),
                "eta30": ours30(rd, 5.0), "eta33": ours33(rd, 5.0, math.inf),
                "nu": ours33(rd, 5.0, 10.0, calibrate=False)}
        for name, out in runs.items():
            widths[name].append(_half(out["interval"]))
            if name in ("eta30", "eta33", "nu"):
                assert 4.5 < out["interval"][0] < out["interval"][1] < 6.5
    for name in ("gaussian", "gram"):
        assert widths[name][0] / widths[name][1] > 3.5              # sigma fell by 4; so did the width
    for name, lo, hi in (("eta30", 0.68, 0.71), ("eta33", 0.68, 0.73), ("nu", 0.50, 0.57)):
        assert all(lo < w < hi for w in widths[name])               # the gap is 1: the width tracks it


def test_a_shared_weight_and_per_reading_weights_part_ways_on_three_readings():
    # readings at 5, 5.05 and 6: one weight for the series widens the answer and keeps the point;
    # weights per reading (here from a handed-over nu) mostly discount the 6 and move the point.
    # Which the definition does is a fork owed a sort (problems/O14.md, 17 Nov).
    cov = np.diag([0.1 ** 2, 0.002 ** 2])
    rd = [Reading(np.array([f]), np.array([1.0]), cov) for f in (5.0, 5.05, 6.0)]
    plain = ours33(rd, 5.0, math.inf, calibrate=False)
    shared = ours30(rd, 5.0)
    assigned = ours33(rd, 5.0, 10.0, calibrate=False)
    assert plain["point"] == pytest.approx(5.35, abs=0.01)
    assert shared["point"] == pytest.approx(plain["point"], abs=0.03) and _half(shared["interval"]) > 0.45
    assert assigned["point"] == pytest.approx(5.08, abs=0.02) and assigned["interval"][1] < 5.3


@pytest.mark.xfail(strict=True, raises=ValueError,
                   reason="ours33 raises on a measured input (L2): at sigma = 0.02 the pair is 50 sigma apart, "
                          "eta = 0.0016, and the warp's rank Phi(-1.96/sqrt(eta)) underflows to 0. ours30 returns "
                          "[4.82, 6.21] on the same input. The numerical face of the break at the report (O18).")
def test_ours33_returns_a_width_on_a_pair_fifty_sigma_apart():
    b = ours33(_pair(B, 0.02), 5.0, math.inf)
    assert _half(b["interval"]) == pytest.approx(0.70, abs=0.03)

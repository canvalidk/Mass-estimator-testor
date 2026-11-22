"""The evidence curve with no starting law (T73, T74; problems/O20.md). Harvested 22 Nov from
explore/e13-unnormalized-evidence (insights 1, 3, 4) and checked here against met.full, against
closed forms, and on the scenes. Nothing in this file is a law of mass."""
import math
import os

import numpy as np
import pytest

from met.conditional_score import gaussian_score_test
from met.evidence import Curve, common_force, leftover, lump
from met.full import Reading, _h_terms, isotropic


def _any_reading(rng, d, k):
    a = rng.normal(size=(2 * d, 2 * d))
    cov = a @ a.T + 0.3 * np.eye(2 * d)
    q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    return Reading(2 * rng.normal(size=d), rng.normal(size=d), 0.5 * (cov + cov.T), subspace=q[:, :k])


def _pair(scene, sigma=0.1):
    cov = np.diag([sigma ** 2, 0.002 ** 2])
    return [Reading(np.array([f]), np.array([1.0]), cov) for f in scene]


@pytest.mark.parametrize("d,k", [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])
def test_the_curve_is_the_leftover_met_full_already_computes(d, k):
    # K = E - H of T36, any covariance, a declared subspace: an independent route to the same number
    rng = np.random.default_rng(100 * d + k)
    readings = [_any_reading(rng, d, k) for _ in range(4)]
    c, u = 1.7, np.linspace(0.01, 1.55, 9)
    reference = sum(r.pre()["energy"] - _h_terms(r.pre(), c, u) for r in readings)
    assert Curve(readings, scale=c).at_angle(u) == pytest.approx(reference, abs=1e-10)


def test_on_one_instrument_the_least_is_the_smaller_eigenvalue_of_the_matrix_of_dot_products():
    # total least squares through the origin; the 2x2 matrix of claude/vector_division_common_line
    # (22 Nov): its smaller eigenvalue is the leftover, its larger eigenvector the mass
    rng = np.random.default_rng(3)
    sf, sa = 0.3, 0.7
    a = rng.normal(size=(5, 3))
    f = 2.0 * a + 0.4 * rng.normal(size=(5, 3))
    readings = [isotropic(x, y, sf, sa) for x, y in zip(f, a)]
    gram = np.array([[np.sum(a * a) / sa ** 2, np.sum(a * f) / (sa * sf)],
                     [np.sum(a * f) / (sa * sf), np.sum(f * f) / sf ** 2]])
    values, vectors = np.linalg.eigh(gram)
    found = leftover(readings)
    assert Curve(readings).scale == pytest.approx(sf / sa)
    assert found["least"] == pytest.approx(values[0], rel=1e-9)
    assert found["mass"] == pytest.approx(vectors[1, 1] / vectors[0, 1] * sf / sa, rel=1e-6)


def test_the_all_zero_experiment_is_the_flat_curve_and_a_zero_reading_adds_exactly_nothing():
    # insight 1: with no starting law, 0/0 is the identity of combination, K = 0 at every mass
    zero = isotropic([0, 0, 0], [0, 0, 0], 2.0, 1.0)
    found = leftover([zero])
    assert found["flat"] and found["least"] == 0.0 and found["largest"] == 0.0 and found["mass"] is None
    rng = np.random.default_rng(8)
    series = [_any_reading(rng, 3, 2) for _ in range(3)]
    masses = np.array([0.0, 0.2, 1.0, 7.0, math.inf])
    scale = Curve(series).scale
    assert np.array_equal(Curve(series + [zero], scale=scale)(masses), Curve(series, scale=scale)(masses))


def test_both_ends_of_the_positive_masses_stay_alive_at_every_finite_sigma():
    # one ordinary reading, F = 10 +- 1, a = 1 +- 1: K(0) = 100, K(infinity) = 1. The curve never
    # rules out a large mass, so its set at the chi-square threshold is [3.28, infinity)
    # (E13 writes this scene as "(0, inf)"; it is unbounded above only).
    curve = Curve([isotropic([10.0], [1.0], 1.0, 1.0)])
    assert curve(0.0) == pytest.approx(100.0) and curve(math.inf) == pytest.approx(1.0)
    assert curve.least()[0] == pytest.approx(0.0, abs=1e-12)
    low = (-20 + math.sqrt(400 + 4 * 2.8415 * 96.1585)) / (2 * 2.8415)
    assert low == pytest.approx(3.28, abs=0.01)
    assert curve(low - 0.01) > 3.8415 > curve(low + 0.01) and curve(1e9) < 3.8415


def test_a_zero_reading_and_an_infinitely_noisy_reading_are_one_reading_in_noise_units():
    # insight 3. A reading whose sigma grows without bound goes to the zero reading in its own
    # noise metric. The curve and T45's score follow; the count beside the leftover does not.
    rng = np.random.default_rng(11)
    series = [_any_reading(rng, 2, 2) for _ in range(3)]
    extra = _any_reading(rng, 2, 2)
    zero = Reading(np.zeros(2), np.zeros(2), extra.cov)
    noisy = Reading(extra.force, extra.acceleration, 1e16 * extra.cov)
    masses = np.array([0.3, 1.0, 4.0])
    scale = Curve(series).scale
    with_zero, with_noisy = Curve(series + [zero], scale)(masses), Curve(series + [noisy], scale)(masses)
    assert with_noisy == pytest.approx(with_zero, abs=1e-12)
    for m in masses:
        a, b = gaussian_score_test(series + [zero], m), gaussian_score_test(series + [noisy], m)
        assert b["statistic"] == pytest.approx(a["statistic"], rel=1e-9)
        assert a["statistic"] == pytest.approx(gaussian_score_test(series, m)["statistic"], rel=1e-12)
    # the least is zero-safe; the probability bound beside it counts the zero reading's components
    plain, padded = leftover(series, scale), leftover(series + [zero], scale)
    assert padded["least"] == plain["least"] and padded["count"] == plain["count"] + 2
    assert padded["bound"] > plain["bound"]


def test_the_two_scenes_differ_in_the_leftover_and_not_in_the_curve_around_it():
    # Can's scenes (R-2026-11-17-3). The far pair's leftover is 49.4 against 0.49, and the curve's
    # rise from its least is the same in both: the pair's improbability does not depend on m.
    near, far = leftover(_pair((5.5, 5.6))), leftover(_pair((5.0, 6.0)))
    assert near["least"] == pytest.approx(0.494, abs=1e-3) and far["least"] == pytest.approx(49.40, abs=1e-2)
    assert near["bound"] == pytest.approx(0.781, abs=1e-3) and far["bound"] == pytest.approx(1.87e-11, rel=0.01)
    a, b = Curve(_pair((5.5, 5.6))), Curve(_pair((5.0, 6.0)))
    for step in (0.05, 0.1386, 0.3):
        rise_a = a(near["mass"] + step) - near["least"]
        rise_b = b(far["mass"] + step) - far["least"]
        assert rise_b == pytest.approx(rise_a, rel=5e-3)
    # two readings, each known to sqrt(0.1^2 + (5.55 * 0.002)^2): at 1.96 combined sigmas the rise is 3.84
    each = math.hypot(0.1, 5.55 * 0.002)
    assert a(near["mass"] + 1.96 * each / math.sqrt(2)) - near["least"] == pytest.approx(3.8416, rel=0.01)


# ---- composition of curves (T74) ------------------------------------------------------------

def _quadratic(centre, sd):
    return lambda m: ((np.asarray(m, dtype=float) - centre) / sd) ** 2


def test_lumping_two_quadratic_curves_adds_the_centres_and_the_variances():
    totals = np.array([4.0, 4.5, 5.0, 5.7, 6.5])
    got = lump(_quadratic(2.0, 0.3), _quadratic(3.0, 0.4), totals)
    assert got == pytest.approx((totals - 5.0) ** 2 / (0.3 ** 2 + 0.4 ** 2), abs=1e-9)


def test_a_body_with_no_evidence_beside_one_at_5_makes_the_total_one_sided():
    # the unmeasured partner's curve is flat, so every total above 5 is fully compatible, and a
    # total below 5 costs what the measured body pays to be that light: [4.02, infinity)
    flat = lambda m: np.zeros_like(np.asarray(m, dtype=float))
    totals = np.array([3.5, 4.02, 4.5, 5.0, 8.0, 50.0])
    got = lump(_quadratic(5.0, 0.5), flat, totals)
    assert got == pytest.approx(np.where(totals < 5, ((totals - 5) / 0.5) ** 2, 0.0), abs=1e-9)
    assert lump(_quadratic(5.0, 0.5), flat, [5 - 1.96 * 0.5])[0] == pytest.approx(3.8416, rel=1e-6)
    assert lump(_quadratic(5.0, 0.5), flat, [math.inf])[0] == pytest.approx(0.0, abs=1e-9)


def _body(rng, m, n, r, sf=0.2, sa=0.1):
    out = []
    for _ in range(n):
        a = rng.normal(size=3)
        a *= r * sa / np.linalg.norm(a)
        out.append(isotropic(m * a + sf * rng.normal(size=3), a + sa * rng.normal(size=3), sf, sa))
    return out


def _swap(r):
    d = r.dimension
    idx = np.r_[d:2 * d, 0:d]
    return Reading(r.acceleration, r.force, r.cov[np.ix_(idx, idx)], subspace=r.subspace)


def test_on_real_series_the_composite_curve_bottoms_where_the_parts_do_and_swaps_exactly():
    rng = np.random.default_rng(5)
    one, two = _body(rng, 2.0, 6, 40.0), _body(rng, 3.0, 6, 40.0)
    c1, c2 = Curve(one), Curve(two)
    totals = np.linspace(4.5, 5.5, 201)
    joined = lump(c1, c2, totals)
    assert joined.min() == pytest.approx(c1.least()[0] + c2.least()[0], rel=1e-3)
    assert abs(totals[np.argmin(joined)] - (leftover(one)["mass"] + leftover(two)["mass"])) < 0.006
    # the common-force composition is the lumping of the swapped readings on inverse masses
    s1, s2 = Curve([_swap(r) for r in one]), Curve([_swap(r) for r in two])
    reduced = np.array([1.0, 1.2, 1.3])
    assert common_force(c1, c2, reduced) == pytest.approx(lump(s1, s2, 1 / reduced), rel=1e-9)


def test_as_the_sigmas_shrink_the_composite_curve_admits_only_the_sum():
    one = Curve([isotropic([2.0], [1.0], 1e-6, 1e-6)])
    two = Curve([isotropic([3.0], [1.0], 1e-6, 1e-6)])
    got = lump(one, two, [4.99999, 5.0, 5.00001])
    assert got[1] == pytest.approx(0.0, abs=1e-6) and got[0] > 6 and got[2] > 6


@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="300 composites per signal size")
@pytest.mark.parametrize("signal,low,high", [(1.0, 0.78, 0.90), (6.0, 0.91, 0.975)])
def test_the_raw_composite_set_under_covers_at_weak_signal(signal, low, high):
    # the composite curve cut at the chi-square threshold inherits the raw curve's weak-signal
    # shortfall (measured 22 Nov: 0.840, 0.910, 0.943 at signal 1, 3, 6; E13 reports 89.3, 90.7, 94.0
    # on 150). E13's conditional calibration of the composite (insight 7) is not harvested.
    rng = np.random.default_rng(7)

    def body(m, n):
        out = []
        for _ in range(n):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            th = math.atan(m)
            out.append(isotropic(signal * math.sin(th) * d + rng.normal(size=3),
                                 signal * math.cos(th) * d + rng.normal(size=3), 1.0, 1.0))
        return out
    hit, reps = 0, 300
    for _ in range(reps):
        c1, c2 = Curve(body(2.0, 5)), Curve(body(3.0, 5))
        hit += lump(c1, c2, [5.0])[0] - c1.least()[0] - c2.least()[0] <= 3.841
    assert low < hit / reps < high

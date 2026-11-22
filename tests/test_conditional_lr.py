"""Exact conditional tests of a candidate mass with no starting law (T75; problems/O20.md).
Harvested 22 Nov from explore/e13-unnormalized-evidence (insights 2, 6, 8) and checked here on an
independent implementation and an independent world. Coverage figures are seeded; the long runs
need MET_SLOW=1. Nothing in this file is a law of mass."""
import math
import os

import numpy as np
import pytest
from scipy import stats

from met.conditional_lr import (conditional_test, drowning_count, floor_test, lump_kept, lump_test,
                                score_test, weak_signal_efficiency, zero_safe)
from met.evidence import Curve, leftover, lump
from met.full import Reading, isotropic

SLOW = pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="hundreds of simulated series")


def _p(readings, mass, **options):
    return conditional_test(readings, mass, **options)["p_value"]


def _scalar(force, sigma=0.1):
    return isotropic([force], [1.0], sigma, 0.002)


def _vector(force, d, sigma=0.1):
    return isotropic([force] + [0.0] * (d - 1), [1.0] + [0.0] * (d - 1), sigma, 0.002)


# ---- the two scenes -------------------------------------------------------------------------

@pytest.mark.parametrize("scene", [(5.5, 5.6), (5.0, 6.0)])
def test_under_the_declared_gaussian_the_far_pair_is_kept_as_narrowly_as_the_near_pair(scene):
    # insight 2 (Kac-Bernstein): the average is independent of the gap, so the exact Gaussian
    # test keeps the same width around either pair: about 1.96 * 0.1 / sqrt(2) = 0.139
    readings = [_scalar(f) for f in scene]
    centre = leftover(readings)["mass"]
    for weight in (None, zero_safe):
        assert _p(readings, centre - 0.10, weight=weight, draws=999) > 0.05
        assert _p(readings, centre + 0.10, weight=weight, draws=999) > 0.05
        assert _p(readings, centre - 0.19, weight=weight, draws=999) < 0.05
        assert _p(readings, centre + 0.19, weight=weight, draws=999) < 0.05


def test_two_scalar_readings_cannot_be_judged_by_direction_alone():
    # two signs have four patterns, so no candidate reaches 5%: near or far, everything is kept
    for scene in ((5.5, 5.6), (5.0, 6.0)):
        readings = [_scalar(f) for f in scene]
        for mass in (0.01, 3.0, 5.5, 9.0, 500.0):
            assert _p(readings, mass, null="spherical") >= 0.25
            assert floor_test(readings, mass, draws=199)["p_value"] >= 0.25


@pytest.mark.parametrize("d", [2, 3])
def test_on_vector_readings_the_far_pair_is_kept_more_widely_with_no_sigma_changed(d):
    # insight 8: the larger of the two p-values (declared Gaussian with the zero-safe weight; any
    # size of noise per reading). Near pair: the Gaussian's width. Far pair: everything between
    # the two readings, 5.5 +- about 0.47, where the Gaussian keeps 5.5 +- 0.14. No sigma is
    # revised and no number is handed over.
    near, far = [_vector(5.5, d), _vector(5.6, d)], [_vector(5.0, d), _vector(6.0, d)]
    keep = lambda readings, m: floor_test(readings, m, draws=399)["p_value"] > 0.05
    assert keep(near, 5.45) and keep(near, 5.65) and not keep(near, 5.36) and not keep(near, 5.74)
    assert keep(far, 5.1) and keep(far, 5.9) and not keep(far, 4.9) and not keep(far, 6.1)
    assert _p(far, 5.1, weight=zero_safe, draws=399) < 0.05 and _p(far, 5.9, weight=zero_safe, draws=399) < 0.05


def test_readings_that_agree_too_well_shrink_the_scale_free_set_below_the_sigma():
    # the cost of judging by direction alone, and why the two tests are joined: ten readings
    # within 0.01 of each other at sigma = 0.1 leave a scale-free set narrower than 0.02, where
    # the declared sigma allows 5.5 +- 0.06; one aligned vector reading leaves a single point
    tight = [_scalar(5.5 + 0.002 * (i - 4.5)) for i in range(10)]
    assert _p(tight, 5.52, null="spherical") < 0.05 and _p(tight, 5.48, null="spherical") < 0.05
    assert floor_test(tight, 5.52, draws=399)["p_value"] > 0.05
    single = [_vector(5.0, 2)]
    assert _p(single, 5.02, null="spherical", draws=399) < 0.05
    assert floor_test(single, 5.1, draws=399)["p_value"] > 0.05


@pytest.mark.parametrize("n,smallest", [(4, 2 / 16), (5, 2 / 32), (6, 2 / 64)])
def test_scalar_readings_on_one_instrument_need_six_before_direction_alone_can_reject(n, smallest):
    # each sign pattern ties with its mirror image about the candidate line, so the smallest
    # p-value is 2 / 2^n: 0.0625 at five readings, 0.031 at six. (E13's insight 8 says five; its
    # fixed grid breaks the tie numerically.)
    readings = [isotropic([2.0 + 0.001 * i], [1.0], 0.01, 0.00001) for i in range(n)]
    for mass in (1.5, 1.99, 2.05, 2.5):
        assert _p(readings, mass, null="spherical") == pytest.approx(smallest)


# ---- zero readings, right angles, opposed readings ------------------------------------------

def test_a_zero_reading_changes_nothing_under_a_zero_safe_weight_or_the_scale_free_null():
    strong = [isotropic([4.0, 0.2, -0.1], [2.0, 0.0, 0.1], 0.2, 0.1)]
    zeros = [isotropic([0, 0, 0], [0, 0, 0], 0.2, 0.1)] * 30
    for mass in (1.7, 2.0, 2.4):
        assert _p(strong + zeros, mass, weight=zero_safe, draws=499) == _p(strong, mass, weight=zero_safe, draws=499)
        assert _p(strong + zeros, mass, null="spherical", draws=499) == _p(strong, mass, null="spherical", draws=499)


def test_zero_readings_move_the_unweighted_gaussian_test():
    # its conditional law gives every reading k random leftover components, a zero reading too:
    # 150 zero readings take a candidate from rejected to kept (measured 22 Nov: 0.04 to 0.12; 1000 take it to 1)
    strong = [isotropic([4.0, 0.2, -0.1], [2.0, 0.0, 0.1], 0.2, 0.1)]
    zeros = [isotropic([0, 0, 0], [0, 0, 0], 0.2, 0.1)] * 150
    assert _p(strong, 1.7, draws=399) < 0.05 < 0.1 < _p(strong + zeros, 1.7, draws=399)


def test_readings_at_right_angles_keep_every_mass():
    # scene C (proof.md G29): one reading says force only, the other acceleration only. On one
    # instrument the curve is flat, so the unweighted test keeps everything with p = 1
    scene = [isotropic([1.0], [0.0], 0.05, 0.05), isotropic([0.0], [1.0], 0.05, 0.05)]
    assert leftover(scene)["flat"] and leftover(scene)["least"] == pytest.approx(400.0)
    for mass in (0.001, 0.5, 1.0, 20.0, 1000.0):
        assert _p(scene, mass, draws=199) == 1.0
        assert _p(scene, mass, weight=zero_safe, draws=199) > 0.05


def test_opposed_readings_are_rejected_on_the_whole_line_and_pushed_to_the_ends_on_positive_masses():
    pair = [isotropic([2.0], [-1.0], 0.1, 0.1), isotropic([1.0], [-2.0], 0.1, 0.1)]
    for mass in (0.001, 0.05, 1.0, 20.0, 1000.0):
        assert _p(pair, mass, domain="line", draws=399) < 0.05          # no positive mass is kept
    assert _p(pair, 0.001, draws=399) > 0.05 and _p(pair, 1000.0, draws=399) > 0.05
    assert _p(pair, 1.0, draws=399) < 0.05


def test_a_zero_safe_weight_keeps_the_one_mass_a_single_opposed_reading_fits_worst():
    # at m = 0.5 the reading (F, a) = (2, -1) has nothing along the candidate line: its weight is
    # zero and the test is empty, p = 1, although the leftover there is 500. E13's "pinhole";
    # the same blindness at right angles as problems/O14.md, 18 Nov night, item 2.
    single = [isotropic([2.0], [-1.0], 0.1, 0.1)]
    blind = conditional_test(single, 0.5, weight=zero_safe, draws=199)
    assert blind["empty"] and blind["p_value"] == 1.0 and blind["weights"] == [0.0]
    seen = conditional_test(single, 0.5, draws=399)
    assert seen["at_candidate"] == pytest.approx(500.0) and seen["p_value"] < 0.05


# ---- coverage -------------------------------------------------------------------------------

def _world(rng, n, m, r, d=3, k=2, rho=0.6, nu=None):
    """Correlated channels, a different noise ratio per reading, a declared 2D subspace of 3D,
    pushes of length r in each reading's own noise units; nu: an unknown noise size per reading."""
    out = []
    for _ in range(n):
        q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        sub = q[:, :k]
        sf, sa = math.exp(rng.uniform(-1, 1)), math.exp(rng.uniform(-1, 1))
        base = np.block([[sf * sf * np.eye(d), rho * sf * sa * np.eye(d)],
                         [rho * sf * sa * np.eye(d), sa * sa * np.eye(d)]])
        mix = np.eye(2 * d) + 0.2 * rng.normal(size=(2 * d, 2 * d))
        cov = mix @ base @ mix.T
        cov = 0.5 * (cov + cov.T)
        direction = sub @ rng.normal(size=k)
        mean = np.concatenate([m * direction, direction])
        chol = np.linalg.cholesky(cov)
        mean = mean * (r / np.linalg.norm(np.linalg.solve(chol, mean)))
        size = 1.0 if nu is None else nu / rng.chisquare(nu)
        y = mean + math.sqrt(size) * (chol @ rng.normal(size=2 * d))
        out.append(Reading(y[:d], y[d:], cov, subspace=sub))
    return out


def _coverage(reps, n, m, r, nu=None, seed=1, draws=199):
    rng = np.random.default_rng(seed)
    hit = dict(gaussian=0, zero_safe=0, spherical=0, either=0, raw=0)
    for rep in range(reps):
        series = _world(rng, n, m, r, nu=nu)
        a = conditional_test(series, m, draws=draws, seed=rep)
        b = _p(series, m, draws=draws, seed=rep + 10 ** 6, weight=zero_safe)
        c = _p(series, m, draws=draws, seed=rep + 2 * 10 ** 6, null="spherical")
        hit["gaussian"] += a["p_value"] > 0.05
        hit["raw"] += a["statistic"] <= 3.841
        hit["zero_safe"] += b > 0.05
        hit["spherical"] += c > 0.05
        hit["either"] += max(b, c) > 0.05
    return {name: count / reps for name, count in hit.items()}


def test_on_the_full_type_the_conditional_tests_cover_where_the_raw_curve_does_not():
    # a short run for the default suite: weak signal, correlated, mixed ratios, declared subspace
    got = _coverage(100, 10, 2.5, 1.0, draws=99)
    assert got["gaussian"] > 0.89 and got["zero_safe"] > 0.89 and got["spherical"] > 0.89
    assert got["raw"] < 0.80


@SLOW
@pytest.mark.parametrize("r,raw_below", [(0.0, 0.56), (1.0, 0.72), (5.0, 0.97)])
def test_coverage_on_the_declared_gaussian_world(r, raw_below):
    # measured 22 Nov, 400 series each: gaussian 0.955, 0.963, 0.948; zero-safe 0.953, 0.960,
    # 0.940; spherical 0.948, 0.945, 0.945; either 0.995, 0.990, 0.970; raw curve at the
    # chi-square threshold 0.477, 0.645, 0.930 (E13 table B2: 95-97%, raw 44%, 67%, 95%)
    got = _coverage(400, 10, 2.5, r)
    for name in ("gaussian", "zero_safe", "spherical"):
        assert 0.915 < got[name] < 0.985
    assert got["either"] >= max(got["zero_safe"], got["spherical"])
    assert got["raw"] < raw_below


@SLOW
def test_only_the_scale_free_null_covers_when_the_noise_size_varies_from_reading_to_reading():
    # measured 22 Nov, 400 series: nu = 2, r = 3: gaussian 0.427, spherical 0.940; nu = 10,
    # r = 1: gaussian 0.797, spherical 0.958 (E13 table B4: 40.7% and 96.0%; 80.0% and 94.0%)
    heavy = _coverage(400, 15, 2.0, 3.0, nu=2)
    assert heavy["gaussian"] < 0.55 and 0.91 < heavy["spherical"] < 0.98
    mild = _coverage(400, 15, 2.0, 1.0, nu=10)
    assert mild["gaussian"] < 0.87 and 0.92 < mild["spherical"] < 0.99


@SLOW
def test_under_the_declared_gaussian_only_the_joined_test_is_right_at_every_level_of_disagreement():
    # ten scalar readings of one mass, the declared Gaussian world, coverage of the true mass by
    # fifths of the leftover (least disagreement first). Measured 22 Nov on 3000 series:
    #   gaussian  0.950 0.957 0.947 0.942 0.957   (the same at every level: Kac-Bernstein)
    #   spherical 0.848 0.957 0.972 0.985 0.998   (short when the readings agree, long when not)
    #   either    0.950 0.965 0.972 0.985 0.998   (never short; wider than the inputs call for)
    rng = np.random.default_rng(22)
    m, cov, rows = 5.5, np.diag([0.1 ** 2, 0.002 ** 2]), []
    for rep in range(3000):
        series = [Reading(np.array([m + 0.1 * rng.normal()]), np.array([1.0 + 0.002 * rng.normal()]), cov)
                  for _ in range(10)]
        rows.append((leftover(series)["least"], _p(series, m, draws=199, seed=rep, weight=zero_safe),
                     _p(series, m, null="spherical")))
    rows = np.array(rows)
    fifths = np.array_split(np.argsort(rows[:, 0]), 5)
    gaussian = [np.mean(rows[part, 1] > 0.05) for part in fifths]
    spherical = [np.mean(rows[part, 2] > 0.05) for part in fifths]
    either = [np.mean(np.maximum(rows[part, 1], rows[part, 2]) > 0.05) for part in fifths]
    assert all(0.91 < g < 0.985 for g in gaussian)
    assert spherical[0] < 0.90 and spherical[-1] > 0.98
    assert all(e >= g for e, g in zip(either, gaussian)) and either[-1] > 0.98


# ---- two bodies: the composite with a kept promise (T74 (6); E13 insight 7) ------------------

def _body(rng, m, n, r):
    out = []
    for _ in range(n):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        th = math.atan(m)
        out.append(isotropic(r * math.sin(th) * d + rng.normal(size=3), r * math.cos(th) * d + rng.normal(size=3),
                             1.0, 1.0))
    return out


def test_a_measured_body_beside_an_unmeasured_one_gives_the_one_sided_normal_bound():
    # 5 +- 0.5 beside a 0/0 body: the total is tested one-sidedly, p = 1 - Phi((5 - M) / 0.5), so the
    # kept totals are 4.18 upward (5 - 1.645 * 0.5), where the raw curve of T74 (4) gives 4.02 upward
    body = [isotropic([5.0], [1.0], 0.5, 1e-6)]
    nothing = [isotropic([0.0], [0.0], 0.5, 1e-6)]
    for total in (4.0, 4.18, 4.35):
        got = lump_test(body, nothing, total, total * (1 - 1e-9), weight=zero_safe, draws=4999)["p_value"]
        assert got == pytest.approx(stats.norm.sf((5 - total) / 0.5), abs=0.012)
    assert lump_test(body, nothing, 6.0, 5.0, weight=zero_safe, draws=199)["p_value"] == 1.0
    assert lump_test(nothing, nothing, 6.0, 5.0, weight=zero_safe, draws=199)["empty"]


def test_at_strong_signal_a_total_is_kept_as_one_parameter_not_two():
    # the largest p-value over the splits of a total. 2.3 standard deviations below the best total
    # a one-parameter test rejects (p = 0.021) and a two-parameter one keeps (0.071): this rejects.
    # Measured 22 Nov with 41 splits and 1999 draws: 0.070, 0.049, 0.036, 0.022, 0.016 at
    # 1.8, 1.96, 2.1, 2.3, 2.45 sd below; 0.098, 0.074, 0.059, 0.040, 0.030 above.
    rng = np.random.default_rng(3)
    one, two = _body(rng, 2.0, 5, 40.0), _body(rng, 3.0, 5, 40.0)
    grid = np.linspace(4.6, 5.4, 161)
    curve = lump(Curve(one), Curve(two), grid)
    best = grid[np.argmin(curve)]
    inside = grid[curve - curve.min() <= 3.841]
    sd = (inside.max() - inside.min()) / (2 * 1.96)
    assert lump_kept(one, two, best - 1.6 * sd, candidates=5, draws=299)[0]
    assert not lump_kept(one, two, best - 2.3 * sd, candidates=5, draws=299)[0]


@SLOW
@pytest.mark.parametrize("signal", [1.0, 3.0])
def test_at_the_true_split_the_composite_test_covers_where_the_raw_composite_does_not(signal):
    # measured 22 Nov, 150 pairs of bodies (masses 2 and 3, five 3D readings each): kept at the true
    # split 0.973 and 0.947 at signal 1 and 3, the raw composite at 3.84 0.900 and 0.920 on the same
    # series (0.840 and 0.910 on T74 (5)'s 300). E13 table B3: 97.3%, 93.3%; raw 89.3%, 90.7%.
    rng = np.random.default_rng(3)
    _body(rng, 2.0, 5, 40.0), _body(rng, 3.0, 5, 40.0)        # the stream position of the 22 Nov run
    if signal == 3.0:
        for _ in range(150):
            _body(rng, 2.0, 5, 1.0), _body(rng, 3.0, 5, 1.0)
    hit = raw = 0
    for rep in range(150):
        one, two = _body(rng, 2.0, 5, signal), _body(rng, 3.0, 5, signal)
        got = lump_test(one, two, 5.0, 2.0, draws=199, seed=rep)
        hit += got["p_value"] > 0.05
        raw += got["statistic"] <= 3.841
    assert hit / 150 > 0.92 and hit > raw


# ---- the weighted conditional score: one family, two closed forms (T76; E13 insight 9) -------

WEIGHTS = dict(one=lambda h: 1.0 + 0.0 * h, smooth=lambda h: h / (1 + h), cap10=lambda h: h * h / (100 + h * h),
               linear=lambda h: h, square=lambda h: h * h)


def test_the_weighted_score_is_exact_under_every_weight_on_the_full_type():
    rng = np.random.default_rng(76)
    reject = dict.fromkeys(("one", "smooth", "linear"), 0)
    reps = 1500
    for _ in range(reps):
        series = _world(rng, 10, 2.5, 1.0)
        for name in reject:
            reject[name] += score_test(series, 2.5, WEIGHTS[name])["p_value"] < 0.05
    for name, count in reject.items():
        assert 0.033 < count / reps < 0.067, name


def test_a_zero_reading_adds_nothing_to_the_score_under_any_weight_the_constant_one_included():
    rng = np.random.default_rng(77)
    series = _world(rng, 4, 2.5, 3.0)
    zero = Reading(np.zeros(3), np.zeros(3), series[0].cov, subspace=series[0].subspace)
    for weight in WEIGHTS.values():
        assert score_test(series + [zero], 1.9, weight)["normal_score"] == score_test(series, 1.9, weight)["normal_score"]


def test_the_two_closed_forms_and_that_no_weight_is_best_at_both():
    # weak signal everywhere: the constant weight is the only one at 100%. One strong reading among
    # noise-only readings: a rising weight survives thousands of times more of them. E13 table C1.
    assert weak_signal_efficiency(WEIGHTS["one"], 1) == pytest.approx(1.0)
    assert weak_signal_efficiency(WEIGHTS["linear"], 1) == pytest.approx(9 / 15)       # E[H^2]^2 / (E[H] E[H^3])
    assert weak_signal_efficiency(WEIGHTS["linear"], 3) == pytest.approx(225 / 315)
    assert weak_signal_efficiency(WEIGHTS["smooth"], 1) == pytest.approx(0.9205, abs=1e-4)
    assert weak_signal_efficiency(WEIGHTS["cap10"], 3) == pytest.approx(0.5907, abs=1e-4)
    assert drowning_count(WEIGHTS["one"], 1, 20.0) == pytest.approx(400.0)
    assert drowning_count(WEIGHTS["one"], 3, 20.0) == pytest.approx(400.0 / 3)
    assert drowning_count(WEIGHTS["linear"], 1, 20.0) == pytest.approx(20.0 ** 6 / 15)
    assert drowning_count(WEIGHTS["smooth"], 1, 20.0) == pytest.approx(852.2, rel=1e-3)
    for k in (1, 2, 3):
        for name, weight in WEIGHTS.items():
            efficiency = weak_signal_efficiency(weight, k)
            assert efficiency <= 1.0 + 1e-9 and (name == "one" or efficiency < 0.99)
            assert math.isfinite(drowning_count(weight, k, 20.0))      # every weight is drowned in the end
        assert drowning_count(WEIGHTS["linear"], k, 20.0) > 1000 * drowning_count(WEIGHTS["one"], k, 20.0)


def test_what_the_two_closed_forms_mean_in_a_simulated_series():
    # scalar readings in their own noise coordinates: along the candidate line U = r cos(d) + noise,
    # across it V = r sin(d) + noise, for a candidate off by the noise angle d
    rng = np.random.default_rng(9)

    def mean_score(weight, signals, shift, reps):
        signals = np.asarray(signals, dtype=float)
        u = signals * math.cos(shift) + rng.normal(size=(reps, signals.size))
        v = signals * math.sin(shift) + rng.normal(size=(reps, signals.size))
        w = weight(u * u)
        return float(np.mean(np.sum(w * u * v, axis=1) / np.sqrt(np.sum(w * w * u * u, axis=1))))
    # 400 noise-only readings halve the squared mean score of one reading at signal 20, weight 1
    alone = mean_score(WEIGHTS["one"], [20.0], 0.15, 20000)
    buried = mean_score(WEIGHTS["one"], [20.0] + [0.0] * 400, 0.15, 4000)
    assert (buried / alone) ** 2 == pytest.approx(0.5, abs=0.04)
    # the same 400 leave the rising weight's score where it was
    alone = mean_score(WEIGHTS["linear"], [20.0], 0.15, 20000)
    buried = mean_score(WEIGHTS["linear"], [20.0] + [0.0] * 400, 0.15, 4000)
    assert (buried / alone) ** 2 == pytest.approx(1.0, abs=0.04)
    # 4000 readings each at signal 0.3: the rising weight has about 0.55 of the constant weight's
    # squared mean score (0.60 in the limit of no signal; (r^2 + 3)^2 (1 + r^2) / E[U^6] at r = 0.26; measured 0.574)
    flat = mean_score(WEIGHTS["one"], [0.3] * 4000, 0.5, 3000)
    rising = mean_score(WEIGHTS["linear"], [0.3] * 4000, 0.5, 3000)
    assert (rising / flat) ** 2 == pytest.approx(0.555, abs=0.05)

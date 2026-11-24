"""The discount under one unknown common scale (T82; problems/O14.md, 24 Nov).

Harvested 24 Nov 2026 from the exploration branch `explore/e14-ours34-walls` (Antigravity with Can,
22 Nov; its Theorems E14-5, E14-6 and E14-7) and checked here on main's own functions
(met.conditional_lr, met.conditional_score, met.evidence); the probe met/ours34.py is called only to
compare. The costs, the comparison with the pair weight of 18 Nov and the three-way conflict are
Claude's. Measurements and checks; no rule is adopted.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from contradiction_in_pairs import pair_weight          # noqa: E402

from met.conditional_lr import score_test                # noqa: E402
from met.conditional_score import gaussian_score_test    # noqa: E402
from met.evidence import Curve, leftover                 # noqa: E402
from met.full import Reading, isotropic                  # noqa: E402
from met.ours34 import State                             # noqa: E402

Z2 = stats.chi2.ppf(0.95, 1)
SA = 1e-4                                                # acceleration nearly exact: direct readings of a mass


def _make(raw, sf, sa=SA):
    return [isotropic(np.atleast_1d(f), np.atleast_1d(a), sf, sa) for f, a in raw]


def _weight(raw, sf, sa=SA):
    return pair_weight([(np.atleast_1d(np.asarray(f, float)) / sf, np.atleast_1d(np.asarray(a, float)) / sa) for f, a in raw])


def _keeps(raw, sf, mass, rule, sa=SA):
    p = State(_make(raw, sf, sa)).p_value(mass)
    if rule == "bell":
        return p["p_t45"] > 0.05
    if rule == "probe":
        return p["p_value"] > 0.05
    return p["z_t45"] ** 2 <= Z2 / _weight(raw, sf, sa)


# --- E14-5: a weight on each reading's score -----------------------------------------------------------

def test_a_weight_on_each_readings_score_rejects_the_mass_where_the_curve_is_least():
    """One reading (F = 20, a = 2) beside a hundred (F = 0.2, a = 2), unit noise: the curve is least at
    mass 1. T45's unweighted score is zero there. T76's weighted scores are not: 3.003 (p = 0.0027) under
    H / (1 + H) and 12.6 under H. So with T81's point, only the constant weight keeps the point it returns."""
    readings = [isotropic([20.0], [2.0], 1.0, 1.0)] + [isotropic([0.2], [2.0], 1.0, 1.0)] * 100
    least = leftover(readings)["mass"]
    assert least == pytest.approx(1.0, abs=1e-6)
    assert score_test(readings, least)["normal_score"] == pytest.approx(0.0, abs=1e-6)
    saturating = score_test(readings, least, weight=lambda h: h / (1.0 + h))
    assert saturating["normal_score"] == pytest.approx(3.0034, abs=1e-3) and saturating["p_value"] == pytest.approx(0.00267, abs=2e-5)
    assert score_test(readings, least, weight=lambda h: h)["normal_score"] == pytest.approx(12.6, abs=1e-3)


def test_a_reading_and_its_equal_parts_have_one_curve_and_only_the_constant_weight_gives_one_score():
    whole = [isotropic([20.0 * math.sqrt(3)], [2.0 * math.sqrt(3)], 1.0, 1.0)]
    parts = [isotropic([20.0], [2.0], 1.0, 1.0)] * 3
    other = [isotropic([1.0], [2.0], 1.0, 1.0)]
    masses = np.array([0.3, 1.0, 4.0])
    assert Curve(whole + other, scale=1.0)(masses) == pytest.approx(Curve(parts + other, scale=1.0)(masses), rel=1e-12)
    assert score_test(whole + other, 4.0)["normal_score"] == pytest.approx(score_test(parts + other, 4.0)["normal_score"], rel=1e-12)
    w = lambda h: h / (1.0 + h)
    assert abs(score_test(whole + other, 4.0, weight=w)["normal_score"] - score_test(parts + other, 4.0, weight=w)["normal_score"]) > 5e-5


# --- E14-6: the curve and the score's variance cannot carry a discount ----------------------------------

def _two_series(n=50, gap=1.0):
    s_a, s_b = gap / math.sqrt(2 * (n - 1)), gap / 2 * math.sqrt(n / (n - 1))
    z = np.random.default_rng(2).normal(size=n)
    honest = 5.5 + s_b * (z - z.mean()) / z.std(ddof=1)
    return [(5.0, 1.0), (6.0, 1.0)], s_a, [(float(x), 1.0) for x in honest], s_b


def test_two_readings_far_apart_and_fifty_honest_readings_have_one_curve_and_one_score_variance():
    far, s_a, honest, s_b = _two_series()
    masses = np.array([5.0, 5.3, 5.5, 5.62, 6.1])
    k_a, k_b = Curve(_make(far, s_a), scale=1.0)(masses), Curve(_make(honest, s_b), scale=1.0)(masses)
    c_a = np.array([gaussian_score_test(_make(far, s_a), m)["variance"] for m in masses])
    c_b = np.array([gaussian_score_test(_make(honest, s_b), m)["variance"] for m in masses])
    assert k_a == pytest.approx(k_b, rel=1e-4) and c_a == pytest.approx(c_b, rel=1e-4)      # equal in the exact-acceleration limit
    # so the bell keeps the same set on both
    for m, kept in ((5.37, True), (5.35, False), (5.63, True), (5.65, False)):
        assert _keeps(far, s_a, m, "bell") is kept and _keeps(honest, s_b, m, "bell") is kept


def test_both_discounts_tell_the_two_series_apart_and_only_one_of_them_reads_a_count():
    far, s_a, honest, s_b = _two_series()
    assert _keeps(far, s_a, 5.0, "probe") and not _keeps(honest, s_b, 5.0, "probe")        # reads M_nz: 2 against 50
    assert _weight(far, s_a) == pytest.approx(0.0204, abs=3e-4) and _weight(honest, s_b) == 1.0
    assert _keeps(far, s_a, 5.0, "pairs") and not _keeps(honest, s_b, 5.0, "pairs")        # no count in it


# --- E14-7: the Student test under one unknown common scale ----------------------------------------------

def _full_type_series(rng, mass, v, sizes=(3, 2, 3, 1, 2), d=3):
    out = []
    for k in sizes:
        q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        n = q[:, :k]
        a = n @ rng.normal(size=k)
        a *= 1.5 / np.linalg.norm(a)
        mix = rng.normal(size=(2 * d, 2 * d))
        cov = 0.3 * (mix @ mix.T) / (2 * d) + np.diag([0.04] * d + [0.01] * d)
        cov = 0.5 * (cov + cov.T)
        y = np.concatenate([mass * a, a]) + math.sqrt(v) * np.linalg.cholesky(cov) @ rng.normal(size=2 * d)
        out.append(Reading(y[:d], y[d:], cov, subspace=n))
    return out


def _student(readings, mass):
    k = float(Curve(readings)(mass))
    z = gaussian_score_test(readings, mass)["normal_score"]
    dof = sum(r.subspace.shape[1] for r in readings) - 1
    return z * math.sqrt(dof / (k - z * z)), dof, z


@pytest.mark.parametrize("v,bell", [(1.0, 0.95), (7.0, 0.55)])
def test_the_student_test_holds_its_level_whatever_the_common_scale_on_the_full_type(v, bell):
    """Full correlated covariances, declared subspaces of dimension 1 to 3, five readings; the true noise
    is the declared one times sqrt(v). Written from main's curve and score: t = Z sqrt((M - 1) / (K - Z^2)).
    Recorded with 4000 series: 0.949 at v = 1, 0.952 at v = 7 (the bell: 0.953 and 0.551)."""
    rng = np.random.default_rng(int(10 * v) + 1)
    kept_t = kept_z = 0
    for i in range(1200):
        readings = _full_type_series(rng, 2.0, v)
        t, dof, z = _student(readings, 2.0)
        kept_t += abs(t) <= stats.t.ppf(0.975, dof)
        kept_z += abs(z) <= 1.96
        if i < 5:
            assert State(readings).p_value(2.0)["p_tjoint"] == pytest.approx(2 * stats.t.sf(abs(t), dof), abs=1e-10)
    assert kept_t / 1200 == pytest.approx(0.95, abs=0.02) and kept_z / 1200 == pytest.approx(bell, abs=0.04)


def test_the_scenes_under_the_joined_test_and_under_the_pair_weight():
    near, far = [(5.5, 1.0), (5.6, 1.0)], [(5.0, 1.0), (6.0, 1.0)]
    half = 0.5 * stats.t.ppf(0.975, 1)
    assert half == pytest.approx(6.3531, abs=1e-4)                       # two scalar readings: mean +- 6.353 x the gap
    # scene A: the joined test keeps 4.915 to 6.185 (E14's own table: 4.917 to 6.181); the pair weight leaves the bell
    assert _keeps(near, 0.1, 4.93, "probe") and not _keeps(near, 0.1, 4.90, "probe")
    assert _keeps(near, 0.1, 6.17, "probe") and not _keeps(near, 0.1, 6.20, "probe")
    assert _weight(near, 0.1) == 1.0 and not _keeps(near, 0.1, 5.40, "pairs") and _keeps(near, 0.1, 5.42, "pairs")
    # scene B: 0 to 11.85 against 4.52 to 6.48
    assert _keeps(far, 0.1, 0.5, "probe") and _keeps(far, 0.1, 11.8, "probe") and not _keeps(far, 0.1, 11.9, "probe")
    assert _keeps(far, 0.1, 4.53, "pairs") and not _keeps(far, 0.1, 4.50, "pairs")
    # in 3D (five leftover directions, four of them quiet) the two nearly coincide: 4.926 to 6.074, 4.935 to 6.065
    e = np.array([1.0, 0.0, 0.0])
    far3 = [(5.0 * e, e), (6.0 * e, e)]
    assert _keeps(far3, 0.1, 4.93, "probe") and not _keeps(far3, 0.1, 4.92, "probe")
    assert _keeps(far3, 0.1, 4.94, "pairs") and not _keeps(far3, 0.1, 4.93, "pairs")


def test_exact_zeros_change_neither_and_tiny_readings_undo_only_the_joined_test():
    e = np.array([1.0, 0.0, 0.0])
    far3 = [(5.0 * e, e), (6.0 * e, e)]
    for rule in ("probe", "pairs"):
        assert _keeps(far3 + [(0 * e, 0 * e)] * 20, 0.1, 5.0, rule)
    tiny = far3 + [(1e-9 * e, 1e-9 * e)] * 20
    assert not _keeps(tiny, 0.1, 5.0, "probe") and not _keeps(tiny, 0.1, 5.35, "probe")     # back to the bell's 5.362 to 5.638
    assert _keeps(tiny, 0.1, 5.0, "pairs") and _weight(tiny, 0.1) == pytest.approx(_weight(far3, 0.1), rel=1e-6)


def test_recording_each_reading_as_ten_equal_parts_narrows_both():
    parts = [(5.0 / math.sqrt(10), 1 / math.sqrt(10))] * 10 + [(6.0 / math.sqrt(10), 1 / math.sqrt(10))] * 10
    for rule, inside, outside in (("probe", 5.27, 5.25), ("pairs", 5.29, 5.26)):
        assert _keeps(parts, 0.1, inside, rule) and not _keeps(parts, 0.1, outside, rule)


# --- what the joined test costs where the declared sigma is right, and what the plug-in costs where it is not ---

def test_how_often_the_joined_test_widens_honest_direct_readings():
    """Under the declared sigma the leftover over its count is chi-square over its degrees; the joined test
    is wider than the bell whenever that exceeds z^2 / t^2. Shares: 0.88, 0.81, 0.74, 0.61 for 2, 3, 5, 20
    readings. Simulated for two readings (3000 series, recorded 0.88; median factor then 5.1)."""
    share = lambda n: stats.chi2.sf((n - 1) * Z2 / stats.t.ppf(0.975, n - 1) ** 2, n - 1)
    assert [round(float(share(n)), 2) for n in (2, 3, 5, 20)] == [0.88, 0.81, 0.74, 0.61]
    rng = np.random.default_rng(201)
    factors = []
    for _ in range(800):
        raw = [(5.0 + 0.1 * rng.normal(), 1.0 + SA * rng.normal()) for _ in range(2)]
        p = State(_make(raw, 0.1)).p_value(5.0)
        factors.append(math.sqrt(max(1.0, stats.t.ppf(0.975, 1) ** 2 * p["v_hat"] / Z2)))
    factors = np.array(factors)
    assert np.mean(factors > 1) == pytest.approx(0.877, abs=0.03) and np.median(factors[factors > 1]) == pytest.approx(5.0, abs=0.6)


def test_under_a_common_inflation_the_plug_in_falls_short_and_keeping_the_bell_on_quiet_pairs_fails():
    # the pair weight uses the normal quantile on an estimated scale: P(|t| <= 1.96) in the limit of a large inflation
    assert [round(float(1 - 2 * stats.t.sf(1.96, n - 1)), 3) for n in (2, 3, 5, 20)] == [0.7, 0.811, 0.878, 0.935]
    # measured, two readings, the noise twice the declared (3000 series recorded: bell 0.680, pairs 0.829, joined 0.970)
    rng = np.random.default_rng(204)
    kept = dict(bell=0, pairs=0, probe=0)
    for _ in range(800):
        raw = [(5.0 + 0.2 * rng.normal(), 1.0 + 2 * SA * rng.normal()) for _ in range(2)]
        for rule in kept:
            kept[rule] += _keeps(raw, 0.1, 5.0, rule)
    assert kept["bell"] / 800 == pytest.approx(0.68, abs=0.05) and kept["pairs"] / 800 == pytest.approx(0.83, abs=0.05)
    assert kept["probe"] / 800 >= 0.94
    # "the bell while a pair's scatter is below its expectation, Student above it": 0.966, 0.914, 0.875 at inflation 1, 2, 4
    z1, z2 = np.random.default_rng(8).normal(size=(2, 400000))
    for v, recorded in ((1.0, 0.966), (2.0, 0.914), (4.0, 0.875)):
        quiet = v * z2 ** 2 <= 1.0
        held = np.where(quiet, np.abs(z1) <= 1.96 / math.sqrt(v), np.abs(z1) <= stats.t.ppf(0.975, 1) * np.abs(z2))
        assert held.mean() == pytest.approx(recorded, abs=4e-3)

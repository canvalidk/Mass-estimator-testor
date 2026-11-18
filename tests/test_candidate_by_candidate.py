"""Judging each candidate mass on its own (problems/O14.md, 18 Nov, evening; tools/candidate_by_candidate.py).

Measurements, one instrument. They pin: the conditional likelihood ratio's tail probability; which forms judge the true
mass exactly; that the afternoon's signal-energy form is not exact even for one reading (a correction to that entry);
that the exact textbook form is moved by exactly zero readings (L20); and that a no-fit form carries no discount."""
import math
import os
import sys

import numpy as np
import pytest
from scipy import integrate, stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import candidate_by_candidate as cc                             # noqa: E402

NO_FIT = ("score", "signal energy >= 0", "conditional LR")


@pytest.mark.parametrize("x,T,n", [(2.0, 5.0, 4), (6.0, 0.5, 10), (3.0, 400.0, 50), (25.0, 3.0, 8), (1.0, 0.0, 3)])
def test_the_conditional_tail_probability_against_an_independent_integral(x, T, n):
    f = lambda j: stats.chi2.pdf(j, n - 1) * stats.chi2.sf(max(0.0, x * (1 - j / (x + T))), 1)
    exact = integrate.quad(f, 0, x + T, limit=400, epsabs=0, epsrel=1e-10)[0] + stats.chi2.sf(x + T, n - 1)
    assert float(cc.clr_p(x, T, n)) == pytest.approx(exact, rel=0.02, abs=3e-4)


def test_the_conditional_tail_probability_is_the_law_of_lr_given_the_along_line_parts():
    """Draw the leftover as N(0, I) beside fixed along-line parts and compute LR = K - K_min from the Gram itself."""
    rng = np.random.default_rng(7)
    t = np.array([1.5, -0.3, 0.0, 2.2, 0.7, -1.1])
    r = rng.normal(size=(400000, t.size))
    T, K, S = float(t @ t), (r * r).sum(1), r @ t
    lr = 0.5 * (K - T + np.sqrt((K - T) ** 2 + 4 * S * S))
    for x in (0.5, 2.0, 6.0):
        assert float((lr >= x).mean()) == pytest.approx(float(cc.clr_p(x, T, t.size)), abs=0.004)


def test_the_rest_is_the_grams_determinant_over_the_along_line_energy():
    """K = S^2/T + lam+ lam- / T at every candidate: the Gaussian deviance is the score part plus the rest."""
    rng = np.random.default_rng(8)
    Y = rng.normal(size=(30, 2)) + 0.8
    for phi in (0.2, 0.9, 1.4):
        u, v = np.array([math.sin(phi), math.cos(phi)]), np.array([math.cos(phi), -math.sin(phi)])
        t, r = Y @ u, Y @ v
        lam = np.linalg.eigvalsh(Y.T @ Y)
        assert (r @ r) == pytest.approx((t @ r) ** 2 / (t @ t) + lam[0] * lam[1] / (t @ t), rel=1e-10)


def test_which_forms_judge_the_true_mass_exactly():
    noise20 = cc.size_at_truth(20, 0.0, reps=40000)
    weak5 = cc.size_at_truth(5, 0.5, reps=40000)
    one = cc.size_at_truth(1, 0.0, reps=40000)
    for cell in (noise20, weak5, one):
        assert cell["score"] == pytest.approx(0.05, abs=0.006)
        assert cell["conditional LR"] == pytest.approx(0.05, abs=0.006)
    assert noise20["gaussian"] > 0.45 and noise20["gaussian, Gram share"] > 0.2          # the fit, at weak signal
    assert weak5["gaussian, Gram share"] > 0.09
    assert one["signal energy >= 0"] < 0.02 and noise20["signal energy >= 0"] > 0.065    # not exact, either way


def test_for_one_reading_the_afternoon_form_is_not_the_gaussian_factor():
    """Corrects the afternoon entry's '(for one scalar reading it does reduce to the Gaussian factor)': the score law
    and the conditional LR do; the signal-energy form gives 2c/(1 + c) of it at half-angle cosine c."""
    r = cc.one_reading(energy=9.0, degrees=(10, 30, 60, 85))
    for name in ("gaussian, Gram share", "score", "conditional LR"):
        assert np.allclose(r[name], 1.0, atol=2e-3), name
    c = np.cos(np.radians([10, 30, 60]))
    assert np.allclose(r["signal energy >= 0"][:3], 2 * c / (1 + c), atol=1e-9)
    assert r["signal energy >= 0"][2] == pytest.approx(2 / 3, abs=1e-9)


def test_on_noise_the_gaussian_law_narrows_with_the_count_and_the_no_fit_forms_keep_a_bounded_tilt():
    """O15 (9 Nov): the tilt on pure noise cannot vanish beside L11's rate; here it stays bounded, whatever the count."""
    r = cc.noise_only(500, n_series=60)
    assert r["gaussian"][1] < 0.4 and r["gaussian"][0] > 25
    for name in NO_FIT:
        assert r[name][1] > 0.97 and r[name][0] < 3.5, name


def test_the_fading_history_that_broke_the_fixed_anchor():
    r = cc.history(10000, beta=0.75, n_series=60)
    assert r["gaussian"]["inside"] < 0.78
    assert r["score"]["false_side"] > 0.6
    for name in ("gaussian, Gram share", "signal energy >= 0", "conditional LR"):
        assert r[name]["inside"] > 0.86 and r[name]["false_side"] < 1e-3, name
    assert r["conditional LR"]["log_width"] == pytest.approx(r["gaussian, Gram share"]["log_width"], rel=0.05)


def test_exactly_zero_readings_move_the_conditional_lr_and_nothing_else():
    """L20 (Z). The other four are functions of the Gram alone; the conditional LR's reference law has the count in it."""
    base, more = cc.history(20, 1.0, n_series=40), cc.history(20, 1.0, n_series=40, zeros=20)
    for name in cc.LAWS:
        if name != "conditional LR":
            assert more[name] == base[name], name
    assert more["conditional LR"]["log_width"] > 3 * base["conditional LR"]["log_width"]


def test_an_exact_statistic_no_zero_reading_enters_sees_the_false_peak():
    r = cc.weighted_leftover(200, n_series=80)
    assert abs(r["truth"][0]) < 0.5 and r["false peak"][0] > 6
    assert r["false peak"][1] > r["false peak"][0]                 # the unweighted leftover sees it better, with the count
    # exact conditional moments of Q = sum t^2 (r^2 - 1) given t, and a zero along-line part is not in it
    rng = np.random.default_rng(9)
    t = np.array([1.2, 0.0, -0.4, 2.0, 0.0])
    rr = rng.normal(size=(400000, t.size))
    q = ((t * t) * (rr * rr - 1)).sum(1)
    assert q.mean() == pytest.approx(0.0, abs=0.03) and q.var() == pytest.approx(2 * (t ** 4).sum(), rel=0.02)
    assert np.array_equal(((t * t)[[0, 2, 3]] * (rr[:, [0, 2, 3]] ** 2 - 1)).sum(1), q)


def test_a_no_fit_form_carries_no_discount_and_holds_the_right_angle_scene_with_the_one_weight():
    _, plain = cc.right_angle_scene(0.01, False)
    w_a, a = cc.right_angle_scene(0.05, True)
    w_b, b = cc.right_angle_scene(0.01, True)
    assert w_b == pytest.approx(w_a / 25, rel=0.02)                 # the weight goes as sigma^2
    for name in cc.LAWS:
        assert plain[name] > 100, name
        assert 1.0 < b[name] < 3.0 and b[name] == pytest.approx(a[name], rel=0.02), name


def test_read_at_its_least_the_along_line_weighted_leftover_reports_nothing_at_right_angles():
    """Is the discount the exact, zero-safe weighted leftover read at its least? Not in scene C: it is zero there for
    every sigma, while the energy-weighted leftover grows as 1/sigma^2. As one 2D reading the same numbers are seen."""
    far = cc.leftover_at_its_least([([5.0], [1.0]), ([6.0], [1.0])], 0.1, 0.002, grid=20001)
    assert far["squared scores"] == pytest.approx(far["along-line weighted"], rel=1e-12)
    assert far["squared scores"] == pytest.approx(24.7, abs=0.2) and far["energy weighted"] == pytest.approx(24.7, abs=0.2)
    two = cc.leftover_at_its_least([([2.0], [0.0]), ([0.0], [1.0])], 0.05, 0.05, grid=20001)
    assert two["squared scores"] < 1e-4 and two["along-line weighted"] < 1e-4 and two["energy weighted"] > 70
    one = cc.leftover_at_its_least([([2.0, 0.0], [0.0, 1.0])], 0.05, 0.05, grid=20001)
    assert one["squared scores"] < 1e-4 and one["along-line weighted"] == pytest.approx(200, rel=0.01)


def test_leaving_exact_zeros_out_of_the_count_restores_z_only_at_exact_zero():
    r = cc.nearly_zero_readings(n_series=20)
    assert r["exact zeros, left out of the count"] == r["nothing added"]
    assert r["nearly zero, counted"] > 3 * r["nothing added"]

"""met/forward.py: the forward model (truth held still) against met.full and the exact laws.

The engines must agree with each other and with the code the equations run; the draws must follow
the standard laws wherever one exists (noncentral chi-square, the Wishart smallest eigenvalue, the
chi-bar-squared boundary mixture, F and t under the inverse-gamma size). Seeds are fixed, so each
KS test is a deterministic check, set at p > 0.001.
"""
import json
import math

import numpy as np
import pytest
from scipy import stats

from met.forward import (boundary_positive_misfit, forward, from_invariants, invariants,
                         null_misfit, null_positive_misfit, on_line, reference, run, summarize,
                         true_mass, validate)
from met.full import Reading
from tools.positivity_supremum import expectations

P_MIN = 1e-3
NULL_Q_PLUS = {1: 0.5 - 1 / math.pi, 2: 1.5 - math.pi / 4, 3: 2 - 2 / math.pi}   # B02, T38
NULL_Q = {1: 0.0, 2: 2 - math.pi / 2, 3: 1.0}


def _close(a, b, tol=1e-9):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    both = np.isfinite(a) & np.isfinite(b)
    assert np.array_equal(np.isfinite(a), np.isfinite(b))
    assert np.all(np.abs(a[both] - b[both]) <= tol * np.maximum(1.0, np.abs(b[both])))


@pytest.mark.parametrize("k,push,angle_deg,c", [(1, 1.5, 30, 1.0), (2, 0.7, 70, 2.5), (3, 2.0, 45, 0.4),
                                                 (3, 0.0, 45, 1.0)])
def test_isotropic_engine_is_met_full(k, push, angle_deg, c):
    """Every measure from the closed forms equals met.full's (reading_state, _h_terms, bartlett_eta),
    including curv and J converted to met.full's parametrisation m = c tan(phi)."""
    truths = from_invariants(k, push, math.radians(angle_deg), readings=4)
    mass = math.tan(math.radians(angle_deg))
    fast = forward(truths, draws=25, seed=11, mass=mass, c=c, engine="isotropic")
    slow = forward(truths, draws=25, seed=11, mass=mass, c=c, engine="general")
    for key in ("E", "H_true", "K_true", "Q", "Q_plus", "g29", "g30", "m_own", "K_total",
                "eta_flat", "eta29", "eta30", "curv30", "J30", "m_series", "latent", "score", "misalignment"):
        _close(fast[key], slow[key])


def test_general_engine_on_correlated_covariance_and_subspace():
    """The general engine runs the full type; on the true line the residual is chi2_k and splits into
    a unit score and k - 1 misalignment dimensions, whatever the correlation."""
    rng = np.random.default_rng(3)
    a = rng.normal(size=(3, 3))
    cov = a @ a.T + 3 * np.eye(3)
    big = np.block([[cov, 0.4 * np.eye(3)], [0.4 * np.eye(3), 0.5 * cov]])
    n = np.linalg.qr(rng.normal(size=(3, 2)))[0]
    acc = n @ np.array([1.0, -0.5])
    truth = Reading(2.0 * acc, acc, big, n)
    res = forward([truth], draws=400, seed=5, mass=2.0, series=False)
    assert res["truth"]["engine"] == "general"
    _close(res["K_true"], res["score"] ** 2 + res["misalignment"], 1e-9)
    assert stats.kstest(res["K_true"].ravel(), stats.chi2(2).cdf).pvalue > P_MIN
    assert stats.kstest(res["score"].ravel(), "norm").pvalue > P_MIN


@pytest.mark.parametrize("k", (1, 2, 3))
@pytest.mark.parametrize("push", (0.5, 30.0))
def test_true_mass_residual_and_its_channels_are_exact_at_every_push(k, push):
    """T36 (4) refined: at the true mass the residual is chi2_k at every push, and it is the sum of a
    standard normal score (1 dim) and a chi2_(k-1) misalignment, independent; latent is chi2_k."""
    truths = from_invariants(k, push, math.radians(35))
    res = forward(truths, draws=20000, seed=20 + k, series=False)
    ref = reference(truths[0])
    for key in ("E", "H_true", "K_true", "score", "latent") + (("misalignment",) if k > 1 else ()):
        assert stats.kstest(res[key].ravel(), ref[key].cdf).pvalue > P_MIN, key
    _close(res["K_true"], res["score"] ** 2 + res["misalignment"], 1e-9)
    assert abs(np.corrcoef(res["score"].ravel() ** 2, res["latent"].ravel())[0, 1]) < 0.04
    if k == 1:
        assert np.all(res["misalignment"] == 0.0)


@pytest.mark.parametrize("k", (2, 3))
def test_zero_push_laws(k):
    """At zero push: Q is the smaller eigenvalue of W_2(k, I); Q_plus is half that and half the smaller
    of two chi2_k; their means are B02's exact values."""
    assert null_misfit(k).mean() == pytest.approx(NULL_Q[k], abs=1e-9)
    assert null_positive_misfit(k).mean() == pytest.approx(NULL_Q_PLUS[k], abs=1e-9)
    res = forward(from_invariants(k, 0.0, math.radians(45)), draws=20000, seed=30 + k, mass=1.0,
                  series=False)
    assert stats.kstest(res["Q"].ravel(), null_misfit(k).cdf).pvalue > P_MIN
    assert stats.kstest(res["Q_plus"].ravel(), null_positive_misfit(k).cdf).pvalue > P_MIN


def test_zero_push_one_dimension():
    assert null_positive_misfit(1).mean() == pytest.approx(NULL_Q_PLUS[1], abs=1e-9)
    res = forward(from_invariants(1, 0.0, math.radians(45)), draws=20000, seed=34, mass=1.0, series=False)
    assert np.all(res["Q"] < 1e-12)                       # zero up to rounding: a rank-one Gram
    assert res["Q_plus"].mean() == pytest.approx(NULL_Q_PLUS[1], abs=4 * res["Q_plus"].std() / math.sqrt(20000))


@pytest.mark.parametrize("k,push,angle_deg", [(2, 1.0, 45), (3, 2.0, 60), (3, 4.0, 80), (1, 1.0, 70)])
def test_means_against_b02_quadrature(k, push, angle_deg):
    """E[Q] and E[Q_plus] at intermediate push against B02's independent quadrature (T38's integral):
    R = push sin(theta0) in the force channel, b = push cos(theta0) in the acceleration channel."""
    th = math.radians(angle_deg)
    exact = expectations(k, push * math.sin(th), push * math.cos(th), tolerance=1e-10)
    res = forward(from_invariants(k, push, th), draws=40000, seed=40 + k, series=False)
    for key, value in (("Q_plus", exact["positive"]), ("Q", exact["unsigned"])):
        draws = res[key].ravel()
        assert abs(draws.mean() - value) < 4 * draws.std() / math.sqrt(draws.size) + 1e-12, key


@pytest.mark.parametrize("k", (2, 3))
def test_boundary_strong_push_is_the_chi_bar_square(k):
    """One channel empty, the other strong: Q_plus tends to 1/2 chi2_(k-1) + 1/2 chi2_k, mean k - 1/2
    (Chernoff's boundary mixture; T38's supremum, approached, not attained)."""
    res = forward(from_invariants(k, 300.0, 0.0), draws=20000, seed=50 + k, series=False)
    law = boundary_positive_misfit(k)
    assert law.mean() == k - 0.5
    assert stats.kstest(res["Q_plus"].ravel(), law.cdf).pvalue > P_MIN


def test_boundary_strong_push_one_dimension_has_an_atom_at_zero():
    res = forward(from_invariants(1, 300.0, 0.0), draws=20000, seed=53, series=False)
    q = res["Q_plus"].ravel()
    atom = q < 1e-9                                       # C > 0: Q_plus = Q = 0 up to rounding
    assert np.mean(atom) == pytest.approx(0.5, abs=0.015)
    assert stats.kstest(q[~atom], stats.chi2(1).cdf).pvalue > P_MIN


@pytest.mark.parametrize("k", (2, 3))
def test_inside_the_half_line_q_plus_tends_to_k_minus_one_and_the_cap_tempers(k):
    """At m0 = s with strong push Q_plus is chi2_(k-1): kappa = k - 1/2 is the supremum over worlds,
    not the typical value, and ours30's cap acts on P(chi2_(k-1) > k - 1/2) of readings whose noise is
    exactly the declared noise (B04's 22-29%)."""
    res = forward(from_invariants(k, 300.0, math.radians(45)), draws=20000, seed=60 + k, series=False)
    assert stats.kstest(res["Q_plus"].ravel(), stats.chi2(k - 1).cdf).pvalue > P_MIN
    expected = stats.chi2(k - 1).sf(k - 0.5)
    assert summarize(res)["tempered30"] == pytest.approx(expected, abs=0.012)


def test_inverse_gamma_size_gives_f_and_t_laws():
    """v ~ InvGamma(nu/2, nu/2): E[1/v] = 1 (G20), K_true/k ~ F(k, nu), score ~ t_nu,
    misalignment/(k-1) ~ F(k-1, nu); the same z as the fixed-size run with the same seed."""
    truths = from_invariants(3, 3.0, math.radians(50))
    res = forward(truths, draws=20000, seed=70, nu_eff=8.0, series=False, keep_noise=True)
    fixed = forward(truths, draws=20000, seed=70, series=False, keep_noise=True)
    assert np.array_equal(res["z"], fixed["z"])
    assert (1 / res["v"]).mean() == pytest.approx(1.0, abs=0.02)
    ref = reference(truths[0], nu_eff=8.0)
    for key in ("K_true", "score", "misalignment", "latent"):
        assert stats.kstest(res[key].ravel(), ref[key].cdf).pvalue > P_MIN, key


def test_series_residual_and_units():
    """K_total ~ chi2 with sum k_i degrees (T36 (4)); scaling the force channel's truth and sigma
    together (L3) leaves every measure but the masses unchanged, draw for draw, and scales m."""
    acc = [[1.0, 0.0, 0.0], [0.3, 0.2, 0.0], [0.0, 0.0, 0.0], [2.0, -1.0, 0.5]]
    truths = on_line(1.7, acc, 0.4, 0.2)
    assert true_mass(truths) == pytest.approx(1.7)
    res = forward(truths, draws=10000, seed=80)
    assert stats.kstest(res["K_total"], stats.chi2(12).cdf).pvalue > P_MIN
    scaled = forward(on_line(1.7 * 9.0, acc, 0.4 * 9.0, 0.2), draws=10000, seed=80)
    for key in ("E", "K_true", "Q", "Q_plus", "g30", "eta30", "score", "misalignment"):
        _close(res[key], scaled[key], 1e-10)
    _close(9.0 * res["m_series"], scaled["m_series"], 1e-9)


def test_zero_signal_reading_has_no_score_direction():
    truths = on_line(2.0, [[1.0, 0.0], [0.0, 0.0]], 1.0, 1.0)
    res = forward(truths, draws=50, seed=1)
    assert np.all(np.isfinite(res["score"][:, 0])) and np.all(np.isnan(res["score"][:, 1]))
    assert res["truth"]["information"][1] == 0.0
    assert res["truth"]["information"][0] == pytest.approx(1.0 / (1.0 + 4.0))   # |a|^2/(sF^2 + m^2 sa^2)


def test_invariants_of_a_truth():
    inv = invariants(from_invariants(3, 2.0, math.radians(60), force_sd=3.0, acceleration_sd=0.5)[0])
    assert inv["push"] == pytest.approx(2.0) and inv["angle_deg"] == pytest.approx(60.0)
    assert inv["R"] == pytest.approx(2.0 * math.sin(math.radians(60)))
    assert inv["s"] == pytest.approx(6.0) and inv["cos_psi"] == pytest.approx(1.0)


def test_preset_validates_and_runs(tmp_path):
    import yaml
    data = yaml.safe_load(open("presets/forward_noise_measures.yaml", encoding="utf-8"))
    spec = validate(data)
    assert spec["kind"] == "forward" and len(spec["cells"]) == 3 * 7 * 2 * 3 * 2
    small = dict(data, grid=dict(dimension=[3], push=[0.0, 30.0], angle_deg=[45.0], readings=[1, 5],
                                 nu_eff=[float("inf")]))
    path = tmp_path / "small.yaml"
    path.write_text(yaml.safe_dump(small), encoding="utf-8")
    out, first = run(path, tmp_path, replicates=3000, progress=None)
    _, second = run(path, tmp_path, replicates=3000, progress=None)
    assert json.dumps(first["rows"]) == json.dumps(second["rows"])      # NaN-safe, bit for bit
    rows = {(r["push"], r["readings"]): r for r in first["rows"]}
    null = rows[(0.0, 1)]
    assert null["Q_plus_null_exact"] == pytest.approx(NULL_Q_PLUS[3], abs=1e-9)
    assert abs(null["Q_plus_mean"] - null["Q_plus_null_exact"]) < 4 * null["Q_plus_se"]
    strong = rows[(30.0, 5)]
    assert abs(strong["K_true_mean"] - 3) < 4 * strong["K_true_se"]
    assert (out / "summary.csv").exists()

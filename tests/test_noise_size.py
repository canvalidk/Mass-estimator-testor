"""T37: the per-reading noise size (definition/theorems/T37.md).

Sigma is the covariance of a reading's noise (L18, L13); given its own noise size v the reading
is Gaussian (L8's conditional form), with E[v] = 1. Under L21's cut the reading's precision read
from its own misfit is g(Q) = E[1/v | Q], with Q | v ~ v chi2_nu at strong signal, nu = k - 1.
"""
import math
import os

import numpy as np
import pytest
from scipy import stats

from met.full import _h_terms, isotropic, ours30, reading_state

_T = np.linspace(-30.0, 40.0, 200001)            # t = ln v


def _cut_precision(logpdf, q, nu):
    lw = logpdf(np.exp(_T)) - nu / 2 * _T - q / (2 * np.exp(_T)) + _T
    w = np.exp(lw - lw.max())
    return float((w * np.exp(-_T)).sum() / w.sum())


def _mean(logpdf):
    lw = logpdf(np.exp(_T)) + _T
    w = np.exp(lw - lw.max())
    return float((w * np.exp(_T)).sum() / w.sum())


def _invgamma(a):
    return lambda v: stats.invgamma.logpdf(v, a, scale=a - 1)


def _pareto(a):
    xm = (a - 1) / a
    return lambda v: np.where(v >= xm, np.log(a) + a * np.log(xm) - (a + 1) * np.log(np.maximum(v, xm)), -np.inf)


FAMILIES = {"invgamma 1.5": _invgamma(1.5), "invgamma 3": _invgamma(3.0), "pareto 1.5": _pareto(1.5),
            "pareto 3": _pareto(3.0),
            "lognormal": lambda v: stats.lognorm.logpdf(v, 1.0, scale=math.exp(-0.5)),
            "uniform": lambda v: np.where((v > 0.5) & (v < 1.5), 0.0, -np.inf)}


def test_the_families_have_mean_one():
    for name, lp in FAMILIES.items():
        assert _mean(lp) == pytest.approx(1.0, abs=2e-3), name


@pytest.mark.parametrize("nu", [0, 1, 2])
def test_cut_precision_inverse_gamma_closed_form(nu):
    a = 3.0
    for q in (0.0, 0.5, 4.0, 50.0):
        assert _cut_precision(_invgamma(a), q, nu) == pytest.approx((2 * a + nu) / (2 * a - 2 + q), rel=1e-4)


@pytest.mark.parametrize("nu", [0, 1, 2])
def test_t37_tail_coefficient_is_nu_plus_twice_the_tail_index(nu):
    for a in (1.5, 3.0):
        for lp in (_invgamma(a), _pareto(a)):
            assert 1e6 * _cut_precision(lp, 1e6, nu) == pytest.approx(nu + 2 * a, rel=1e-3)


@pytest.mark.parametrize("nu", [0, 2])
def test_t37_tails_lighter_than_a_power_fail_t30(nu):
    for name in ("lognormal", "uniform"):
        lp = FAMILIES[name]
        assert 1e6 * _cut_precision(lp, 1e6, nu) > 1.3 * 1e4 * _cut_precision(lp, 1e4, nu), name


@pytest.mark.parametrize("nu", [0, 1, 2])
def test_t37_the_clean_end_is_above_one(nu):
    for name, lp in FAMILIES.items():
        assert _cut_precision(lp, 1e-9, nu) > 1.05, name


@pytest.mark.parametrize("k", [1, 2, 3])
def test_t37_k_minus_half_needs_a_noise_size_with_infinite_mean(k):
    """The tail coefficient nu + 2 alpha equals ours30's kappa = k - 1/2 only at alpha = 1/4; a Pareto
    tail of index 1/4 has no mean (its truncated mean grows without bound)."""
    nu, kappa = k - 1, k - 0.5
    alpha = (kappa - nu) / 2
    assert alpha == 0.25
    x = np.logspace(0, 12, 4001)
    trunc = [float(np.trapezoid(x[x <= c] * alpha * x[x <= c] ** (-alpha - 1), x[x <= c])) for c in (1e4, 1e8, 1e12)]
    assert trunc[2] > 10 * trunc[1] > 100 * trunc[0] > 0


@pytest.mark.parametrize("k", [2, 3])
def test_t37_cut_factor_is_scale_free_and_keeps_the_direction(k):
    """A fixed 30 degree pair, declared noise shrinking: the cut factor for an inverse-gamma noise size
    converges to a sigma-free curve, peaked at the reading's own best line."""
    a, th = 3.0, math.radians(30)
    f = np.zeros(k); f[0], f[1] = 2 * math.cos(th), 2 * math.sin(th)
    acc = np.zeros(k); acc[0] = 1.0
    phi = np.linspace(0.3, 1.4, 401)
    curves = []
    for sd in (1e-3, 1e-5, 1e-7):
        r = isotropic(f, acc, sd, sd)
        st = reading_state(r, 1.0)
        along = np.maximum(r.pre()["energy"] - st["Q"] - _h_terms(r.pre(), 1.0, phi), 0.0)
        curves.append(-(a + (k - 1) / 2) * np.log1p(along / (2 * (a - 1) + st["Q"])))
    assert np.allclose(curves[1], curves[2], atol=1e-6)
    best = phi[int(np.argmax(curves[2]))]
    assert best == pytest.approx(st["phi"], abs=0.005)   # the reading's own best line (noise units)


def _mixture_interval(readings, a):
    lm = np.linspace(math.log(1e-2), math.log(1e2), 6001)
    m = np.exp(lm)
    lp = np.log(1 / (1 + m * m)) + lm
    for r in readings:
        st = reading_state(r, 1.0)
        along = np.maximum(r.pre()["energy"] - st["Q"] - _h_terms(r.pre(), 1.0, np.arctan(m)), 0.0)
        lp = lp - (a + (r.k - 1) / 2) * np.log1p(along / (2 * (a - 1) + st["Q"]))
    p = np.exp(lp - lp.max())
    cdf = np.cumsum(p) / p.sum()
    return m[np.searchsorted(cdf, 0.025)], m[np.searchsorted(cdf, 0.975)]


def test_t37_world_check_sigma_is_the_covariance():
    """3D, strong signal, N = 20, noise size InvGamma(3) (mean 1): the mixture law (a = 3 declared)
    is narrower than ours30 at no loss of coverage. MET_SLOW=1: 300 series, the numbers in T37."""
    n = 300 if os.environ.get("MET_SLOW") else 40
    rng = np.random.default_rng(3011)
    cov_o, cov_m, w_o, w_m = [], [], [], []
    for _ in range(n):
        rd = []
        for _ in range(20):
            u = rng.normal(size=3); u /= np.linalg.norm(u)
            x = 20.0 * u / math.sqrt(2)
            v = 2.0 / rng.gamma(3.0)
            rd.append(isotropic(x + math.sqrt(v) * rng.normal(size=3), x + math.sqrt(v) * rng.normal(size=3), 1.0, 1.0))
        lo, hi = ours30(rd, 1.0, readout=False)["interval"]
        cov_o.append(lo <= 1 <= hi); w_o.append(math.log(hi / lo))
        lo, hi = _mixture_interval(rd, 3.0)
        cov_m.append(lo <= 1 <= hi); w_m.append(math.log(hi / lo))
    assert np.median(w_m) < 0.95 * np.median(w_o)
    assert np.mean(cov_m) > 0.88
    if os.environ.get("MET_SLOW"):
        print(f"\nours30: cover {np.mean(cov_o):.2f} width {np.median(w_o):.3f} | "
              f"mixture: cover {np.mean(cov_m):.2f} width {np.median(w_m):.3f}")

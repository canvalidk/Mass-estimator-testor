"""ours33: a probe. ours32's reading level with the normalised residual density (T48) and eta applied as a
calibration of the output law (T47), not as a power on the likelihood.

NOT the definition. ours30 is the working equation (R-2026-10-29-1). Assembled by Claude (4 Nov) from
two results Codex set down on 3 Nov and did not build on; problems/O18.md states its holes before any
result. Built to break (GOAL.md); each change answers a break T39 proved in ours32.

Per reading (T48): the declared noise's density of the reading's k-dimensional residual at m, with the
noise size v_i ~ IG(nu_eff_i/2, nu_eff_i/2) (T37, G20) integrated out,

  l_i(phi) = -(nu_eff_i + k_i)/2 * log(1 + K_i(phi)/nu_eff_i),     K_i = E_i - H_i  (T36's residual),

on the whole projective line, up to an m-independent constant. ours31/32 had exponent (nu + k - 1)/2,
the fitted count; T48's integration of the normalised Gaussian over the candidate plane gives k. This
is the change that makes a single reading's raw law concentrate for every nu > 0 (T39 part 4 fails
for nu < 1 under the old exponent; T48 section 3).

Raw law (eta = 1): p_0(m) ~ pi(m) prod_i exp(l_i(m)) on m > 0, proper (every factor in (0, 1]).

Series (T31, as ours32): phi_hat = argmax sum_i l_i on the whole line, curv = -sum l_i'', J = sum l_i'^2,
eta = min(1, curv / J) (1 if J = 0 or curv <= 0).

Calibration (T47): the quoted law is the normal-score warp of the raw law's CDF,
F_eta = Phi(sqrt(eta) Phi^{-1}(F_0)). Its median is the raw median exactly; its central 95% interval is
the raw law's quantiles at Phi(-+ z_0.975 / sqrt(eta)). The raw law is never raised to a power, so
the law stays proper and concentrating whenever the raw law does (T39 part 5 fails for the powered
law at nu = 10; T47 section 3). eta's job is unchanged: it widens the report; it never moves the point.

Point: the raw law's median. All-zero series: eta = 1 and the law is pi, so the point is c (T24).
"""
import math

import numpy as np
from scipy.stats import norm

from met.full import geometric_mean, _h_terms, _log_prior, _quantiles
from met.ours31 import COARSE_HALF, COARSE_NODES, FINE_NODES, SUPPORT_NATS, _per_reading
from met.ours32 import _residual_energy, bartlett_eta

Z95 = float(norm.ppf(0.975))


def _l_terms(states, c, phi, order):
    """Per-reading log factors with T48's exponent (nu_eff + k)/2, up to m-independent constants, and
    their first two phi-derivatives. The same residual as met.ours32, one more residual dimension."""
    per = []
    for r, ne, _, _ in states:
        a = 0.5 * (ne + r.k)
        residual = _residual_energy(r, c, phi)
        l0 = -a * np.log1p(residual / ne)
        if order == 0:
            per.append((l0,))
            continue
        _, h1, h2 = _h_terms(r.pre(), c, phi, 2)
        denominator = ne + residual
        l1 = a * h1 / denominator
        l2 = a * (h2 / denominator + (h1 / denominator) ** 2)
        per.append((l0, l1, l2))
    return per


def warp_probabilities(eta, probs):
    """T47: the raw-law probability ranks whose quantiles are the calibrated law's at `probs`.
    Phi(Phi^{-1}(p) / sqrt(eta)); p = 1/2 is fixed for every eta."""
    eta = float(eta)
    return [float(norm.cdf(norm.ppf(p) / math.sqrt(eta))) for p in probs]


def _tail_quantiles(lm, lp, p):
    """The quantiles of m at ranks p and 1 - p of the law exp(lp) on the grid lm, for p that may be far
    below floating-point resolution of 1 - p (the warp asks for Phi(-z/sqrt(eta))). Both tails are
    accumulated from their own end in log space (trapezoid masses, logaddexp), so a rank of 1e-80 is
    resolved on either side; the crossing is found by linear interpolation of the log-CDF in ln m.
    Returns (m_lo, m_hi)."""
    lp = np.asarray(lp, dtype=float)
    step = lm[1] - lm[0]
    cell = np.logaddexp(lp[:-1], lp[1:]) + math.log(0.5 * step)      # log mass of each cell
    total = np.logaddexp.reduce(cell)
    target = math.log(p) + total

    def crossing(cells, nodes):
        cum = np.logaddexp.accumulate(cells)
        j = int(np.searchsorted(cum, target))
        if j >= len(cum):
            return float(nodes[-1])
        if j == 0:
            return float(nodes[1]) if cum[0] >= target else float(nodes[0])
        a, b = cum[j - 1], cum[j]
        frac = (target - a) / (b - a) if b > a else 1.0
        return float(nodes[j] + frac * (nodes[j + 1] - nodes[j]))

    lo = crossing(cell, lm)
    hi = -crossing(cell[::-1], -lm[::-1])
    return math.exp(lo), math.exp(hi)


def ours33(readings, c, nu_eff, *, calibrate=True, grid=None, eta_method="fitted"):
    """The ours33 probe on one series. Same arguments as ours32. calibrate=False gives eta = 1 (the raw
    law; the quoted interval is then the raw interval). eta_method="fitted" retains the original
    H/J rule. eta_method="conditional" explicitly selects B15's study-only model variance H/C;
    its fixed-candidate exactness is scoped in met.eta_conditional, not a sampling coverage claim.

    Returns a dict: point, median (the raw median), interval (95%, warped by eta), raw_interval (the raw
    law's 95%), eta, law (the raw law's weights on lm), lm, Q_plus, nu_eff."""
    if not c > 0:
        raise ValueError("the starting law's centre c must be positive")
    if eta_method not in ("fitted", "conditional"):
        raise ValueError("eta_method must be 'fitted' or 'conditional'")
    readings = list(readings)
    if np.ndim(nu_eff) == 0:
        nu_eff = [float(nu_eff)] * len(readings)
    nu_eff = [float(x) for x in nu_eff]
    if len(nu_eff) != len(readings) or not all(x > 0 and math.isfinite(x) for x in nu_eff):
        raise ValueError("nu_eff: one positive finite value per reading")
    states = _per_reading(readings, c, nu_eff)
    if calibrate and eta_method == "conditional":
        from met.eta_conditional import eta_c
        eta = eta_c(readings, c, nu_eff, factors="ours33")[0]
    else:
        eta = bartlett_eta(states, c, terms=_l_terms)[0] if calibrate else 1.0

    def log_law(lm):
        if not states:
            return _log_prior(lm, c)
        phi = np.arctan(np.exp(lm - math.log(c)))
        return _log_prior(lm, c) + sum(p[0] for p in _l_terms(states, c, phi, 0))

    if grid is not None:
        half, nodes = grid
        lm = np.linspace(math.log(c) - half, math.log(c) + half, int(nodes))
    else:
        coarse = np.linspace(math.log(c) - COARSE_HALF, math.log(c) + COARSE_HALF, COARSE_NODES)
        lp = log_law(coarse)
        keep = np.nonzero(lp > lp.max() - SUPPORT_NATS)[0]
        step = coarse[1] - coarse[0]
        lm = np.linspace(coarse[keep[0]] - 2 * step, coarse[keep[-1]] + 2 * step, FINE_NODES)
    lp = log_law(lm)
    (rlo, med, rhi), p = _quantiles(lm, lp, [0.025, 0.5, 0.975])
    # The warped ranks can be far into the raw law's tails (eta = 1/90 asks for Phi(-18.6) ~ 1e-77), where
    # the adaptive grid, cut at e^-40 of the peak, has no nodes. The prior's tails are Cauchy, so the
    # rank p sits near |ln(m/c)| ~ -ln p: widen a dedicated grid to reach it.
    plo = warp_probabilities(eta, [0.025])[0]
    if eta >= 1.0:
        lo, hi = rlo, rhi
    else:
        half = max(COARSE_HALF, 12.0 - math.log(plo))
        if grid is None:
            lm_w = np.linspace(math.log(c) - half, math.log(c) + half, int(FINE_NODES * half / COARSE_HALF))
            lp_w = log_law(lm_w)
        else:
            lm_w, lp_w = lm, lp
        lo, hi = _tail_quantiles(lm_w, lp_w, plo)
    return dict(point=float(med), median=float(med), geometric_mean=geometric_mean(lm, p),
                interval=(float(lo), float(hi)), raw_interval=(float(rlo), float(rhi)), eta=eta, law=p, lm=lm,
                Q_plus=[s[2]["Q_plus"] for s in states], nu_eff=nu_eff)

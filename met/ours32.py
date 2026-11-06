"""ours32: a probe. ours31's reading level with eta put back (Can, 1 Nov: "put it back in ... test it").

NOT the definition. ours30 is the working equation (R-2026-10-29-1). ours31 (problems/O18.md) broke at
weak signal in 3D with many readings; the break was re-attributed to the flat latent measure, which
eta had been patching. ours32 puts eta back on ours31's reading-level factor so the reading level can
be tested with the series-level patch in place; the signal-size question (what eta stands in for) is
left for the next equation (O15, B14).

Per reading (as ours31):

  l_i(phi) = -(nu_eff_i + k_i - 1)/2 * log(1 + (H+_i - H_i(phi)) / (nu_eff_i + Q+_i))

  on the whole projective line. H_i <= E_i gives H+_i - H_i >= -Q+_i, so the argument of the log is
  at least nu_eff/(nu_eff + Q+) > 0: l_i is smooth on the whole line, and on m > 0 it is log f_i.

Series (as ours30, T31): phi_hat = argmax sum_i l_i on the whole line, curv = -sum l_i'', J = sum l_i'^2,
eta = min(1, curv / J) (eta = 1 if J = 0 or curv <= 0).

Law: p(m) ~ pi(m) prod_i f_i(m)^eta on m > 0. Point: the law's median.
"""
import math

import numpy as np

from met.full import geometric_mean, _best_line, _h_terms, _log_prior, _quantiles
from met.ours31 import COARSE_HALF, COARSE_NODES, FINE_NODES, SUPPORT_NATS, _per_reading


def _l_terms(states, c, phi, order):
    """Per-reading l_i and its first two phi-derivatives on an array phi."""
    per = []
    for r, ne, st, h_pos_max in states:
        a = 0.5 * (ne + r.k - 1)
        b = ne + st["Q_plus"]
        h, h1, h2 = _h_terms(r.pre(), c, phi, 2)
        u = 1.0 + (h_pos_max - h) / b
        l0 = -a * np.log(u)
        if order == 0:
            per.append((l0,))
            continue
        l1 = a * h1 / (b * u)
        l2 = a * (h2 / (b * u) + h1 ** 2 / (b * u) ** 2)
        per.append((l0, l1, l2))
    return per


def bartlett_eta(states, c):
    if not states:
        return 1.0, 0.0, 0.0, 0.0

    def fun(phi, order):
        per = _l_terms(states, c, phi, 2)
        tot = [sum(p[j] for p in per) for j in range(3)]
        return tuple(tot) if order else tot[0]

    phi_hat = _best_line(fun)
    per = _l_terms(states, c, np.array([phi_hat]), 2)
    d1 = np.array([p[1][0] for p in per])
    d2 = np.array([p[2][0] for p in per])
    curv, j = float(-d2.sum()), float((d1 ** 2).sum())
    if j <= 1e-300 or curv <= 0:
        return 1.0, curv, j, phi_hat
    return float(min(1.0, curv / j)), curv, j, phi_hat


def ours32(readings, c, nu_eff, *, calibrate=True, grid=None):
    """The ours32 probe on one series. As ours31 (same arguments), with eta (T31) on the cut factors.
    calibrate=False gives eta = 1, which is ours31.

    Returns a dict: point, median, interval (95%, eta applied), eta, law, lm, Q_plus, nu_eff."""
    if not c > 0:
        raise ValueError("the starting law's centre c must be positive")
    readings = list(readings)
    if np.ndim(nu_eff) == 0:
        nu_eff = [float(nu_eff)] * len(readings)
    nu_eff = [float(x) for x in nu_eff]
    if len(nu_eff) != len(readings) or not all(x > 0 and math.isfinite(x) for x in nu_eff):
        raise ValueError("nu_eff: one positive finite value per reading")
    states = _per_reading(readings, c, nu_eff)
    eta = bartlett_eta(states, c)[0] if calibrate else 1.0

    def log_law(lm):
        if not states:
            return _log_prior(lm, c)
        phi = np.arctan(np.exp(lm - math.log(c)))
        per = _l_terms(states, c, phi, 0)
        total = sum(p[0] for p in per)
        return _log_prior(lm, c) + eta * np.minimum(total, 0.0)

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
    (lo, med, hi), p = _quantiles(lm, lp, [0.025, 0.5, 0.975])
    return dict(point=float(med), median=float(med), geometric_mean=geometric_mean(lm, p), interval=(float(lo), float(hi)), eta=eta,
                law=p, lm=lm, Q_plus=[s[2]["Q_plus"] for s in states], nu_eff=nu_eff)

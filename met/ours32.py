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

Numerically the within-span residual K = E - H is evaluated directly, and l_i is stored as
-(nu_eff_i + k_i - 1)/2 * log1p(K_i/nu_eff_i), dropping an m-independent constant. The law,
scores and curvature are unchanged. The series maximum uses refined finite search seeds;
this is a numerical search, not a proof of the global maximum for every multimodal input.
"""
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from met.full import SCAN, geometric_mean, _h_terms, _log_prior, _quantiles
from met.ours31 import COARSE_HALF, COARSE_NODES, FINE_NODES, SUPPORT_NATS, _per_reading


def _residual_energy(reading, c, phi):
    """Within-span residual, computed directly rather than as E - H.

    The fixed declared span is whitened once. QR of its mass subspace then
    gives the residual without subtracting two signal-sized squared norms.
    """
    pre = reading.pre()
    if "ours32_residual" not in pre:
        n = reading.subspace
        zero = np.zeros_like(n)
        span = np.block([[n, zero], [zero, n]])
        chol = np.linalg.cholesky(reading.cov)
        white_span = np.linalg.solve(chol, span)
        basis, _ = np.linalg.qr(white_span, mode="reduced")
        y = np.linalg.solve(chol, np.r_[reading.force, reading.acceleration])
        pre["ours32_residual"] = (basis.T @ y, basis.T @ white_span)
    y, span = pre["ours32_residual"]
    k = reading.k
    force, acceleration = c * np.sin(phi), np.cos(phi)
    scale = np.hypot(force, acceleration)
    candidate = ((force / scale)[:, None, None] * span[None, :, :k]
                 + (acceleration / scale)[:, None, None] * span[None, :, k:])
    basis, _ = np.linalg.qr(candidate, mode="reduced")
    coordinates = np.einsum("nij,i->nj", basis, y)
    residual = y[None, :] - np.einsum("nij,nj->ni", basis, coordinates)
    return np.einsum("ni,ni->n", residual, residual)


def _l_terms(states, c, phi, order):
    """Per-reading log factors, up to m-independent constants, and derivatives."""
    per = []
    for r, ne, _, _ in states:
        a = 0.5 * (ne + r.k - 1)
        residual = _residual_energy(r, c, phi)
        # l differs from the peak-normalised cut by an m-independent constant:
        # -a log((ne + K)/(ne + Q+)) = const - a log1p(K/ne).
        # Keeping that constant out avoids cancellation on negative best lines.
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


def bartlett_eta(states, c, terms=None):
    """eta (T31) on the per-reading terms. `terms(states, c, phi, order)` defaults to this probe's
    _l_terms; met.ours33 passes its own (a different exponent). Default behaviour unchanged."""
    if not states:
        return 1.0, 0.0, 0.0, 0.0
    if terms is None:
        terms = _l_terms

    def fun(phi, order):
        per = terms(states, c, phi, 2 if order else 0)
        tot = [sum(p[j] for p in per) for j in range(3 if order else 1)]
        return tuple(tot) if order else tot[0]

    # Student factors can have maxima much narrower than the initial scan.
    # A fixed scan-sized gradient step can oscillate in their convex tails.
    # Add own lines to the scan before locating local maxima. This captures
    # narrow peaks without doing an entire series fit for every reading.
    step = math.pi / SCAN
    grid = -math.pi / 2 + step * (np.arange(SCAN) + 0.5)
    own_lines = [float(s[2]["phi"]) for s in states if s[2]["E"] > 0]
    own_lines = (np.asarray(own_lines) + math.pi / 2) % math.pi - math.pi / 2
    grid = np.unique(np.r_[grid, own_lines])
    values = fun(grid, 0)
    left, right = np.roll(values, 1), np.roll(values, -1)
    maxima = np.flatnonzero((values >= left) & (values >= right)
                           & ((values > left) | (values > right)))
    if not maxima.size:
        maxima = np.array([np.argmax(values)])

    def value(x):
        return float(fun(np.array([x]), 0)[0])

    def score(x):
        return float(fun(np.array([x]), 2)[1][0])

    candidates = []
    for index in maxima:
        seed = float(grid[index])
        low = float(grid[index - 1]) if index else float(grid[-1] - math.pi)
        high = float(grid[index + 1]) if index + 1 < len(grid) else float(grid[0] + math.pi)
        fitted = minimize_scalar(lambda x: -value(x), bounds=(low, high), method="bounded",
                                 options={"xatol": 1e-14})
        candidates.extend((seed, float(fitted.x)))
        if score(low) > 0 and score(high) < 0:
            candidates.append(brentq(score, low, high, xtol=1e-14, rtol=1e-14))
    phi_hat = max(candidates, key=value)
    # Values at nearby candidates can tie at floating-point precision while
    # their scores still differ enough to perturb H/J. Refine stationarity
    # inside the selected peak, where bounded search has already found it.
    for _ in range(8):
        _, first, second = (float(x[0]) for x in fun(np.array([phi_hat]), 2))
        if second >= 0:
            break
        move = -first / second
        if abs(move) > step:
            break
        phi_hat += move
        if abs(move) < 1e-15:
            break
    phi_hat = (phi_hat + math.pi / 2) % math.pi - math.pi / 2
    per = terms(states, c, np.array([phi_hat]), 2)
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
        return _log_prior(lm, c) + eta * total

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

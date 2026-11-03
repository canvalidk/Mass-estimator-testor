"""The forward model: hold the truth still and ask what the noise does to every measure.

The equations (ours29 to ours32) map readings to a law of mass. This maps the other way. Given each
reading's true pair (F*_i, a*_i), its declared covariance Sigma_i and subspace N_i (L18: the declared
Sigma is the true covariance), and optionally the law of the noise size (T37, G20: nu_eff), it draws

    y_i = (F*_i; a*_i) + sqrt(v_i) L_i z_i,     z_i ~ N(0, I_2d),  L_i L_i^T = Sigma_i,
    v_i = 1                                     (L13: the fixed-size Gaussian), or
    v_i ~ InvGamma(nu_eff/2, nu_eff/2)          (E[1/v] = 1: practice's normalisation, T37 G20)

and returns, draw by draw, the noise measures the equations read, computed as met.full computes them.

  per reading (draws, N)
    E             whitened energy of the reading in the span (T36)
    H_true        explained energy at the true mass m0
    K_true        residual at the true mass, E - H(m0) (T36 (4))
    Q             misfit on the whole projective line (ours29; nu = k - 1)
    Q_plus        misfit at the best positive mass (ours30 to ours32; kappa = k - 1/2)
    g29, g30      min(1, nu/Q), min(1, kappa/Q_plus) (L21, T30)
    w31           ours31/32's local weight (nu_d + k - 1)/(nu_d + Q_plus): its cut factor's curvature at
                  the reading's best line over the Gaussian's, = E[1/v | Q_plus] under ours31's size
                  posterior (T37 (4)); nu_d is the equation's declared nu_eff (`nu_declared`)
    m_own         the reading's own best line, signed (the whole projective line)
    v             the drawn noise size
    latent, score, misalignment     the noise vector's three channels (below)
  per series (draws,)
    K_total       sum of K_true (T36 (4))
    eta_flat, eta29, eta30          eta = min(1, curv/J) with g = 1, g29, g30 (T31)
    eta32         ours32's eta, on the cut factors (nu_declared given)
    curv30, J30   eta30's curvature and squared-score sum
    m_series      the series' best line under ours30's weights, signed

The noise vector's channels. In a reading's noise metric (T36 (1)) its noise, inside the span, is a
2k-vector. With the truth on the line L(m0) it splits orthogonally into

    latent        k dims along L(m0): what re-estimates the true acceleration; H moves with it
    score         1 dim, signed, along the direction in which L(m) moves the signal as m grows: the
                  efficient score for m with the latent pair profiled out, in units of
                  sqrt(information); positive means the noise says "heavier"
    misalignment  k - 1 dims: the codirectionality conditions (L8); what Q reads at strong signal

Given v the three are exactly N(0, v I), independent, at every signal size except zero (where the
score direction does not exist), and K_true = score^2 + misalignment. They are the noise alone; the
fitted measures (Q, Q_plus, eta) are signal plus noise refitted, and weak signal enters only there.
`latent` and `misalignment` are returned as squared lengths, `score` as the signed component.

What the distributions depend on (L9, L16). For isotropic uncorrelated channels on the full subspace,
every measure's law depends on the truth only through, per reading: k, the push (the signal's length
in noise units, sqrt(|F*|^2/sigma_F^2 + |a*|^2/sigma_a^2)), the angle theta0 = atan(m0 / s) of the
true line in the whitened plane (s = sigma_F / sigma_a; theta0 = 45 deg is m0 = s), and nu_eff; and
for the series measures, N. The direction of a* drops out, also across readings. `invariants` returns
these. With correlated or anisotropic Sigma the whitened geometry carries more (O17, G5).

Standard mathematics (GOAL.md rule 6). E is noncentral chi^2_2k(push^2); H_true is noncentral
chi^2_k(push^2) and K_true central chi^2_k, independent (Cochran: complementary projections); with the
inverse-gamma size, K_true/k ~ F(k, nu_eff) and score ~ Student t_nu_eff (the quadratic forms of a
multivariate t). Q is the smallest eigenvalue of the 2x2 Gram of a k x 2 Gaussian matrix with a
rank-one mean: the orthogonal-regression (total least squares) residual, a noncentral Wishart
eigenvalue (Anderson 1951; Gleser 1981). At zero push its law is the smallest eigenvalue of
W_2(k, I); as the push grows it tends to chi^2_(k-1). Q_plus is the same residual with the slope
restricted to a cone: at zero push it is the half-and-half mixture of Q and min(chi^2_k, chi^2_k), and
at strong push on the boundary (one channel without signal) the chi-bar-squared mixture
1/2 chi^2_(k-1) + 1/2 chi^2_k (Chernoff 1954; Self and Liang 1987): T38's k - 1/2. eta's curv/J is the
ratio of the observed information to the outer product of the scores at the fitted line: the
information-matrix equality (Bartlett's second identity), White's (1982) test statistic.

Two engines, the same numbers: `isotropic` (closed forms on the 2x2 Gram, vectorised, for series of
isotropic uncorrelated readings on the full subspace sharing one s) and `general` (met.full per
draw, any Sigma and subspace). tests/test_forward.py checks they agree.

This is the world, not the definition: it holds the true accelerations fixed, which the definition
may not (L10: no size law). Its use is to say what the measures do when the truth is known: whether
a numerator such as kappa is the expectation it claims, how often a cap tempers a reading the declared
noise explains, what eta does at finite N.
"""
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.special import gamma as gamma_fn
from scipy.special import gammaincc

from .full import (Reading, _h_terms, bartlett_eta, exact_mass, instrument_ratio, isotropic,
                   misfit_weight, positivity_numerator, reading_state)

PER_READING = ("E", "H_true", "K_true", "Q", "Q_plus", "g29", "g30", "w31", "m_own", "v",
               "latent", "score", "misalignment")
PER_SERIES = ("K_total", "eta_flat", "eta29", "eta30", "curv30", "J30", "m_series", "eta32")
SCAN32 = 1024
NEWTON32 = 30


# ---------------------------------------------------------------- truths


def on_line(mass, accelerations, force_sd, acceleration_sd):
    """Truths on one mass line: F*_i = mass a*_i, isotropic uncorrelated channels.
    accelerations: a list of length-d vectors (zero allowed)."""
    if not mass > 0:
        raise ValueError("the true mass must be positive (L1)")
    return [isotropic(mass * np.asarray(a, dtype=float), a, force_sd, acceleration_sd)
            for a in accelerations]


def from_invariants(k, push, angle, *, readings=1, force_sd=1.0, acceleration_sd=1.0):
    """N identical truths with the given push (signal length in noise units) and angle theta0 (radians,
    in [0, pi/2]; m0 = s tan theta0), acceleration along the first axis. Angle 0 is m0 = 0 and pi/2
    is m0 = infinity: the positive half-line's two boundaries, allowed here as worlds (T38's suprema)
    though not as masses (L1)."""
    if k not in (1, 2, 3) or not push >= 0 or not 0 <= angle <= math.pi / 2:
        raise ValueError("k in 1..3, push >= 0, angle in [0, pi/2]")
    e1 = np.zeros(k)
    e1[0] = 1.0
    force = push * math.sin(angle) * force_sd * e1
    acc = push * math.cos(angle) * acceleration_sd * e1
    return [isotropic(force, acc, force_sd, acceleration_sd) for _ in range(int(readings))]


def true_mass(truths):
    """The common mass of the truths if every true pair is exactly codirectional with one positive
    ratio (zero pairs fit every mass); None otherwise."""
    masses = []
    for t in truths:
        kind, m = exact_mass(t.force, t.acceleration)
        if kind == "none":
            return None
        if kind == "one":
            masses.append(m)
    if not masses:
        return None
    m0 = masses[0]
    return m0 if all(abs(m - m0) <= 1e-12 * m0 for m in masses) else None


def invariants(truth, mass=None):
    """What a reading's measures depend on (isotropic uncorrelated, full subspace): k, the channel
    strengths R = |F*|/sigma_F and b = |a*|/sigma_a, the push sqrt(R^2 + b^2), cos psi between F* and a*
    (1 on the line), the angle theta0 = atan(m0/s) of the true line if m0 is given or on the line.
    None if the reading is not isotropic on the full subspace."""
    s = instrument_ratio(truth)
    if s is None or truth.k != truth.dimension:
        return None
    d = truth.dimension
    sf, sa = math.sqrt(truth.cov[0, 0]), math.sqrt(truth.cov[d, d])
    r, b = float(np.linalg.norm(truth.force)) / sf, float(np.linalg.norm(truth.acceleration)) / sa
    cos_psi = (float(truth.force @ truth.acceleration) / (r * sf * b * sa)) if r > 0 and b > 0 else float("nan")
    m0 = mass
    if m0 is None:
        kind, value = exact_mass(truth.force, truth.acceleration)
        m0 = value if kind == "one" else None
    angle = math.atan(m0 / s) if m0 is not None else (math.atan2(r, b) if r > 0 or b > 0 else float("nan"))
    return dict(k=truth.k, R=r, b=b, push=math.hypot(r, b), cos_psi=cos_psi, angle=angle,
                angle_deg=math.degrees(angle) if angle == angle else float("nan"), s=s)


# ---------------------------------------------------------------- noise


def draw_noise(truths, draws, rng, nu_eff=None):
    """eps (draws, N, 2d) with eps = sqrt(v) L z, and v (draws, N). z is drawn first, all at once, then
    v: so the same seed gives the same z whether or not a size law is declared."""
    n, d = len(truths), truths[0].dimension
    if any(t.dimension != d for t in truths):
        raise ValueError("every truth must have the same dimension")
    z = rng.standard_normal((draws, n, 2 * d))
    if nu_eff is None or math.isinf(nu_eff):
        v = np.ones((draws, n))
    else:
        if not nu_eff > 0:
            raise ValueError("nu_eff must be positive")
        v = stats.invgamma.rvs(nu_eff / 2, scale=nu_eff / 2, size=(draws, n), random_state=rng)
    chol = np.stack([np.linalg.cholesky(t.cov) for t in truths])          # (N, 2d, 2d)
    eps = np.einsum("nij,rnj->rni", chol, z) * np.sqrt(v)[..., None]
    return eps, v, z


# ---------------------------------------------------------------- the noise vector's channels


def _channel_maps(truth, m0):
    """For one truth on the line at m0: the map from a noise vector (2d) to its whitened coordinates
    in the span (2k), the orthonormal bases of L(m0) and its complement there, and the unit score
    direction in the complement (None if the truth has no signal), with the reading's information
    about m at the truth (|P_perp C'(m0) alpha*|^2: the efficient information, latent profiled out)."""
    d, k, n = truth.dimension, truth.k, truth.subspace
    si = np.linalg.inv(truth.cov)
    span = np.block([[n, np.zeros_like(n)], [np.zeros_like(n), n]])           # (2d, 2k)
    gram = span.T @ si @ span
    upper = np.linalg.cholesky(gram).T                                        # gram = upper^T upper
    to_beta = np.linalg.solve(gram, span.T @ si)                              # (2k, 2d)
    whiten = upper @ to_beta
    line = upper @ np.vstack([m0 * np.eye(k), np.eye(k)])                     # (2k, k)
    full_q, _ = np.linalg.qr(line, mode="complete")
    basis_l, basis_perp = full_q[:, :k], full_q[:, k:]
    beta = to_beta @ np.concatenate([truth.force, truth.acceleration])
    alpha = beta[k:]
    t = basis_perp.T @ upper @ np.concatenate([alpha, np.zeros(k)])
    info = float(t @ t)
    direction = t / math.sqrt(info) if info > 0 else None
    return whiten, basis_l, basis_perp, direction, info


def channels(truths, eps, m0):
    """latent (squared), score (signed), misalignment (squared), each (draws, N), and the information
    about m of each reading at the truth (N,). NaN where the reading has no signal (score undefined)."""
    draws, n = eps.shape[:2]
    latent = np.full((draws, n), np.nan)
    score = np.full((draws, n), np.nan)
    mis = np.full((draws, n), np.nan)
    info = np.zeros(n)
    for i, t in enumerate(truths):
        whiten, basis_l, basis_perp, direction, info[i] = _channel_maps(t, m0)
        w = eps[:, i, :] @ whiten.T
        latent[:, i] = ((w @ basis_l) ** 2).sum(axis=1)
        perp = w @ basis_perp
        if direction is not None:
            score[:, i] = perp @ direction
            mis[:, i] = np.maximum((perp ** 2).sum(axis=1) - score[:, i] ** 2, 0.0)
    return latent, score, mis, info


# ---------------------------------------------------------------- engines


def _common_isotropic(truths):
    """s if every truth is isotropic uncorrelated on the full subspace with one common s; else None."""
    ss = [instrument_ratio(t) for t in truths]
    if any(s is None for s in ss) or any(t.k != t.dimension for t in truths):
        return None
    return ss[0] if all(abs(s - ss[0]) <= 1e-12 * ss[0] for s in ss) else None


def _isotropic_engine(truths, eps, m0, s):
    """Closed forms on each reading's 2x2 Gram in the whitened plane (T32, T38's setup)."""
    d = truths[0].dimension
    sf = np.array([math.sqrt(t.cov[0, 0]) for t in truths])
    sa = np.array([math.sqrt(t.cov[d, d]) for t in truths])
    fstar = np.stack([t.force for t in truths]) / sf[:, None]                 # (N, d)
    astar = np.stack([t.acceleration for t in truths]) / sa[:, None]
    x = fstar[None] + eps[..., :d] / sf[None, :, None]                        # (draws, N, d)
    y = astar[None] + eps[..., d:] / sa[None, :, None]
    a_, b_, c_ = (x * x).sum(-1), (y * y).sum(-1), (x * y).sum(-1)
    e = a_ + b_
    lmax = 0.5 * e + np.sqrt((0.5 * (a_ - b_)) ** 2 + c_ ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(lmax > 0, (a_ * b_ - c_ ** 2) / lmax, 0.0)
    q = np.maximum(q, 0.0)
    q_plus = np.where(c_ > 0, q, np.minimum(a_, b_))
    q_plus = np.maximum(q_plus, q)
    out = dict(E=e, Q=q, Q_plus=q_plus,
               m_own=s * np.tan(0.5 * np.arctan2(2 * c_, b_ - a_)))
    if m0 is not None:
        th = math.atan(m0 / s)
        h = a_ * math.sin(th) ** 2 + 2 * c_ * math.sin(th) * math.cos(th) + b_ * math.cos(th) ** 2
        out["H_true"] = h
        out["K_true"] = np.maximum(e - h, 0.0)
    out["_gram"] = (a_, b_, c_)
    out["_observed"] = (x, y)
    return out


def _isotropic_eta(gram, g, s, c):
    """eta = min(1, curv/J) at the series' best line on the whole projective line (T31), from the
    g-weighted Gram, in the whitened angle theta (m = s tan theta). eta is the same in every
    parametrisation (T31); curv and J are returned in met.full's, m = c tan phi, by the chain rule
    at the stationary point (the first-derivative term vanishes there)."""
    a_, b_, c_ = gram
    ag, bg, cg = (g * a_).sum(-1), (g * b_).sum(-1), (g * c_).sum(-1)
    theta = 0.5 * np.arctan2(2 * cg, bg - ag)                                 # (draws,)
    s2, c2 = np.sin(2 * theta)[:, None], np.cos(2 * theta)[:, None]
    d1 = 0.5 * g * ((a_ - b_) * s2 + 2 * c_ * c2)
    d2 = 0.5 * g * (2 * (a_ - b_) * c2 - 4 * c_ * s2)
    curv, j = -d2.sum(-1), (d1 ** 2).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        eta = np.where((j <= 1e-300) | (curv <= 0), 1.0, np.minimum(1.0, curv / j))
    phi = np.arctan2(s * np.sin(theta), c * np.cos(theta))
    r = c / s
    dtheta_dphi = r / (np.cos(phi) ** 2 + r * r * np.sin(phi) ** 2)
    return eta, curv * dtheta_dphi ** 2, j * dtheta_dphi ** 2, s * np.tan(theta)


def _isotropic_eta32(gram, q_plus, k, nu_d, *, observations=None):
    """Vectorised ours32 eta, with a robust fallback for unresolved narrow peaks.

    The scan uses Grams; Newton/final residuals use the observed whitened vectors
    when supplied. Q+ contributes only a mass-independent factor and is omitted.
    Failed stationarity or sub-grid reading peaks use the production finite search
    on those observations, without passing any investigator truth to the fit.
    """
    a_, b_, c_ = gram
    e = a_ + b_
    aa = 0.5 * (nu_d + k - 1)

    def terms(theta, order):                     # theta (draws, M) -> per-reading l (draws, N, M)
        s2, c2 = np.sin(2 * theta)[:, None, :], np.cos(2 * theta)[:, None, :]
        A, B, C = a_[..., None], b_[..., None], c_[..., None]
        residual = np.maximum(0.5 * (A + B) + 0.5 * (A - B) * c2 - C * s2, 0.0)
        if observations is not None and order:
            x, y = observations
            sine, cosine = np.sin(theta)[:, None, :, None], np.cos(theta)[:, None, :, None]
            off = x[:, :, None, :] * cosine - y[:, :, None, :] * sine
            along = x[:, :, None, :] * sine + y[:, :, None, :] * cosine
            residual = np.sum(off**2, axis=-1)
            h1 = 2 * np.sum(off * along, axis=-1)
            h2 = 2 * (residual - np.sum(along**2, axis=-1))
        l0 = -aa * np.log1p(residual / nu_d)
        if order == 0:
            return l0
        if observations is None:
            h1 = (A - B) * s2 + 2 * C * c2
            h2 = 2 * (A - B) * c2 - 4 * C * s2
        denominator = nu_d + residual
        return l0, aa * h1 / denominator, aa * (h2 / denominator + (h1 / denominator) ** 2)

    draws = a_.shape[0]
    grid = -math.pi / 2 + math.pi * (np.arange(SCAN32) + 0.5) / SCAN32
    total = terms(np.broadcast_to(grid, (draws, SCAN32)), 0).sum(axis=1)
    theta = grid[np.argmax(total, axis=1)].astype(float)
    step = math.pi / SCAN32
    for _ in range(NEWTON32):
        _, d1, d2 = (x.sum(axis=1)[:, 0] for x in terms(theta[:, None], 2))
        move = np.where(d2 < 0, np.clip(-d1 / np.where(d2 < 0, d2, -1.0), -step, step),
                        step * np.sign(d1))
        theta = theta + move
        if np.all(np.abs(move) < 1e-15):
            break
    _, d1, d2 = (x[..., 0] for x in terms(theta[:, None], 2))
    curv, j = -d2.sum(-1), (d1 ** 2).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        eta = np.where((j <= 1e-300) | (curv <= 0), 1.0, np.minimum(1.0, curv / j))
    if observations is not None:
        narrow = np.any(e * step**2 > 0.25 * nu_d, axis=1)
        score_scale = 1 + np.sqrt(j) + np.sqrt(np.abs(curv))
        unresolved = (np.abs(d1.sum(-1)) > 1e-8 * score_scale)
        unresolved |= ~np.isfinite(eta) | ~np.isfinite(curv) | ~np.isfinite(j)
        unresolved |= (curv <= 0) & (j > 1e-20)
        indices = np.flatnonzero(narrow | unresolved)
        if indices.size:
            from .ours31 import _per_reading
            from .ours32 import bartlett_eta as stable_eta

            x, y = observations
            for index in indices:
                readings = [isotropic(f, a, 1.0, 1.0) for f, a in zip(x[index], y[index])]
                eta[index] = stable_eta(_per_reading(readings, 1.0, [nu_d] * len(readings)), 1.0)[0]
    return eta


def _general_engine(truths, eps, m0, c, series=True, nu_d=None):
    """met.full per draw: any covariance, any subspace."""
    draws, n = eps.shape[:2]
    d = truths[0].dimension
    keys = ("E", "Q", "Q_plus", "m_own") + (("H_true", "K_true") if m0 is not None else ())
    out = {key: np.empty((draws, n)) for key in keys}
    out.update(g29=np.empty((draws, n)), g30=np.empty((draws, n)))
    ser = {key: np.empty(draws) for key in ("eta_flat", "eta29", "eta30", "curv30", "J30", "m_series")
           + (("eta32",) if nu_d is not None else ())}
    phi0 = math.atan(m0 / c) if m0 is not None else None
    for r in range(draws):
        readings = [Reading(t.force + eps[r, i, :d], t.acceleration + eps[r, i, d:], t.cov, t.subspace)
                    for i, t in enumerate(truths)]
        g29, g30 = [], []
        for i, rd in enumerate(readings):
            st = reading_state(rd, c)
            out["E"][r, i], out["Q"][r, i], out["Q_plus"][r, i] = st["E"], st["Q"], st["Q_plus"]
            out["m_own"][r, i] = c * math.tan(st["phi"])
            g29.append(misfit_weight(st["Q"], st["nu"]))
            g30.append(misfit_weight(st["Q_plus"], st["kappa"]))
            if phi0 is not None:
                h = float(_h_terms(rd.pre(), c, np.array([phi0]))[0])
                out["H_true"][r, i], out["K_true"][r, i] = h, max(st["E"] - h, 0.0)
        out["g29"][r], out["g30"][r] = g29, g30
        if series:
            e30, curv, j, phi = bartlett_eta(readings, c, g30)
            ser["eta30"][r], ser["curv30"][r], ser["J30"][r] = e30, curv, j
            ser["m_series"][r] = c * math.tan(phi)
            ser["eta29"][r] = bartlett_eta(readings, c, g29)[0]
            ser["eta_flat"][r] = bartlett_eta(readings, c, [1.0] * n)[0]
            if nu_d is not None:
                from .ours31 import _per_reading
                from .ours32 import bartlett_eta as eta32
                ser["eta32"][r] = eta32(_per_reading(readings, c, [nu_d] * n), c)[0]
    if not series:
        ser = {}
    return out, ser


def forward(truths, *, draws, seed, mass=None, nu_eff=None, nu_declared=None, c=None, engine="auto",
            series=True, keep_noise=False):
    """Hold the truths still, draw the noise `draws` times, and return every measure's draws.

    truths: a list of met.full.Reading holding the TRUE force and acceleration and the covariance
      (the declared one, believed: L18) and subspace. `on_line` and `from_invariants` build them.
    mass: where K, H and the channels are read; default the truths' common exact mass (None if they
      have none: then only E, Q, Q_plus, g, eta are returned). Give it for an all-zero truth.
    nu_eff: the WORLD's noise size: None or inf for the fixed-size Gaussian (L13); a number for the
      inverse-gamma size (T37).
    nu_declared: the EQUATION's nu_eff, for ours31/32's measures (w31, eta32); default the world's when
      finite, else none (w31 and eta32 are then not computed). Give it to run ours31/32 in a world
      whose noise differs from their declaration (e.g. the fixed-size Gaussian).
    c: the parametrisation of the line, m = c tan(phi); default s. Only curv30 and J30 depend on it
      (as the square of dphi); every other measure, eta included, does not (tested).
    engine: "auto" (isotropic when the series allows it), "isotropic" or "general" (met.full).
    series: compute the series measures (eta, m_series). keep_noise: return eps and z too.

    Returns a dict: per-reading arrays (draws, N), per-series arrays (draws,), and `truth`: m0, nu_eff,
    the information about m of each reading at the truth, its invariants, the engine used."""
    truths = list(truths)
    if not truths:
        raise ValueError("forward needs at least one truth")
    if type(draws) is not int or draws < 1:
        raise ValueError("draws must be a positive integer")
    rng = np.random.default_rng(seed)
    m0 = true_mass(truths) if mass is None else float(mass)
    if m0 is not None and not m0 > 0:
        raise ValueError("mass must be positive")
    eps, v, z = draw_noise(truths, draws, rng, nu_eff)
    nu_d = nu_declared
    if nu_d is None and nu_eff is not None and not math.isinf(nu_eff):
        nu_d = float(nu_eff)
    if nu_d is not None and not (nu_d > 0 and math.isfinite(nu_d)):
        raise ValueError("nu_declared must be positive and finite (ours31/32 take no infinite nu_eff)")
    s = _common_isotropic(truths)
    if c is None:
        c = s if s is not None else 1.0
    c = float(c)
    if not c > 0:
        raise ValueError("c must be positive")
    if engine == "auto":
        engine = "isotropic" if s is not None else "general"
    if engine == "isotropic":
        if s is None:
            raise ValueError("the isotropic engine needs isotropic uncorrelated readings on the full "
                             "subspace with one common s")
        out = _isotropic_engine(truths, eps, m0, s)
        gram = out.pop("_gram")
        observations = out.pop("_observed")
        k = truths[0].k
        g29 = np.where(out["Q"] > 0, np.minimum(1.0, (k - 1) / np.where(out["Q"] > 0, out["Q"], 1.0)), 1.0)
        g29 = g29 if k > 1 else np.ones_like(g29)
        kappa = positivity_numerator(k)
        g30 = np.where(out["Q_plus"] > 0, np.minimum(1.0, kappa / np.where(out["Q_plus"] > 0, out["Q_plus"], 1.0)), 1.0)
        out.update(g29=g29, g30=g30)
        ser = {}
        if series:
            e30, curv, j, m_series = _isotropic_eta(gram, g30, s, c)
            ser = dict(eta30=e30, curv30=curv, J30=j, m_series=m_series,
                       eta29=_isotropic_eta(gram, g29, s, c)[0],
                       eta_flat=_isotropic_eta(gram, np.ones_like(g30), s, c)[0])
            if nu_d is not None:
                ser["eta32"] = _isotropic_eta32(gram, out["Q_plus"], k, nu_d,
                                               observations=observations)
    elif engine == "general":
        out, ser = _general_engine(truths, eps, m0, c, series, nu_d)
    else:
        raise ValueError(f"unknown engine {engine!r}")
    out["v"] = v
    if nu_d is not None:
        kk = np.array([tr.k for tr in truths], dtype=float)
        out["w31"] = (nu_d + kk - 1) / (nu_d + out["Q_plus"])
    info = np.full(len(truths), np.nan)
    if m0 is not None:
        out["latent"], out["score"], out["misalignment"], info = channels(truths, eps, m0)
        out["K_total"] = out["K_true"].sum(axis=1)
    out.update(ser)
    out["truth"] = dict(mass=m0, nu_eff=nu_eff, nu_declared=nu_d, engine=engine, draws=draws, seed=seed, c=c,
                        information=info, k=[t.k for t in truths],
                        invariants=[invariants(t, m0) for t in truths])
    if keep_noise:
        out["eps"], out["z"] = eps, z
    return out


# ---------------------------------------------------------------- reference laws


class _Law:
    """A one-dimensional law given by a cdf on [0, inf) and a mean (for the laws scipy lacks)."""

    def __init__(self, name, cdf, mean):
        self.name, self._cdf, self._mean = name, cdf, mean

    def cdf(self, x):
        return self._cdf(np.asarray(x, dtype=float))

    def mean(self):
        return self._mean

    def __repr__(self):
        return f"<{self.name}, mean {self._mean:.6g}>"


def _wishart_min_density(x, k):
    """Unnormalised density of the smaller eigenvalue of W_2(k, I), k >= 2: the joint eigenvalue
    density (l1 l2)^a e^{-(l1+l2)/2} (l1 - l2), a = (k-3)/2, integrated over l1 > x in closed form."""
    a = 0.5 * (k - 3)
    x = np.asarray(x, dtype=float)
    up2 = 2 ** (a + 2) * gammaincc(a + 2, x / 2) * gamma_fn(a + 2)
    up1 = 2 ** (a + 1) * gammaincc(a + 1, x / 2) * gamma_fn(a + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(x > 0, x ** a * np.exp(-x / 2) * (up2 - x * up1), 0.0)


def _tabulated(name, density, upper=200.0, nodes=4000):
    """A law from a density on (0, inf): cdf tabulated by quadrature between nodes, mean by quadrature."""
    total = quad(density, 0, np.inf, limit=200)[0]
    mean = quad(lambda t: t * density(t), 0, np.inf, limit=200)[0] / total
    grid = np.concatenate([[0.0], np.geomspace(1e-8, upper, nodes)])
    pieces = [quad(density, grid[i], grid[i + 1], limit=100)[0] for i in range(len(grid) - 1)]
    cdf_grid = np.concatenate([[0.0], np.cumsum(pieces)]) / total
    return _Law(name, lambda x: np.interp(x, grid, cdf_grid, right=1.0), mean)


@lru_cache(maxsize=None)
def null_misfit(k):
    """The law of Q at zero push: the smaller eigenvalue of W_2(k, I) (k >= 2); Q = 0 for k = 1.
    Exact means: 2 - pi/2 at k = 2, 1 at k = 3 (Exp(1))."""
    if k == 1:
        return _Law("Q at zero push, k=1: 0", lambda x: (x >= 0).astype(float), 0.0)
    return _tabulated(f"smallest eigenvalue of W_2({k}, I)",
                      lambda t: float(_wishart_min_density(np.array([t]), k)[0]))


@lru_cache(maxsize=None)
def null_positive_misfit(k):
    """The law of Q_plus at zero push: (1/2) Q's null law + (1/2) min(chi^2_k, chi^2_k). Exact because
    flipping the acceleration noise's sign keeps A, B and Q and swaps C > 0 with C < 0.
    Means: 1/2 - 1/pi, 3/2 - pi/4, 2 - 2/pi at k = 1, 2, 3 (tests/test_positivity_supremum.py NULL)."""
    q = null_misfit(k)
    chi = stats.chi2(k)
    mean_min = quad(lambda t: chi.sf(t) ** 2, 0, np.inf)[0]
    return _Law(f"Q_plus at zero push, k={k}",
                lambda x: 0.5 * q.cdf(x) + 0.5 * (1 - chi.sf(x) ** 2), 0.5 * q.mean() + 0.5 * mean_min)


def boundary_positive_misfit(k):
    """The limit law of Q_plus as the push grows with one channel empty (theta0 = 0 or pi/2): the
    chi-bar-squared mixture (1/2) chi^2_(k-1) + (1/2) chi^2_k, mean k - 1/2 (T38's supremum)."""
    lo = stats.chi2(k - 1) if k > 1 else None
    hi = stats.chi2(k)
    return _Law(f"1/2 chi2_{k - 1} + 1/2 chi2_{k}",
                lambda x: 0.5 * (lo.cdf(x) if lo is not None else (x >= 0).astype(float)) + 0.5 * hi.cdf(x),
                k - 0.5)


def reference(truth, *, mass=None, nu_eff=None):
    """The exact laws of a reading's measures that standard results give, for this truth.
    Keys present only where exact: E, H_true (fixed size only), K_true, score, misalignment, latent;
    K_true is noncentral chi^2_k(delta) when the truth is off the line at `mass` (delta returned).
    Limits for Q and Q_plus are in null_misfit, null_positive_misfit, boundary_positive_misfit."""
    k = truth.k
    pre = truth.pre()
    push2 = pre["energy"]
    m0 = mass
    if m0 is None:
        kind, value = exact_mass(truth.force, truth.acceleration)
        m0 = value if kind == "one" else None
    out = dict(push=math.sqrt(push2))
    fixed = nu_eff is None or math.isinf(nu_eff)
    if fixed:
        out["E"] = stats.ncx2(2 * k, push2) if push2 > 0 else stats.chi2(2 * k)
    if m0 is None:
        return out
    s = 1.0
    h = float(_h_terms(pre, s, np.array([math.atan(m0 / s)]))[0])
    delta = max(push2 - h, 0.0)                                               # the signal's own residual
    out["delta"] = delta
    if fixed:
        out["K_true"] = stats.ncx2(k, delta) if delta > 1e-12 * max(push2, 1.0) else stats.chi2(k)
        out["H_true"] = stats.ncx2(k, h) if h > 0 else stats.chi2(k)
        out["latent"] = stats.chi2(k)
        if delta <= 1e-12 * max(push2, 1.0) and push2 > 0:
            out["score"] = stats.norm()
            if k > 1:
                out["misalignment"] = stats.chi2(k - 1)
        return out
    if delta <= 1e-12 * max(push2, 1.0):
        out["K_true"] = stats.f(k, nu_eff, scale=k)                          # K = v chi2_k, 1/v ~ chi2_nu/nu
        out["latent"] = stats.f(k, nu_eff, scale=k)
        if push2 > 0:
            out["score"] = stats.t(nu_eff)
            if k > 1:
                out["misalignment"] = stats.f(k - 1, nu_eff, scale=k - 1)
    return out


# ---------------------------------------------------------------- calibration of a weight


def calibration(result, key, bins=10):
    """Is a per-reading weight the precision the reading actually had? The world knows v; L21's cut
    reads g as E[1/v | the reading] (T37). Returns, over equal-count bins of the weight, the mean
    weight and the mean true precision 1/v in each, the overall means, and the ratio E[w]/E[1/v].
    A calibrated weight has the two columns equal."""
    w = np.asarray(result[key], dtype=float).ravel()
    p = 1.0 / np.asarray(result["v"], dtype=float).ravel()
    order = np.argsort(w, kind="stable")
    parts = np.array_split(order, bins)
    table = [(float(w[i].mean()), float(p[i].mean())) for i in parts if i.size]
    return dict(weight=float(w.mean()), precision=float(p.mean()), ratio=float(w.mean() / p.mean()),
                bins=table, worst_bin_gap=max(abs(a - b) for a, b in table))


# ---------------------------------------------------------------- summary


def summarize(result):
    """Means, SDs, quartiles and the tail fractions the equations care about, per measure, pooled over
    draws and readings (per-reading measures) or over draws (series measures)."""
    out = {}
    for key in PER_READING + PER_SERIES:
        if key not in result:
            continue
        a = np.asarray(result[key], dtype=float).ravel()
        a = a[np.isfinite(a)]
        if not a.size:
            continue
        q25, q50, q75 = np.quantile(a, [0.25, 0.5, 0.75])
        out[key] = dict(mean=float(a.mean()), sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                        se=float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else 0.0,
                        q25=float(q25), median=float(q50), q75=float(q75), n=int(a.size))
    k = result["truth"]["k"][0]
    out["tempered30"] = float(np.mean(result["g30"] < 1))      # P(Q_plus > kappa): the cap acts
    out["tempered29"] = float(np.mean(result["g29"] < 1)) if k > 1 else 0.0
    if "eta30" in result:
        out["eta30_below_1"] = float(np.mean(result["eta30"] < 1))
    if "m_series" in result:
        out["m_series_negative"] = float(np.mean(result["m_series"] < 0))
    return out


# ---------------------------------------------------------------- the preset kind `forward`

GRID_KEYS = ("dimension", "push", "angle_deg", "readings", "nu_eff")


def validate(data):
    """A `forward` preset: name, question, seed, draws, and a grid of the invariants."""
    if not isinstance(data, dict) or data.get("kind") != "forward":
        raise ValueError("a forward preset needs kind: forward")
    allowed = {"kind", "name", "question", "notes", "seed", "draws", "grid"}
    extra = set(data) - allowed
    if extra:
        raise ValueError(f"unknown forward setting(s): {sorted(extra)}")
    for key in ("name", "question"):
        if not isinstance(data.get(key), str) or not data[key].strip():
            raise ValueError(f"forward.{key} must be a nonempty string")
    if type(data.get("seed")) is not int or data["seed"] < 0:
        raise ValueError("forward.seed must be a nonnegative integer")
    if type(data.get("draws")) is not int or data["draws"] < 2:
        raise ValueError("forward.draws must be an integer >= 2")
    grid = data.get("grid")
    if not isinstance(grid, dict) or not set(GRID_KEYS) <= set(grid) or set(grid) - set(GRID_KEYS) - {"nu_declared"}:
        raise ValueError(f"forward.grid needs {list(GRID_KEYS)}, each a list, and optionally nu_declared")
    if "nu_declared" in grid and (not isinstance(grid["nu_declared"], list) or not grid["nu_declared"] or any(
            not (isinstance(x, (int, float)) and 0 < x < float("inf")) for x in grid["nu_declared"])):
        raise ValueError("nu_declared: positive finite numbers (ours31/32's declared nu_eff)")
    for key in GRID_KEYS:
        if not isinstance(grid[key], list) or not grid[key]:
            raise ValueError(f"forward.grid.{key} must be a nonempty list")
    if any(k not in (1, 2, 3) for k in grid["dimension"]):
        raise ValueError("dimension: 1, 2 or 3")
    if any(not (isinstance(x, (int, float)) and x >= 0) for x in grid["push"]):
        raise ValueError("push: nonnegative numbers")
    if any(not (isinstance(x, (int, float)) and 0 <= x <= 90) for x in grid["angle_deg"]):
        raise ValueError("angle_deg: in [0, 90]")
    if any(type(n) is not int or n < 1 for n in grid["readings"]):
        raise ValueError("readings: positive integers")
    if any(not (isinstance(x, (int, float)) and x > 0) for x in grid["nu_eff"]):
        raise ValueError("nu_eff: positive numbers (.inf for the fixed-size Gaussian)")
    cells = []
    keys = GRID_KEYS + (("nu_declared",) if "nu_declared" in grid else ())
    for index, combo in enumerate(itertools.product(*(grid[key] for key in keys))):
        cells.append(dict(id=f"f{index:03d}", **dict(zip(keys, combo))))
    return dict(kind="forward", raw=data, cells=cells)


def _cell_seed(seed, index):
    return int(np.random.SeedSequence([seed, index]).generate_state(1)[0])


def run_cell(cell, seed, draws):
    """One cell: N identical truths at the cell's invariants (sigma_F = sigma_a = 1, so s = 1 and
    m0 = tan theta0). Returns the summary row with the exact reference means beside the draws'."""
    k, push, angle = cell["dimension"], float(cell["push"]), math.radians(cell["angle_deg"])
    nu = float(cell["nu_eff"])
    truths = from_invariants(k, push, angle, readings=cell["readings"])
    on = 0 < angle < math.pi / 2
    m0 = math.tan(angle) if on else None
    if push == 0:
        m0 = math.tan(angle) if on else 1.0
    res = forward(truths, draws=draws, seed=seed, mass=m0, nu_eff=None if math.isinf(nu) else nu,
                  nu_declared=cell.get("nu_declared"))
    summ = summarize(res)
    row = dict(cell)
    row["mass"] = m0
    row["kappa"] = positivity_numerator(k)
    for key in ("E", "K_true", "score", "misalignment", "Q", "Q_plus", "g30", "g29", "w31",
                "eta30", "eta29", "eta_flat", "K_total"):
        if key in summ:
            row[f"{key}_mean"] = summ[key]["mean"]
            row[f"{key}_se"] = summ[key]["se"]
    row["score_var"] = summ["score"]["sd"] ** 2 if "score" in summ else float("nan")
    for key in ("tempered30", "tempered29", "eta30_below_1", "m_series_negative"):
        row[key] = summ.get(key, float("nan"))
    if "eta30" in summ:
        row["eta30_sd"] = summ["eta30"]["sd"]
    fixed = math.isinf(nu)
    for key in ("g30", "w31"):
        if key in res:
            cal = calibration(res, key)
            row[f"{key}_over_precision"] = cal["ratio"]           # E[w] / E[1/v]
            row[f"{key}_worst_bin_gap"] = cal["worst_bin_gap"]    # max over deciles |E[w] - E[1/v]|
    if "eta32" in res:
        row["eta32_mean"] = float(res["eta32"].mean())
        row["eta32_sd"] = float(res["eta32"].std(ddof=1))
        row["eta32_below_1"] = float(np.mean(res["eta32"] < 1))
    row["E_exact"] = 2 * k + push ** 2 if fixed else float("nan")
    row["K_true_exact"] = k if fixed else (k * nu / (nu - 2) if nu > 2 else float("inf"))
    row["Q_null_exact"] = null_misfit(k).mean() if (push == 0 and fixed) else float("nan")
    row["Q_plus_null_exact"] = null_positive_misfit(k).mean() if (push == 0 and fixed) else float("nan")
    row["Q_strong_limit"] = k - 1 if fixed else float("nan")
    return row


def run(path, out_root="results", replicates=None, progress=print):
    """Run a forward preset: one summary row per cell, written as JSON and CSV beside a copy of the
    preset. `replicates` overrides draws (recorded)."""
    import yaml
    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    raw = spec["raw"]
    draws = raw["draws"] if replicates is None else int(replicates)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base = Path(out_root) / raw["name"]
    out, suffix = base / stamp, 1
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    (out / "preset.yaml").write_bytes(source)
    rows = []
    for index, cell in enumerate(spec["cells"]):
        if progress:
            progress(f"cell {cell['id']} k={cell['dimension']} push={cell['push']} "
                     f"angle={cell['angle_deg']} N={cell['readings']} nu_eff={cell['nu_eff']}")
        rows.append(run_cell(cell, _cell_seed(raw["seed"], index), draws))
    record = dict(name=raw["name"], seed=raw["seed"], draws=draws, rows=rows)
    (out / "summary.json").write_text(json.dumps(record, indent=1, default=float), encoding="utf-8")
    keys = list(rows[0])
    for row in rows[1:]:
        keys += [key for key in row if key not in keys]
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    if progress:
        progress(f"wrote {out}")
    return out, record

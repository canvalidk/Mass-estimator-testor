"""ours29 in the full input type: any error covariance per reading, a declared subspace.

The rest of met/ is a batch pipeline built on closed forms for isotropic and per-axis
noise. This module is the definition's full type (definition/type.md) for one series:

  a reading declares its force and acceleration (any values, including zero), an error
  covariance Sigma (2d x 2d, positive definite, any correlation, including between the
  channels), a candidate subspace N (d x k, orthonormal columns; default the identity),
  and optionally its window (bookkeeping: L18, T33; the equation does not use it).

  the starting law pi is an input: half-Cauchy with centre c (definition/equation.md).
  There is no default centre. c = s = sigma_F / sigma_a is the principal branch for one
  noise pair (L3', O11); `instrument_ratio` computes it where it is defined.

The equation (definition/equation.md; the ours29 record, section 3):

  per reading, on the whole projective line m = c tan(phi):
     C(phi) = (c sin(phi) N; cos(phi) N),  G = C^T Sigma^-1 C,  b = C^T Sigma^-1 y,
     H(phi) = b^T G^-1 b, E = whitened energy of y in the span of (N 0; 0 N),
     Q = E - max H,  nu = k - 1,  g = min(1, nu / Q),  l(phi) = g H(phi) / 2
  series: phi_hat = argmax sum l;  curv = -sum l'';  J = sum l'^2  (at phi_hat)
     eta = min(1, curv / J), eta = 1 if J = 0 or curv <= 0
  law:  p(m) ~ pi(m) exp(eta sum l(m)) on m > 0
  readout: m_hat = int m A p / int A p,  A(m) = E[ |stacked a*| | m ],
     each reading's candidate alpha | m ~ N(G_m^-1 b_m, G_m^-1) with C_m = (m N; N).

Differences from reference/noise_vs_signal/ours29_reference.py, all numerical:
  - l' and l'' are analytic (the reference uses finite differences with h = 1e-4, which
    is what puts its rotation checks at 3-7e-8);
  - the best lines are found by Newton steps from a scan, to machine precision;
  - the law's grid adapts to the law's support, so narrow laws are resolved;
  - the readout's E|X| integral carries its two tails analytically.
Pass `grid=(half_width, nodes)` to use the reference's fixed grid (log c +- half_width).

Exact inputs (sigma = 0) do not enter `ours29`; `exact_mass` is the exact branch (T1).
"""
from dataclasses import dataclass, field
import math

import numpy as np

SCAN = 2048            # points on the projective line for the first scan of a best line
NEWTON_STEPS = 30
SUPPORT_NATS = 40.0    # the law's support: log density within this of its maximum
FINE_NODES = 2001      # nodes on the law's support
COARSE_HALF = 25.0     # first look at the law: log c +- this, in log m
COARSE_NODES = 4001
T_HALF = 40.0          # E|X| integral: log t within +- this of 1 / E|X|^2
T_NODES = 801


@dataclass
class Reading:
    """One reading in the full type. force, acceleration: length-d vectors. cov: (2d, 2d),
    the error covariance of (force; acceleration). subspace: (d, k) with orthonormal
    columns, or None for the identity. window: optional (start, end), carried, not used."""
    force: np.ndarray
    acceleration: np.ndarray
    cov: np.ndarray
    subspace: np.ndarray = None
    window: tuple = None
    _pre: dict = field(default=None, repr=False)

    def __post_init__(self):
        self.force = np.atleast_1d(np.asarray(self.force, dtype=float))
        self.acceleration = np.atleast_1d(np.asarray(self.acceleration, dtype=float))
        d = self.force.shape[0]
        if self.force.shape != (d,) or self.acceleration.shape != (d,) or d not in (1, 2, 3):
            raise ValueError("force and acceleration must be vectors of the same length, 1 to 3")
        self.cov = np.asarray(self.cov, dtype=float)
        if self.cov.shape != (2 * d, 2 * d):
            raise ValueError(f"cov must be {2 * d} x {2 * d} for a {d}-dimensional reading")
        if not np.allclose(self.cov, self.cov.T, rtol=1e-12, atol=0.0):
            raise ValueError("cov must be symmetric")
        if np.linalg.eigvalsh(self.cov).min() <= 0.0:
            raise ValueError("cov must be positive definite: every declared sigma > 0 (L2)")
        n = np.eye(d) if self.subspace is None else np.asarray(self.subspace, dtype=float)
        if n.ndim == 1:
            n = n[:, None]
        if n.shape[0] != d or not 1 <= n.shape[1] <= d:
            raise ValueError("subspace must be d x k with 1 <= k <= d")
        if not np.allclose(n.T @ n, np.eye(n.shape[1]), atol=1e-10):
            raise ValueError("subspace columns must be orthonormal")
        self.subspace = n

    @property
    def dimension(self):
        return self.force.shape[0]

    @property
    def k(self):
        return self.subspace.shape[1]

    def pre(self):
        """Blocks of the whitened problem projected on the subspace (cached)."""
        if self._pre is None:
            d, n = self.dimension, self.subspace
            si = np.linalg.inv(self.cov)
            y = np.concatenate([self.force, self.acceleration])
            k_ff = n.T @ si[:d, :d] @ n
            k_fa = n.T @ si[:d, d:] @ n
            k_aa = n.T @ si[d:, d:] @ n
            u = n.T @ (si[:d, :] @ y)
            v = n.T @ (si[d:, :] @ y)
            span = np.block([[n, np.zeros_like(n)], [np.zeros_like(n), n]])
            gp, bp = span.T @ si @ span, span.T @ si @ y
            energy = float(bp @ np.linalg.solve(gp, bp))
            self._pre = dict(k_ff=k_ff, k_s=k_fa + k_fa.T, k_aa=k_aa, u=u, v=v, energy=energy)
        return self._pre


def isotropic(force, acceleration, force_sd, acceleration_sd, window=None):
    """A reading with one SD per channel on every axis and no correlation."""
    force = np.atleast_1d(np.asarray(force, dtype=float))
    d = force.shape[0]
    cov = np.diag(np.r_[np.full(d, float(force_sd) ** 2), np.full(d, float(acceleration_sd) ** 2)])
    return Reading(force, acceleration, cov, window=window)


def instrument_ratio(reading):
    """s = sigma_F / sigma_a where the reading has one noise pair (isotropic, uncorrelated
    channels); None otherwise, since then no single ratio is singled out (O11)."""
    d = reading.dimension
    c = reading.cov
    sf, sa = np.diag(c)[:d], np.diag(c)[d:]
    off = c - np.diag(np.diag(c))
    if np.any(off != 0.0) or np.ptp(sf) > 0 or np.ptp(sa) > 0:
        return None
    return float(math.sqrt(sf[0] / sa[0]))


def _h_terms(pre, c, phi, order=0):
    """H(phi) and, for order >= 1 and 2, its first and second derivatives, for an array phi."""
    phi = np.asarray(phi, dtype=float)
    s, co = np.sin(phi)[:, None, None], np.cos(phi)[:, None, None]
    kff, ks, kaa = pre["k_ff"][None], pre["k_s"][None], pre["k_aa"][None]
    g = c * c * s * s * kff + c * s * co * ks + co * co * kaa
    b = c * s[..., 0] * pre["u"][None] + co[..., 0] * pre["v"][None]
    w = np.linalg.solve(g, b[..., None])[..., 0]
    h = np.einsum("ni,ni->n", b, w)
    if order == 0:
        return h
    g1 = 2 * c * c * s * co * kff + c * (co * co - s * s) * ks - 2 * s * co * kaa
    b1 = c * co[..., 0] * pre["u"][None] - s[..., 0] * pre["v"][None]
    h1 = 2 * np.einsum("ni,ni->n", b1, w) - np.einsum("ni,nij,nj->n", w, g1, w)
    if order == 1:
        return h, h1
    cp_g_cp = c * c * co * co * kff - c * co * s * ks + s * s * kaa       # C'^T Si C'
    g2 = -2 * g + 2 * cp_g_cp
    w1 = np.linalg.solve(g, (b1 - np.einsum("nij,nj->ni", g1, w))[..., None])[..., 0]
    # H'' = 2 b''.w + 2 b'.w' - 2 w'.G'.w - w.G''.w, with b'' = -b
    h2 = (-2 * h + 2 * np.einsum("ni,ni->n", b1, w1)
          - 2 * np.einsum("ni,nij,nj->n", w1, g1, w) - np.einsum("ni,nij,nj->n", w, g2, w))
    return h, h1, h2


def _best_line(fun):
    """argmax on the projective line (-pi/2, pi/2] of a pi-periodic smooth function given as
    fun(phi, order) -> (value, d1, d2): a scan, then Newton steps from the best scan point."""
    grid = -math.pi / 2 + math.pi * (np.arange(SCAN) + 0.5) / SCAN
    values = fun(grid, 0)
    phi = float(grid[int(np.argmax(values))])
    step = math.pi / SCAN
    for _ in range(NEWTON_STEPS):
        _, d1, d2 = (float(x[0]) for x in fun(np.array([phi]), 2))
        if d2 < 0:
            move = -d1 / d2
            move = max(-step, min(step, move))
        else:                                   # not concave here: climb by the gradient
            move = step * math.copysign(1.0, d1) if d1 != 0 else 0.0
        phi += move
        if abs(move) < 1e-15:
            break
    phi = (phi + math.pi / 2) % math.pi - math.pi / 2
    return phi


def reading_state(reading, c):
    """The reading's best line on the whole projective line, its misfit Q, k and nu."""
    pre = reading.pre()

    def fun(phi, order):
        return _h_terms(pre, c, phi, order) if order else _h_terms(pre, c, phi)

    phi = _best_line(fun)
    h_max = float(_h_terms(pre, c, np.array([phi]))[0])
    return dict(Q=max(pre["energy"] - h_max, 0.0), E=pre["energy"], k=reading.k, nu=reading.k - 1, phi=phi)


def misfit_weight(q, nu, form="kink"):
    """g(Q): 1 for no misfit or nu = 0; kink min(1, nu/Q) (T30, O12 closed); smooth 1/(1+Q/nu)."""
    if nu <= 0 or q <= 0:
        return 1.0
    if form == "kink":
        return min(1.0, nu / q)
    if form == "smooth":
        return 1.0 / (1.0 + q / nu)
    raise ValueError(f"unknown misfit form {form!r}")


def _series_terms(readings, c, g, phi, order):
    """Sum and per-reading l = g H / 2 (and derivatives) on an array phi."""
    per = []
    for r, gi in zip(readings, g):
        terms = _h_terms(r.pre(), c, phi, order) if order else (_h_terms(r.pre(), c, phi),)
        per.append([0.5 * gi * t for t in terms])
    total = [sum(p[j] for p in per) for j in range(len(per[0]))]
    return total, per


def bartlett_eta(readings, c, g):
    """eta = min(1, curv / J) at the series' best line on the whole projective line (T31).
    Returns (eta, curv, J, phi_hat)."""
    if not readings:
        return 1.0, 0.0, 0.0, 0.0

    def fun(phi, order):
        total, _ = _series_terms(readings, c, g, phi, order)
        return tuple(total) if order else total[0]

    phi_hat = _best_line(fun)
    _, per = _series_terms(readings, c, g, np.array([phi_hat]), 2)
    d1 = np.array([p[1][0] for p in per])
    d2 = np.array([p[2][0] for p in per])
    curv, j = float(-d2.sum()), float((d1 ** 2).sum())
    if j <= 1e-300 or curv <= 0:
        return 1.0, curv, j, phi_hat
    return float(min(1.0, curv / j)), curv, j, phi_hat


def _log_prior(lm, c):
    """log of the half-Cauchy as a density in log m."""
    return math.log(c) + lm - np.logaddexp(2 * math.log(c), 2 * lm)


def _log_law(readings, c, g, eta, lm):
    if not readings:
        return _log_prior(lm, c)
    phi = np.arctan(np.exp(lm - math.log(c)))
    total, _ = _series_terms(readings, c, g, phi, 0)
    return _log_prior(lm, c) + eta * total[0]


def _grid(readings, c, g, eta, grid):
    if grid is not None:
        half, nodes = grid
        return np.linspace(math.log(c) - half, math.log(c) + half, int(nodes))
    lm = np.linspace(math.log(c) - COARSE_HALF, math.log(c) + COARSE_HALF, COARSE_NODES)
    lp = _log_law(readings, c, g, eta, lm)
    keep = np.nonzero(lp > lp.max() - SUPPORT_NATS)[0]
    step = lm[1] - lm[0]
    lo, hi = lm[keep[0]] - 2 * step, lm[keep[-1]] + 2 * step
    return np.linspace(lo, hi, FINE_NODES)


def _quantiles(lm, lp, probs):
    """Quantiles of m and the law's weights on lm (summing to 1), by met.analyst's
    fourth-order CDF and cubic-Hermite inversion (the same numerics as the batch pipeline)."""
    from met.analyst import _cdf, _quantiles as _invert
    dens = np.exp(lp - lp.max())[None, :]
    h = np.array([lm[1] - lm[0]])
    cdf_raw = _cdf(dens, h)                       # normalised to 1 at the last node
    area = h[0] * (dens[0].sum() - 0.5 * (dens[0, 0] + dens[0, -1]))
    q = _invert(lm[None, :], dens / area, cdf_raw, list(probs))[0]
    p = dens[0] / dens[0].sum()
    return np.exp(q), p


def _stacked_mean_length(m, readings):
    """A(m) = E|X| for X the stacked candidates alpha_i | m ~ N(mu_i, G_i^-1), at each m.

    E|X| = (1 / (2 sqrt(pi))) int_0^inf (1 - E exp(-t |X|^2)) t^(-3/2) dt, on a grid in
    log t centred on 1 / E|X|^2, with the two tails carried analytically:
    below t0 the integrand is E|X|^2 t^(-1/2) to first order (2 E|X|^2 sqrt(t0)); above t1
    it is t^(-3/2) to first order (2 / sqrt(t1))."""
    m = np.asarray(m, dtype=float)
    mus, lams, rots = [], [], []
    for r in readings:
        pre = r.pre()
        g = (m[:, None, None] ** 2 * pre["k_ff"][None] + m[:, None, None] * pre["k_s"][None]
             + pre["k_aa"][None])
        b = m[:, None] * pre["u"][None] + pre["v"][None]
        mu = np.linalg.solve(g, b[..., None])[..., 0]
        lam, rot = np.linalg.eigh(np.linalg.inv(g))
        mus.append(np.einsum("mji,mj->mi", rot, mu) ** 2)
        lams.append(lam)
    ex2 = sum((mt + lam).sum(axis=1) for mt, lam in zip(mus, lams))            # E|X|^2, (M,)
    lt = np.linspace(-T_HALF, T_HALF, T_NODES)
    t = np.exp(lt)[None, :] / np.maximum(ex2, 1e-300)[:, None]                 # (M, T)
    log_l = np.zeros_like(t)
    t2 = 2 * t
    for mt, lam in zip(mus, lams):
        # One component at a time in 2D, in place, with the operations and the order of
        # addition the (M, T, k) form had: bit for bit the same, about 2.5 times faster
        # (numpy's reduction over a last axis of length <= 3 was most of the time).
        part = None
        for j in range(lam.shape[1]):
            tl = np.multiply(t2, lam[:, j:j + 1])
            term = np.log1p(tl)
            term *= -0.5
            tl += 1
            push = np.multiply(t, mt[:, j:j + 1])
            push /= tl
            term -= push
            if part is None:
                part = term
            else:
                part += term
        log_l += part
    integrand = -np.expm1(log_l) * t ** -0.5                                    # d(log t) form
    body = np.trapezoid(integrand, lt, axis=1)
    t0, t1 = t[:, 0], t[:, -1]
    tails = 2 * ex2 * np.sqrt(t0) + 2 / np.sqrt(t1)
    return (body + tails) / (2 * math.sqrt(math.pi))


def ours29(readings, c, *, misfit=True, calibrate=True, misfit_form="kink", readout=True, grid=None):
    """The definition of mass on one series, in the full type.

    readings: a list of Reading. c: the starting law's centre (the input pi; no default).
    misfit=False gives g_i = 1; calibrate=False gives eta = 1 (the raw law); both off is ours27.
    grid: None for the adaptive grid, or (half_width, nodes) for log c +- half_width.

    Returns a dict: readout, median, interval (95%, quoted: eta applied), raw_interval (eta = 1),
    eta, g, Q, phi_hat, law (normalised weights on lm), lm, A (on lm, where computed)."""
    if not c > 0:
        raise ValueError("the starting law's centre c must be positive")
    readings = list(readings)
    states = [reading_state(r, c) for r in readings]
    g = [misfit_weight(s["Q"], s["nu"], misfit_form) if misfit else 1.0 for s in states]
    if calibrate:
        eta, curv, j, phi_hat = bartlett_eta(readings, c, g)
    else:
        eta, curv, j, phi_hat = 1.0, float("nan"), float("nan"), float("nan")
    lm = _grid(readings, c, g, eta, grid)
    lp = _log_law(readings, c, g, eta, lm)
    (lo, med, hi), p = _quantiles(lm, lp, [0.025, 0.5, 0.975])
    out = dict(median=float(med), interval=(float(lo), float(hi)), eta=eta, g=g,
               Q=[s["Q"] for s in states], curv=curv, J=j, phi_hat=phi_hat, law=p, lm=lm)
    if eta < 1.0:
        lm_raw = _grid(readings, c, g, 1.0, grid)
        (rlo, _, rhi), _ = _quantiles(lm_raw, _log_law(readings, c, g, 1.0, lm_raw), [0.025, 0.5, 0.975])
        out["raw_interval"] = (float(rlo), float(rhi))
    else:
        out["raw_interval"] = out["interval"]
    if readout:
        if not readings:
            out["readout"] = float(c)          # no candidates: the input's centre (T24, TASK.md 2)
            return out
        keep = p > 1e-14 * p.max()
        a = np.zeros_like(lm)
        a[keep] = _stacked_mean_length(np.exp(lm[keep]), readings)
        m = np.exp(lm)
        out["readout"] = float((m * a * p).sum() / (a * p).sum())
        out["A"] = a
    return out


def exact_mass(force, acceleration, rtol=1e-12):
    """The exact branch (T1), for sigma = 0: one mass, no mass, or every mass.

    Pairing form (10 Oct): m1 = F.a / |a|^2, m2 = |F|^2 / (F.a). One mass iff the pair is
    codirectional with positive sizes. Returns ("one", m), ("none", reason) or ("every", None)."""
    f = np.atleast_1d(np.asarray(force, dtype=float))
    a = np.atleast_1d(np.asarray(acceleration, dtype=float))
    ff, aa, fa = float(f @ f), float(a @ a), float(f @ a)
    if ff == 0.0 and aa == 0.0:
        return ("every", None)
    if aa == 0.0:
        return ("none", "acceleration is zero and force is not: m2 is infinite")
    if ff == 0.0:
        return ("none", "force is zero and acceleration is not: m1 is zero")
    if fa * fa < ff * aa * (1 - rtol):
        return ("none", "not codirectional: the two pairings disagree (P < 1)")
    if fa <= 0.0:
        return ("none", "anti-aligned: the pairings agree on a negative value")
    return ("one", fa / aa)

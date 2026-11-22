"""A family of pooled laws of mass on one scalar instrument: the witness of exploration E12.

Harvested 22 Nov 2026 from the exploration branch `explore/e12-scalar-mass-proof` (Codex, 18 Nov;
`e12-scalar-proof.md` section 2, `e12-pooled-witness.md`), written again from the formulas and
checked by Claude (T78). It is a witness: a function that meets four stated requirements on this
one instrument, in a family with free numbers `alpha` and `b`. Not an equation of the definition,
not a probe, and nothing selects a member.

One instrument, scalar readings, independent channels. In noise units Y_i = (F_i / sigma_F,
a_i / sigma_a), a candidate mass is the angle theta with m = s tan(theta), s = sigma_F / sigma_a.

  B = sum Y_i Y_i'            L, l its eigenvalues, H = L - l, theta_hat its principal line
  Q = the least of v'Bv over the closed arc of positive masses     (T73's leftover)
  C = sum h_b(|Y_i|^2),  h_b(x) = (1 - z)(1 + (2 + 1/b) z),  z = exp(-b x / 2)
                              (a count of readings that is 0 at a zero reading; mean >= 1 each)
  d = min(1, (C / Q) (L / Q)^alpha)          the discount: count over leftover, times a factor
  K = d H^2 / L                              the Gaussian's H, times the weak-signal share H / L,
                                             times the discount
  K' = min(K, 1 + min(B_11, B_22)^2)         (matters for a reading on an axis)
  U(theta) = - sum log(1 + (v'Y_i)^2 / (1 + |Y_i|^2)) / (1 + sum |Y_i|^2 / (1 + |Y_i|^2))
  log factor = -K' sin^2(theta - theta_hat) / 2 + exp(-H^2 / L) (Q / (1 + Q)) U(theta)
                                             (the second term matters when the eigenvalues tie)
  law = prior times exp(log factor) on 0 < theta < pi/2;  point = its median
  report = the hull of the law's central interval and T45's score set at the same level
"""
import math

import numpy as np
from scipy import stats

NODES = 20001
LOCAL = 6001


def _rows(force, acceleration, force_sd, acceleration_sd):
    rows = np.column_stack([np.asarray(force, dtype=float) / force_sd,
                            np.asarray(acceleration, dtype=float) / acceleration_sd])
    return rows.reshape(-1, 2)


def count(energies, b=1.0):
    """C: the smooth, zero-safe count of readings from their squared lengths in noise units."""
    z = np.exp(-0.5 * b * np.asarray(energies, dtype=float))
    return float(np.sum((1.0 - z) * (1.0 + (2.0 + 1.0 / b) * z)))


def summary(rows, alpha, b=1.0):
    """The pooled quantities of the docstring for whitened rows (N, 2)."""
    rows = np.asarray(rows, dtype=float).reshape(-1, 2)
    scatter = rows.T @ rows
    values, vectors = np.linalg.eigh(scatter)
    largest, smallest = float(values[1]), max(float(values[0]), 0.0)
    gap = largest - smallest
    f, a = vectors[:, 1]
    angle = math.atan2(f, a) % math.pi if gap > 0 else 0.0          # the principal line, mod pi
    if scatter[0, 1] > 0:
        least = smallest                                            # the principal line has positive slope
    else:
        least = float(min(scatter[0, 0], scatter[1, 1]))            # else the better of the two axes
    energies = np.einsum("ij,ij->i", rows, rows)
    counted = count(energies, b)
    discount = 1.0
    if least > 0 and largest > 0:
        discount = min(1.0, (counted / least) * (largest / least) ** alpha)
    precision = discount * gap * gap / largest if largest > 0 else 0.0
    capped = min(precision, 1.0 + float(min(scatter[0, 0], scatter[1, 1])) ** 2)
    return dict(scatter=scatter, largest=largest, gap=gap, angle=angle, least=least, count=counted,
                discount=discount, precision=precision, capped=capped, energies=energies)


def log_factor(rows, angles, alpha, b=1.0):
    """The law's log density relative to the prior, up to a constant, at the given angles."""
    rows = np.asarray(rows, dtype=float).reshape(-1, 2)
    angles = np.asarray(angles, dtype=float)
    s = summary(rows, alpha, b)
    if s["largest"] == 0.0:
        return np.zeros_like(angles)
    out = -0.5 * s["capped"] * np.sin(angles - s["angle"]) ** 2
    gate = math.exp(-s["gap"] ** 2 / s["largest"]) * s["least"] / (1.0 + s["least"])
    if gate > 1e-17:                                             # the term is at most `gate` in size
        across = np.cos(angles)[:, None] * rows[:, 0][None, :] - np.sin(angles)[:, None] * rows[:, 1][None, :]
        shape = np.sum(np.log1p(across ** 2 / (1.0 + s["energies"])[None, :]), axis=1)
        out = out - gate * shape / (1.0 + float(np.sum(s["energies"] / (1.0 + s["energies"]))))
    return out


def _log_prior(angles, ratio, centre):
    """Half-Cauchy in mass with the given centre, as a density in the noise angle. centre None is
    the centre at the instrument's ratio: uniform in the angle."""
    if centre is None:
        return np.zeros_like(angles)
    c = centre / ratio
    return np.log(c) - np.log(c * c * np.cos(angles) ** 2 + np.sin(angles) ** 2)


def witness(force, acceleration, force_sd, acceleration_sd, *, alpha, b=1.0, centre=None, level=0.95):
    """The witness law's point, its central interval, and the hull report at `level`.

    Returns point (the median, in mass), interval, report (the hull with the score set; an end may
    be 0.0 or inf), discount, precision, least, count."""
    ratio = force_sd / acceleration_sd
    rows = _rows(force, acceleration, force_sd, acceleration_sd)
    s = summary(rows, alpha, b)
    grid = np.linspace(0.0, 0.5 * math.pi, NODES)[1:-1]
    if s["capped"] > 0:
        # refine where the law can be narrow: at the principal line if it has positive slope,
        # otherwise at both ends of the arc
        half = min(14.0 / math.sqrt(s["capped"]), 0.5 * math.pi)
        centres = [s["angle"]] if s["angle"] <= 0.5 * math.pi else [0.0, 0.5 * math.pi]
        local = np.concatenate([np.linspace(c - half, c + half, LOCAL) for c in centres])
        grid = np.unique(np.concatenate([grid, local[(local > 0) & (local < 0.5 * math.pi)]]))
    log_density = log_factor(rows, grid, alpha, b) + _log_prior(grid, ratio, centre)
    density = np.exp(log_density - log_density.max())
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    tail = 0.5 * (1.0 - level)
    lo, mid, hi = (ratio * math.tan(float(np.interp(q, cdf, grid))) for q in (tail, 0.5, 1.0 - tail))

    # T45's score set on this slice: (v'Bu)^2 <= z^2 u'Bu, the ends of the arc included
    z2 = float(stats.chi2.ppf(level, 1))
    full = np.concatenate([[0.0], grid, [0.5 * math.pi]])
    u = np.stack([np.sin(full), np.cos(full)])
    v = np.stack([np.cos(full), -np.sin(full)])
    score = np.einsum("ia,ij,ja->a", v, s["scatter"], u)
    variance = np.einsum("ia,ij,ja->a", u, s["scatter"], u)
    kept = score ** 2 <= z2 * variance
    report_lo, report_hi = lo, hi
    if kept.any():
        first, last = full[kept][0], full[kept][-1]
        report_lo = min(lo, ratio * math.tan(first))
        report_hi = math.inf if last >= 0.5 * math.pi else max(hi, ratio * math.tan(last))
    return dict(point=mid, interval=(lo, hi), report=(report_lo, report_hi), discount=s["discount"],
                precision=s["capped"], least=s["least"], count=s["count"])

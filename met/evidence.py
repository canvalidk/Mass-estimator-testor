"""The evidence curve of a series: T36's leftover at every candidate mass, with no starting law.

Harvested 22 Nov 2026 from the exploration branch `explore/e13-unnormalized-evidence` (Can, 21 Nov;
insights 1, 3 and 4), rewritten on met.full's reading type and checked by Claude (T73, T74;
`problems/O20.md`). Nothing here is a law of mass: there is no pi, no normalisation and no point.

  K(m) = sum_i K_i(m),  K_i(m) = the squared length, in the reading's noise metric, of the part of
         the reading that the candidate mass does not explain (T36), inside the declared span.
  exp(-K(m) / 2) is the declared Gaussian's probability of the readings at the best compatible
         pair for m, up to a factor that does not depend on m.

Candidates are carried as an angle u in [0, pi/2] with m = scale * tan(u), so both ends of the
positive masses (m -> 0 and m -> infinity) are ordinary points of the curve. `scale` must come
from the declarations, never from the data; the default is the geometric mean over readings of
sqrt(trace Sigma_FF / trace Sigma_aa), which is sigma_F / sigma_a on one isotropic instrument.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy import stats

from met.full import Reading

GRID = 721                 # angles on the arc for the first scan of a least value
GOLDEN_STEPS = 48
_GR = 0.5 * (math.sqrt(5.0) - 1.0)


@dataclass
class Span:
    """One reading in the coordinates of its own whitened span W (dimension 2k): the data z,
    and the two k-column blocks whose combination scale*sin(u)*bf + cos(u)*ba spans the
    candidate's mean space at m = scale*tan(u)."""
    z: np.ndarray        # (2k,)
    bf: np.ndarray       # (2k, k): a unit force push along each subspace direction
    ba: np.ndarray       # (2k, k): a unit acceleration along each subspace direction
    k: int

    @property
    def energy(self):
        return float(self.z @ self.z)


def span_of(reading):
    if not isinstance(reading, Reading):
        raise ValueError("each reading must be a met.full.Reading")
    d, n = reading.dimension, reading.subspace
    k = n.shape[1]
    chol = np.linalg.cholesky(reading.cov)
    y = np.concatenate([reading.force, reading.acceleration])
    zero = np.zeros_like(n)
    push_f, push_a = np.vstack([n, zero]), np.vstack([zero, n])
    white = np.linalg.solve(chol, np.column_stack([push_f, push_a, y]))
    basis, _ = np.linalg.qr(white[:, :2 * k])                     # (2d, 2k), orthonormal
    z = basis.T @ white[:, -1]
    if not np.any(y):
        z = np.zeros(2 * k)                                       # an exact zero stays exact
    return Span(z=z, bf=basis.T @ white[:, :k], ba=basis.T @ white[:, k:2 * k], k=k)


def declared_scale(readings):
    """A mass scale from the declared covariances alone (never the data)."""
    logs = []
    for r in readings:
        d = r.dimension
        logs.append(0.5 * math.log(np.trace(r.cov[:d, :d]) / np.trace(r.cov[d:, d:])))
    return math.exp(sum(logs) / len(logs)) if logs else 1.0


def angle_of(mass, scale):
    mass = np.asarray(mass, dtype=float)
    if np.any(mass < 0) or np.any(np.isnan(mass)):
        raise ValueError("a candidate mass is in [0, infinity]")
    return np.where(np.isinf(mass), 0.5 * math.pi, np.arctan(mass / scale))


def unexplained(span, u, z, scale):
    """K: the squared length of the part of z outside the candidate's mean space. `u` is an array
    of angles (n,). `z` is (2k,) for one data vector, or (n, 2k) for one data vector per angle.
    Returns (n,). Computed from the leftover vector itself, not as energy minus the explained
    part, so a leftover of order 1 beside a signal of 1e5 noise units keeps its digits."""
    u = np.atleast_1d(np.asarray(u, dtype=float))
    b = (scale * np.sin(u))[:, None, None] * span.bf[None] + np.cos(u)[:, None, None] * span.ba[None]
    gram = np.einsum("nij,nil->njl", b, b)
    data = np.broadcast_to(z, (u.shape[0], z.shape[-1])) if z.ndim == 1 else z
    rhs = np.einsum("nij,ni->nj", b, data)
    rest = data - np.einsum("nij,nj->ni", b, np.linalg.solve(gram, rhs[..., None])[..., 0])
    return np.einsum("ni,ni->n", rest, rest)


def split(span, u, scale):
    """At one candidate angle: orthonormal bases (2k, k) of the candidate's mean space and of
    the rest of the span, and the data's coordinates U (along) and V (leftover) in them."""
    b = scale * math.sin(u) * span.bf + math.cos(u) * span.ba
    q, _ = np.linalg.qr(b, mode="complete")
    par, perp = q[:, :span.k], q[:, span.k:]
    return par, perp, par.T @ span.z, perp.T @ span.z


class Curve:
    """K(m) for one series. Call it with masses (0 and inf allowed) or use `at_angle`."""

    def __init__(self, readings, scale=None):
        self.readings = tuple(readings)
        self.spans = [span_of(r) for r in self.readings]
        self.scale = float(declared_scale(self.readings) if scale is None else scale)
        if not self.scale > 0 or not math.isfinite(self.scale):
            raise ValueError("scale must be finite and positive")
        self.count = sum(s.k for s in self.spans)                 # sum of k_i: noise components

    def at_angle(self, u):
        u = np.asarray(u, dtype=float)
        flat = np.atleast_1d(u).ravel()
        total = np.zeros(flat.shape)
        for s in self.spans:
            total += unexplained(s, flat, s.z, self.scale)
        return total.reshape(u.shape) if u.ndim else float(total[0])

    def __call__(self, mass):
        return self.at_angle(angle_of(mass, self.scale))

    def mass_at(self, u):
        return math.inf if u >= 0.5 * math.pi else self.scale * math.tan(u)

    def least(self, lo=0.0, hi=0.5 * math.pi):
        """The least value of K on the closed arc [lo, hi] and an angle where it is taken."""
        return _least_on(self.at_angle, lo, hi)


def _least_on(value_at, lo, hi, fine=GRID):
    """The least of a function of one variable on [lo, hi]: a scan, then golden section in the
    bracket around the best scan point. `value_at` takes and returns arrays."""
    grid = np.linspace(lo, hi, fine)
    values = np.asarray(value_at(grid), dtype=float)
    j = int(np.argmin(values))
    a, b = float(grid[max(j - 1, 0)]), float(grid[min(j + 1, fine - 1)])
    one = lambda x: float(np.asarray(value_at(np.array([x])), dtype=float)[0])
    c, d = b - _GR * (b - a), a + _GR * (b - a)
    fc, fd = one(c), one(d)
    for _ in range(GOLDEN_STEPS):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - _GR * (b - a)
            fc = one(c)
        else:
            a, c, fc = c, d, fd
            d = a + _GR * (b - a)
            fd = one(d)
    best_x, best = (c, fc) if fc < fd else (d, fd)
    if values[j] <= best:
        best_x, best = float(grid[j]), float(values[j])
    return float(best), float(best_x)


def leftover(readings, scale=None):
    """What no single positive mass explains: the least of K over the closed arc of positive
    masses, where it is taken, and the largest value the curve reaches.

    `flat` is True when the curve is constant (the all-zero experiment, and any series whose
    readings carry nothing along any candidate line). `bound` is P(chi2_count >= least): at the
    true mass K is chi-square with `count` = sum k_i components and the least is no larger, so
    this is an upper bound on how probable a leftover this large was. The count takes in every
    reading, zero readings included, so the bound is moved by a zero reading although the least
    itself is not (T73)."""
    curve = Curve(readings, scale)
    least, u = curve.least()
    top = float(np.max(curve.at_angle(np.linspace(0.0, 0.5 * math.pi, GRID))))
    flat = top - least <= 1e-12 * max(1.0, top)
    return dict(least=least, mass=None if flat else curve.mass_at(u), angle=u, flat=bool(flat),
                at_end=bool(u <= 0.0 or u >= 0.5 * math.pi), largest=top, count=curve.count,
                bound=float(stats.chi2.sf(least, curve.count)) if curve.count else 1.0)


def lump(first, second, total):
    """Two bodies joined at a common acceleration, M = m1 + m2 (L5), with no starting law:
    K(M) = inf over m1 + m2 = M of first(m1) + second(m2). `first`, `second`: curves (callables
    of mass, 0 and inf allowed). Returns the composite's curve at each total."""
    total = np.atleast_1d(np.asarray(total, dtype=float))
    out = np.empty(total.shape)
    for i, big in enumerate(total):
        if math.isinf(big):
            # an unbounded total: one part is unbounded and the other is anything
            out[i] = min(float(first(math.inf)) + _least_of(second), _least_of(first) + float(second(math.inf)))
            continue
        out[i] = _least_on(lambda t, big=big: first(t * big) + second((1.0 - t) * big), 0.0, 1.0)[0]
    return out


def common_force(first, second, total):
    """Two bodies under a common force, 1/M = 1/m1 + 1/m2 (L14): the same operation on inverse
    masses."""
    total = np.atleast_1d(np.asarray(total, dtype=float))
    inverse = np.where(total > 0, 1.0 / np.where(total > 0, total, 1.0), math.inf)
    return lump(lambda w: first(_reciprocal(w)), lambda w: second(_reciprocal(w)), inverse)


def _reciprocal(w):
    w = np.asarray(w, dtype=float)
    with np.errstate(divide="ignore"):
        return np.where(w > 0, 1.0 / np.where(w > 0, w, 1.0), math.inf) if w.ndim else (
            1.0 / float(w) if float(w) > 0 else math.inf)


def _least_of(curve):
    if isinstance(curve, Curve):
        return curve.least()[0]

    def at_angle(u):
        with np.errstate(divide="ignore", over="ignore"):
            return curve(np.where(u >= 0.5 * math.pi, math.inf, np.tan(u)))
    return _least_on(at_angle, 0.0, 0.5 * math.pi)[0]

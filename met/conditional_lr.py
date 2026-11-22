"""Exact conditional tests of one candidate mass on the full reading type, with no starting law.

Harvested 22 Nov 2026 from the exploration branch `explore/e13-unnormalized-evidence` (Can, 21 Nov;
insights 6 to 9), rewritten on met.full's reading type and checked by Claude (T75; the two-body test
`lump_test` is T74 (6); the weighted score `score_test` and its two closed forms are T76;
`problems/O20.md`). A research reference beside `met/conditional_score.py` (T45): it returns tests
of candidate masses, not a law of mass and not a point.

At a candidate m0 each reading splits, in its own noise metric, into the part along the
candidate's mean space, U_i, and the leftover, V_i (T36, T45). If m0 is the true mass and the
world's errors have the declared Gaussian law, then given every U_i the leftovers are independent
standard normal vectors, whatever the pushes were. So the law of ANY function of the readings
is known given the U_i, and comparing its observed value with that law is an exact test of m0.

  statistic:  LR(m0) = W(m0) - inf_m W(m),  W(m) = sum_i w_i K_i(m),  w_i = weight(|U_i(m0)|^2)
              (the drop of the weighted leftover from the candidate to its least; weight None
              is w_i = 1, the Gaussian likelihood-ratio statistic). The weights are fixed by
              the conditioning. A weight with weight(0) = 0 makes an exactly zero reading drop out.
  null "gaussian":   V_i ~ N(0, I) given U_i. Exact when the world is the declared Gaussian.
  null "spherical":  V_i = |V_i| times a uniform direction, given U_i and |V_i|. Exact whenever
              each reading's error, in its noise metric, has a law that does not depend on
              direction: the declared Gaussian, and the same with any unknown size per reading.

The p-value is the share of the conditional law at or above the observed statistic. It is
computed from `draws` random leftovers (seeded; (1 + count) / (draws + 1)), or, for the spherical
null on scalar readings, from every one of the 2^n sign patterns when n <= `enumerate_up_to`.
The infimum is taken on the closed arc of positive masses (domain "positive") or on the whole
projective line of slopes (domain "line"), by one fixed rule applied alike to the observed
readings and to every draw: a scan of `grid` angles, then golden section around the best.
"""
import math

import numpy as np

from met.evidence import angle_of, declared_scale, span_of, split, unexplained, _GR

GRID = 361
GOLDEN_STEPS = 36


def _projectors(span, angles, scale):
    """The projector on the candidate's mean space at each grid angle: (G, 2k, 2k)."""
    b = (scale * np.sin(angles))[:, None, None] * span.bf[None] + np.cos(angles)[:, None, None] * span.ba[None]
    gram = np.einsum("gij,gil->gjl", b, b)
    return np.einsum("gij,gjl,gml->gim", b, np.linalg.inv(gram), b)


def _golden(value_at, a, b):
    """The least of value_at on [a, b], one bracket per draw (arrays), by golden section."""
    c, d = b - _GR * (b - a), a + _GR * (b - a)
    fc, fd = value_at(c), value_at(d)
    for _ in range(GOLDEN_STEPS):
        left = fc < fd                                       # the least lies in [a, d]; else in [c, b]
        a, b = np.where(left, a, c), np.where(left, d, b)
        c_new, d_new = b - _GR * (b - a), a + _GR * (b - a)
        fresh = value_at(np.where(left, c_new, d_new))
        fc, fd = np.where(left, fresh, fd), np.where(left, fc, fresh)
        c, d = c_new, d_new
    return np.minimum(fc, fd)


def _scan(weights, data, projectors, draws):
    """The weighted leftover of every draw at every angle of a stack of projectors: (D, G).
    The fast form (energy less the explained part); it only chooses a bracket."""
    out = np.zeros((draws, projectors[0].shape[0])) if projectors else np.zeros((draws, 1))
    for w, z, p in zip(weights, data, projectors):
        outer = (z[:, :, None] * z[:, None, :]).reshape(draws, -1)
        out += w * (np.einsum("di,di->d", z, z)[:, None] - outer @ p.reshape(p.shape[0], -1).T)
    return out


class _Statistic:
    def __init__(self, spans, weights, scale, u0, domain, grid):
        self.spans, self.weights, self.scale, self.u0 = spans, weights, scale, u0
        if domain == "positive":
            self.lo, self.hi, self.clip = 0.0, 0.5 * math.pi, True
            self.angles = np.linspace(0.0, 0.5 * math.pi, grid)
        elif domain == "line":
            self.lo, self.hi, self.clip = -0.5 * math.pi, 0.5 * math.pi, False
            self.angles = -0.5 * math.pi + math.pi * np.arange(grid) / grid
        else:
            raise ValueError("domain is 'positive' or 'line'")
        self.step = float(self.angles[1] - self.angles[0])
        self.proj = [_projectors(s, self.angles, scale) for s in spans]

    def _at(self, data, u):
        """W at one angle per draw. data: list of (D, 2k); u: (D,)."""
        total = np.zeros(u.shape)
        for s, w, z in zip(self.spans, self.weights, data):
            total += w * unexplained(s, u, z, self.scale)
        return total

    def least(self, data, draws):
        """inf of W over the domain, for each draw."""
        if not self.spans:
            return np.zeros(draws)
        j = np.argmin(_scan(self.weights, data, self.proj, draws), axis=1)
        a, b = self.angles[j] - self.step, self.angles[j] + self.step
        if self.clip:
            a, b = np.maximum(a, self.lo), np.minimum(b, self.hi)
        return np.minimum(self._at(data, self.angles[j]), _golden(lambda u: self._at(data, u), a, b))

    def __call__(self, data):
        """LR for each draw. data: list over readings of arrays (D, 2k)."""
        draws = data[0].shape[0]
        w0 = self._at(data, np.full(draws, self.u0))
        least = np.minimum(self.least(data, draws), w0)
        return np.maximum(w0 - least, 0.0), w0


def conditional_test(readings, mass, *, null="gaussian", weight=None, domain="positive",
                     draws=999, seed=0, scale=None, grid=GRID, enumerate_up_to=12):
    """Exact conditional test of one candidate mass (see the module docstring).

    Returns p_value, statistic (the observed LR), at_candidate (the weighted leftover at the
    candidate), weights, and `empty` (True when no reading has a positive weight, or the series
    is empty: the test then accepts, p_value = 1)."""
    readings = tuple(readings)
    if null not in ("gaussian", "spherical"):
        raise ValueError("null is 'gaussian' or 'spherical'")
    scale = float(declared_scale(readings) if scale is None else scale) if readings else 1.0
    u0 = float(angle_of(mass, scale))
    spans, weights, bases, perps, norms, all_weights = [], [], [], [], [], []
    for r in readings:
        s = span_of(r)
        par, perp, along, rest = split(s, u0, scale)
        h = float(along @ along)
        w = 1.0 if weight is None else float(weight(h))
        if w < 0 or not math.isfinite(w):
            raise ValueError("a weight is finite and not negative")
        all_weights.append(w)
        if w == 0.0 or (null == "spherical" and not np.any(s.z)):
            continue                      # no weight; or nothing to turn (it would add zero anyway)
        spans.append(s)
        weights.append(w)
        bases.append(par @ along)
        perps.append(perp)
        norms.append(float(np.linalg.norm(rest)))
    if not spans:
        return dict(p_value=1.0, statistic=0.0, at_candidate=0.0, weights=all_weights, empty=True)
    stat = _Statistic(spans, weights, scale, u0, domain, grid)
    observed, w0 = stat([s.z[None, :] for s in spans])
    observed, w0 = float(observed[0]), float(w0[0])
    out = dict(statistic=observed, at_candidate=w0, weights=all_weights, empty=False)
    if observed <= 1e-13 * max(1.0, w0):
        return dict(out, p_value=1.0)
    tol = 1e-9 * max(1.0, observed)

    if null == "spherical" and all(s.k == 1 for s in spans):
        moving = [i for i, n in enumerate(norms) if n > 0.0]
        if len(moving) <= enumerate_up_to:
            patterns = 1 << len(moving)
            signs = np.ones((patterns, len(spans)))
            for bit, i in enumerate(moving):
                signs[:, i] = np.where((np.arange(patterns) >> bit) & 1, 1.0, -1.0)
            data = [bases[i][None, :] + (signs[:, i] * norms[i])[:, None] * perps[i][:, 0][None, :]
                    for i in range(len(spans))]
            values, _ = stat(data)
            return dict(out, p_value=float(np.mean(values >= observed - tol)), enumerated=patterns)

    rng = np.random.default_rng(seed)
    data = []
    for s, base, perp, norm in zip(spans, bases, perps, norms):
        v = rng.normal(size=(draws, s.k))
        if null == "spherical":
            if s.k == 1:
                v = np.where(v >= 0, 1.0, -1.0) * norm
            else:
                v *= norm / np.linalg.norm(v, axis=1, keepdims=True)
        data.append(base[None, :] + v @ perp.T)
    values, _ = stat(data)
    return dict(out, p_value=float((1 + np.sum(values >= observed - tol)) / (draws + 1)))


def zero_safe(h):
    """A weight that vanishes on a reading with nothing along the candidate: h / (1 + h).
    One member of a family; nothing selects it (T75)."""
    return h / (1.0 + h)


def floor_test(readings, mass, **options):
    """The larger of two p-values at one candidate: the Gaussian null with the zero-safe weight,
    and the spherical null unweighted. A candidate is kept when either test keeps it."""
    gaussian = conditional_test(readings, mass, null="gaussian", weight=zero_safe, **options)
    spherical = conditional_test(readings, mass, null="spherical", **options)
    return dict(p_value=max(gaussian["p_value"], spherical["p_value"]),
                gaussian=gaussian["p_value"], spherical=spherical["p_value"])


def kept(readings, masses, *, level=0.95, test=conditional_test, **options):
    """Which of the given candidate masses the test keeps at the stated level, and each p-value.
    Pointwise tests: the kept set may be unbounded, in pieces, or all of the candidates."""
    if not 0.0 < level < 1.0:
        raise ValueError("level lies strictly between zero and one")
    masses = np.atleast_1d(np.asarray(masses, dtype=float))
    seed = options.pop("seed", 0)
    p = np.array([test(readings, m, seed=seed + 7919 * i, **options)["p_value"] for i, m in enumerate(masses)])
    return p > 1.0 - level, p



# ---- two bodies (T74 (6)) ----------------------------------------------------------------------

SPLITS = 181


def _prepare(readings, mass, weight, scale):
    """One body at its candidate mass: the readings with a positive weight, split there."""
    scale = float(declared_scale(readings) if scale is None else scale) if readings else 1.0
    u0 = float(angle_of(mass, scale))
    spans, weights, bases, perps = [], [], [], []
    for r in readings:
        s = span_of(r)
        par, perp, along, _ = split(s, u0, scale)
        w = 1.0 if weight is None else float(weight(float(along @ along)))
        if w > 0.0:
            spans.append(s)
            weights.append(w)
            bases.append(par @ along)
            perps.append(perp)
    return dict(spans=spans, weights=weights, bases=bases, perps=perps, scale=scale, u0=u0)


def lump_test(first, second, total, first_mass, *, weight=None, draws=199, seed=0, scales=(None, None),
              grid=GRID, splits=SPLITS):
    """Exact conditional test that two independently measured bodies have masses `first_mass` and
    `total - first_mass` (harvested from E13 insight 7; T74 (6)).

    statistic: the composite curve at `total`, less the two bodies' own least values:
        inf over m1 in [0, total] of W1(m1) + W2(total - m1), minus inf W1, minus inf W2
    with the weights fixed at the candidate split. Given each reading's along-line part at the
    candidate split, the leftovers are independent standard normal vectors, so the statistic's
    law is known. The larger this p-value is over all splits of one total, the test of that total
    (`lump_kept`). A body with no reading of positive weight contributes nothing."""
    total, first_mass = float(total), float(first_mass)
    if not 0.0 < first_mass < total or not math.isfinite(total):
        raise ValueError("the candidate split has 0 < first_mass < total, both finite")
    one = _prepare(tuple(first), first_mass, weight, scales[0])
    two = _prepare(tuple(second), total - first_mass, weight, scales[1])
    if not one["spans"] and not two["spans"]:
        return dict(p_value=1.0, statistic=0.0, empty=True)
    stats_ = [_Statistic(b["spans"], b["weights"], b["scale"], b["u0"], "positive", grid) for b in (one, two)]
    share = np.linspace(0.0, 1.0, splits)
    step = float(share[1] - share[0])
    angles = [angle_of(share * total, one["scale"]), angle_of((1.0 - share) * total, two["scale"])]
    projectors = [[_projectors(s, a, b["scale"]) for s in b["spans"]] for b, a in zip((one, two), angles)]

    def statistic(data1, data2, n):
        least = stats_[0].least(data1, n) + stats_[1].least(data2, n)
        joined = _scan(one["weights"], data1, projectors[0], n) + _scan(two["weights"], data2, projectors[1], n)

        def at(t):
            return (stats_[0]._at(data1, angle_of(t * total, one["scale"]))
                    + stats_[1]._at(data2, angle_of((1.0 - t) * total, two["scale"])))
        j = np.argmin(joined, axis=1)
        best = np.minimum(at(share[j]), _golden(at, np.maximum(share[j] - step, 0.0), np.minimum(share[j] + step, 1.0)))
        return np.maximum(best - least, 0.0)

    observed = float(statistic([s.z[None, :] for s in one["spans"]], [s.z[None, :] for s in two["spans"]], 1)[0])
    if observed <= 1e-12:
        return dict(p_value=1.0, statistic=0.0, empty=False)
    rng = np.random.default_rng(seed)
    drawn = [[base[None, :] + rng.normal(size=(draws, s.k)) @ perp.T
              for s, base, perp in zip(b["spans"], b["bases"], b["perps"])] for b in (one, two)]
    values = statistic(drawn[0], drawn[1], draws)
    tol = 1e-9 * max(1.0, observed)
    return dict(p_value=float((1 + np.sum(values >= observed - tol)) / (draws + 1)), statistic=observed, empty=False)


def lump_kept(first, second, total, *, level=0.95, candidates=25, **options):
    """Whether a total is kept: the largest p-value over a grid of candidate splits of it. A test
    of the total at the stated level, since the true split is one of the splits (up to the grid)."""
    seed = options.pop("seed", 0)
    shares = (np.arange(candidates) + 0.5) / candidates
    p = max(lump_test(first, second, total, t * total, seed=seed + 7919 * i, **options)["p_value"]
            for i, t in enumerate(shares))
    return p > 1.0 - level, p


# ---- the weighted conditional score (T76) --------------------------------------------------------

def score_test(readings, mass, weight=None):
    """T45's conditional score with a weight on each reading, w_i = weight(|U_i|^2), fixed by the
    conditioning (harvested from E13 insight 9; T76). Given every U_i,
    sum w_i S_i / sqrt(sum w_i^2 C_i) is exactly standard normal at the true mass. weight None is
    T45 itself. An exactly zero reading has S_i = C_i = 0 under every weight."""
    from met.conditional_score import _geometry
    num = den = 0.0
    for r in readings:
        along, rest, tangent = _geometry(r, mass)
        w = 1.0 if weight is None else float(weight(float(along @ along)))
        num += w * float(rest @ tangent)
        den += w * w * float(tangent @ tangent)
    z = 0.0 if den == 0.0 else num / math.sqrt(den)
    return dict(normal_score=z, statistic=z * z, p_value=float(math.erfc(abs(z) / math.sqrt(2.0))),
                degenerate=den == 0.0)


def _chi_moment(function, k):
    from scipy import integrate, stats
    return integrate.quad(lambda h: stats.chi2.pdf(h, k) * function(h), 0.0, np.inf, limit=200)[0]


def weak_signal_efficiency(weight, k):
    """(E[H w])^2 / (E[H] E[H w^2]), H chi-square with k components: the squared noncentrality of
    the weighted score, relative to the unweighted one, when every reading's signal is weak. At
    most 1, with equality only for a constant weight (Cauchy-Schwarz)."""
    hw = _chi_moment(lambda h: h * weight(h), k)
    hw2 = _chi_moment(lambda h: h * weight(h) ** 2, k)
    return hw * hw / (k * hw2)


def drowning_count(weight, k, signal):
    """w(r^2)^2 r^2 / E[H w^2]: the number of noise-only readings whose weighted noise equals one
    reading of signal r, so that the pair's squared noncentrality is halved."""
    r2 = float(signal) ** 2
    return weight(r2) ** 2 * r2 / _chi_moment(lambda h: h * weight(h) ** 2, k)


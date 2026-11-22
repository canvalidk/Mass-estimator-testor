"""Exact conditional tests of one candidate mass on the full reading type, with no starting law.

Harvested 22 Nov 2026 from the exploration branch `explore/e13-unnormalized-evidence` (Can, 21 Nov;
insights 6 and 8), rewritten on met.full's reading type and checked by Claude (T75;
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

    def __call__(self, data):
        """LR for each draw. data: list over readings of arrays (D, 2k)."""
        draws = data[0].shape[0]
        w0 = self._at(data, np.full(draws, self.u0))
        on_grid = np.zeros((draws, self.angles.shape[0]))
        # the scan only chooses the bracket, so the fast form (energy less the explained part) will do
        for w, z, p in zip(self.weights, data, self.proj):
            outer = (z[:, :, None] * z[:, None, :]).reshape(draws, -1)
            on_grid += w * (np.einsum("di,di->d", z, z)[:, None] - outer @ p.reshape(p.shape[0], -1).T)
        j = np.argmin(on_grid, axis=1)
        best = self._at(data, self.angles[j])
        a, b = self.angles[j] - self.step, self.angles[j] + self.step
        if self.clip:
            a, b = np.maximum(a, self.lo), np.minimum(b, self.hi)
        c, d = b - _GR * (b - a), a + _GR * (b - a)
        fc, fd = self._at(data, c), self._at(data, d)
        for _ in range(GOLDEN_STEPS):
            left = fc < fd                                   # the least lies in [a, d]; else in [c, b]
            a, b = np.where(left, a, c), np.where(left, d, b)
            c_new, d_new = b - _GR * (b - a), a + _GR * (b - a)
            fresh = self._at(data, np.where(left, c_new, d_new))
            fc, fd = np.where(left, fresh, fd), np.where(left, fc, fresh)
            c, d = c_new, d_new
        least = np.minimum(np.minimum(best, np.minimum(fc, fd)), w0)
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

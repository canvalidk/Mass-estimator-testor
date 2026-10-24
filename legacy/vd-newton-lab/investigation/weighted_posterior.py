"""Flat-law mass posterior multiplied by a declared factor of latent mass.

The original ``flat_joint`` law is P proportional to exp(-Q/2) df dalpha dOmega.
In standardized units (both channel SDs one, so the mass unit is
sigma_F/sigma_a) write z = log M with M = f/alpha = tan(theta).  A factor
w(M) that depends on latent mass only defines P_w proportional to w(M) P.
The readout is still the composition-forced E_{P_w}[f] / E_{P_w}[alpha];
only the law changes.

Every factor below is a declared modelling choice.  Newton II, the Gaussian
error model and the readout theorem do not select one.

``SechTilt(lam, centre)``
    w = exp(lam * (sech(log(M/centre)) - 1)).  With centre 1 this is the
    fixed-experiment family rho_lambda = exp(lambda sin 2 theta) of the
    investigation README, section 4 (the constant factor cancels).
``LogNormalFactor(centre, width)``
    w = exp(-(log(M/centre))**2 / (2 width**2)): the flat law tempered by a
    stated belief that log mass lies near log(centre).
``LogNormalPrior(centre, width)``
    Replaces the flat law's implied uniform-angle log-mass reference
    sech(z)/2 by a normal density on log mass with the stated centre and
    width.  With exactly zero readings the posterior is that normal law.

Centres are standardized masses: physical mass divided by sigma_F/sigma_a.
They may be scalars or one value per reading.  Readings are standardized as
in ``comparison_posterior.infer_batch``; ``standardize`` converts physical
readings with isotropic supplied SDs.  Integration reuses that module's
panels, radial kernels and continuous summaries unchanged, adding panel
edges around each factor's centre.  This is investigation code for the
isotropic, known-noise, single-reading model only.
"""

import math

import numpy as np

from comparison_posterior import _LIMIT, _quadrature, _radial, _summary


_FIXED_EDGES = np.array([-_LIMIT, -20., -12., -8., -6., -4., -2., 0.,
                         2., 4., 6., 8., 12., 20., _LIMIT])
_ADAPTIVE = np.array([-8., -4., -2., 0., 2., 4., 8.])


def _log_sech_half(z):
    """log(sech(z)/2) = log(1/(e^z + e^-z)), stable for large |z|."""
    a = np.abs(z)
    return -(a + np.log1p(np.exp(-2 * a)))


def _positive_centres(centre):
    centre = np.asarray(centre, dtype=float)
    if centre.ndim > 1 or not np.all(np.isfinite(centre)) or np.any(centre <= 0):
        raise ValueError("centre must be a positive finite scalar or one value per reading")
    return centre


class _Factor:
    """Common handling of scalar or per-reading centres."""

    def __init__(self, centre):
        self.centre = _positive_centres(centre)
        self.log_centre = np.log(self.centre)

    def _shift(self, z):
        c = self.log_centre
        if c.ndim == 0:
            return z - c
        if z.shape[0] != c.shape[0]:
            raise ValueError("per-reading centres must match the number of readings")
        return z - c.reshape((-1,) + (1,) * (z.ndim - 1))

    @property
    def exchange_symmetric(self):
        # Force/acceleration exchange maps z to -z; only centre 1 keeps the
        # factor symmetric, which justifies the exact median for p == q.
        return self.log_centre.ndim == 0 and self.log_centre == 0

    def edges(self, count):
        return None


class SechTilt(_Factor):
    def __init__(self, lam, centre=1.):
        if isinstance(lam, bool) or not math.isfinite(lam):
            raise ValueError("lam must be finite")
        super().__init__(centre)
        self.lam = float(lam)

    def log_weight(self, z):
        a = np.exp(-np.abs(self._shift(z)))
        return self.lam * (2 * a / (1 + a * a) - 1)

    def edges(self, count):
        # Tilts with lam <= 4 vary on a log-mass scale of about one or more;
        # the fixed panels resolve them (checked by order refinement in tests).
        if self.lam <= 4:
            return None
        width = 2 / math.sqrt(self.lam)
        offsets = width * np.array([-4., -2., -1., 0., 1., 2., 4.])
        c = np.broadcast_to(self.log_centre, (count,))
        return c[:, None] + offsets

    def describe(self):
        return {"family": "sech_tilt", "lam": self.lam, "centre": self.centre.tolist()}


class LogNormalFactor(_Factor):
    def __init__(self, centre, width):
        if isinstance(width, bool) or not math.isfinite(width) or width <= 0:
            raise ValueError("width must be positive and finite")
        super().__init__(centre)
        self.width = float(width)

    def log_weight(self, z):
        x = self._shift(z)
        return -x * x / (2 * self.width**2)

    def edges(self, count):
        offsets = self.width * np.array([-8., -4., -2., -1., 0., 1., 2., 4., 8.])
        c = np.broadcast_to(self.log_centre, (count,))
        return c[:, None] + offsets

    def describe(self):
        return {"family": "lognormal_factor", "width": self.width, "centre": self.centre.tolist()}


class LogNormalPrior(LogNormalFactor):
    def log_weight(self, z):
        return super().log_weight(z) - _log_sech_half(z)

    def describe(self):
        return {"family": "lognormal_prior", "width": self.width, "centre": self.centre.tolist()}


def standardize(force, acceleration, force_sd, acceleration_sd, rtol=1e-9):
    """Divide physical readings by isotropic supplied SDs; return (x, y, s)."""
    force = np.atleast_2d(np.asarray(force, dtype=float))
    acceleration = np.atleast_2d(np.asarray(acceleration, dtype=float))
    if force.shape != acceleration.shape or force.shape[1] != 3:
        raise ValueError("force and acceleration must have identical (n,3) shapes")
    sd = []
    for value in (force_sd, acceleration_sd):
        value = np.broadcast_to(np.asarray(value, dtype=float), force.shape) if np.ndim(value) != 1 \
            else np.broadcast_to(np.asarray(value, dtype=float)[:, None], force.shape)
        if not np.all(np.isfinite(value)) or np.any(value <= 0):
            raise ValueError("supplied SDs must be positive and finite")
        if np.any(value.max(axis=1) > value.min(axis=1) * (1 + rtol)):
            raise ValueError("this reduction requires isotropic supplied SDs within each reading")
        sd.append(value[:, 0])
    return force / sd[0][:, None], acceleration / sd[1][:, None], sd[0] / sd[1]


def _setup(force, acceleration, factors, order):
    force = np.atleast_2d(np.asarray(force, dtype=float))
    acceleration = np.atleast_2d(np.asarray(acceleration, dtype=float))
    if force.shape != acceleration.shape or force.ndim != 2 or force.shape[1] != 3:
        raise ValueError("force and acceleration must have identical (n,3) shapes")
    if not len(force) or not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration)):
        raise ValueError("at least one finite reading is required")
    if order < 24:
        raise ValueError("order must be >= 24")
    p, q = np.sum(force**2, axis=1), np.sum(acceleration**2, axis=1)
    dot = np.sum(force * acceleration, axis=1)
    peak = np.where(dot >= 0, (p + q + np.hypot(p - q, 2 * dot)) / 2, np.maximum(p, q))
    centre = .5 * np.log((p + 1) / (q + 1))
    scale = np.sqrt(1 / (p + 1) + 1 / (q + 1))
    # Identical to comparison_posterior.infer_batch, plus factor-centred edges.
    parts = [np.broadcast_to(_FIXED_EDGES, (len(p), len(_FIXED_EDGES))),
             np.clip(centre[:, None] + scale[:, None] * _ADAPTIVE, -_LIMIT + 1e-8, _LIMIT - 1e-8)]
    for factor in factors:
        if factor is not None:
            extra = factor.edges(len(p))
            if extra is not None:
                parts.append(np.clip(extra, -_LIMIT + 1e-8, _LIMIT - 1e-8))
    edges = np.sort(np.concatenate(parts, axis=1), axis=1)
    half = np.diff(edges, axis=1) / 2
    midpoint = (edges[:, 1:] + edges[:, :-1]) / 2
    x, weights, vander, projection = _quadrature(order)
    z = midpoint[:, :, None] + half[:, :, None] * x
    theta = np.arctan(np.exp(z))
    sn, cs = np.sin(theta), np.cos(theta)
    jacobian = 1 / (2 * np.cosh(z))
    h2 = np.maximum(0, p[:, None, None] * sn**2 + q[:, None, None] * cs**2
                    + 2 * dot[:, None, None] * sn * cs)
    return dict(p=p, q=q, dot=dot, peak=peak, edges=edges, half=half, z=z, sn=sn, cs=cs,
                jacobian=jacobian, h2=h2, weights=weights, vander=vander, projection=projection)


def _factor_values(factor, z, count):
    if factor is None:
        return np.ones(z.shape), np.ones(count), np.ones(count)
    ends = np.broadcast_to(np.array([-_LIMIT, _LIMIT]), (count, 2))
    ends_w = np.exp(factor.log_weight(ends))
    return np.exp(factor.log_weight(z)), ends_w[:, 0], ends_w[:, 1]


def weighted_points(force, acceleration, factors, order=96):
    """Ratio-of-means points, shape (len(factors), n); ``None`` means flat.

    Only first radial moments are needed, so no error functions are
    evaluated.  All factors share one panel set per reading.
    """
    factors = list(factors)
    g = _setup(force, acceleration, factors, order)
    first = np.exp((g["h2"] - g["peak"][:, None, None]) / 2)
    base = first * g["jacobian"] * g["weights"] * g["half"][:, :, None]
    num_base, den_base = base * g["sn"], base * g["cs"]
    tail_left = np.exp((g["q"] - g["peak"]) / 2) * math.exp(-_LIMIT)
    tail_right = np.exp((g["p"] - g["peak"]) / 2) * math.exp(-_LIMIT)
    count = len(g["p"])
    result = np.empty((len(factors), count))
    tilt_cache = {}
    for k, factor in enumerate(factors):
        if isinstance(factor, SechTilt) and factor.log_centre.ndim == 0:
            # Tilts sharing a centre share sech(z - c) - 1; only exp(lam * .) differs.
            c = float(factor.log_centre)
            if c not in tilt_cache:
                a = np.exp(-np.abs(g["z"] - c))
                ends = np.exp(-np.abs(np.array([-_LIMIT, _LIMIT]) - c))
                tilt_cache[c] = (2 * a / (1 + a * a) - 1, 2 * ends / (1 + ends * ends) - 1)
            shape, ends = tilt_cache[c]
            w = np.exp(factor.lam * shape)
            w_left, w_right = np.exp(factor.lam * ends)
        else:
            w, w_left, w_right = _factor_values(factor, g["z"], count)
        numerator = np.sum(num_base * w, axis=(1, 2)) + tail_right * w_right
        denominator = np.sum(den_base * w, axis=(1, 2)) + tail_left * w_left
        result[k] = numerator / denominator
    return result


def infer_weighted(force, acceleration, factor=None, truth=1., order=96):
    """Point and continuous log-mass summaries for one declared factor.

    Returns ``{'point', 'reciprocal_root', 'distribution'}``; the
    distribution keys match ``comparison_posterior.infer_batch``.  With
    ``factor=None`` every value equals that function's ``flat_joint`` output.
    Endpoint tails beyond |z| = 32 hold each factor at its endpoint value.
    """
    g = _setup(force, acceleration, [factor], order)
    count = len(g["p"])
    truth = np.asarray(truth, dtype=float)
    if truth.shape not in ((), (count,)) or not np.all(np.isfinite(truth)) or np.any(truth <= 0):
        raise ValueError("truth must be finite and positive")
    radial = _radial(g["h2"], g["peak"][:, None, None])
    endpoints = _radial(np.stack((g["q"], g["p"]), axis=1), g["peak"][:, None])
    w, w_left, w_right = _factor_values(factor, g["z"], count)
    first = radial[2]
    numerator = np.sum(np.sum(first * g["sn"] * g["jacobian"] * w * g["weights"], axis=2) * g["half"], axis=1)
    denominator = np.sum(np.sum(first * g["cs"] * g["jacobian"] * w * g["weights"], axis=2) * g["half"], axis=1)
    numerator += endpoints[2][:, 1] * math.exp(-_LIMIT) * w_right
    denominator += endpoints[2][:, 0] * math.exp(-_LIMIT) * w_left
    y = np.broadcast_to(np.log(truth), (count,))
    ytheta = np.arctan(np.broadcast_to(truth, (count,)))
    ys, yc = np.sin(ytheta), np.cos(ytheta)
    yh2 = np.maximum(0, g["p"] * ys**2 + g["q"] * yc**2 + 2 * g["dot"] * ys * yc)
    yradial = _radial(yh2, g["peak"])
    log_jacobian = -abs(y) - np.log1p(np.exp(-2 * abs(y)))
    truth_log_kernel = (yh2 - g["peak"]) / 2 + np.log(yradial[4]) + log_jacobian
    if factor is not None:
        truth_log_kernel = truth_log_kernel + factor.log_weight(y.copy())
    summary, root_point = _summary(radial[0] * g["jacobian"] * w, g["edges"], g["z"], g["half"],
                                   g["weights"], g["vander"], g["projection"],
                                   endpoints[0][:, 0] * w_left, endpoints[0][:, 1] * w_right,
                                   y, truth_log_kernel)
    if factor is None or factor.exchange_symmetric:
        summary["log_median"] = np.where(g["p"] == g["q"], 0., summary["log_median"])
    return {"point": numerator / denominator, "reciprocal_root": root_point, "distribution": summary}

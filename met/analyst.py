"""The analyst: stages 2-5 of the pipeline.

An estimator is a declared model plus a readout. It sees only the readings and
the standard deviations it was supplied, never the world's settings.

    [reduce] -> per-reading law -> combine (rule + prior) -> readouts

or, for a direct rule, readings -> pooled pair -> one number.
"""

import math

import numpy as np

from .law import NUISANCES, REFERENCES, ReadingSet

REDUCES = ("none", "pool_pair")
RULES = ("single", "likelihood_product", "posterior_product")
COORDINATES = ("mass", "log_mass")
POINT_READOUTS = ("ratio_of_means", "median", "geometric", "reciprocal_root")
DIRECT_RULES = ("norm_ratio", "dot_acceleration", "dot_force")
PRIOR_FLAT = ("flat_mass", "flat_log_mass", "flat_inverse_mass")
PRIOR_PARAMETRIC = ("uniform_angle", "sech_tilt", "lognormal")

# Numerical settings: the only settings with defaults. Recorded in every output.
DEFAULT_NUMERICS = {"span": 60.0, "coarse_step": 0.05, "cutoff": 40.0, "grid_points": 2001}


# ---------------------------------------------------------------- priors

def _center(value, readings):
    """log of a prior's centre: a positive number, or 'first_reading' = log s of reading 1."""
    if value == "first_reading":
        return np.log(readings.s[:, 0])[:, None]
    return math.log(float(value))


def prior_log_density(factors, u, readings):
    """Sum of declared prior factors, as a log density in u = log m. u: (B, G)."""
    total = np.zeros_like(u)
    for factor in factors:
        if isinstance(factor, str):
            kind, params = factor, {}
        else:
            (kind, params), = factor.items()
        if kind == "flat_mass":
            total = total + u
        elif kind == "flat_log_mass":
            pass
        elif kind == "flat_inverse_mass":
            total = total - u
        elif kind == "uniform_angle":
            z = u - _center(params["center"], readings)
            total = total - np.logaddexp(z, -z)
        elif kind == "sech_tilt":
            z = u - _center(params["center"], readings)
            total = total + float(params["lambda"]) * (1.0 / np.cosh(np.clip(z, -700, 700)) - 1.0)
        elif kind == "lognormal":
            z = u - _center(params["center"], readings)
            width = float(params["width"])
            total = total - z * z / (2.0 * width * width)
        else:
            raise ValueError(f"unknown prior factor {kind!r}")
    return total


def validate_prior(factors, rule):
    if not isinstance(factors, list):
        raise ValueError("combine.prior must be a list of factors (use [] only with posterior_product or single)")
    if rule == "likelihood_product" and not factors:
        raise ValueError("likelihood_product needs at least one prior factor: the prior is never implicit")
    for factor in factors:
        if isinstance(factor, str):
            if factor not in PRIOR_FLAT:
                raise ValueError(f"unknown prior factor {factor!r}")
            continue
        if not isinstance(factor, dict) or len(factor) != 1:
            raise ValueError(f"prior factor must be a name or a one-key mapping, got {factor!r}")
        (kind, params), = factor.items()
        required = {"uniform_angle": {"center"}, "sech_tilt": {"lambda", "center"},
                    "lognormal": {"center", "width"}}.get(kind)
        if required is None:
            raise ValueError(f"unknown prior factor {kind!r}")
        if not isinstance(params, dict) or set(params) != required:
            raise ValueError(f"prior factor {kind} needs exactly: {', '.join(sorted(required))}")
        center = params["center"]
        if center != "first_reading" and not (isinstance(center, (int, float)) and center > 0):
            raise ValueError(f"prior centre must be a positive number or 'first_reading', got {center!r}")


# ---------------------------------------------------------------- joint law

def joint_log_density(readings, u, law, combine):
    """Joint log density in u (up to a constant) and A(u) = sum_i E[alpha_i | m]."""
    log_like, log_ref, cond_alpha = readings.reading_terms(u, law["nuisance"])
    n = readings.readings
    rule = combine["rule"]
    if rule == "single":
        if n != 1:
            raise ValueError("rule 'single' needs exactly one reading per series after reduce")
        lp = log_like[:, 0] + log_ref[:, 0]
    elif rule == "likelihood_product":
        lp = np.sum(log_like, axis=1)
    elif rule == "posterior_product":
        lp = np.sum(log_like + log_ref, axis=1)
        if combine["coordinate"] == "mass":
            lp = lp + (1 - n) * u
    else:
        raise ValueError(f"unknown rule {rule!r}")
    lp = lp + prior_log_density(combine["prior"], u, readings)
    return lp, np.sum(cond_alpha, axis=1)


def _trapezoid_weights(n):
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    return w


def _cdf(density, h):
    """CDF at the nodes of a uniform grid, fourth-order accurate.

    Cumulative trapezoid with the Euler-Maclaurin endpoint correction
    -h^2/12 (f'(u_k) - f'(u_0)); f' from second-order differences.
    Normalised so the last node is exactly 1.
    """
    cells = 0.5 * (density[:, 1:] + density[:, :-1]) * h[:, None]
    trap = np.concatenate([np.zeros((density.shape[0], 1)), np.cumsum(cells, axis=1)], axis=1)
    slope = np.gradient(density, axis=1, edge_order=2) / h[:, None]
    cdf = trap - (h[:, None] ** 2 / 12.0) * (slope - slope[:, :1])
    return cdf / cdf[:, -1:]


def _quantiles(u, density, cdf, probabilities, newton_steps=4):
    """Invert the CDF by cubic Hermite interpolation (values and slopes) in each cell."""
    h = u[:, 1] - u[:, 0]
    out = np.empty((u.shape[0], len(probabilities)))
    rows = np.arange(u.shape[0])
    for j, p in enumerate(probabilities):
        k = np.clip(np.sum(cdf < p, axis=1) - 1, 0, u.shape[1] - 2)
        c0, c1 = cdf[rows, k], cdf[rows, k + 1]
        m0, m1 = density[rows, k] * h, density[rows, k + 1] * h
        gap = c1 - c0
        s = np.clip(np.where(gap > 0, (p - c0) / np.where(gap > 0, gap, 1.0), 0.5), 0.0, 1.0)
        for _ in range(newton_steps):
            s2, s3 = s * s, s * s * s
            value = ((2 * s3 - 3 * s2 + 1) * c0 + (s3 - 2 * s2 + s) * m0
                     + (-2 * s3 + 3 * s2) * c1 + (s3 - s2) * m1)
            slope = ((6 * s2 - 6 * s) * c0 + (3 * s2 - 4 * s + 1) * m0
                     + (-6 * s2 + 6 * s) * c1 + (3 * s2 - 2 * s) * m1)
            step = np.where(slope > 0, (value - p) / np.where(slope > 0, slope, 1.0), 0.0)
            s = np.clip(s - step, 0.0, 1.0)
        out[:, j] = u[rows, k] + s * h
    return out


def posterior_summaries(readings, law, combine, readouts, numerics=None):
    """Evaluate the joint law on an adaptive log-mass grid and compute readouts.

    Returns {readout: array}. Points are (B,), intervals (B, 2); a series whose
    law is improper or runs off the grid gets NaN in every readout, and
    'unresolved' is True for it.
    """
    num = dict(DEFAULT_NUMERICS, **(numerics or {}))
    b = readings.series
    center = np.log(readings.s[:, 0])
    offsets = np.arange(-num["span"], num["span"] + num["coarse_step"] / 2, num["coarse_step"])
    coarse = center[:, None] + offsets[None, :]
    lp, _ = joint_log_density(readings, coarse, law, combine)
    top = np.max(lp, axis=1, keepdims=True)
    keep = lp > top - num["cutoff"]
    first = np.argmax(keep, axis=1)
    last = keep.shape[1] - 1 - np.argmax(keep[:, ::-1], axis=1)
    unresolved = (~np.isfinite(top[:, 0])) | (first == 0) | (last == keep.shape[1] - 1)
    rows = np.arange(b)
    lo = coarse[rows, np.maximum(first - 2, 0)]
    hi = coarse[rows, np.minimum(last + 2, keep.shape[1] - 1)]

    g = int(num["grid_points"])
    t = np.linspace(0.0, 1.0, g)
    u = lo[:, None] + (hi - lo)[:, None] * t[None, :]
    lp, acc = joint_log_density(readings, u, law, combine)
    dens = np.exp(lp - np.max(lp, axis=1, keepdims=True))
    h = (hi - lo) / (g - 1)
    w = dens * _trapezoid_weights(g)[None, :] * h[:, None]
    z = np.sum(w, axis=1, keepdims=True)
    w = w / z
    dens = dens / z
    mass = np.exp(u)

    out = {"unresolved": unresolved}
    wanted = set(readouts)
    if "ratio_of_means" in wanted:
        out["ratio_of_means"] = np.sum(w * mass * acc, axis=1) / np.sum(w * acc, axis=1)
    if "geometric" in wanted or "log_sd" in wanted:
        mean_u = np.sum(w * u, axis=1)
        if "geometric" in wanted:
            out["geometric"] = np.exp(mean_u)
        if "log_sd" in wanted:
            out["log_sd"] = np.sqrt(np.maximum(np.sum(w * u * u, axis=1) - mean_u**2, 0.0))
    if "reciprocal_root" in wanted:
        out["reciprocal_root"] = np.sum(w * np.exp(0.5 * u), axis=1) / np.sum(w * np.exp(-0.5 * u), axis=1)
    probabilities, slots = [], []
    if "median" in wanted:
        probabilities.append(0.5)
        slots.append(("median", None))
    for name in wanted:
        if name.startswith("interval_"):
            content = interval_content(name)
            probabilities += [(1 - content) / 2, (1 + content) / 2]
            slots.append((name, content))
    if probabilities:
        q = np.exp(_quantiles(u, dens, _cdf(dens, h), probabilities))
        j = 0
        for name, content in slots:
            if content is None:
                out[name] = q[:, j]
                j += 1
            else:
                out[name] = q[:, j:j + 2]
                j += 2
    for name, value in out.items():
        if name != "unresolved":
            value[unresolved] = np.nan
    return out


def interval_content(name):
    try:
        pct = float(name[len("interval_"):])
    except ValueError:
        raise ValueError(f"bad interval readout {name!r}; use e.g. interval_95") from None
    if not 0 < pct < 100:
        raise ValueError(f"interval content must be between 0 and 100: {name!r}")
    return pct / 100.0


# ---------------------------------------------------------------- direct rules

def direct_rule(readings, rule):
    """A rule that skips the law, applied to the inverse-variance pooled pair."""
    pair = readings.pooled()
    f, a = pair.force[:, 0, :], pair.acceleration[:, 0, :]
    if rule == "norm_ratio":
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.linalg.norm(f, axis=1) / np.linalg.norm(a, axis=1)
    fa = np.sum(f * a, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        if rule == "dot_acceleration":
            value = fa / np.sum(a * a, axis=1)
        elif rule == "dot_force":
            value = np.sum(f * f, axis=1) / fa
        else:
            raise ValueError(f"unknown direct rule {rule!r}")
    return np.where(value > 0, value, np.nan)


# ---------------------------------------------------------------- estimators

class Estimator:
    """A validated estimator specification from a preset."""

    def __init__(self, spec):
        if not isinstance(spec, dict) or "name" not in spec:
            raise ValueError("each estimator needs a name")
        self.name = str(spec["name"])
        self.spec = spec
        where = f"estimator {self.name!r}"
        if "direct" in spec:
            extra = set(spec) - {"name", "direct"}
            if extra:
                raise ValueError(f"{where}: a direct rule takes only 'name' and 'direct' (got {sorted(extra)})")
            if spec["direct"] not in DIRECT_RULES:
                raise ValueError(f"{where}: direct rule must be one of {DIRECT_RULES}")
            self.direct = spec["direct"]
            self.readouts = [self.direct]
            return
        self.direct = None
        required = {"name", "reduce", "law", "combine", "readouts"}
        missing = required - set(spec)
        if missing:
            raise ValueError(f"{where}: missing required setting(s) {sorted(missing)}")
        extra = set(spec) - required - {"numerics"}
        if extra:
            raise ValueError(f"{where}: unknown setting(s) {sorted(extra)}")
        if spec["reduce"] not in REDUCES:
            raise ValueError(f"{where}: reduce must be one of {REDUCES}")
        law = spec["law"]
        if not isinstance(law, dict) or set(law) != {"reference", "nuisance"}:
            raise ValueError(f"{where}: law needs exactly reference and nuisance")
        if law["reference"] not in REFERENCES:
            raise ValueError(f"{where}: reference must be one of {REFERENCES} (others not built yet)")
        if law["nuisance"] not in NUISANCES:
            raise ValueError(f"{where}: nuisance must be one of {NUISANCES}")
        combine = spec["combine"]
        if not isinstance(combine, dict) or "rule" not in combine or "prior" not in combine:
            raise ValueError(f"{where}: combine needs rule and prior")
        rule = combine["rule"]
        if rule not in RULES:
            raise ValueError(f"{where}: rule must be one of {RULES}")
        allowed = {"rule", "prior"} | ({"coordinate"} if rule == "posterior_product" else set())
        if set(combine) != allowed:
            raise ValueError(f"{where}: combine for {rule} needs exactly {sorted(allowed)}")
        if rule == "posterior_product" and combine["coordinate"] not in COORDINATES:
            raise ValueError(f"{where}: coordinate must be one of {COORDINATES}")
        try:
            validate_prior(combine["prior"], rule)
        except ValueError as error:
            raise ValueError(f"{where}: {error}") from None
        readouts = spec["readouts"]
        if not isinstance(readouts, list) or not readouts:
            raise ValueError(f"{where}: readouts must be a nonempty list")
        for name in readouts:
            if name in POINT_READOUTS or name == "log_sd":
                continue
            if isinstance(name, str) and name.startswith("interval_"):
                interval_content(name)
                continue
            raise ValueError(f"{where}: unknown readout {name!r}")
        self.readouts = list(readouts)
        self.numerics = dict(DEFAULT_NUMERICS, **spec.get("numerics", {}))

    def evaluate(self, readings):
        if self.direct is not None:
            return {self.direct: direct_rule(readings, self.direct),
                    "unresolved": np.zeros(readings.series, dtype=bool)}
        if self.spec["reduce"] == "pool_pair":
            readings = readings.pooled()
        return posterior_summaries(readings, self.spec["law"], self.spec["combine"],
                                   self.readouts, self.numerics)

    def point_readouts(self):
        return [r for r in self.readouts if r in POINT_READOUTS or r in DIRECT_RULES]

    def interval_readouts(self):
        return [r for r in self.readouts if isinstance(r, str) and r.startswith("interval_")]


def estimate(force, acceleration, force_sd, acceleration_sd, readouts=("ratio_of_means",),
             law=None, combine=None, numerics=None):
    """Convenience call for one series of readings (N x d arrays, or d for one reading)."""
    force = np.atleast_2d(np.asarray(force, dtype=float))
    acceleration = np.atleast_2d(np.asarray(acceleration, dtype=float))
    n = force.shape[0]
    readings = ReadingSet(force[None], acceleration[None], np.full((1, n), force_sd),
                          np.full((1, n), acceleration_sd))
    law = law or {"reference": "flat", "nuisance": "radius"}
    combine = combine or {"rule": "single", "prior": []}
    result = posterior_summaries(readings, law, combine, list(readouts), numerics)
    return {k: (v[0] if v.ndim else v) for k, v in result.items()}

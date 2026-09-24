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
RULES = ("single", "likelihood_product", "posterior_product", "sequential")
CARRIES = ("exact", "curve_only", "lognormal_fit", "tilt_fit", "vonmises_state")
COORDINATES = ("mass", "log_mass")
POINT_READOUTS = ("ratio_of_means", "median", "geometric", "reciprocal_root")
DIRECT_RULES = ("norm_ratio", "dot_acceleration", "dot_force")
PRIOR_FLAT = ("flat_mass", "flat_log_mass", "flat_inverse_mass")
PRIOR_PARAMETRIC = ("uniform_angle", "sech_tilt", "lognormal")

# Numerical settings: the only settings with defaults. Recorded in every output.
#   span, coarse_step   first search grid: log s of reading 1 +- span, this step
#   cutoff              the law is cut where it falls exp(-cutoff) below its peak
#   min_support_steps   refine the search until the support spans this many steps
#   refine_points       points per refinement pass
#   grid_points         points of the final uniform grid
DEFAULT_NUMERICS = {"span": 60.0, "coarse_step": 0.05, "cutoff": 40.0, "grid_points": 2001,
                    "min_support_steps": 100, "refine_points": 401}
_INTEGER_NUMERICS = {"grid_points": 201, "min_support_steps": 10, "refine_points": 101}


def _number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number, got {value!r}")
    return float(value)


def validate_numerics(numerics):
    if not isinstance(numerics, dict):
        raise ValueError("numerics must be a mapping")
    unknown = set(numerics) - set(DEFAULT_NUMERICS)
    if unknown:
        raise ValueError(f"unknown numerics setting(s) {sorted(unknown)}; allowed: {sorted(DEFAULT_NUMERICS)}")
    for key, value in numerics.items():
        if key in _INTEGER_NUMERICS:
            if type(value) is not int or value < _INTEGER_NUMERICS[key]:
                raise ValueError(f"numerics.{key} must be an integer >= {_INTEGER_NUMERICS[key]}")
        elif _number(value, f"numerics.{key}") <= 0:
            raise ValueError(f"numerics.{key} must be positive")


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
        if center != "first_reading" and not (
                not isinstance(center, bool) and isinstance(center, (int, float))
                and math.isfinite(center) and center > 0):
            raise ValueError(f"prior centre must be a positive number or 'first_reading', got {center!r}")
        if kind == "sech_tilt":
            _number(params["lambda"], "sech_tilt lambda")
        if kind == "lognormal" and _number(params["width"], "lognormal width") <= 0:
            raise ValueError("lognormal width must be positive")


# Tail slopes of the joint log density in u = log m, as u -> +inf and u -> -inf.
# The radial kernel tends to a constant at both ends whatever the data, so
# whether a law can be normalised is decided by the settings and N alone.
_LIKE_SLOPES = {"radius": (0, 0), "acceleration": (-2, 0), "force": (0, 2)}
_REF_SLOPES = {"radius": (-1, 1), "acceleration": (1, 1), "force": (-1, -1)}
_PRIOR_SLOPES = {"flat_mass": (1, 1), "flat_log_mass": (0, 0), "flat_inverse_mass": (-1, -1),
                 "uniform_angle": (-1, 1), "sech_tilt": (0, 0)}


def tail_slopes(law, combine, n):
    """(slope as u -> +inf, slope as u -> -inf, confined) for n readings after reduce.

    `confined` is True when a log-normal factor makes both tails decay regardless.
    The law is proper iff confined, or the first slope is < 0 and the second > 0.
    """
    like = _LIKE_SLOPES[law["nuisance"]]
    ref = _REF_SLOPES[law["nuisance"]]
    rule = combine["rule"]
    if rule == "sequential":
        carry = combine["carry"]
        if carry == "lognormal_fit":
            return like[0], like[1], True
        if carry == "tilt_fit":
            return like[0] - 1, like[1] + 1, False
        if carry == "vonmises_state":
            plus, minus, confined = tail_slopes(law, combine["old"], 0)
            return plus + like[0], minus + like[1], confined
        plus, minus, confined = tail_slopes(law, combine["old"], n - 1)
        return plus + like[0], minus + like[1], confined
    if rule == "single":
        plus, minus = like[0] + ref[0], like[1] + ref[1]
    elif rule == "likelihood_product":
        plus, minus = n * like[0], n * like[1]
    elif combine["coordinate"] == "mass":
        plus, minus = n * (like[0] + ref[0]) + (1 - n), n * (like[1] + ref[1]) + (1 - n)
    else:
        plus, minus = n * (like[0] + ref[0]), n * (like[1] + ref[1])
    confined = False
    for factor in combine["prior"]:
        kind = factor if isinstance(factor, str) else next(iter(factor))
        if kind == "lognormal":
            confined = True
            continue
        plus += _PRIOR_SLOPES[kind][0]
        minus += _PRIOR_SLOPES[kind][1]
    return plus, minus, confined


def check_proper(law, combine, n):
    """Raise if the declared law cannot be normalised for n readings."""
    if combine["rule"] == "sequential":
        if n < 2:
            raise ValueError(f"rule 'sequential' needs at least two readings (old and new), but this world gives {n}")
        try:
            check_proper(law, combine["old"], n - 1)
        except ValueError as error:
            raise ValueError(f"the old readings' law: {error}") from None
    if combine["rule"] == "single" and n != 1:
        raise ValueError(f"rule 'single' needs exactly one reading after reduce, but this world gives {n}")
    plus, minus, confined = tail_slopes(law, combine, n)
    if not confined and not (plus < 0 and minus > 0):
        raise ValueError(
            f"this law cannot be normalised for N = {n}: its log density has slope {plus} as m -> infinity "
            f"and {minus} as m -> 0 (in log m); it needs < 0 and > 0. Change the prior or the rule.")


# ---------------------------------------------------------------- joint law

def joint_log_density(readings, u, law, combine):
    """Joint log density in u (up to a constant) and A(u) = sum_i E[alpha_i | m]."""
    rule = combine["rule"]
    if rule == "sequential":
        return _sequential_log_density(readings, u, law, combine)
    log_like, log_ref, cond_alpha = readings.reading_terms(u, law["nuisance"], law["reference"])
    n = readings.readings
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


# ---------------------------------------------------------------- the sequential rule
#
# The last reading of a series is "new"; the ones before it are "old". Under a
# rule that counts the prior once, the exact joint law factorises as
#
#     p_all(m) = p_old(m) * L_new(m),        A_all(m) = A_old(m) + E[alpha_new | m],
#
# so the old information enters exactly where a prior does, plus the old
# conditional acceleration A_old that the ratio-of-means readout needs.
# `carry` says what of the old information is kept:
#
#   exact          p_old and A_old, evaluated exactly (equals the rule applied to all readings)
#   curve_only     p_old exactly, A_old dropped
#   lognormal_fit  p_old replaced by a normal law in log m with the same mean and SD
#   tilt_fit       p_old replaced by the tilted half-Cauchy  sech(z) exp(lambda (sech z - 1)),
#                  z = log(m / m0), with m0 and lambda set so its mean and SD in log m match
#
#   vonmises_state the old readings' likelihoods replaced by their von Mises parts in the
#                  doubled angle 2 theta_i (m = s_i tan theta_i):
#                      log L_i ~ ((Q_i - P_i)/4) cos 2theta_i + (D_i/2) sin 2theta_i,
#                  i.e. the complex number c_i = (1/4) sum_k (y_k + i x_k)^2, times the old prior.
#                  With a common s the old state is the single vector sum_i c_i. Exact
#                  under the cartesian reference in every d, and under the flat
#                  reference in d = 2; with flat in d = 1, 3 it drops the factor
#                  ~ h^(2-d). Needs the old rule to be likelihood_product, radius.
#
# The fits keep two numbers of the old information, in the slots the prior
# parameters occupy. A_old is dropped by every carry except `exact`.

def tilt_sd(lam, points=40001):
    """SD in z of sech(z) exp(lam (sech z - 1)), on a grid scaled to the law's width."""
    width = 80.0 if lam < 50 else 15.0 / math.sqrt(1.0 + lam)
    z = np.linspace(-width, width, points)
    sech = 1.0 / np.cosh(z)
    logd = np.log(sech) + lam * (sech - 1.0)
    d = np.exp(logd - logd.max())
    return math.sqrt(np.sum(d * z * z) / np.sum(d))


def _tilt_table():
    lams = np.concatenate([np.linspace(-40, 0, 401)[:-1], np.geomspace(1e-3, 2e5, 1200)])
    lams = np.unique(np.concatenate([lams, [0.0]]))
    return lams, np.array([tilt_sd(lam) for lam in lams])


_TILT = None


def tilt_lambda_for_sd(sd):
    """lambda such that sech(z) exp(lambda (sech z - 1)) has SD `sd` in z (vectorised).

    Beyond the table, lambda ~ 1/sd^2 - 1 (the law is then nearly normal).
    Wider than lambda = -40 allows is clipped to -40.
    """
    global _TILT
    if _TILT is None:
        _TILT = _tilt_table()
    lams, sds = _TILT
    sd = np.asarray(sd, dtype=float)
    inside = np.interp(-sd, -sds, lams)
    return np.where(sd < sds[-1], 1.0 / np.maximum(sd, 1e-300) ** 2 - 1.0, inside)


def _old_fit(readings, law, combine):
    """Mean and SD of log m under the old readings' law, per series (cached on the readings)."""
    key = ("old_fit", id(combine))
    if key not in readings.series_data:
        old = readings.take_readings(slice(0, readings.readings - 1))
        old.series_data = {}
        out = posterior_summaries(old, law, combine["old"], ["geometric", "log_sd"])
        readings.series_data[key] = np.stack([np.log(out["geometric"]), out["log_sd"]], axis=1)
    return readings.series_data[key]


def _sequential_log_density(readings, u, law, combine):
    carry = combine["carry"]
    n = readings.readings
    new = readings.take_readings(slice(n - 1, n))
    like_new, _, cond_new = new.reading_terms(u, law["nuisance"], law["reference"])
    lp = like_new[:, 0]
    acc = cond_new[:, 0]
    if carry in ("exact", "curve_only"):
        old = readings.take_readings(slice(0, n - 1))
        lp_old, acc_old = joint_log_density(old, u, law, combine["old"])
        lp = lp + lp_old
        if carry == "exact":
            acc = acc + acc_old
        return lp, acc
    if carry == "vonmises_state":
        old = readings.take_readings(slice(0, n - 1))
        z_old = u[:, None, :] - np.log(old.s)[:, :, None]
        theta = np.arctan(np.exp(np.clip(z_old, -700, 700)))
        vm = (((old.Q - old.P) / 4.0)[:, :, None] * np.cos(2 * theta)
              + (old.D / 2.0)[:, :, None] * np.sin(2 * theta))
        return lp + np.sum(vm, axis=1) + prior_log_density(combine["old"]["prior"], u, readings), acc
    fit = _old_fit(readings, law, combine)
    mu, sd = fit[:, 0:1], fit[:, 1:2]
    z = u - mu
    if carry == "lognormal_fit":
        return lp - z * z / (2.0 * sd * sd), acc
    lam = tilt_lambda_for_sd(sd)
    return lp - np.logaddexp(z, -z) + lam * (1.0 / np.cosh(np.clip(z, -700, 700)) - 1.0), acc


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


def _support(readings, grid, law, combine, cutoff):
    """Where the law is within exp(-cutoff) of its peak on `grid` (B, G).

    Returns lo, hi (two grid steps outside the support), the support's width in
    grid steps, and whether the support touches either end of the grid.
    """
    lp, _ = joint_log_density(readings, grid, law, combine)
    top = np.max(lp, axis=1, keepdims=True)
    keep = lp > top - cutoff
    first = np.argmax(keep, axis=1)
    last = keep.shape[1] - 1 - np.argmax(keep[:, ::-1], axis=1)
    touches = (~np.isfinite(top[:, 0])) | (first == 0) | (last == keep.shape[1] - 1)
    rows = np.arange(grid.shape[0])
    lo = grid[rows, np.maximum(first - 2, 0)]
    hi = grid[rows, np.minimum(last + 2, keep.shape[1] - 1)]
    return lo, hi, (last - first).astype(float), touches


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
    lo, hi, span_steps, unresolved = _support(readings, coarse, law, combine, num["cutoff"])
    # A safety net only: laws that cannot be normalised are refused at validation
    # (check_proper). Off-grid here means a proper law sits outside the search span.
    for _ in range(12):
        narrow = np.flatnonzero(~unresolved & (span_steps < num["min_support_steps"]))
        if narrow.size == 0:
            break
        grid = lo[narrow, None] + (hi - lo)[narrow, None] * np.linspace(0, 1, num["refine_points"])[None, :]
        lo_n, hi_n, steps_n, bad_n = _support(readings.subset(narrow), grid, law, combine, num["cutoff"])
        lo[narrow], hi[narrow], span_steps[narrow] = lo_n, hi_n, steps_n
        unresolved[narrow] |= bad_n

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
            out["log_sd"] = np.sqrt(np.sum(w * (u - mean_u[:, None]) ** 2, axis=1))
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
        self.declared_mismatch = spec.get("declared_mismatch")
        if self.declared_mismatch is not None and (
                not isinstance(self.declared_mismatch, str) or not self.declared_mismatch.strip()):
            raise ValueError(f"{where}: declared_mismatch must be a sentence saying why")
        if "direct" in spec:
            extra = set(spec) - {"name", "direct", "reduce", "declared_mismatch"}
            if extra:
                raise ValueError(f"{where}: a direct rule takes name, direct, reduce (got {sorted(extra)})")
            if spec["direct"] not in DIRECT_RULES:
                raise ValueError(f"{where}: direct rule must be one of {DIRECT_RULES}")
            if spec.get("reduce") != "pool_pair":
                raise ValueError(f"{where}: a direct rule acts on one pair, so it must state reduce: pool_pair")
            self.direct = spec["direct"]
            self.readouts = [self.direct]
            self.numerics = None
            return
        self.direct = None
        required = {"name", "reduce", "law", "combine", "readouts"}
        missing = required - set(spec)
        if missing:
            raise ValueError(f"{where}: missing required setting(s) {sorted(missing)}")
        extra = set(spec) - required - {"numerics", "declared_mismatch"}
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
        rule = combine.get("rule") if isinstance(combine, dict) else None
        if rule == "sequential":
            if set(combine) != {"rule", "old", "carry"}:
                raise ValueError(f"{where}: combine for sequential needs exactly old, carry and rule "
                                 "(the prior belongs to the old readings' law)")
            if combine["carry"] not in CARRIES:
                raise ValueError(f"{where}: carry must be one of {CARRIES}")
            old = combine["old"]
            if not isinstance(old, dict) or old.get("rule") not in ("single", "likelihood_product", "posterior_product"):
                raise ValueError(f"{where}: combine.old must be a non-sequential combine block")
            self._validate_combine(old, f"{where} (old readings)")
            if combine["carry"] == "vonmises_state" and (
                    old["rule"] != "likelihood_product" or law["nuisance"] != "radius"):
                raise ValueError(f"{where}: vonmises_state approximates the radius-nuisance likelihood "
                                 "product, so it needs old rule likelihood_product and nuisance radius")
        else:
            self._validate_combine(combine, where)
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
        try:
            validate_numerics(spec.get("numerics", {}))
        except ValueError as error:
            raise ValueError(f"{where}: {error}") from None
        self.numerics = dict(DEFAULT_NUMERICS, **spec.get("numerics", {}))

    @staticmethod
    def _validate_combine(combine, where):
        if not isinstance(combine, dict) or "rule" not in combine or "prior" not in combine:
            raise ValueError(f"{where}: combine needs rule and prior")
        rule = combine["rule"]
        if rule not in RULES or rule == "sequential":
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

    @property
    def pools(self):
        """Whether this estimator averages a series' readings into one pair."""
        return self.spec.get("reduce") == "pool_pair"

    def check_world(self, world):
        """Refuse settings that don't fit this world: improper laws, 'single' with
        several readings, and undeclared model/world mismatches."""
        where = f"estimator {self.name!r}"
        n = 1 if self.pools else world["readings"]
        if self.pools and world["design"] == "new_excitation" and world["readings"] > 1 \
                and self.declared_mismatch is None:
            raise ValueError(
                f"{where} pools readings, which assumes they share one latent pair, but the world's design is "
                f"new_excitation. If that mismatch is intended, add `declared_mismatch: <why>` to the estimator.")
        if self.direct is None:
            try:
                check_proper(self.spec["law"], self.spec["combine"], n)
            except ValueError as error:
                raise ValueError(f"{where}: {error}") from None

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

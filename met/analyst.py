"""The analyst: stages 2-5 of the pipeline.

An estimator is a declared model plus a readout. It sees only the readings and
the standard deviations it was supplied, never the world's settings.

    [reduce] -> per-reading law -> combine (rule + prior) -> readouts

or, for a direct rule, readings -> pooled pair -> one number.
"""

import math

import numpy as np

from .law import NUISANCES, REFERENCES, ReadingSet, nuisance_power

REDUCES = ("none", "pool_pair")
RULES = ("single", "likelihood_product", "posterior_product", "sequential", "hierarchical")
CARRIES = ("exact", "curve_only", "lognormal_fit", "tilt_fit", "vonmises_state")
COORDINATES = ("mass", "log_mass")
CALIBRATIONS = ("sandwich",)
POINT_READOUTS = ("ratio_of_means", "median", "geometric", "reciprocal_root")
DIRECT_RULES = ("norm_ratio", "dot_acceleration", "dot_force")
PRIOR_FLAT = ("flat_mass", "flat_log_mass", "flat_inverse_mass")
PRIOR_PARAMETRIC = ("uniform_angle", "sech_tilt", "lognormal", "angle_power")

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
        elif kind == "angle_power":
            # sin(theta)^p cos(theta)^q with tan(theta) = m / C. {sin: 1, cos: 1} is uniform_angle.
            z = u - _center(params["center"], readings)
            total = total - 0.5 * float(params["sin"]) * np.logaddexp(0.0, -2.0 * z) \
                          - 0.5 * float(params["cos"]) * np.logaddexp(0.0, 2.0 * z)
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
                    "lognormal": {"center", "width"}, "angle_power": {"sin", "cos", "center"}}.get(kind)
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
        if kind == "angle_power":
            _number(params["sin"], "angle_power sin")
            _number(params["cos"], "angle_power cos")
        if kind == "lognormal" and _number(params["width"], "lognormal width") <= 0:
            raise ValueError("lognormal width must be positive")


# Tail slopes of the joint log density in u = log m, as u -> +inf and u -> -inf.
# The radial kernel tends to a constant at both ends whatever the data, so
# whether a law can be normalised is decided by the settings, N and d alone.
# k is the nuisance power: 2 for the flat reference, d for the cartesian one
# (so the table depends on the dimension only under the cartesian reference).
def _like_slopes(nuisance, k):
    return {"radius": (0, 0), "acceleration": (-k, 0), "force": (0, k)}[nuisance]


def _ref_slopes(nuisance, k):
    return {"radius": (-1, 1), "acceleration": (k - 1, 1), "force": (-1, 1 - k)}[nuisance]


_PRIOR_SLOPES = {"flat_mass": (1, 1), "flat_log_mass": (0, 0), "flat_inverse_mass": (-1, -1),
                 "uniform_angle": (-1, 1), "sech_tilt": (0, 0)}


def tail_slopes(law, combine, n, d):
    """(slope as u -> +inf, slope as u -> -inf, confined) for n readings of dimension d after reduce.

    `confined` is True when a log-normal factor makes both tails decay regardless.
    The law is proper iff confined, or the first slope is < 0 and the second > 0.
    """
    k = nuisance_power(d, law["reference"])
    like = _like_slopes(law["nuisance"], k)
    ref = _ref_slopes(law["nuisance"], k)
    rule = combine["rule"]
    if rule == "sequential":
        carry = combine["carry"]
        if carry == "lognormal_fit":
            return like[0], like[1], True
        if carry == "tilt_fit":
            return like[0] - 1, like[1] + 1, False
        if carry == "vonmises_state":
            plus, minus, confined = tail_slopes(law, combine["old"], 0, d)
            return plus + like[0], minus + like[1], confined
        plus, minus, confined = tail_slopes(law, combine["old"], n - 1, d)
        return plus + like[0], minus + like[1], confined
    if rule == "hierarchical":
        plus, minus = 0, 0
    elif rule == "single":
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
        if kind == "angle_power":
            params = factor["angle_power"]
            plus -= float(params["cos"])
            minus += float(params["sin"])
            continue
        plus += _PRIOR_SLOPES[kind][0]
        minus += _PRIOR_SLOPES[kind][1]
    return plus, minus, confined


def check_proper(law, combine, n, d):
    """Raise if the declared law cannot be normalised for n readings of dimension d."""
    if combine["rule"] == "sequential":
        if n < 2:
            raise ValueError(f"rule 'sequential' needs at least two readings (old and new), but this world gives {n}")
        try:
            check_proper(law, combine["old"], n - 1, d)
        except ValueError as error:
            raise ValueError(f"the old readings' law: {error}") from None
    if combine["rule"] == "single" and n != 1:
        raise ValueError(f"rule 'single' needs exactly one reading after reduce, but this world gives {n}")
    plus, minus, confined = tail_slopes(law, combine, n, d)
    if not confined and not (plus < 0 and minus > 0):
        raise ValueError(
            f"this law cannot be normalised for N = {n} in {d}D: its log density has slope {plus} as m -> infinity "
            f"and {minus} as m -> 0 (in log m); it needs < 0 and > 0. Change the prior or the rule.")


# ---------------------------------------------------------------- joint law

def joint_log_density(readings, u, law, combine):
    """Joint log density in u (up to a constant) and A(u) = sum_i E[alpha_i | m]."""
    rule = combine["rule"]
    if rule == "sequential":
        return _sequential_log_density(readings, u, law, combine)
    if rule == "hierarchical":
        return _hierarchical_log_density(readings, u, law, combine)
    log_like, log_ref, cond_alpha = readings.reading_terms(u, law["nuisance"], reference=law["reference"])
    n = readings.readings
    if rule == "single":
        if n != 1:
            raise ValueError("rule 'single' needs exactly one reading per series after reduce")
        lp = log_like[:, 0] + log_ref[:, 0]
    elif rule == "likelihood_product":
        lp = np.sum(log_like, axis=1)
        if combine.get("calibrate") == "sandwich":
            lp = sandwich_weight(readings, law, combine)[:, None] * lp
    elif rule == "posterior_product":
        lp = np.sum(log_like + log_ref, axis=1)
        if combine["coordinate"] == "mass":
            lp = lp + (1 - n) * u
    else:
        raise ValueError(f"unknown rule {rule!r}")
    lp = lp + prior_log_density(combine["prior"], u, readings)
    return lp, np.sum(cond_alpha, axis=1)


# ---------------------------------------------------------------- sandwich calibration (our25)
#
# With many readings sharing only the mass, the declared law's curvature H
# understates the spread of its own peak: the peak's variance is J / H^2, the
# law's is 1 / H, where J is the variance of the summed per-reading score. For
# the cartesian law, per reading, J = d + r*^2 and H = r*^2 exactly (r* the true
# pair's noise-unit length). The calibrated law raises the likelihood to the
# power w = H / J, so its width matches the peak's spread; the prior is not
# tempered. J and H are estimated from the readings at the likelihood's peak:
#
#     H = -sum_i l_i''(u*),   J = N/(N-1) sum_i (l_i'(u*) - mean l')^2,
#
# (centred scores, so a grid-sized miss of the peak does not bias J). The ratio
# does not depend on the coordinate. w is clipped to [W_FLOOR, 1]: never narrower
# than the declared law, and if the readings show no curvature (H <= 0: the
# mass is not identified) the law falls back to the prior alone.

W_FLOOR = 1e-6


def _likelihood_total(readings, u, law):
    log_like, _, _ = readings.reading_terms(u, law["nuisance"], reference=law["reference"])
    return log_like


def sandwich_weight(readings, law, combine):
    """Per-series tempering power w = H/J (cached on the readings)."""
    key = ("sandwich", id(combine))
    if key in readings.series_data:
        return readings.series_data[key]
    b = readings.series
    center = np.log(readings.s[:, 0])
    coarse = center[:, None] + np.arange(-30.0, 30.0001, 0.05)[None, :]
    total = np.sum(_likelihood_total(readings, coarse, law), axis=1)
    peak = coarse[np.arange(b), np.argmax(total, axis=1)]
    for width in (0.1, 0.002):
        local = peak[:, None] + np.linspace(-width, width, 201)[None, :]
        total = np.sum(_likelihood_total(readings, local, law), axis=1)
        k = np.clip(np.argmax(total, axis=1), 1, 199)
        rows = np.arange(b)
        y0, y1, y2 = total[rows, k - 1], total[rows, k], total[rows, k + 1]
        step = local[0, 1] - local[0, 0]
        curv = y0 - 2 * y1 + y2
        shift = np.where(curv < 0, 0.5 * (y0 - y2) / np.where(curv < 0, curv, -1.0), 0.0)
        peak = local[rows, k] + np.clip(shift, -1, 1) * step
    h = 1e-3
    stencil = peak[:, None] + np.array([-h, 0.0, h])[None, :]
    ll = _likelihood_total(readings, stencil, law)                    # (B, N, 3)
    grad = (ll[:, :, 2] - ll[:, :, 0]) / (2 * h)
    hess = (ll[:, :, 2] - 2 * ll[:, :, 1] + ll[:, :, 0]) / h**2
    n = readings.readings
    H = -np.sum(hess, axis=1)
    J = n / (n - 1) * np.sum((grad - grad.mean(axis=1, keepdims=True)) ** 2, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where((H > 0) & (J > 0), H / J, W_FLOOR)
    w = np.clip(w, W_FLOOR, 1.0)
    readings.series_data[key] = w
    return w


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
    like_new, _, cond_new = new.reading_terms(u, law["nuisance"], reference=law["reference"])
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


# ---------------------------------------------------------------- the hierarchical rule
#
# Each reading's standardised latent vector v_i (f/sigma_F = r sin theta, alpha/sigma_a =
# r cos theta, v = r u) gets a shared prior N(0, omega^2 I) instead of a flat one.
# Integrating every v_i out exactly leaves, with n = N d and beta = omega^2/(1+omega^2),
#
#     L(theta, beta) = (1 - beta)^(n/2) exp(beta H(theta) / 2),   H = sum_i h_i(theta)^2.
#
# The flat (cartesian) excitation prior is the beta -> 1 limit. The excitation scale is
# learned: beta ~ Beta(1, b) is integrated out exactly,
#
#     int_0^1 (1-beta)^(n/2+b-1) e^(beta z) d(beta) = e^z int_0^1 g^(c-1) e^(-g z) dg,
#     c = n/2 + b,  z = H/2,
#
# a lower incomplete gamma function. Given theta, v_i ~ N(beta w_i, beta I), so
# E[alpha_i | theta, beta] = sigma_a cos(theta) sqrt(beta) chi_d(sqrt(beta) h_i); the beta
# average uses 24 quantile midpoints of beta's conditional law (1 - beta is a Gamma(c, z)
# truncated to (0, 1)).

_HIER_NODES = (np.arange(24) + 0.5) / 24


def _log_lower_gamma_integral(c, z):
    """log int_0^1 g^(c-1) e^(-g z) dg for c > 0, z >= 0 (vectorised, stable)."""
    from scipy.special import gammainc, gammaln
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    big = z > c
    if np.any(big):
        zb = z[big]
        out[big] = gammaln(c) - c * np.log(zb) + np.log(gammainc(c, zb))
    small = ~big
    if np.any(small):
        zs = z[small]
        # int_0^1 g^(c-1) e^(-gz) dg = e^(-z) sum_k z^k / (c (c+1) ... (c+k)), a positive series
        term = np.full_like(zs, 1.0 / c)
        total = term.copy()
        for k in range(1, 2000):
            term = term * zs / (c + k)
            total += term
            if np.all(term < 1e-17 * total):
                break
        out[small] = -zs + np.log(total)
    return out


def _hierarchical_log_density(readings, u, law, combine):
    from scipy.special import gammainc, gammaincinv
    b = float(combine["excitation"]["beta_b"])
    like, _, _ = readings.reading_terms(u, "radius", reference="cartesian")
    h2 = 2.0 * like                                  # (B, N, G): h_i(theta)^2
    z = 0.5 * np.sum(h2, axis=1)                     # (B, G): H/2
    n = readings.readings * readings.dimension
    c = 0.5 * n + b
    fixed = combine.get("_fixed_beta")
    if fixed is not None:
        # declared oracle: the true excitation scale is given, not learned
        lp = 0.5 * n * math.log1p(-fixed) + fixed * z + prior_log_density(combine["prior"], u, readings)
    else:
        lp = z + _log_lower_gamma_integral(c, z) + prior_log_density(combine["prior"], u, readings)
    # A(theta): average over beta | theta, on a sub-grid, interpolated in u (it is smooth)
    g = u.shape[1]
    idx = np.unique(np.linspace(0, g - 1, min(g, 201)).round().astype(int))
    zs = z[:, idx]
    if fixed is not None:
        beta = np.full(zs.shape + (1,), float(fixed))
    else:
        frac = gammainc(c, np.maximum(zs, 1e-300))[..., None] * _HIER_NODES
        with np.errstate(divide="ignore", invalid="ignore"):
            gam = np.where(frac > 1e-280, gammaincinv(c, frac) / np.maximum(zs, 1e-300)[..., None],
                           _HIER_NODES ** (1.0 / c))
        beta = np.clip(1.0 - gam, 0.0, 1.0)          # (B, S, 24)
    from .law import cartesian_mean_radius
    h = np.sqrt(np.maximum(h2[:, :, idx], 0.0))      # (B, N, S)
    sq = np.sqrt(beta)[:, None, :, :]                # (B, 1, S, 24)
    radius = np.mean(sq * cartesian_mean_radius(sq * h[..., None], readings.dimension), axis=3)
    zrel = u[:, idx][:, None, :] - np.log(readings.s)[:, :, None]
    cos = np.exp(-0.5 * np.logaddexp(0.0, 2.0 * zrel))
    acc_sub = np.sum(readings.acceleration_sd[:, :, None] * cos * radius, axis=1)   # (B, S)
    if len(idx) == g:
        acc = acc_sub
    else:
        acc = np.stack([np.interp(u[i], u[i, idx], acc_sub[i]) for i in range(u.shape[0])])
    return lp, acc


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


# ---------------------------------------------------------------- declared oracles
#
# An oracle is the one place the analyst is shown part of the truth. It is only
# ever declared by name in a preset (`oracle: ...`), so it cannot happen by accident.
#
#   known_excitation  told every reading's true acceleration vector a*_i. Only the
#                     force readings are then uncertain: F_i = m a*_i + noise, so with a
#                     flat prior on m > 0 the law is a normal (mean sum F.a*/sigma^2 over
#                     sum |a*|^2/sigma^2) truncated to m > 0. A ceiling: no method that
#                     sees only the readings can know more.
#   known_beta        the hierarchical rule with the true excitation scale plugged in
#                     instead of learned (beta from the world's settings). In the Gaussian-
#                     excitation world it is the best any method can do with the same
#                     readings and the same mass prior.
ORACLES = ("known_excitation", "known_beta")


def world_beta(world):
    """The true beta = omega^2 / (1 + omega^2) of a world, from its settings.

    omega^2 is the per-component variance of the standardised latent vector:
    R^2 (1 + (m/s)^2) / d, with R the RMS acceleration SNR of the pushes.
    """
    spec = world["acceleration_snr"]
    if isinstance(spec, dict) and "normal_rms" in spec:
        r2 = float(spec["normal_rms"]) ** 2
    elif isinstance(spec, dict):
        lo, hi = spec["uniform"]
        r2 = (lo * lo + lo * hi + hi * hi) / 3.0
    else:
        r2 = float(spec) ** 2
    s = world["noise"]["force_sd"] / world["noise"]["acceleration_sd"]
    omega2 = r2 * (1.0 + (float(world["mass"]) / s) ** 2) / world["dimension"]
    return omega2 / (1.0 + omega2)


def known_excitation_summaries(readings, true_acceleration, readouts):
    from scipy.stats import truncnorm
    w = 1.0 / readings.force_sd ** 2                                    # (B, N)
    prec = np.sum(w * np.sum(true_acceleration ** 2, axis=2), axis=1)   # (B,)
    num = np.sum(w * np.sum(readings.force * true_acceleration, axis=2), axis=1)
    ok = prec > 0
    mu = np.where(ok, num / np.where(ok, prec, 1.0), 0.0)
    sd = np.where(ok, 1.0 / np.sqrt(np.where(ok, prec, 1.0)), 1.0)
    lo = -mu / sd
    out = {"unresolved": ~ok}
    for name in readouts:
        if name == "ratio_of_means":      # alpha is known exactly, so E[f]/E[alpha] = E[m]
            out[name] = truncnorm.mean(lo, np.inf, loc=mu, scale=sd)
        elif name == "median":
            out[name] = truncnorm.ppf(0.5, lo, np.inf, loc=mu, scale=sd)
        else:
            content = interval_content(name)
            out[name] = np.stack([truncnorm.ppf((1 - content) / 2, lo, np.inf, loc=mu, scale=sd),
                                  truncnorm.ppf((1 + content) / 2, lo, np.inf, loc=mu, scale=sd)], axis=1)
    for name, value in out.items():
        if name != "unresolved":
            value[~ok] = np.nan
    return out


# ---------------------------------------------------------------- estimators

class Estimator:
    """A validated estimator specification from a preset."""

    def __init__(self, spec):
        if not isinstance(spec, dict) or "name" not in spec:
            raise ValueError("each estimator needs a name")
        self.name = str(spec["name"])
        self.spec = spec
        where = f"estimator {self.name!r}"
        self.oracle = spec.get("oracle")
        if self.oracle is not None and self.oracle not in ORACLES:
            raise ValueError(f"{where}: oracle must be one of {ORACLES}")
        if self.oracle == "known_excitation":
            extra = set(spec) - {"name", "oracle", "readouts"}
            if extra or "readouts" not in spec:
                raise ValueError(f"{where}: a known_excitation oracle takes exactly name, oracle, readouts")
            for r in spec["readouts"]:
                if r not in ("ratio_of_means", "median") and not str(r).startswith("interval_"):
                    raise ValueError(f"{where}: known_excitation offers ratio_of_means, median and intervals")
                if str(r).startswith("interval_"):
                    interval_content(r)
            self.direct, self.readouts, self.numerics, self.declared_mismatch = None, list(spec["readouts"]), None, None
            return
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
        extra = set(spec) - required - {"numerics", "declared_mismatch", "oracle"}
        if extra:
            raise ValueError(f"{where}: unknown setting(s) {sorted(extra)}")
        if spec["reduce"] not in REDUCES:
            raise ValueError(f"{where}: reduce must be one of {REDUCES}")
        law = spec["law"]
        if not isinstance(law, dict) or set(law) != {"reference", "nuisance"}:
            raise ValueError(f"{where}: law needs exactly reference and nuisance")
        if law["reference"] not in REFERENCES:
            raise ValueError(f"{where}: reference must be one of {REFERENCES}")
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
            if "calibrate" in old:
                raise ValueError(f"{where}: calibrate is not supported inside a sequential rule")
            if combine["carry"] == "vonmises_state" and (
                    old["rule"] != "likelihood_product" or law["nuisance"] != "radius"):
                raise ValueError(f"{where}: vonmises_state approximates the radius-nuisance likelihood "
                                 "product, so it needs old rule likelihood_product and nuisance radius")
        else:
            self._validate_combine(combine, where)
            if self.oracle == "known_beta" and rule != "hierarchical":
                raise ValueError(f"{where}: the known_beta oracle is the hierarchical rule with beta given")
            if rule == "hierarchical" and (law["reference"] != "cartesian" or law["nuisance"] != "radius"):
                raise ValueError(f"{where}: the hierarchical rule is built on the cartesian reference "
                                 "(it is its proper version), so it needs reference cartesian and nuisance radius")
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
        if rule == "hierarchical":
            if set(combine) != {"rule", "prior", "excitation"}:
                raise ValueError(f"{where}: combine for hierarchical needs exactly rule, prior and excitation")
            exc = combine["excitation"]
            if not isinstance(exc, dict) or set(exc) != {"beta_b"} or isinstance(exc["beta_b"], bool) \
                    or not isinstance(exc["beta_b"], (int, float)) or not exc["beta_b"] > 0:
                raise ValueError(f"{where}: excitation must be {{beta_b: b > 0}} (prior Beta(1, b) on "
                                 "beta = omega^2/(1+omega^2))")
            try:
                validate_prior(combine["prior"], "likelihood_product")
            except ValueError as error:
                raise ValueError(f"{where}: {error}") from None
            return
        allowed = {"rule", "prior"} | ({"coordinate"} if rule == "posterior_product" else set())
        optional = {"calibrate"} if rule == "likelihood_product" else set()
        if not (allowed <= set(combine) <= allowed | optional):
            raise ValueError(f"{where}: combine for {rule} needs exactly {sorted(allowed)}"
                             + (f" (optional: {sorted(optional)})" if optional else ""))
        if "calibrate" in combine and combine["calibrate"] not in CALIBRATIONS:
            raise ValueError(f"{where}: calibrate must be one of {CALIBRATIONS}")
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
        if self.oracle == "known_excitation":
            return
        n = 1 if self.pools else world["readings"]
        if self.pools and world["design"] == "new_excitation" and world["readings"] > 1 \
                and self.declared_mismatch is None:
            raise ValueError(
                f"{where} pools readings, which assumes they share one latent pair, but the world's design is "
                f"new_excitation. If that mismatch is intended, add `declared_mismatch: <why>` to the estimator.")
        if self.direct is None:
            try:
                check_proper(self.spec["law"], self.spec["combine"], n, world["dimension"])
            except ValueError as error:
                raise ValueError(f"{where}: {error}") from None
            if self.spec["combine"].get("calibrate"):
                if n < 2:
                    raise ValueError(f"{where}: calibrate estimates a spread across readings, so it needs "
                                     f"at least two readings after reduce (this world gives {n})")
                if self.spec["law"]["nuisance"] != "radius":
                    raise ValueError(f"{where}: calibrate needs nuisance radius (the per-reading "
                                     "likelihood must carry no mass factor)")

    def evaluate(self, readings, truth=None, world=None):
        if self.oracle == "known_excitation":
            return known_excitation_summaries(readings, truth["true_acceleration"], self.readouts)
        if self.oracle == "known_beta":
            combine = dict(self.spec["combine"], _fixed_beta=world_beta(world))
            return posterior_summaries(readings, self.spec["law"], combine, self.readouts, self.numerics)
        if self.direct is not None:
            return {self.direct: direct_rule(readings, self.direct),
                    "unresolved": np.zeros(readings.series, dtype=bool)}
        if self.spec["reduce"] == "pool_pair":
            readings = readings.pooled()
        return posterior_summaries(readings, self.spec["law"], self.spec["combine"],
                                   self.readouts, self.numerics)

    def point_readouts(self):
        return [r for r in self.readouts if r in POINT_READOUTS or r in DIRECT_RULES]

    def check_world_oracle(self, world):
        if self.oracle == "known_beta":
            world_beta(world)

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

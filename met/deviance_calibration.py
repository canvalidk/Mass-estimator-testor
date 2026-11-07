"""Research-only deviance operators and deterministic references for T40.

No estimator uses these functions. Family, scale and scale unit are explicit,
unforced choices, not defaults or an adopted correction. Oscillation-derived
scales respect exact Z and continuity, but their L19 admissibility is unknown.

Both maps have slope eta at zero, slope one at infinity, and a multiplicative
eta semigroup at fixed scale. Fixed scale does NOT give the full eta-powered
Gaussian local limit. Growing scale can recover that reference; it is a
sandwich reference when the cap leaves C/J unchanged, including C/J = 1.

Tail preservation repairs T39's nu=10 pair when eta stays away from zero and
the raw law concentrates. It does not repair raw nu<1 or arbitrary eta->0
paths. These are compatibility diagnostics, not coverage rankings.
eta=0 is a boundary extension only; tail-preservation statements require eta>0.
"""
import math

import numpy as np
from scipy import integrate, optimize, special

from .cut_consistency import opposed_pair_eta


def _nonnegative(values, name):
    values = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return values


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _eta(value):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("eta must be finite and lie in [0, 1]")
    return value


def deviance_loss(deviance, eta, *, family, scale):
    """Transform D=Lhat-L >= 0, retaining array shape (or returning a scalar).

    logistic: a*log(1 + eta*expm1(D/a)).
    lambert:  a*W(eta*(D/a)*exp(D/a)), on the nonnegative real branch.

    The log domain and a tail-offset solve avoid exponent overflow. A zero
    scale is accepted only when D is identically zero: every positive-scale
    operator already has the same forced zero value there.
    """
    d, eta = _nonnegative(deviance, "deviance"), _eta(eta)
    if family not in ("logistic", "lambert"):
        raise ValueError("family must be 'logistic' or 'lambert'")
    a = float(scale)
    if not math.isfinite(a) or a < 0 or (a == 0 and np.any(d > 0)):
        raise ValueError("scale must be positive and finite unless deviance is identically zero")
    result = np.zeros_like(d)
    if eta == 1:
        result = d.copy()
    elif eta and a:
        nonzero = d > 0
        values = d[nonzero]
        with np.errstate(over="ignore", under="ignore"):
            x = values / a
        log_eta = math.log(eta)
        transformed = np.empty_like(x)
        tiny = x < 1e-8
        if family == "logistic":
            small = (x <= math.log(2)) & ~tiny
            transformed[small] = a * np.log1p(eta * np.expm1(x[small]))
            other = ~small & ~tiny
            xo, vo = x[other], values[other]
            shift = xo + log_eta
            dominant = shift >= 0
            rest = np.empty_like(xo)
            # Write the large branch in physical D units to retain its offset.
            with np.errstate(under="ignore"):
                rest[dominant] = (vo[dominant] + a * log_eta
                                  + a * np.log1p(np.exp(math.log1p(-eta) - shift[dominant])))
                log_expm1 = xo[~dominant] + np.log1p(-np.exp(-xo[~dominant]))
                z = log_eta + log_expm1
                rest[~dominant] = a * np.logaddexp(0, z)
            transformed[other] = rest
        else:
            # For large x solve k-log(1-k/x)=-log(eta), k=(D-F)/a.
            # Computing a*W(exp(log_eta+log(x)+x)) would lose the offset first.
            tail = x > max(50.0, -4 * log_eta)
            xt = x[tail]
            target = -log_eta
            k = np.full_like(xt, target) / (1 + 1 / xt)
            for _ in range(6):
                residual = k - np.log1p(-k / xt) - target
                k -= residual / (1 + 1 / (xt - k))
            transformed[tail] = values[tail] - a * k
            middle = ~tail & ~tiny
            xm = x[middle]
            omega = special.wrightomega(log_eta + np.log(xm) + xm)
            transformed[middle] = a * omega
        # Evaluate the near-zero series in physical D units: eta*x or W's
        # argument may underflow before multiplication by a huge finite scale.
        coefficient = 0.5 if family == "logistic" else 1.0
        cubic = ((1-eta)*(1-2*eta)/6 if family == "logistic"
                 else (1-eta)*(1-3*eta)/2)
        transformed[tiny] = eta * values[tiny] * (1 + coefficient*(1-eta)*x[tiny] + cubic*x[tiny]**2)
        result[nonzero] = transformed
    return float(result) if result.ndim == 0 else result


def effective_scale(oscillations, *, unit, family):
    """Sum continuous bounded weights of each factor's nonnegative oscillation.

    x = oscillation/unit. rational uses x/(1+x), exponential uses 1-exp(-x).
    Each weight is in [0,1], is zero exactly for a flat factor, and vanishes
    continuously there. The unit and weight family are explicit, unforced.
    """
    values, unit = _nonnegative(oscillations, "oscillations"), _positive(unit, "unit")
    if family == "rational":
        small = values <= unit
        weights = np.empty_like(values)
        ratio = values[small] / unit
        weights[small] = ratio / (1 + ratio)
        weights[~small] = 1 / (1 + unit / values[~small])
    elif family == "exponential":
        with np.errstate(over="ignore"):
            weights = -np.expm1(-values / unit)
    else:
        raise ValueError("family must be 'rational' or 'exponential'")
    return float(np.sum(weights))


def local_gaussian_reference(eta, *, family, scale, probability):
    """Continuous local limit in z with eta-powered Gaussian variance one.

    Its density is proportional to exp[-F(z^2/(2*eta))]. Return variance and
    the central probability interval's positive endpoint. These deterministic
    integrals do not evaluate repeated-world estimator coverage.
    """
    eta = _eta(eta)
    if eta == 0:
        raise ValueError("the standardized local reference requires eta > 0")
    _positive(scale, "scale")
    probability = float(probability)
    if not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    deviance_loss(0.0, eta, family=family, scale=scale)

    def kernel(z):
        return math.exp(-deviance_loss(z*z / (2 * eta), eta, family=family, scale=scale))

    normalizer = 2 * integrate.quad(kernel, 0, np.inf, epsabs=2e-12, epsrel=2e-12)[0]
    variance = 2 * integrate.quad(lambda z: z*z * kernel(z), 0, np.inf,
                                  epsabs=2e-12, epsrel=2e-12)[0] / normalizer

    def centered_probability(stop):
        return 2 * integrate.quad(kernel, 0, stop, epsabs=2e-12, epsrel=2e-12)[0] / normalizer

    upper = 1.0
    while centered_probability(upper) < probability:
        upper *= 2
    half_width = optimize.brentq(lambda z: centered_probability(z) - probability, 0, upper,
                                 xtol=2e-12)
    return dict(variance=variance, half_width=half_width, normalizer=normalizer,
                probability=probability)


def pair_angle_cdf(angle, *, truth_angle, push, nu_eff, error_size, family, scale):
    """Independent transformed-law CDF for T39's strict-maximum opposed pair.

    c=s=1, so the angular prior is uniform. Deviance is the explicit polynomial
    denominator ratio, eta is T39's closed form, and no estimator grid or fitted
    values enter. Return (CDF, quadrature error estimate), not a certified bound.
    """
    theta0, radius = float(truth_angle), _positive(push, "push")
    nu, error = _positive(nu_eff, "nu_eff"), _positive(error_size, "error_size")
    if not math.isfinite(theta0) or not 0 < theta0 < math.pi / 2:
        raise ValueError("truth_angle must lie strictly between zero and pi/2")
    eta = opposed_pair_eta(push=radius, nu_eff=nu, error_size=error)
    deviance_loss(0.0, eta, family=family, scale=scale)
    _positive(scale, "scale")
    theta = float(angle)
    if not math.isfinite(theta):
        raise ValueError("angle must be finite")
    if theta <= 0:
        return 0.0, 0.0
    if theta >= math.pi / 2:
        return 1.0, 0.0
    error2, radius2 = error*error, radius*radius
    baseline = nu + error2
    linear = 2 * (radius2 * (nu-error2) - error2*baseline) / baseline**2
    quadratic_root = (radius2 + error2) / baseline

    def kernel(t):
        z = math.sin(t-theta0)**2
        d = 0.5 * nu * math.log1p(linear*z + (quadratic_root*z)**2)
        return math.exp(-deviance_loss(d, eta, family=family, scale=scale))

    width = math.sqrt(nu) / radius
    points = [theta0 + factor*width for factor in (-100, -10, -1, 0, 1, 10, 100)]

    def integral(stop):
        interior = sorted(x for x in points if 0 < x < stop)
        return integrate.quad(kernel, 0, stop, points=interior, limit=200,
                              epsabs=2e-12, epsrel=2e-11)

    numerator, ne = integral(theta)
    denominator, de = integral(math.pi / 2)
    cdf = numerator / denominator
    return cdf, (ne + cdf*de) / denominator


def pair_angle_quantile(probability, *, truth_angle, push, nu_eff, error_size, family, scale):
    """Quantile of pair_angle_cdf; mass is tan(angle) when c=s=1."""
    probability = float(probability)
    if not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    arguments = dict(truth_angle=truth_angle, push=push, nu_eff=nu_eff,
                     error_size=error_size, family=family, scale=scale)
    pair_angle_cdf(truth_angle, **arguments)
    return optimize.brentq(lambda theta: pair_angle_cdf(theta, **arguments)[0] - probability,
                           0, math.pi / 2, xtol=2e-12)

"""Research-only normalization references for T44, not an estimator adapter.

The single-reading reference has positive observed angle, unit channel noise
and uniform angular prior c=s=1. Variance laws/weights are supplied and frozen.
Normalizing a factor average reweights q(V) by its component evidence Z(V);
averaging normalized conditional laws keeps q(V) unchanged.

component_temperature is INSIDE each Gaussian factor before normalization.
At temperature 1 the factor average is the raw ours31/ours32 operation. At
other temperatures it is NOT ours32's outer power (E_q L_V)^eta. Nor does a
product of single-reading cut marginals define the genuine whole-series cut:
that construction needs one conditional posterior normalized across the series.
No physical signal-size population or production equation is adopted here.
"""
import math

import numpy as np
from scipy import integrate, optimize, special


HALF_PI = math.pi/2


def _finite(value, name):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _controls(epsabs, epsrel, limit):
    absolute, relative = _finite(epsabs, "epsabs"), _finite(epsrel, "epsrel")
    if absolute <= 0 or relative <= 0 or isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("positive quadrature tolerances and integer limit are required")
    return dict(epsabs=absolute, epsrel=relative, limit=limit)


def _reading(observed_angle, radius, component_temperature):
    center = _finite(observed_angle, "observed_angle")
    radius = _finite(radius, "radius")
    temperature = _finite(component_temperature, "component_temperature")
    if not 0 < center < HALF_PI or radius < 0 or not 0 <= temperature <= 1:
        raise ValueError("positive interior observed angle, nonnegative radius and temperature in [0,1] required")
    return center, radius, temperature


def _component(angle, center, precision, controls):
    theta = _finite(angle, "angle")
    stop = min(HALF_PI, max(0., theta))
    if precision == 0:
        return dict(cdf=stop/HALF_PI, evidence=1., cdf_error=0., evidence_error=0.)
    if math.isinf(precision):
        cdf = 0. if stop < center else (1. if stop > center else .5)
        return dict(cdf=cdf, evidence=0., cdf_error=0., evidence_error=0.)

    def kernel(t):
        return math.exp(-.5*precision*math.sin(t-center)**2)

    width = 1/math.sqrt(precision)
    points = [center+x*width for x in (-10, -1, 0, 1, 10)]

    def integral(upper):
        interior = sorted(set(x for x in points if 0 < x < upper))
        return integrate.quad(kernel, 0, upper, points=interior, **controls)

    total, te = integral(HALF_PI)
    if stop == 0:
        numerator, ne = 0., 0.
    elif stop == HALF_PI:
        numerator, ne = total, 0.
    else:
        numerator, ne = integral(stop)
    cdf = numerator/total
    return dict(cdf=cdf, evidence=total/HALF_PI,
                cdf_error=(ne+cdf*te)/total, evidence_error=te/HALF_PI)


def gaussian_component(angle, *, observed_angle, radius, variance,
                       component_temperature, epsabs, epsrel, limit):
    """CDF and evidence Z(V)=int pi exp[-temperature*A/(2V)].

    A(theta)=radius^2*sin^2(theta-observed_angle). Error fields are quadrature
    estimates, not rigorous bounds. All numerical controls are explicit.
    """
    center, r, temperature = _reading(observed_angle, radius, component_temperature)
    variance = _finite(variance, "variance")
    if variance <= 0:
        raise ValueError("variance must be positive")
    return _component(angle, center, temperature*r*r/variance, _controls(epsabs, epsrel, limit))


def discrete_mixture_cdf(angle, *, observed_angle, radius, variances, weights,
                         kind, component_temperature, epsabs, epsrel, limit):
    """Frozen proper-q mixture; kind='factor_average' or 'normalized_cut'.

    Return the CDF, component evidences and effective variance weights. The
    factor average's q* is proportional to q*Z, while the normalized cut's
    variance marginal is the supplied q. This is a single-reading identity.
    """
    vs, q = np.asarray(variances, dtype=float), np.asarray(weights, dtype=float)
    if (vs.ndim != 1 or not vs.size or q.shape != vs.shape or np.any(~np.isfinite(vs))
            or np.any(vs <= 0) or np.any(~np.isfinite(q)) or np.any(q < 0)
            or not math.isclose(float(q.sum()), 1., rel_tol=1e-12, abs_tol=1e-12)):
        raise ValueError("variances and weights must specify a proper finite positive-variance law")
    if kind not in ("factor_average", "normalized_cut"):
        raise ValueError("kind must be 'factor_average' or 'normalized_cut'")
    q = q/q.sum()
    components = [gaussian_component(angle, observed_angle=observed_angle, radius=radius,
                  variance=v, component_temperature=component_temperature,
                  epsabs=epsabs, epsrel=epsrel, limit=limit) for v in vs]
    evidence = np.array([x["evidence"] for x in components])
    effective = q if kind == "normalized_cut" else q*evidence/float(q@evidence)
    cdf = float(effective @ np.array([x["cdf"] for x in components]))
    error = float(effective @ np.array([x["cdf_error"] for x in components]))
    return dict(cdf=cdf, error=error, component_evidence=evidence,
                effective_weights=effective.copy())


def discrete_mixture_quantile(probability, *, observed_angle, radius, variances,
                              weights, kind, component_temperature, epsabs, epsrel, limit):
    """Angular quantile of discrete_mixture_cdf; mass is tan(angle)."""
    probability = _finite(probability, "probability")
    if not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    arguments = dict(observed_angle=observed_angle, radius=radius, variances=variances,
                     weights=weights, kind=kind, component_temperature=component_temperature,
                     epsabs=epsabs, epsrel=epsrel, limit=limit)
    discrete_mixture_cdf(observed_angle, **arguments)
    return optimize.brentq(lambda t: discrete_mixture_cdf(t, **arguments)["cdf"]-probability,
                           0, HALF_PI, xtol=2e-12)


def ig_normalized_cut_cdf(angle, *, observed_angle, radius, nu_eff,
                          component_temperature, epsabs, epsrel, limit):
    """N1 normalized cut with q(V)=IG(nu/2,nu/2), any nu>0.

    Precision's Gamma quantile maps its probability integral to [0,1], avoiding
    the endpoint singularity at nu<2 and resolving large-V tails. Only a single
    positive-angle reading is modeled, so Q+=0 and q is fixed. Temperature
    stays inside each normalized component; no outer powered law is computed.
    """
    center, r, temperature = _reading(observed_angle, radius, component_temperature)
    nu = _finite(nu_eff, "nu_eff")
    if nu <= 0:
        raise ValueError("nu_eff must be positive")
    theta = _finite(angle, "angle")
    controls = _controls(epsabs, epsrel, limit)
    if theta <= 0:
        return 0., 0.
    if theta >= HALF_PI:
        return 1., 0.
    if r == 0 or temperature == 0:
        return theta/HALF_PI, 0.
    alpha, inner_errors = nu/2, []

    def integrand(probability):
        gamma = special.gammaincinv(alpha, probability)
        component = _component(theta, center, temperature*r*r*gamma/alpha, controls)
        inner_errors.append(component["cdf_error"])
        return component["cdf"]

    cdf, error = integrate.quad(integrand, 0, 1, **controls)
    return cdf, error+max(inner_errors, default=0.)


def positive_q_bias(nu_eff, *, epsabs, epsrel, limit):
    """Exact null crossmoment of w(Q+)XY with w=nu/(nu+Q+).

    B_nu=(1/pi)*int exp(-t)*t/(nu+t)dt, equal to
    [1-nu*exp(nu)*E1(nu)]/pi. The integral avoids large-nu cancellation and
    exponential overflow. This is a proof diagnostic, not a fitted weight.
    """
    nu = _finite(nu_eff, "nu_eff")
    if nu <= 0:
        raise ValueError("nu_eff must be positive")
    value, error = integrate.quad(lambda t: math.exp(-t)*t/(nu+t), 0, np.inf,
                                   **_controls(epsabs, epsrel, limit))
    return value/math.pi, error/math.pi

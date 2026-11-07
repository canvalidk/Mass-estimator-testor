"""Independent deterministic references for T39; these do not fit or replace a law.

The population gap is a Laplace-transform integral for a noncentral chi-square
residual. The single-reading references integrate in the physical noise angle,
where the half-Cauchy with c = sigma_F/sigma_a is uniform. No estimator grid,
fitted noise, latent-size distribution, or Monte Carlo draws enter these checks.
"""
import math

import numpy as np
from scipy import integrate, optimize, special


HALF_PI = math.pi / 2


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def population_loss_gap(k, nu_eff, displacement_squared, *, noise_scale):
    """E[-ell(m) + ell(m0)] for K = v*chi2_k(delta/v).

    delta is the squared distance of the true mean to the candidate plane;
    noise_scale is v, not its SD. Return (gap, quadrature_error_estimate).
    Integrating this conditional gap over a declared scale mixture preserves
    positivity. The error estimate is numerical, not a rigorous bound.
    """
    if k not in (1, 2, 3):
        raise ValueError("k must be 1, 2 or 3")
    nu = _positive(nu_eff, "nu_eff")
    scale = _positive(noise_scale, "noise_scale")
    delta = float(displacement_squared)
    if not math.isfinite(delta) or delta < 0:
        raise ValueError("displacement_squared must be nonnegative and finite")
    if delta == 0:
        return 0.0, 0.0

    def integrand(t):
        if t == 0:
            return delta / nu
        denominator = 1 + 2 * scale * t / nu
        central = math.exp(-t - 0.5 * k * math.log(denominator))
        return central * -math.expm1(-delta * t / (nu * denominator)) / t

    value, error = integrate.quad(integrand, 0, np.inf, epsabs=2e-12, epsrel=2e-11)
    coefficient = 0.5 * (nu + k - 1)
    return coefficient * value, coefficient * error


def _angle_inputs(truth_angle, nu_eff):
    theta0 = float(truth_angle)
    if not math.isfinite(theta0) or not 0 < theta0 < HALF_PI:
        raise ValueError("truth_angle must lie strictly between 0 and pi/2")
    return theta0, _positive(nu_eff, "nu_eff")


def aligned_angle_cdf(angle, *, truth_angle, push, nu_eff):
    """One aligned 1D reading, c = s; independent finite-push CDF in theta.

    The whitened data are push*(sin(truth_angle), cos(truth_angle)). eta is
    mathematically 1 for one reading. Adaptive quadrature resolves the peak
    separately from its tails. Return (CDF, propagated error estimate).
    """
    theta0, nu = _angle_inputs(truth_angle, nu_eff)
    radius = _positive(push, "push")
    theta = float(angle)
    if not math.isfinite(theta):
        raise ValueError("angle must be finite")
    if theta <= 0:
        return 0.0, 0.0
    if theta >= HALF_PI:
        return 1.0, 0.0

    def density(x):
        residual = radius * math.sin(x - theta0) / math.sqrt(nu)
        return math.exp(-0.5 * nu * math.log1p(residual * residual))

    width = math.sqrt(nu) / radius
    points = [theta0 + multiple * width for multiple in (-10, -1, 0, 1, 10)]

    def integral(stop):
        interior = sorted(set(x for x in points if 0 < x < stop))
        return integrate.quad(density, 0, stop, points=interior, limit=200,
                              epsabs=2e-12, epsrel=2e-11)

    numerator, numerator_error = integral(theta)
    total, total_error = integral(HALF_PI)
    cdf = numerator / total
    return cdf, (numerator_error + cdf * total_error) / total


def aligned_angle_quantile(probability, *, truth_angle, push, nu_eff):
    """Quantile in theta of aligned_angle_cdf, not a mass-grid readout."""
    p = float(probability)
    if not math.isfinite(p) or not 0 < p < 1:
        raise ValueError("probability must lie strictly between 0 and 1")
    _angle_inputs(truth_angle, nu_eff)
    _positive(push, "push")
    return optimize.brentq(
        lambda theta: aligned_angle_cdf(theta, truth_angle=truth_angle,
                                       push=push, nu_eff=nu_eff)[0] - p,
        0, HALF_PI, xtol=2e-13)


def _sine_power_integral(length, nu):
    # Integral_0^length sin(x)^(-nu) dx; length <= pi/2 and 0 < nu < 1.
    shape = 0.5 * (1 - nu)
    return 0.5 * special.beta(shape, 0.5) * special.betainc(
        shape, 0.5, math.sin(length) ** 2)


def heavy_tail_limit_cdf(angle, *, truth_angle, tail_index):
    """T39's nonconcentrating sine-power limit, 0 < tail_index < 1.

    The single-reading index is nu_eff; the tempered two-reading index is
    2*nu_eff*eta_limit. Density is proportional to |sin(theta-theta0)|**(-index).
    Incomplete beta functions integrate the singularity without sampling it.
    """
    theta0, nu = _angle_inputs(truth_angle, tail_index)
    if nu >= 1:
        raise ValueError("the nonconcentrating limit requires tail_index < 1")
    theta = float(angle)
    if not math.isfinite(theta):
        raise ValueError("angle must be finite")
    if theta <= 0:
        return 0.0
    if theta >= HALF_PI:
        return 1.0
    left = _sine_power_integral(theta0, nu)
    total = left + _sine_power_integral(HALF_PI - theta0, nu)
    if theta <= theta0:
        return (left - _sine_power_integral(theta0 - theta, nu)) / total
    return (left + _sine_power_integral(theta - theta0, nu)) / total


def heavy_tail_limit_quantile(probability, *, truth_angle, tail_index):
    """Quantile in theta of the nondegenerate limit; mass is tan(theta) at s=1."""
    p = float(probability)
    if not math.isfinite(p) or not 0 < p < 1:
        raise ValueError("probability must lie strictly between 0 and 1")
    # Validate the family even if a root algorithm would not visit its interior.
    heavy_tail_limit_cdf(truth_angle, truth_angle=truth_angle, tail_index=tail_index)
    return optimize.brentq(
        lambda theta: heavy_tail_limit_cdf(theta, truth_angle=truth_angle,
                                           tail_index=tail_index) - p,
        0, HALF_PI, xtol=2e-13)


def opposed_pair_eta(*, push, nu_eff, error_size):
    """Exact eta for y_plus/minus = push*u0 +/- error_size*v0 (1D each).

    Only the regime in which the central line is the strict whole-line
    maximum is admitted. The condition is explicit; no numerical fit enters.
    """
    radius = _positive(push, "push")
    nu = _positive(nu_eff, "nu_eff")
    error = _positive(error_size, "error_size")
    error2 = error * error
    numerator = radius**2 * (nu - error2) - error2 * (nu + error2)
    if numerator <= 0:
        raise ValueError("the central line must be the strict whole-line maximum")
    return min(1.0, numerator / (nu * radius**2 * error2))


def opposed_pair_angle_cdf(angle, *, truth_angle, push, nu_eff, error_size, calibrate):
    """Independent continuous CDF of the T39 two-reading counterexample.

    c = s, so the angular prior is uniform. Both the calibrated and raw law
    use exactly the same two readings. Return (CDF, quadrature error estimate).
    """
    theta0, nu = _angle_inputs(truth_angle, nu_eff)
    radius = _positive(push, "push")
    error = _positive(error_size, "error_size")
    eta = opposed_pair_eta(push=radius, nu_eff=nu, error_size=error) if calibrate else 1.0
    theta = float(angle)
    if not math.isfinite(theta):
        raise ValueError("angle must be finite")
    if theta <= 0:
        return 0.0, 0.0
    if theta >= HALF_PI:
        return 1.0, 0.0

    def density(x):
        sine, cosine = math.sin(x - theta0), math.cos(x - theta0)
        residuals = (radius * sine - error * cosine, radius * sine + error * cosine)
        return math.exp(-0.5 * nu * eta * sum(math.log1p(r*r / nu) for r in residuals))

    width, displacement = math.sqrt(nu) / radius, math.atan(error / radius)
    points = [center + multiple * width
              for center in (theta0 - displacement, theta0, theta0 + displacement)
              for multiple in (-10, -1, 0, 1, 10)]

    def integral(stop):
        interior = sorted(set(x for x in points if 0 < x < stop))
        return integrate.quad(density, 0, stop, points=interior, limit=200,
                              epsabs=2e-12, epsrel=2e-11)

    numerator, numerator_error = integral(theta)
    total, total_error = integral(HALF_PI)
    cdf = numerator / total
    return cdf, (numerator_error + cdf * total_error) / total


def opposed_pair_angle_quantile(probability, *, truth_angle, push, nu_eff, error_size, calibrate):
    """Angular quantile of opposed_pair_angle_cdf, independent of native grids."""
    p = float(probability)
    if not math.isfinite(p) or not 0 < p < 1:
        raise ValueError("probability must lie strictly between 0 and 1")
    _angle_inputs(truth_angle, nu_eff)
    _positive(push, "push")
    _positive(error_size, "error_size")
    if calibrate:
        opposed_pair_eta(push=push, nu_eff=nu_eff, error_size=error_size)
    return optimize.brentq(
        lambda theta: opposed_pair_angle_cdf(theta, truth_angle=truth_angle,
                                            push=push, nu_eff=nu_eff, error_size=error_size,
                                            calibrate=calibrate)[0] - p,
        0, HALF_PI, xtol=2e-13)

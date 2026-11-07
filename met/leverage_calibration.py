"""Independent research references for flat-Gaussian leverage calibration.

These are exact score identities, limiting coverage integrals, and an unadopted
spectral calibration control. Observations are standardized scalar (F, a)
pairs with identity Gaussian covariance and a common noise angle. No helper
is an estimator adapter or changes the production law. Count and energy are
proof inputs to the score identities only; the spectral control accepts data
alone.
"""
import math

import numpy as np
from scipy import integrate, special, stats

from met.gaussian_information import flat_gaussian_fit


def _count(value):
    if type(value) is not int or value < 1:
        raise ValueError("count must be a positive integer")
    return value


def _energy(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError("energy must be finite and nonnegative")
    return value


def _probability(value):
    value = float(value)
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError("probability must lie strictly between zero and one")
    return value, math.sqrt(2) * float(special.erfinv(value))


def score_log_cf(count, energy, t):
    """Exact log characteristic function of sum_i (r_i + xi_i) W_i.

    xi_i and W_i are independent standard Gaussians; energy=sum_i r_i**2.
    The answer depends on the count and total energy, not their allocation.
    t may be a finite scalar or array. The logarithm is real since this CF
    is strictly positive. logaddexp avoids forming t**2 at large t.
    """
    count, energy = _count(count), _energy(energy)
    t = np.asarray(t, dtype=float)
    if np.any(~np.isfinite(t)):
        raise ValueError("t must be finite")
    with np.errstate(divide="ignore"):
        log_one_plus_square = np.logaddexp(0., 2 * np.log(np.abs(t)))
    fraction = (t / np.hypot(1., t))**2
    result = -0.5 * count * log_one_plus_square - 0.5 * energy * fraction
    return float(result) if result.ndim == 0 else result


def score_moments(count, energy):
    """Exact first, second and fourth moments of the same aggregate score."""
    count, energy = _count(count), _energy(energy)
    variance = count + energy
    return dict(mean=0., variance=variance,
                fourth=3 * variance**2 + 6 * count + 12 * energy)


def one_dominant_coverage(tau, probability):
    """Limiting central-interval coverage with one dominant reading.

    tau is the limiting ratio N/Delta, not an estimator input. Conditional
    on G~N(0,1), the standardized fitted error is W+sqrt(tau)*G and the
    quoted radius is z*sqrt(max(1, tau*(1+G**2))), with W independent of G.
    Return (coverage, quadrature error estimate). The integral splits at
    the cap boundary. Standard-normal tails beyond 12 contribute at most
    2*sf(12), which is included in the error estimate. It is not a rigorous
    bound on quadrature or floating-point error.
    """
    tau = float(tau)
    if not math.isfinite(tau) or tau < 0:
        raise ValueError("tau must be finite and nonnegative")
    probability, z = _probability(probability)
    if tau == 0:
        return probability, 0.
    root_tau = math.sqrt(tau)

    def integrand(g):
        shift = root_tau * g
        radius = z * max(1., root_tau * math.hypot(1., g))
        conditional = stats.norm.cdf(radius-shift) - stats.norm.cdf(-radius-shift)
        return 2 * stats.norm.pdf(g) * conditional

    points = [0., 12.]
    if tau < 1:
        boundary = math.sqrt(1-tau) / root_tau
        if boundary < 12:
            points.append(boundary)
    points.sort()
    value, error = 0., 2 * stats.norm.sf(12.)
    for low, high in zip(points[:-1], points[1:]):
        piece, estimate = integrate.quad(integrand, low, high,
                                         epsabs=2e-12, epsrel=2e-11)
        value += piece
        error += estimate
    return float(value), float(error)


def equal_dominant_coverage(count, probability):
    """Strong-signal limiting coverage for equal dominant signal shares.

    For count=M, V=chi2_(M-1)/M is independent of the standard-normal
    fitted error. The cap uses max(1,V). The one-dimensional tail integral
    adds the nonnegative excess coverage to the requested nominal level.
    No simulated observations or production law enter this reference.
    """
    count = _count(count)
    probability, z = _probability(probability)
    if count == 1:
        return probability, 0.

    def excess(v):
        return z * stats.norm.pdf(z*math.sqrt(v)) / math.sqrt(v) * stats.chi2.sf(count*v, count-1)

    value, error = integrate.quad(excess, 1., np.inf,
                                  epsabs=2e-12, epsrel=2e-11)
    return float(probability+value), float(error)


def spectral_calibration(observations):
    """Unadopted data-only control eta_sp=H/lambda_max and K_sp=H**2/lambda_max.

    The existing flat-Gaussian reference supplies the eigensystem. Tiny
    negative eigenvalues from roundoff are clipped to the PSD domain. A
    zero gap gives a flat angular law (prior_only=True); at zero scatter
    eta_sp=1 is a harmless fallback because the precision is zero. Neither
    true signal energy nor a theoretical observation count is accepted.
    """
    fitted = flat_gaussian_fit(observations)
    eigenvalues = np.maximum(fitted["eigenvalues"], 0.)
    largest = float(eigenvalues[1])
    gap = float(largest-eigenvalues[0])
    eta = gap/largest if largest > 0 else 1.
    return dict(angle=fitted["angle"], curv=gap, eta_sp=eta,
                precision_sp=eta*gap, eigenvalues=eigenvalues,
                scatter=fitted["scatter"], prior_only=gap == 0.)

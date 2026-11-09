"""Research controls for E10's null-washout/consistency compatibility question.

None of these helpers is a production estimator or an adopted correction.
The logarithmic gate has an unforced schedule. The ordered scatter changes
under permutations and exact-zero padding; it is deliberately a countercontrol.
"""
import math

import numpy as np
from scipy.special import logsumexp

from met.gaussian_information import flat_gaussian_fit


def mixture_log_second_moment(radii):
    """Log E0[L^2] for the proof-only independent fair-sign mixture.

    Each whitened Gaussian mean is +/-r_i*u for one unit direction u.
    This is a likelihood comparison, not a declared signal population.
    """
    r = np.asarray(radii, dtype=float)
    if r.ndim != 1 or not np.all(np.isfinite(r)):
        raise ValueError("radii must be a finite vector")
    squared = r*r
    if not np.all(np.isfinite(squared)):
        raise ValueError("squared radii exceed floating-point range")
    # log(cosh(x)), keeping the small-x cancellation out of the calculation.
    small = squared < 20
    terms = np.empty_like(squared)
    terms[small] = np.log1p(2*np.sinh(squared[small]/2)**2)
    z = squared[~small]
    terms[~small] = z-math.log(2)+np.log1p(np.exp(-2*z))
    return float(terms.sum())


def _log_alpha(value):
    value = float(value)
    if math.isnan(value) or value > 0:
        raise ValueError("log null probability must lie in [-infinity, 0]")
    return value


def fixed_signal_log_transfer_bound(log_null_probability, radii):
    """Log capped Cauchy--Schwarz bound on a proof-mixture event probability."""
    log_alpha = _log_alpha(log_null_probability)
    return min(0., .5*(log_alpha+mixture_log_second_moment(radii)))


def homogeneous_log_transfer_bound(log_null_probability, *, count, energy):
    """Log bound .5*log(alpha)+Delta^2/(4*N), capped at zero.

    Requires N equal signal radii with total squared radius Delta. It is
    not an upper bound for arbitrary allocations of the same total energy.
    """
    log_alpha, energy = _log_alpha(log_null_probability), float(energy)
    if count < 1 or int(count) != count or not math.isfinite(energy) or energy < 0:
        raise ValueError("positive integer count and finite nonnegative energy required")
    return min(0., .5*log_alpha+.25*(energy/math.sqrt(count))**2)


def symmetric_mixture_log_product_bound(squared_radii):
    """Log product-of-row-means bound for a sign-and-permutation mixture.

    The exact second moment is the permutation average of products
    cosh(sqrt(a_i*a_perm(i))). Bounding cosh(z)<=exp(z²/2), then negative
    association of sampling without replacement, gives this upper bound.
    """
    a = np.asarray(squared_radii, dtype=float)
    if a.ndim != 1 or len(a) == 0 or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError("squared radii must be a finite nonempty nonnegative vector")
    return float(np.sum(logsumexp(.5*np.outer(a, a), axis=1)-math.log(len(a))))


def gated_precision(precision, energy):
    """K*=K^2/[K+log(e+E)], a continuous unforced comparison schedule."""
    k, e = float(precision), float(energy)
    if not math.isfinite(k) or not math.isfinite(e) or min(k, e) < 0:
        raise ValueError("precision and energy must be finite and nonnegative")
    return k*(k/(k+math.log(math.e+e)))


def washout_fit(observations):
    """T42 fit plus the research gate; inherited fit fields stay unchanged."""
    fit = flat_gaussian_fit(observations)
    energy = float(np.trace(fit['scatter']))
    return dict(fit, energy=energy,
                washout_precision=gated_precision(fit['precision'], energy))


def ordered_scatter(observations):
    """Centered scatter with row weight 1/sqrt(i), divided by log(N+1).

    The returned matrix defines the angular tilt u_theta' matrix u_theta.
    Row position and count are intentional, physically unlicensed features.
    """
    y = np.asarray(observations, dtype=float)
    if y.ndim != 2 or y.shape[1] != 2 or len(y) == 0 or not np.all(np.isfinite(y)):
        raise ValueError("observations must be a nonempty finite N by 2 array")
    weights = 1/np.sqrt(np.arange(1, len(y)+1, dtype=float))
    centered = (y.T*weights)@y-weights.sum()*np.eye(2)
    return dict(scatter=centered, normalizer=math.log(len(y)+1),
                matrix=centered/math.log(len(y)+1))

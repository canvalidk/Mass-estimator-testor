"""Research-only Gaussian testing certificates for T41.

The auxiliary independent signs average allowed deterministic signal worlds
for a testing impossibility argument. They are NOT a declared physical signal
population, a nuisance prior supplied to an estimator, or a mass estimator.
All TV quantities returned below evaluate analytic upper bounds, not exact TV
or interval-certified floating-point enclosures.
The binary error bounds use equal testing weights on two hypotheses only.

Total signal energy alone can diverge while these testing errors stay bounded
away from zero. The separate scatter identity illustrates a stronger sufficient
information condition; no estimator or oracle adapter is implemented here.
"""
import math

import numpy as np


def _nonnegative(value, name):
    result = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def _radii(values):
    result = _nonnegative(values, "radii")
    if result.ndim != 1:
        raise ValueError("radii must be a one-dimensional sequence")
    return result


def _finite(value, name):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def log_cosh(value):
    """Stable log(cosh(x)); extended-real infinities are allowed, NaN is not."""
    x = np.asarray(value, dtype=float)
    if np.any(np.isnan(x)):
        raise ValueError("log_cosh does not accept NaN")
    absolute = np.abs(x)
    result = np.empty_like(absolute)
    small = absolute <= 20
    result[small] = np.log1p(2*np.sinh(absolute[small]/2)**2)
    with np.errstate(under="ignore"):
        result[~small] = absolute[~small] - math.log(2) + np.log1p(np.exp(-absolute[~small])**2)
    return float(result) if result.ndim == 0 else result


def sign_mixture_log_ratio(points, *, radius, angle):
    """Log density ratio to N(0,I2) for means +/-radius*u_angle, equal signs."""
    r, theta = float(_nonnegative(radius, "radius")), _finite(angle, "angle")
    y = np.asarray(points, dtype=float)
    if y.ndim == 0 or y.shape[-1] != 2 or np.any(~np.isfinite(y)):
        raise ValueError("points must be finite vectors with final dimension two")
    direction = np.array([math.cos(theta), math.sin(theta)])
    result = log_cosh(r*(y @ direction)) - 0.5*r*r
    return float(result) if np.ndim(result) == 0 else result


def log_cross_moment(radius, *, angle_difference):
    """Exact log E0[L_u L_v] = log cosh(radius^2*cos(angle_difference))."""
    r = float(_nonnegative(radius, "radius"))
    difference = _finite(angle_difference, "angle_difference")
    with np.errstate(over="ignore"):
        squared = np.float64(r)*r
    return float(log_cosh(squared*math.cos(difference)))


def _log_expm1(nonnegative):
    if nonnegative == 0:
        return -math.inf
    if nonnegative <= math.log(2):
        return math.log(math.expm1(nonnegative))
    return nonnegative + math.log1p(-math.exp(-nonnegative))


def _chi_square(log_second):
    if log_second > math.log(np.finfo(float).max):
        return math.inf
    return math.expm1(log_second)


def _tv_from_log(log_bound):
    return 1.0 if log_bound >= 0 else math.exp(log_bound)


def product_testing_certificate(radii, *, angle_difference):
    """Finite independent-reading testing bounds against null and each other.

    S=sum log cosh(r_i^2), C=sum log cosh(r_i^2*cos(delta)). Then chi2_null
    =exp(S)-1. The triangle bound on two-angle TV is sqrt(chi2_null);
    direct Cauchy gives sqrt((exp(S)-exp(C))/2), usually sharper.
    """
    radii = _radii(radii)
    difference = _finite(angle_difference, "angle_difference")
    with np.errstate(over="ignore"):
        squared = radii*radii
    same = float(np.sum(log_cosh(squared)))
    cross = float(np.sum(log_cosh(squared*math.cos(difference))))
    chi_square = _chi_square(same)
    log_chi = _log_expm1(same)
    each_null = _tv_from_log(0.5*log_chi - math.log(2))
    triangle = _tv_from_log(0.5*log_chi)
    if difference == 0:
        direct = 0.0
    elif math.isinf(same):
        direct = 1.0  # Moment overflow leaves only the trivial upper bound.
    else:
        gap = max(0.0, same-cross)
        if gap <= 8*np.finfo(float).eps*max(1., same):
            # Do not turn an unresolved subtraction into a false claim of
            # identical distributions. E0[L_u L_v]>=1 still gives this bound.
            direct = _tv_from_log(0.5*(log_chi-math.log(2)))
        else:
            direct = _tv_from_log(0.5*(cross + _log_expm1(gap) - math.log(2)))
    tv = min(triangle, direct)
    return dict(total_energy=float(np.sum(squared)), log_second_moment=same,
                log_cross_moment=cross, chi_square_null=chi_square,
                tv_each_to_null_upper=each_null, tv_triangle_upper=triangle,
                tv_direct_upper=direct, tv_upper=tv, binary_error_lower=0.5*(1-tv))


def harmonic_radii(*, scale, count):
    """The explicit proof schedule r_i=scale/sqrt(i), i=1,...,count."""
    b = _finite(scale, "scale")
    if b <= 0:
        raise ValueError("scale must be positive")
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError("count must be a positive integer")
    return b/np.sqrt(np.arange(1, count+1, dtype=float))


def harmonic_testing_certificate(*, scale):
    """Uniform-in-N bounds for r_i=scale/sqrt(i), despite divergent energy.

    log cosh(x)<=x^2/2 and sum i^-2=pi^2/6 give S_N<=A=b^4*pi^2/12.
    The triangle TV bound is sqrt(expm1(A)); the angle-free direct-moment
    bound uses E0[L_u L_v]>=1 and improves this by sqrt(2).
    """
    b = _finite(scale, "scale")
    if b <= 0:
        raise ValueError("scale must be positive")
    with np.errstate(over="ignore"):
        upper = float(np.float64(b)**4 * (math.pi**2/12))
    chi_square = _chi_square(upper)
    log_chi = _log_expm1(upper)
    triangle = _tv_from_log(0.5*log_chi)
    direct = _tv_from_log(0.5*(log_chi-math.log(2)))
    return dict(log_second_moment_upper=upper, chi_square_null_upper=chi_square,
                tv_triangle_upper=triangle, binary_error_triangle_lower=0.5*(1-triangle),
                tv_direct_upper=direct, binary_error_direct_lower=0.5*(1-direct))


def harmonic_energy_bounds(*, scale, count):
    """b^2 log(N+1)<=sum r_i^2<=b^2(1+log N); the lower bound diverges."""
    b = _finite(scale, "scale")
    # Validate the schedule without allocating it.
    if b <= 0 or isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError("scale must be positive and count a positive integer")
    return dict(lower=b*b*math.log1p(count), upper=b*b*(1+math.log(count)))


def scatter_fluctuation_second_moment(radii):
    """Exact E||sum YY^T-E(sum YY^T)||_F^2=6(N+sum r_i^2).

    Independent 2D Gaussians have means of the stated lengths. Independent
    sign mixing gives the same scatter law. The mean eigengap for a common
    direction is Delta=sum r_i^2; Delta^2/(N+Delta)->infinity controls relative
    scatter fluctuation. This function implements no direction estimator.
    """
    radii = _radii(radii)
    with np.errstate(over="ignore"):
        energy = float(np.sum(radii*radii))
    return 6*(len(radii)+energy)

"""B21 research likelihoods on the common-instrument slice (T69).

These functions do not change the estimator. ``d`` is the number of scalar
rows N*k, not the number of nonzero observations. The integrated comparison
uses the unforced Lebesgue measure dDelta on the one total signal energy.
"""
import math

import numpy as np
from scipy.special import erf, gammainc, gammaln, hyp0f1, hyp1f1, ive


def _inputs(h, d):
    if not np.isfinite(d) or d < 1 or int(d) != d:
        raise ValueError("d must be a positive integer number of scalar rows")
    h = np.asarray(h, dtype=float)
    if np.any(~np.isfinite(h)) or np.any(h < 0):
        raise ValueError("h must be finite and nonnegative")
    return h, int(d)


def _bessel_ratio_over_x(x, d):
    """I_(d/2)(x) / [x I_(d/2-1)(x)], including its value 1/d at zero."""
    x = np.asarray(x, dtype=float)
    result = np.full_like(x, 1.0/d)
    small = (x > 0) & (x < 100)
    z = x[small]**2/4
    result[small] = hyp0f1(d/2+1, z) / (d*hyp0f1(d/2, z))
    large = x >= 100
    xx = x[large]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = ive(d/2, xx) / (xx*ive(d/2-1, xx))
    bad = ~np.isfinite(ratio)
    if np.any(bad):
        z = xx[bad]**2/4
        ratio[bad] = hyp0f1(d/2+1, z) / (d*hyp0f1(d/2, z))
    result[large] = ratio
    if np.any(~np.isfinite(result)):
        raise FloatingPointError("Bessel ratio outside the supported numerical range")
    return result


def _log_0f1(d, x):
    """log 0F1(d/2; x**2/4), with exponential Bessel scaling."""
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)
    small = (x > 0) & (x < 100)
    result[small] = np.log(hyp0f1(d/2, x[small]**2/4))
    large = x >= 100
    xx = x[large]
    with np.errstate(divide="ignore", invalid="ignore"):
        values = gammaln(d/2) + (1-d/2)*np.log(xx/2) + np.log(ive(d/2-1, xx)) + xx
    bad = ~np.isfinite(values)
    if np.any(bad):
        values[bad] = np.log(hyp0f1(d/2, xx[bad]**2/4))
    result[large] = values
    if np.any(~np.isfinite(result)):
        raise FloatingPointError("hypergeometric kernel outside the supported numerical range")
    return result


def invariant_log_kernel(h, d, delta):
    """log K_d(h, Delta), omitting the central-Wishart factor independent of mass."""
    h, d = _inputs(h, d)
    delta = np.asarray(delta, dtype=float)
    if np.any(~np.isfinite(delta)) or np.any(delta < 0):
        raise ValueError("delta must be finite and nonnegative")
    h, delta = np.broadcast_arrays(h, delta)
    return -.5*delta + _log_0f1(d, np.sqrt(delta)*np.sqrt(h))


def profile_delta(h, d):
    """Maximizing total energy Delta >= 0; exactly zero for h <= d.

    For h > d the Bessel ratio divided by x decreases strictly from 1/d
    to zero. Bisection in x**2/h avoids the stationary root x=0.
    """
    h, d = _inputs(h, d)
    result = np.zeros_like(h)
    active = h > d
    hh = h[active]
    if hh.size:
        lo, hi = np.zeros_like(hh), hh.copy()
        for _ in range(48):
            mid = (lo+hi)/2
            derivative = _bessel_ratio_over_x(np.sqrt(hh)*np.sqrt(mid), d)-1/hh
            positive = derivative > 0
            lo = np.where(positive, mid, lo)
            hi = np.where(positive, hi, mid)
        result[active] = (lo+hi)/2
    return result


def profile_log_kernel(h, d):
    """The profile log likelihood, zero on h <= d, with Delta eliminated."""
    h, d = _inputs(h, d)
    return invariant_log_kernel(h, d, profile_delta(h, d))


def profile_curvature_ratio(h, d):
    """Profile/flat angular curvature ratio at an interior stationary maximum of h."""
    h, d = _inputs(h, d)
    return np.divide(profile_delta(h, d), h, out=np.zeros_like(h), where=h > 0)


def integrated_log_kernel(h, d):
    """log integral K_d(h, Delta) dDelta, less the constant log 2.

    The closed form is log 1F1(1; d/2; h/2). Flat dDelta is an unforced
    improper nuisance measure, not a declared per-reading size population.
    """
    h, d = _inputs(h, d)
    z = h/2
    if d == 2:
        return z.copy()
    result = np.zeros_like(z)
    positive = z > 0
    zz = z[positive]
    if d == 1:
        result[positive] = np.logaddexp(0., zz+.5*np.log(math.pi*zz)+np.log(erf(np.sqrt(zz))))
    else:
        direct = zz < 500
        values = np.empty_like(zz)
        values[direct] = np.log(hyp1f1(1, d/2, zz[direct]))
        big = zz[~direct]
        with np.errstate(divide="ignore"):
            values[~direct] = gammaln(d/2)+big+(1-d/2)*np.log(big)+np.log(gammainc(d/2-1, big))
        result[positive] = values
    if np.any(~np.isfinite(result)):
        raise FloatingPointError("integrated kernel outside the supported numerical range")
    return result


def gram_energy(gram, theta):
    """h(theta) for v=(sin(theta),cos(theta)) and a 2 by 2 Gram matrix."""
    gram, theta = np.asarray(gram, dtype=float), np.asarray(theta, dtype=float)
    if gram.shape != (2, 2) or not np.all(np.isfinite(gram)) or not np.allclose(gram, gram.T):
        raise ValueError("gram must be a finite symmetric 2 by 2 matrix")
    if np.linalg.eigvalsh(gram).min() < -1e-10*max(1., np.max(np.abs(gram))):
        raise ValueError("gram must be positive semidefinite")
    s, c = np.sin(theta), np.cos(theta)
    return np.maximum(0., gram[0, 0]*s*s+2*gram[0, 1]*s*c+gram[1, 1]*c*c)


def angular_law(gram, d, kind="profile", nodes=4001):
    """Research grid on theta in [0, pi/2], with a half-Cauchy prior at c=s.

    Fixed grids must be refined for concentrated laws. Endpoints describe
    limits only; no zero or infinite mass is returned as an admissible point.
    """
    if nodes < 3 or int(nodes) != nodes:
        raise ValueError("nodes must be an integer at least three")
    theta = np.linspace(0., math.pi/2, int(nodes))
    h = gram_energy(gram, theta)
    kernels = {"profile": profile_log_kernel, "integrated": integrated_log_kernel,
               "flat": lambda h, d: h/2}
    if kind not in kernels:
        raise ValueError("kind must be profile, integrated or flat")
    log_kernel = kernels[kind](h, d)
    density = np.exp(log_kernel-log_kernel.max())
    density /= np.trapezoid(density, theta)
    return theta, density, log_kernel

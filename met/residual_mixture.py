"""Research references for T48's normalized Gaussian residual mixture.

The specified flat latent measure is arclength in the declared noise metric,
independent of V. Fully normalized residual Gaussians mixed with the same
V~IG(nu/2,nu/2) give exponent (nu+residual_dimension)/2. This measure choice
and joint integration differ from the fitted-Q cut; no production change,
new interpretation of nu, or adoption is made. A Student density in residual
coordinates is not itself a normalized probability density of mass.
"""
import math

import numpy as np
from scipy import integrate, optimize, special


HALF_PI = math.pi/2


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def student_residual_log_density(residual_squared, *, dimension, nu_eff):
    """Proper whitened Student residual density evaluated at ||residual||^2."""
    q = np.asarray(residual_squared, dtype=float)
    if np.any(~np.isfinite(q)) or np.any(q < 0):
        raise ValueError("residual_squared must be finite and nonnegative")
    if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    nu = _positive(nu_eff, "nu_eff")
    coefficient = special.gammaln((nu+dimension)/2)-special.gammaln(nu/2)-dimension*math.log(math.pi*nu)/2
    result = coefficient-(nu+dimension)*np.log1p(q/nu)/2
    return float(result) if result.ndim == 0 else result


def student_residual_density(residual_squared, *, dimension, nu_eff):
    """Density corresponding to student_residual_log_density, not a mass law."""
    result = np.exp(student_residual_log_density(residual_squared, dimension=dimension, nu_eff=nu_eff))
    return float(result) if np.ndim(result) == 0 else result


def _controls(epsabs, epsrel, limit):
    absolute, relative = _positive(epsabs, "epsabs"), _positive(epsrel, "epsrel")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")
    return dict(epsabs=absolute, epsrel=relative, limit=limit)


def _cdf(angle, center, kernel, width, controls):
    theta, center = float(angle), float(center)
    if not math.isfinite(theta) or not math.isfinite(center) or not 0 < center < HALF_PI:
        raise ValueError("angle must be finite and the center strictly between zero and pi/2")
    if theta <= 0:
        return 0., 0.
    if theta >= HALF_PI:
        return 1., 0.
    points = [center+x*width for x in (-100, -10, -1, 0, 1, 10, 100)]

    def integral(stop):
        interior = sorted(set(x for x in points if 0 < x < stop))
        return integrate.quad(kernel, 0, stop, points=interior, **controls)

    numerator, ne = integral(theta)
    total, te = integral(HALF_PI)
    probability = numerator/total
    return probability, (ne+probability*te)/total


def aligned_angle_cdf(angle, *, observed_angle, radius, nu_eff, epsabs, epsrel, limit):
    """Corrected raw N1 angular law, uniform prior c=s=1, positive fitted line.

    Kernel [1+R^2*sin^2(theta-center)/nu]^(-(nu+1)/2). The observations are
    held fixed; this is continuous quadrature, not an estimator's grid.
    """
    r, nu = _positive(radius, "radius"), _positive(nu_eff, "nu_eff")

    def kernel(theta):
        q = (r*math.sin(theta-observed_angle))**2
        return math.exp(-.5*(nu+1)*math.log1p(q/nu))

    return _cdf(angle, observed_angle, kernel, math.sqrt(nu)/r, _controls(epsabs, epsrel, limit))


def aligned_angle_quantile(probability, *, observed_angle, radius, nu_eff, epsabs, epsrel, limit):
    """Angular quantile of aligned_angle_cdf; physical mass is tan(angle)."""
    probability = float(probability)
    if not math.isfinite(probability) or not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    arguments = dict(observed_angle=observed_angle, radius=radius, nu_eff=nu_eff,
                     epsabs=epsabs, epsrel=epsrel, limit=limit)
    aligned_angle_cdf(observed_angle, **arguments)
    return optimize.brentq(lambda theta: aligned_angle_cdf(theta, **arguments)[0]-probability,
                           0, HALF_PI, xtol=2e-12)


def student_pair_parameters(*, radius, nu_eff, error_size):
    """Exact global-fit Bartlett quantities for R*u0 +/- d*v0, rank one each.

    The strict central-maximum regime is required. Prefactor a=(nu+1)/2 is
    fixed by this reference model; it is not supplied as a free correction.
    """
    r, nu = _positive(radius, "radius"), _positive(nu_eff, "nu_eff")
    d = float(error_size)
    if not math.isfinite(d) or d < 0:
        raise ValueError("error_size must be finite and nonnegative")
    r2, d2, a = r*r, d*d, (nu+1)/2
    numerator = r2*(nu-d2)-d2*(nu+d2)
    if numerator <= 0:
        raise ValueError("the central line must be the strict global projective maximum")
    curv = 4*a*numerator/(nu+d2)**2
    j = 8*a*a*r2*d2/(nu+d2)**2
    eta = 1. if j == 0 else min(1., curv/j)
    return dict(coefficient=a, curv=curv, J=j, eta=eta,
                raw_tail_index=4*a, quoted_tail_index=4*a*eta)


def pair_angle_cdf(angle, *, truth_angle, radius, nu_eff, error_size, calibrate,
                    epsabs, epsrel, limit):
    """Exact-family raw or globally powered corrected-Student angular CDF.

    This investigator reference uses the proved global central fit of the
    symmetric pair, not hidden truth passed to a production estimator.
    """
    parameters = student_pair_parameters(radius=radius, nu_eff=nu_eff, error_size=error_size)
    if not isinstance(calibrate, bool):
        raise ValueError("calibrate must be explicit True or False")
    eta = parameters["eta"] if calibrate else 1.
    a = parameters["coefficient"]
    r2, d2 = radius*radius, error_size*error_size
    baseline = nu_eff+d2
    linear = 2*(r2*(nu_eff-d2)-d2*baseline)/baseline**2
    square_root = (r2+d2)/baseline

    def kernel(theta):
        z = math.sin(theta-truth_angle)**2
        return math.exp(-a*eta*math.log1p(linear*z+(square_root*z)**2))

    return _cdf(angle, truth_angle, kernel, math.sqrt(nu_eff)/radius, _controls(epsabs, epsrel, limit))

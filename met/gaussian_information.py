"""Independent research references for T42's flat-Gaussian information limit.

Scope: independent scalar readings (F,a), unit Gaussian channel noise, common
noise ratio s=1, and the uniform angular prior c=s=1. The eigh fit and the
continuous integrals are deterministic diagnostics, not a new estimator or an
adoption of the flat law. They do not apply the finite-nu ours32 factors.
Known signal radii occur only in moment identities, never in the fitted law.
"""
import math

import numpy as np
from scipy import integrate, optimize


HALF_PI = math.pi/2


def flat_gaussian_fit(observations):
    """Independent 2x2 eigensystem fit to rows (force, acceleration).

    Return the whole-projective-line angle, curvature, raw squared-score sum,
    existing Bartlett cap/fallback, and precision eta*curvature. At an exact
    eigenvalue tie angle=0 is bookkeeping; the angular law is flat regardless.
    No true mass, latent radii or investigator truth enters this calculation.
    """
    y = np.asarray(observations, dtype=float)
    if y.ndim != 2 or y.shape[1] != 2 or np.any(~np.isfinite(y)):
        raise ValueError("observations must be a finite N by 2 array")
    scatter = y.T @ y
    if np.any(~np.isfinite(scatter)):
        raise ValueError("scatter exceeds the floating-point range")
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    gap = float(max(0., eigenvalues[1]-eigenvalues[0]))
    principal = eigenvectors[:, 1]
    angle = ((math.atan2(principal[0], principal[1])+HALF_PI) % math.pi)-HALF_PI if gap else 0.
    along = y @ np.array([math.sin(angle), math.cos(angle)])
    perpendicular = y @ np.array([math.cos(angle), -math.sin(angle)])
    j = float(np.sum((along*perpendicular)**2))
    eta = 1. if j <= 1e-300 or gap <= 0 else min(1., gap/j)
    return dict(angle=angle, curv=gap, J=j, eta=eta, precision=eta*gap,
                eigenvalues=eigenvalues, scatter=scatter)


def gaussian_angle_cdf(angle, *, fitted_angle, precision):
    """Continuous CDF proportional to exp[-precision*sin^2(theta-fit)/2].

    The domain is (0,pi/2), prior uniform, and precision may be zero. A
    mass-independent baseline is removed for negative fitted lines. Peak and
    boundary splits assist quadrature; its returned error is an estimate,
    not a rigorous bound or a certification of extreme floating-point inputs.
    """
    theta, fitted, kappa = float(angle), float(fitted_angle), float(precision)
    if not all(math.isfinite(x) for x in (theta, fitted, kappa)) or kappa < 0:
        raise ValueError("angles and precision must be finite; precision must be nonnegative")
    if theta <= 0:
        return 0., 0.
    if theta >= HALF_PI:
        return 1., 0.
    if kappa == 0:
        return theta/HALF_PI, 0.
    fitted = (fitted+HALF_PI) % math.pi-HALF_PI
    inside = 0 <= fitted <= HALF_PI
    baseline = 0. if inside else min(math.sin(fitted)**2, math.cos(fitted)**2)

    def kernel(t):
        excess = max(0., math.sin(t-fitted)**2-baseline)
        return math.exp(-0.5*kappa*excess)

    root = math.sqrt(kappa)
    points = []
    centers = [0., HALF_PI] + ([fitted] if inside else [])
    for center in centers:
        width = 1/max(root, kappa*abs(math.sin(2*(center-fitted))))
        points.extend(center+multiple*width for multiple in (-10, -1, 0, 1, 10))

    def integral(stop):
        interior = sorted(set(x for x in points if 0 < x < stop))
        return integrate.quad(kernel, 0, stop, points=interior, limit=200,
                              epsabs=2e-12, epsrel=2e-11)

    numerator, ne = integral(theta)
    total, te = integral(HALF_PI)
    cdf = numerator/total
    return cdf, (ne+cdf*te)/total


def gaussian_angle_quantile(probability, *, fitted_angle, precision):
    """Quantile of gaussian_angle_cdf; physical mass is tan(theta) at s=1."""
    probability = float(probability)
    if not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    gaussian_angle_cdf(HALF_PI/2, fitted_angle=fitted_angle, precision=precision)
    return optimize.brentq(
        lambda theta: gaussian_angle_cdf(theta, fitted_angle=fitted_angle,
                                         precision=precision)[0]-probability,
        0, HALF_PI, xtol=2e-12)


def score_polynomial_moments(radii):
    """Exact independent-Gaussian aggregate moments at the true signal basis.

    A_i=r_i+Z_i, B_i=W_i; U_i=(A_i B_i)^2,
    V_i=(A_i^2-B_i^2)^2, W_i=A_i B_i(A_i^2-B_i^2).
    Signs of r_i do not matter. These known radii are proof inputs only.
    """
    r = np.asarray(radii, dtype=float)
    if r.ndim != 1 or np.any(~np.isfinite(r)):
        raise ValueError("radii must be a finite one-dimensional sequence")
    n, delta = len(r), float(np.sum(r*r))
    fourth, sixth = float(np.sum(r**4)), float(np.sum(r**6))
    return dict(mean_U=n+delta, var_U=8*n+16*delta+2*fourth,
                mean_V=4*n+4*delta+fourth,
                var_V=128*n+256*delta+144*fourth+16*sixth,
                mean_W=0., var_W=12*n+24*delta+9*fourth+sixth,
                cov_UV=8*n+16*delta+4*fourth, cov_UW=0., cov_VW=0.)

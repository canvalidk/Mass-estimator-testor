"""Small user-facing readouts of the existing law, with no rejection gate.

Intervals and log spread describe latent M=f/alpha conditional on the declared
Gaussian / flat-positive-magnitude model. Local uncertainty is an asymptotic
approximation, NOT a finite exact standard deviation of that ratio law.
"""
from dataclasses import dataclass
import math
import numpy as np
from estimator import gauss, positive_normal_mean, validate


def known_direction_mass(force_projection, acceleration_projection,
                         sigma_force, sigma_acceleration):
    """Exact scalar reduction for known orientation, independent isotropic noise.

    Inputs are signed projections on the supplied physical orientation, not
    measured vector norms with an orientation inferred from the same data.
    """
    return (positive_normal_mean(force_projection, sigma_force)
            / positive_normal_mean(acceleration_projection, sigma_acceleration))


def aligned_local_uncertainty(force, acceleration, covariance):
    """Leading small-noise mass spread at positive, nonzero aligned readings.

    Fits common direction as a nuisance parameter, retaining full covariance.
    No promise of accuracy at low signal: use the full interval in that case.
    """
    force, acceleration, cov = validate(force, acceleration, covariance)
    f, a = np.linalg.norm(force), np.linalg.norm(acceleration)
    if f == 0 or a == 0:
        raise ValueError("Local approximation requires both signals to be nonzero")
    u = acceleration / a
    if not np.allclose(force / f, u, rtol=0, atol=1e-10):
        raise ValueError("Local aligned formula requires positive alignment")
    d = force.size
    # Orthonormal coordinates in the tangent plane to the common direction.
    basis = np.linalg.svd(u[None, :], full_matrices=True)[2][1:].T
    jac = np.zeros((2*d, d+1))
    jac[:d, 0], jac[d:, 1] = u, u
    jac[:d, 2:], jac[d:, 2:] = f * basis, a * basis
    information = jac.T @ np.linalg.solve(cov, jac)
    gradient = np.r_[1/a, -f/a**2, np.zeros(d-1)]
    return math.sqrt(float(gradient @ np.linalg.solve(information, gradient)))


def log_mass_moments(result, *, order=32):
    """Mean log(M/mass_scale) and SD(log M) of the angle-density interpolant.

    This uses the same piecewise linear density as Estimate.quantile. Fourth-
    power substitutions resolve integrable log singularities at both endpoints.
    Refine the original estimator grid to check accuracy for a particular input.
    """
    if not isinstance(order, int) or order < 8:
        raise ValueError("Use integer quadrature order >= 8")
    theta, density = result.ratio_angle, result.angle_density
    x, w = gauss(order)
    t, w = (x+1)/2, w/2
    widths = np.diff(theta)
    # Ordinary cells, excluding the two singular endpoint cells.
    angle = theta[1:-2, None] + widths[1:-1, None] * t
    weights = (widths[1:-1, None] * w
               * (density[1:-2, None]*(1-t) + density[2:-1, None]*t))
    logs = [np.log(np.tan(angle)).ravel()]
    masses = [weights.ravel()]
    for left in (True, False):
        h = widths[0] if left else widths[-1]
        z = t**4
        endpoint, inner = ((density[0], density[1]) if left
                           else (density[-1], density[-2]))
        masses.append(w * (4*h*t**3) * (endpoint*(1-z) + inner*z))
        log = np.log(np.tan(h*z))
        logs.append(log if left else -log)
    logs, masses = np.concatenate(logs), np.concatenate(masses)
    masses /= masses.sum()
    mean = float(masses @ logs)
    spread = math.sqrt(float(masses @ (logs-mean)**2))
    return mean, spread


@dataclass(frozen=True)
class MassReadout:
    mass: float
    lower: float
    upper: float
    probability: float
    log_standard_deviation: float


def summarize(result, *, probability=.95):
    """Point, equal-tail probability interval, and finite logarithmic spread.

    The interval is conditional model probability, not an automatically
    calibrated repeated-experiment confidence interval. Its centre need not
    equal the ratio-of-means point. Log SD does not assert a lognormal law.
    """
    if not math.isfinite(probability) or not 0 < probability < 1:
        raise ValueError("Probability must lie strictly between zero and one")
    tail = (1-probability)/2
    return MassReadout(result.mass, result.quantile(tail),
                       result.quantile(1-tail), probability,
                       log_mass_moments(result)[1])

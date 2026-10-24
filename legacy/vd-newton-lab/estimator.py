"""Exploratory reference for the VD compatible-pair ratio-of-means estimator.

Flat df d(alpha) dOmega measure; positive definite Gaussian covariance.
This is a numerical laboratory, not a production admission policy.
The radial integral is analytic; direction and ratio angle are quadrature.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np


@lru_cache(None)
def gauss(n):
    return np.polynomial.legendre.leggauss(n)


def directions(dimension, order):
    if dimension == 1:
        return np.array([[1.], [-1.]]), np.array([.5, .5])
    if dimension == 2:
        angle = np.arange(order) * (2 * math.pi / order)
        return np.column_stack((np.cos(angle), np.sin(angle))), np.full(order, 1 / order)
    if dimension == 3:
        z, w = gauss(order)
        phi = np.arange(2 * order) * math.pi / order
        z, phi = np.meshgrid(z, phi, indexing="ij")
        radius = np.sqrt(1 - z * z)
        return np.column_stack((radius.ravel() * np.cos(phi).ravel(),
                                radius.ravel() * np.sin(phi).ravel(), z.ravel())), np.repeat(w / (4 * order), 2 * order)
    raise ValueError("Only one, two, or three spatial dimensions are supported")


def radial_logs(t):
    """log J1, log J2 where Jk = integral_0^inf x^k exp(-x*x/2+t*x) dx.

    A scaled erfc recurrence handles t >= -4. For the negative tail,
    x=v/abs(t) gives a bounded, noncancelling Gauss integral on v in [0,48].
    The omitted tail there is < 4e-17 relative for J2 at t <= -4.
    """
    t = np.asarray(t, dtype=float)
    one, two = np.empty_like(t), np.empty_like(t)
    direct = t >= -4
    x = t[direct]
    phi_cdf = np.fromiter((.5 * math.erfc(-v / math.sqrt(2)) for v in x), float, count=x.size)
    k0 = math.sqrt(2 * math.pi) * phi_cdf
    edge = np.exp(-.5 * x * x)
    one[direct] = .5 * x * x + np.log(edge + x * k0)
    two[direct] = .5 * x * x + np.log((1 + x * x) * k0 + x * edge)
    rate = -t[~direct]
    if rate.size:
        v, w = gauss(64)
        v, w = 24 * (v + 1), 24 * w
        kernel = np.exp(-v[:, None] - .5 * (v[:, None] / rate) ** 2)
        one[~direct] = np.log((w * v) @ kernel) - 2 * np.log(rate)
        two[~direct] = np.log((w * v * v) @ kernel) - 3 * np.log(rate)
    return one, two


@dataclass
class Estimate:
    mass: float
    mean_force: float
    mean_acceleration: float
    ratio_angle: np.ndarray
    angle_density: np.ndarray
    mass_scale: float

    def quantile(self, p):
        if not 0 < p < 1:
            raise ValueError("Quantile probability must lie strictly between zero and one")
        dx = np.diff(self.ratio_angle)
        cdf = np.r_[0., np.cumsum(dx * (self.angle_density[:-1] + self.angle_density[1:]) / 2)]
        theta = np.interp(p, cdf / cdf[-1], self.ratio_angle)
        return float(self.mass_scale * np.tan(theta))

    def log_density(self):
        theta = self.ratio_angle[1:-1]
        return np.log(np.tan(theta)), self.angle_density[1:-1] * np.sin(theta) * np.cos(theta)


def validate(force, acceleration, covariance):
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    if force.ndim != 1 or force.shape != acceleration.shape or force.size not in (1, 2, 3):
        raise ValueError("Force and acceleration must have the same dimension (1, 2, or 3)")
    cov = np.asarray(covariance, float)
    if cov.shape != (2 * force.size, 2 * force.size) or not np.all(np.isfinite(cov)):
        raise ValueError("Supply the full finite covariance in [force, acceleration] order")
    if not np.all(np.isfinite(np.r_[force, acceleration])):
        raise ValueError("Readings must be finite")
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=0):
        raise ValueError("Covariance must be symmetric")
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Empirical covariance must be positive definite; use exact_recovery for exact inputs") from exc
    return force, acceleration, cov


def isotropic_covariance(dimension, sigma_force, sigma_acceleration):
    if not (math.isfinite(sigma_force) and math.isfinite(sigma_acceleration)
            and sigma_force > 0 and sigma_acceleration > 0):
        raise ValueError("Empirical standard uncertainties must be finite and positive")
    return np.diag([sigma_force ** 2] * dimension + [sigma_acceleration ** 2] * dimension)


def estimate(force, acceleration, covariance, *, direction_order=128, ratio_order=512, known_direction=None):
    if not isinstance(ratio_order, int) or ratio_order < 16 or ratio_order % 2 or not isinstance(direction_order, int) or direction_order < 4:
        raise ValueError("Use an even ratio_order >=16 and an integer direction_order >=4")
    force, acceleration, cov = validate(force, acceleration, covariance)
    d = force.size
    # Scaling is numerical only: df d(alpha) remains the specified flat measure.
    sf, sa = np.sqrt(np.trace(cov[:d, :d]) / d), np.sqrt(np.trace(cov[d:, d:]) / d)
    scale = np.r_[np.full(d, sf), np.full(d, sa)]
    y = np.r_[force, acceleration] / scale
    precision = np.linalg.inv(cov / np.outer(scale, scale))
    if known_direction is None:
        u, weights = directions(d, direction_order)
    else:
        u = np.asarray(known_direction, float)
        if u.shape != (d,) or not np.isfinite(u).all() or not np.isclose(u @ u, 1, atol=1e-12):
            raise ValueError("Known direction must be a finite unit vector of the input dimension")
        u, weights = u[None, :], np.ones(1)
    theta = np.linspace(0, math.pi / 2, ratio_order + 1)
    log_z, log_r = [], []
    # Integrate radial r at each ratio angle and common direction. The Jacobian
    # is r, hence J1 for probability and J2 for first magnitude moments.
    for start in range(0, theta.size, 32):
        th = theta[start:start + 32]
        v = np.concatenate((np.sin(th)[:, None, None] * u[None, :, :],
                            np.cos(th)[:, None, None] * u[None, :, :]), axis=2)
        aa = np.einsum("...i,ij,...j->...", v, precision, v)
        bb = np.einsum("...i,i->...", v, precision @ y)
        j1, j2 = radial_logs((bb / np.sqrt(aa)).ravel())
        l1 = j1.reshape(aa.shape) - np.log(aa)
        l2 = j2.reshape(aa.shape) - 1.5 * np.log(aa)
        # Constant -y' precision y/2 cancels in the normalized law.
        for logs, dest in ((l1, log_z), (l2, log_r)):
            peak = logs.max(axis=1)
            dest.extend(peak + np.log(np.exp(logs - peak[:, None]) @ weights))
    log_z, log_r = np.asarray(log_z), np.asarray(log_r)
    peak = max(log_z.max(), log_r.max())
    density, radial = np.exp(log_z - peak), np.exp(log_r - peak)
    def integrate(values):
        return (theta[1] - theta[0]) / 3 * (values[0] + values[-1] + 4 * values[1:-1:2].sum() + 2 * values[2:-1:2].sum())
    integral = integrate(density)
    fbar = sf * integrate(np.sin(theta) * radial) / integral
    abar = sa * integrate(np.cos(theta) * radial) / integral
    return Estimate(float(fbar / abar), float(fbar), float(abar), theta, density / integral, float(sf / sa))


def positive_normal_mean(reading, sigma):
    if sigma <= 0 or not math.isfinite(sigma) or not math.isfinite(reading):
        raise ValueError("Finite reading and positive standard uncertainty required")
    z = reading / sigma
    if z >= -4:
        return sigma * (z + math.exp(-z * z / 2) / (math.sqrt(2 * math.pi) * .5 * math.erfc(-z / math.sqrt(2))))
    # J1/J0 evaluated independently, avoiding cancellation in the inverse Mills ratio.
    v, w = gauss(80)
    v, w = 24 * (v + 1), 24 * w
    kernel = np.exp(-v - .5 * (v / -z) ** 2)
    return float(sigma / -z * (w * v @ kernel) / (w @ kernel))


def exact_recovery(force, acceleration):
    """Exact comparison of supplied floating representations; no tolerance repair."""
    f, a = np.asarray(force, float), np.asarray(acceleration, float)
    if f.ndim != 1 or f.shape != a.shape or not np.isfinite(np.r_[f, a]).all():
        raise ValueError("Matching finite exact vectors required")
    if not np.any(a):
        return ("not_identified", None) if not np.any(f) else ("contradiction", None)
    i = np.flatnonzero(a)[0]
    mass = f[i] / a[i]
    if mass > 0 and np.array_equal(f, mass * a):
        return "mass", float(mass)
    return "contradiction", None


def profile(force, acceleration, covariance, mass):
    f, a, cov = validate(force, acceleration, covariance)
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError("Fixed candidate mass must be finite and positive")
    d = f.size
    # Algebraically scaled for large mass: avoids a needless m*m overflow.
    if mass > 1:
        residual = f / mass - a
        residual_cov = cov[:d, :d] / mass / mass + cov[d:, d:] - (cov[:d, d:] + cov[d:, :d]) / mass
    else:
        residual = f - mass * a
        residual_cov = cov[:d, :d] + mass ** 2 * cov[d:, d:] - mass * (cov[:d, d:] + cov[d:, :d])
    return float(residual @ np.linalg.solve(residual_cov, residual))


def isotropic_qmin(force, acceleration, sigma_force, sigma_acceleration):
    f, a = np.asarray(force) / sigma_force, np.asarray(acceleration) / sigma_acceleration
    p, q, dot = float(f @ f), float(a @ a), float(f @ a)
    if dot <= 0:
        return min(p, q)
    denominator = p + q + math.hypot(p - q, 2 * dot)
    return max(0., 2 * (p * q - dot * dot) / denominator) if denominator else 0.


def compatibility_intervals(force, acceleration, sigma_force, sigma_acceleration, threshold):
    """All m>0 satisfying absolute Q(m)<=threshold, NOT a confidence claim.

    Intervals include finite roots; 0 and infinity are open parameter boundaries.
    Empty, singleton, bounded, one-sided, all-positive and two-ray sets are kept.
    """
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError("Nonnegative finite discrepancy threshold required")
    f, a = np.asarray(force) / sigma_force, np.asarray(acceleration) / sigma_acceleration
    aa, bb, cc = float(a @ a) - threshold, -2 * float(f @ a), float(f @ f) - threshold
    roots = []
    if aa == 0:
        if bb != 0:
            roots = [-cc / bb]
    else:
        disc = bb * bb - 4 * aa * cc
        if disc >= 0:
            root = math.sqrt(disc)
            stable = -.5 * (bb + math.copysign(root, bb))
            roots = [stable / aa, cc / stable] if stable else [0.]
    roots = sorted(set(float(r) for r in roots if r > 0 and math.isfinite(r)))
    edges = [0.] + roots + [math.inf]
    intervals = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        probe = (lo + hi) / 2 if math.isfinite(hi) else max(1., lo * 2)
        if (aa * probe + bb) * probe + cc <= 0:
            intervals.append([lo, hi])
    for root in roots:
        if not any(lo <= root <= hi for lo, hi in intervals):
            intervals.append([root, root])
    intervals.sort()
    merged = []
    for lo, hi in intervals:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [[lo * sigma_force / sigma_acceleration, hi * sigma_force / sigma_acceleration] for lo, hi in merged]


def gate(candidate_mass, qmin, threshold):
    if not all(math.isfinite(x) for x in (candidate_mass, qmin, threshold)) or candidate_mass <= 0 or qmin < 0 or threshold < 0:
        raise ValueError("Finite positive candidate, nonnegative discrepancy and threshold required")
    return {"status": "poor_model_compatibility" if qmin > threshold else "empirical_candidate",
            "returned_mass": None if qmin > threshold else candidate_mass,
            "conditional_mass": candidate_mass, "qmin": qmin, "threshold": threshold}

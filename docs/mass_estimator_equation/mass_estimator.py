"""VD mass estimator: ratio of compatible magnitude means, with uncertainty.

Python 3.10+; NumPy is the only external dependency. No workspace imports.
Model: y | (f, alpha, u) ~ N([f*u, alpha*u], covariance), f, alpha > 0,
with flat df d(alpha) and uniform common direction. Newton II is a premise.
Point = E[f]/E[alpha]. Intervals describe M=f/alpha under that same law.
See README.md for the complete specification and source/version check.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np


@dataclass(frozen=True)
class MassEstimate:
    """Mass-unit outputs use the supplied force-unit / acceleration-unit.

    interval is an equal-tail conditional probability interval, not a
    calibrated confidence interval. log_standard_deviation is dimensionless.
    numerical_change is the largest scaled change across the last two grids;
    it is a refinement diagnostic, not a bound on integration error.
    """
    mass: float
    interval: tuple[float, float]
    probability: float
    log_standard_deviation: float
    mean_force_magnitude: float
    mean_acceleration_magnitude: float
    numerical_change: float
    direction_order: int
    ratio_order: int

    @property
    def inverse_mass(self) -> float:
        return 1.0 / self.mass

    @property
    def inverse_mass_interval(self) -> tuple[float, float]:
        return (1.0 / self.interval[1], 1.0 / self.interval[0])


@lru_cache(maxsize=16)
def _gauss(order):
    return np.polynomial.legendre.leggauss(order)


def _directions(dimension, order):
    if dimension == 1:
        return np.array([[1.0], [-1.0]]), np.array([0.5, 0.5])
    if dimension == 2:
        phi = np.arange(order) * (2 * math.pi / order)
        return np.column_stack((np.cos(phi), np.sin(phi))), np.full(order, 1 / order)
    z, weights = _gauss(order)
    z, phi = np.meshgrid(z, np.arange(2 * order) * math.pi / order, indexing="ij")
    radius = np.sqrt(1 - z * z)
    u = np.column_stack((radius.ravel() * np.cos(phi).ravel(),
                         radius.ravel() * np.sin(phi).ravel(), z.ravel()))
    return u, np.repeat(weights / (4 * order), 2 * order)


def _radial_logs(t):
    """log integrals of x^k exp(-x^2/2+t*x), x>0, for k=1,2.

    Adapted from the original laboratory radial reduction. For t < -4,
    x=v/abs(t) avoids cancellation; the v>48 tail is negligible here.
    """
    t = np.asarray(t, dtype=float)
    one, two = np.empty_like(t), np.empty_like(t)
    direct = t >= -4
    v = t[direct]
    cdf = np.fromiter((0.5 * math.erfc(-x / math.sqrt(2)) for x in v),
                      dtype=float, count=v.size)
    k0 = math.sqrt(2 * math.pi) * cdf
    edge = np.exp(-0.5 * v * v)
    one[direct] = 0.5 * v * v + np.log(edge + v * k0)
    two[direct] = 0.5 * v * v + np.log((1 + v * v) * k0 + v * edge)
    rate = -t[~direct]
    if rate.size:
        v, w = _gauss(64)
        v, w = 24 * (v + 1), 24 * w
        kernel = np.exp(-v[:, None] - 0.5 * (v[:, None] / rate) ** 2)
        one[~direct] = np.log((w * v) @ kernel) - 2 * np.log(rate)
        two[~direct] = np.log((w * v * v) @ kernel) - 3 * np.log(rate)
    return one, two


def _simpson(values, step):
    return step / 3 * (values[0] + values[-1] + 4 * values[1:-1:2].sum()
                       + 2 * values[2:-1:2].sum())


def _quantile(theta, density, mass_scale, probability):
    """Invert the continuous piecewise-linear density, retaining both tails."""
    upper = probability > 0.5
    p = 1 - probability if upper else probability
    values = density[::-1] if upper else density
    widths = np.diff(theta)
    cell_mass = widths * (values[:-1] + values[1:]) / 2
    cumulative = np.r_[0.0, np.cumsum(cell_mass)]
    target = p * cumulative[-1]
    j = min(int(np.searchsorted(cumulative, target, side="right") - 1), len(widths) - 1)
    height_area = (target - cumulative[j]) / widths[j]
    a, delta = values[j], values[j + 1] - values[j]
    denominator = a + math.sqrt(max(0.0, a * a + 2 * delta * height_area))
    fraction = 0.0 if height_area == 0 else 2 * height_area / denominator
    angle = theta[j] + widths[j] * min(1.0, max(0.0, fraction))
    tangent = math.tan(angle)
    if tangent <= 0:
        raise ArithmeticError("Requested tail is not numerically resolved")
    return mass_scale / tangent if upper else mass_scale * tangent


def _log_spread(theta, density):
    """SD(log(M/s)) with endpoint transformations for the integrable logs."""
    x, w = _gauss(32)
    t, w = (x + 1) / 2, w / 2
    h = np.diff(theta)
    angle = theta[1:-2, None] + h[1:-1, None] * t
    weights = h[1:-1, None] * w * (
        density[1:-2, None] * (1 - t) + density[2:-1, None] * t)
    logs, masses = [np.log(np.tan(angle)).ravel()], [weights.ravel()]
    for left in (True, False):
        width = h[0] if left else h[-1]
        endpoint, inner = ((density[0], density[1]) if left
                           else (density[-1], density[-2]))
        z = t ** 4
        masses.append(w * 4 * width * t ** 3 * (endpoint * (1 - z) + inner * z))
        values = np.log(np.tan(width * z))
        logs.append(values if left else -values)
    logs, masses = np.concatenate(logs), np.concatenate(masses)
    masses /= masses.sum()
    mean = float(masses @ logs)
    return math.sqrt(float(masses @ (logs - mean) ** 2))


def _calculate(force, acceleration, cov, probability, direction_order,
               ratio_order, known_direction):
    d = len(force)
    sf = math.sqrt(float(np.mean(np.diag(cov)[:d])))
    sa = math.sqrt(float(np.mean(np.diag(cov)[d:])))
    scale = np.r_[np.full(d, sf), np.full(d, sa)]
    y = np.r_[force, acceleration] / scale
    scaled_cov = cov / np.outer(scale, scale)
    theta = np.linspace(0, math.pi / 2, ratio_order + 1)
    sn, cs = np.sin(theta), np.cos(theta)

    # Exact direction/radius reduction for independent, isotropic 3D errors.
    # Testing exact diagonal equality prevents silently dropping correlation.
    diagonal = np.diag(cov)
    isotropic = (d == 3 and known_direction is None
                 and np.array_equal(cov, np.diag(diagonal))
                 and np.all(diagonal[:d] == diagonal[0])
                 and np.all(diagonal[d:] == diagonal[d]))
    if isotropic:
        h2 = np.sum((sn[:, None] * y[:d] + cs[:, None] * y[d:]) ** 2, axis=1)
        h = np.sqrt(h2)
        erf = np.fromiter((math.erf(float(v) / math.sqrt(2)) for v in h), float)
        factor = np.full_like(h, math.sqrt(2 / math.pi))
        np.divide(erf, h, out=factor, where=h > 1e-12)
        radial = np.exp((h2 - h2.max()) / 2)
        density = radial * factor
        # Common sqrt(pi/2) cancels between probability and first moments.
    else:
        precision = np.linalg.solve(scaled_cov, np.eye(2 * d))
        if known_direction is None:
            u, weights = _directions(d, direction_order)
        else:
            u, weights = known_direction[None, :], np.ones(1)
        log_z, log_r = [], []
        for start in range(0, len(theta), 16):
            th = theta[start:start + 16]
            v = np.concatenate((np.sin(th)[:, None, None] * u[None, :, :],
                                np.cos(th)[:, None, None] * u[None, :, :]), axis=2)
            aa = np.einsum("...i,ij,...j->...", v, precision, v)
            bb = np.einsum("...i,i->...", v, precision @ y)
            j1, j2 = _radial_logs((bb / np.sqrt(aa)).ravel())
            l1 = j1.reshape(aa.shape) - np.log(aa)
            l2 = j2.reshape(aa.shape) - 1.5 * np.log(aa)
            for logs, dest in ((l1, log_z), (l2, log_r)):
                peak = logs.max(axis=1)
                dest.extend(peak + np.log(np.exp(logs - peak[:, None]) @ weights))
        log_z, log_r = np.asarray(log_z), np.asarray(log_r)
        peak = max(log_z.max(), log_r.max())
        density, radial = np.exp(log_z - peak), np.exp(log_r - peak)

    step = theta[1] - theta[0]
    normalizer = _simpson(density, step)
    fbar = sf * _simpson(sn * radial, step) / normalizer
    abar = sa * _simpson(cs * radial, step) / normalizer
    density /= normalizer
    tail = (1 - probability) / 2
    lower = _quantile(theta, density, sf / sa, tail)
    upper = _quantile(theta, density, sf / sa, 1 - tail)
    sd = _log_spread(theta, density)
    return np.array([fbar / abar, lower, upper, sd, fbar, abar], dtype=float)


def estimate_mass(force, acceleration, covariance, *, probability=0.95,
                  known_direction=None, direction_order=24, ratio_order=512,
                  rtol=2e-4, max_refinements=3) -> MassEstimate:
    """Calculate the specified mass point, interval and finite log spread.

    force, acceleration: finite 1D arrays of the same length d=1,2,3.
    covariance: finite symmetric positive-definite (2d,2d) matrix, ordered
      [F1,...,Fd,a1,...,ad]. Includes every relevant error covariance.
    probability: central equal-tail probability, strictly between 0 and 1.
    known_direction: optional externally supplied oriented unit vector.
      None integrates a uniform unknown direction; never estimate this option
      from the same input vectors and treat it as known.
    direction_order, ratio_order: starting quadrature resolution. Ratio order
      is an even integer >=32; direction order is an integer >=8.
    rtol: refinement tolerance (0 < rtol < 1). Compares relative changes in
      point, endpoints and magnitude means; log-SD uses scale max(1, log-SD).
    max_refinements: 1..5 grid doublings, checking every returned quantity.

    Each successful call compares at least two resolutions. Numerical error
    or lack of observed convergence raises ArithmeticError, never a physical
    rejection. The convergence diagnostic is not a certified error bound.
    Unknown calibration variance, exact zero covariance, hard bounds and
    changing-excitation multi-trial inference are outside this interface.
    """
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    cov = np.asarray(covariance, float)
    if (force.ndim != 1 or force.shape != acceleration.shape
            or force.size not in (1, 2, 3)):
        raise ValueError("Force and acceleration must be matching 1D arrays of length 1, 2 or 3")
    d = len(force)
    if (cov.shape != (2 * d, 2 * d) or not np.isfinite(cov).all()
            or not np.isfinite(np.r_[force, acceleration]).all()):
        raise ValueError("Supply finite vectors and a finite (2d, 2d) covariance")
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=0):
        raise ValueError("Covariance must be symmetric")
    cov = (cov + cov.T) / 2  # remove only validated roundoff asymmetry
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Covariance must be positive definite; exact inputs are a separate operation") from exc
    if not math.isfinite(probability) or not 0 < probability < 1:
        raise ValueError("probability must be strictly between zero and one")
    if not math.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must lie strictly between zero and one")
    for value, label, minimum, maximum in (
            (direction_order, "direction_order", 8, 512),
            (ratio_order, "ratio_order", 32, 65536),
            (max_refinements, "max_refinements", 1, 5)):
        if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
            raise ValueError(f"{label} must be an integer in [{minimum}, {maximum}]")
    if ratio_order % 2:
        raise ValueError("ratio_order must be even")
    if known_direction is not None:
        known_direction = np.asarray(known_direction, float)
        if (known_direction.shape != (d,) or not np.isfinite(known_direction).all()
                or not np.isclose(known_direction @ known_direction, 1, rtol=0, atol=1e-12)):
            raise ValueError("known_direction must be an oriented unit vector of length d")
        known_direction = known_direction / np.linalg.norm(known_direction)
    previous = None
    for refinement in range(max_refinements + 1):
        nd, nr = direction_order * 2 ** refinement, ratio_order * 2 ** refinement
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
                current = _calculate(force, acceleration, cov, probability, nd, nr, known_direction)
        except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
            raise ArithmeticError("Numerical range exceeded; rescale units or review the covariance") from exc
        if not np.isfinite(current).all() or np.any(current <= 0):
            raise ArithmeticError("Mass calculation is not numerically resolved")
        if previous is not None:
            scales = np.maximum(np.abs(previous), np.abs(current))
            scales[3] = max(1.0, scales[3])
            change = float(np.max(np.abs(current - previous) / scales))
            if change <= rtol:
                mass, lower, upper, sd, fbar, abar = map(float, current)
                return MassEstimate(mass, (lower, upper), float(probability), sd,
                                    fbar, abar, change, nd, nr)
        previous = current
    raise ArithmeticError(f"Quadrature did not converge: last scaled change {change:.3g} > {rtol:g}. "
                          "Increase starting orders or max_refinements; inspect difficult cases separately.")

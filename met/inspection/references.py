"""Independent compact-angle references for one unrestricted isotropic reading.

These integrate the direct Newton residual cost and a closed-form Rice cloud mean.
Weights are supplied by the caller; these functions do not justify an estimator's
choice of g or eta. QUADPACK errors are estimates, not certified bounds.
"""
import math
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.special import i0e, i1e

from ..full import Reading
from ..misalignment import isotropic_minima


def _problem(observation, prior_center, weight):
    reading = Reading(observation.force, observation.acceleration,
                      observation.covariance, observation.subspace)
    bounds = isotropic_minima(reading)
    if not math.isfinite(weight) or weight < 0 or not math.isfinite(prior_center) or prior_center <= 0:
        raise ValueError("reference needs a nonnegative weight and positive prior center")
    d = reading.dimension
    sf, sa = math.sqrt(reading.cov[0, 0]), math.sqrt(reading.cov[d, d])
    ratio, c = sf / sa, prior_center
    f, a = reading.force / sf, reading.acceleration / sa

    def kernel(phi):
        co, si = math.cos(phi), math.sin(phi)
        residual = co * f - si * a
        excess = max(0.0, float(residual @ residual) - bounds["positive_min_cost"])
        denominator = math.hypot(c * co, ratio * si)
        jacobian = (c / denominator) * (ratio / denominator)
        return jacobian * math.exp(-0.5 * weight * excess)

    energy = float(f @ f + a @ a)
    width = min(0.1, 1 / math.sqrt(1 + weight * energy))
    edge_width = min(0.1, 1 / (1 + weight * energy))
    peak = math.atan(bounds["positive_min_mass"] / ratio) if bounds["positive_min_mass"] else math.pi / 4
    cuts = [0.0, math.pi / 2, math.pi / 4, peak]
    for multiplier in (1, 4, 16):
        cuts.extend([edge_width * multiplier, math.pi / 2 - edge_width * multiplier,
                     peak - width * multiplier, peak + width * multiplier,
                     math.atan2(multiplier * c, ratio), math.atan2(c, multiplier * ratio)])
    cuts = sorted(set(x for x in cuts if 0 <= x <= math.pi / 2))
    return dict(reading=reading, bounds=bounds, sf=sf, sa=sa, ratio=ratio,
                f=f, a=a, kernel=kernel, cuts=cuts, warnings=[])


def _integral(problem, lo, hi, multiplier=None):
    points = [lo, *[x for x in problem["cuts"] if lo < x < hi], hi]
    kernel = problem["kernel"]
    function = kernel if multiplier is None else lambda phi: kernel(phi) * multiplier(phi)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", IntegrationWarning)
        values = [quad(function, x, y, epsabs=1e-13, epsrel=2e-11, limit=100)
                  for x, y in zip(points[:-1], points[1:]) if x < y]
    problem["warnings"].extend(str(w.message) for w in caught)
    return sum(x[0] for x in values), sum(x[1] for x in values)


def angle_reference(observation, prior_center, weight, masses):
    """CDF from direct residual cost in phi=atan(m/(sigma_F/sigma_a)).

    A constant correction cost leaves exactly the supplied half-Cauchy prior and
    uses its analytic CDF. Otherwise adaptive angle integration is independent of
    the estimator's native log-mass grid. Warnings remain part of the result.
    """
    problem = _problem(observation, prior_center, weight)
    masses = np.asarray(masses, dtype=float)
    if masses.ndim != 1 or not np.all(np.isfinite(masses) & (masses > 0)):
        raise ValueError("reference masses must be a finite positive vector")
    if weight == 0 or problem["bounds"]["positive_cost_span"] == 0:
        return dict(cdf=2 / math.pi * np.arctan2(masses, prior_center), normalizer=math.pi / 2,
                    quadrature_error_bound=0.0, quadrature_warnings=[],
                    note="Analytic half-Cauchy CDF: the weighted correction cost is constant in mass.")
    z, z_error = _integral(problem, 0, math.pi / 2)
    if not math.isfinite(z) or z <= 0:
        raise ValueError("angle reference has no resolved normalizer")
    probabilities, errors = [], []
    for mass in masses:
        phi = math.atan(float(mass) / problem["ratio"])
        # Integrate the shorter angle interval, then take the complementary CDF
        # where appropriate. This avoids crossing a nearly empty central region.
        if phi <= math.pi / 4:
            value, error = _integral(problem, 0, phi)
            probabilities.append(value / z)
        else:
            value, error = _integral(problem, phi, math.pi / 2)
            probabilities.append(1 - value / z)
        errors.append(error / z + value * z_error / z ** 2)
    return dict(cdf=np.array(probabilities), normalizer=z,
                quadrature_error_bound=max(errors, default=0.0),
                quadrature_warnings=list(dict.fromkeys(problem["warnings"])),
                note="QUADPACK estimated error, not a rigorous bound; independent of the native log-mass grid.")


def cloud_readout_reference(observation, prior_center, weight):
    """Ratio of two independent Rice/angle integrals for a full 2D reading.

    With x=F/sigma_F, y=a/sigma_a and b=sin(phi)x+cos(phi)y,
    the acceleration cloud is sigma_a*cos(phi)*(b+Z), Z~N(0,I_2).
    Thus A=sigma_a*cos(phi)*RiceMean(b), while m*A uses sigma_F*sin(phi).
    Keeping these bounded factors avoids multiplying a tiny A by a large mass.
    Neither g nor eta changes this conditional cloud.
    """
    problem = _problem(observation, prior_center, weight)
    if problem["reading"].dimension != 2:
        raise ValueError("Rice cloud readout reference requires a full isotropic 2D reading")

    def rice(phi):
        b = math.sin(phi) * problem["f"] + math.cos(phi) * problem["a"]
        z = float(b @ b) / 4
        return math.sqrt(math.pi / 2) * ((1 + 2 * z) * i0e(z) + 2 * z * i1e(z))

    numerator, numerator_error = _integral(
        problem, 0, math.pi / 2, lambda phi: problem["sf"] * math.sin(phi) * rice(phi))
    denominator, denominator_error = _integral(
        problem, 0, math.pi / 2, lambda phi: problem["sa"] * math.cos(phi) * rice(phi))
    if not math.isfinite(numerator) or not math.isfinite(denominator) or numerator <= 0 or denominator <= 0:
        raise ValueError("cloud readout reference has unresolved positive integrals")
    point = numerator / denominator
    error = numerator_error / denominator + numerator * denominator_error / denominator ** 2
    return dict(point=point, numerator=numerator, denominator=denominator,
                numerator_error_estimate=numerator_error, denominator_error_estimate=denominator_error,
                quadrature_error_estimate=error,
                quadrature_warnings=list(dict.fromkeys(problem["warnings"])),
                note="Independent 2D Rice mean and angle integrals; QUADPACK errors are estimates, not rigorous bounds.")


def equal_perpendicular_point(center):
    """Exact cloud point when F and a are perpendicular unit vectors and s=1.

    A(m) is a constant times (1+m²)^-1/2. If
    B(c)=integral_0^1 c/[c²+(1-c²)u²]du, the point is B(1/c)/B(c).
    These equivalent forms avoid squaring a large c or taking acosh(1/c)
    when 1/c would overflow. They also apply to the all-zero isotropic cloud.
    """
    if isinstance(center, bool) or not math.isfinite(center) or center <= 0:
        raise ValueError("center must be finite and positive")
    c = float(center)
    if c == 1:
        return 1.0
    if c < 1:
        root = math.sqrt(1 - c) * math.sqrt(1 + c)
        return c * (math.log1p(root) - math.log(c)) / math.acos(c)
    return c * (math.acos(1 / c) / math.acosh(c))

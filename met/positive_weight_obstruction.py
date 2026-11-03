"""Research diagnostics for T46's scalar positive-misfit weighting obstruction.

The class is a fixed nonincreasing weight w(Q+) multiplying a reading's
Gaussian quadratic criterion. Its boundedness/monotonicity are mathematical
assumptions, not established by finite numerical sampling. These functions
do not change an estimator or cover nonlinear residual-dependent factors.
All quadrature controls and the weight's value at zero are explicit.
"""
import math

import numpy as np
from scipy import integrate


def _controls(epsabs, epsrel, limit):
    if (not math.isfinite(epsabs) or not math.isfinite(epsrel)
            or epsabs <= 0 or epsrel <= 0 or isinstance(limit, bool)
            or not isinstance(limit, int) or limit < 1):
        raise ValueError("positive finite tolerances and an integer limit are required")
    return dict(epsabs=epsabs, epsrel=epsrel, limit=limit)


def _weight(weight, zero_weight):
    w0 = float(zero_weight)
    if not math.isfinite(w0) or w0 <= 0 or not callable(weight):
        raise ValueError("a callable weight and positive finite zero_weight are required")
    if not math.isclose(float(weight(0.)), w0, rel_tol=1e-12, abs_tol=0.):
        raise ValueError("zero_weight must equal the supplied function at zero")

    def evaluated(q):
        result = float(weight(q))
        if not math.isfinite(result) or not 0 <= result <= w0:
            raise ValueError("quadrature encountered a weight outside [0, zero_weight]")
        return result

    return evaluated, w0


def null_cross_moment(weight, *, zero_weight, epsabs, epsrel, limit):
    """B_w=(1/pi)*int exp(-q)[w(0)-w(q)]dq and its error estimate.

    The formula is exact for bounded measurable weights. Nonconstant
    nonincreasing weights make B_w strictly positive. Numerical errors are
    estimates, not outward-rounded certificates of an arbitrary callable.
    """
    evaluated, w0 = _weight(weight, zero_weight)
    value, error = integrate.quad(lambda q: math.exp(-q)*(w0-evaluated(q)), 0, np.inf,
                                   **_controls(epsabs, epsrel, limit))
    return value/math.pi, error/math.pi


def weak_signal_certificate(weight, *, zero_weight, angle, radius, epsabs, epsrel, limit):
    """Evaluate the analytic Gaussian-shift inequalities used by T46.

    The world mean is radius*(sin(angle),cos(angle)); it is not passed to an
    estimator. The score's change is bounded by w(0)*sqrt(expm1(radius^2)).
    A strict score sign and positive crossmoment certify a wrong interior
    population fit when the underlying quadrature value is established.
    """
    theta, r = float(angle), float(radius)
    if not math.isfinite(theta) or not 0 < theta < math.pi/2 or not math.isfinite(r) or r <= 0:
        raise ValueError("angle must be positive and interior, and radius positive and finite")
    bias, error = null_cross_moment(weight, zero_weight=zero_weight,
                                    epsabs=epsabs, epsrel=epsrel, limit=limit)
    null_score = bias*math.cos(2*theta)
    shift = float(zero_weight)*math.sqrt(math.expm1(r*r))
    lower, upper = null_score-shift, null_score+shift
    cross_lower = bias-shift
    return dict(bias=bias, quadrature_error=error, null_score=null_score,
                score_shift_upper=shift, score_lower=lower, score_upper=upper,
                cross_lower=cross_lower,
                wrong_interior_fit=(cross_lower > 0 and (upper < 0 or lower > 0)))


def polar_weighted_scatter(mean, weight, *, zero_weight, angular_order,
                            epsabs, epsrel, limit):
    """Independent 2D polar quadrature of E_mu[w(Q+)YY^T].

    Angular Gauss-Legendre rules split at axes and diagonals. Each radial
    integral is adaptive, so a kink in w need not fall on a fixed radial
    node. radial_error does not include angular discretization error.
    """
    mu = np.asarray(mean, dtype=float)
    if mu.shape != (2,) or np.any(~np.isfinite(mu)):
        raise ValueError("mean must be a finite vector of length two")
    if isinstance(angular_order, bool) or not isinstance(angular_order, int) or angular_order < 2:
        raise ValueError("angular_order must be an integer at least two")
    evaluated, w0 = _weight(weight, zero_weight)
    controls = _controls(epsabs, epsrel, limit)
    nodes, weights = np.polynomial.legendre.leggauss(angular_order)
    result, radial_error = np.zeros((2, 2)), 0.
    norm2 = float(mu@mu)
    for sector in range(8):
        low, high = sector*math.pi/4, (sector+1)*math.pi/4
        angles = (low+high)/2+(high-low)*nodes/2
        angular_weights = (high-low)*weights/2
        for angle, aw in zip(angles, angular_weights):
            direction = np.array([math.cos(angle), math.sin(angle)])
            coefficient = 0. if direction[0]*direction[1] >= 0 else 2*min(direction**2)
            projection = float(mu@direction)

            def radial(t):
                w = w0 if coefficient == 0 else evaluated(coefficient*t)
                exponent = -t-.5*norm2+math.sqrt(2*t)*projection
                return 2*t*w*math.exp(exponent)

            value, error = integrate.quad(radial, 0, np.inf, **controls)
            result += aw*value*np.outer(direction, direction)/(2*math.pi)
            radial_error += aw*error/(2*math.pi)
    return dict(scatter=result, radial_error=radial_error)

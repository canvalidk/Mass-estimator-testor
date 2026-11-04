"""Research-only energy-class consistency certificates for T53.

The model is Y_i=r_i*u+Z_i with independent standard 2D Gaussian Z_i and
arbitrary fixed signed r_i. True energy is used to state/check a theorem,
not supplied to a mass estimator. The lower bound averages signs at one
fixed magnitude while keeping every nuisance vector's total energy exact.
No signal population, production estimator or replacement L11 is selected.
Analytic probability bounds are evaluated in ordinary floating point.
"""
import math

import numpy as np


_SMALL_SIGNAL_CONSTANT = math.exp(.5)*math.cosh(1.)


def _energy(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError("energy must be finite and nonnegative")
    return value


def _count(value):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < 1):
        raise ValueError("count must be a positive integer")
    return int(value)


def _correlation(value):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("correlation must lie in [0,1]")
    return value


def _log_sinh(value):
    if value == 0:
        return -math.inf
    if value <= 20:
        return math.log(math.sinh(value))
    return value-math.log(2)+math.log1p(-math.exp(-2*value))


def _log_sinh_from_log(log_value):
    return log_value if log_value < -36 else _log_sinh(math.exp(log_value))


def _log_softplus(log_value):
    return log_value if log_value < -36 else math.log(float(np.logaddexp(0., log_value)))


def _exp_or_infinity(value):
    return math.exp(value) if value <= math.log(np.finfo(float).max) else math.inf


def _testing_bound(log_kl_upper):
    kl_upper = _exp_or_infinity(log_kl_upper)
    log_tv = (log_kl_upper/2 if log_kl_upper < -36
              else math.log(-math.expm1(-kl_upper))/2)
    tv = math.exp(log_tv)
    # Rationalization avoids subtracting a TV bound rounded to one.
    log_error = -kl_upper-math.log(2)-math.log1p(tv)
    return dict(tv_upper=tv, log_tv_upper=log_tv, log_binary_error_lower=log_error,
                binary_error_lower=math.exp(log_error))


def _one_reading_bound(log_energy, correlation):
    energy = _exp_or_infinity(log_energy)
    if log_energy == -math.inf or correlation == 1:
        log_bound, log_small = -math.inf, -math.inf
    else:
        log_bound = (math.log(4)+energy/2
            +_log_sinh_from_log(log_energy+math.log1p(correlation)-math.log(2))
            +_log_sinh_from_log(log_energy+math.log1p(-correlation)-math.log(2)))
        log_small = (math.log(_SMALL_SIGNAL_CONSTANT)+math.log1p(-correlation)
                     +math.log1p(correlation)+2*log_energy)
    return dict(chi_square_upper=_exp_or_infinity(log_bound),
                log_chi_square_upper=log_bound,
                small_signal_chi_square_upper=(_exp_or_infinity(log_small)
                    if log_energy <= 0 else None),
                log_small_signal_chi_square_upper=(log_small if log_energy <= 0 else None))


def one_reading_relative_chi_square_bound(*, energy, correlation):
    """Bound chi2(q_v || q_u) for equal-sign mixtures at radius sqrt(t).

    With t=energy and rho=u'v, q_u/phi=exp(-t/2)*cosh(sqrt(t)*u'Y).
    The bound is 2*exp(t/2)*(cosh(t)-cosh(rho*t)). For t<=1 it is at
    most exp(.5)*cosh(1)*(1-rho^2)*t^2. The exact chi-square itself is
    not evaluated here; tests integrate it independently.
    """
    energy, correlation = _energy(energy), _correlation(correlation)
    return _one_reading_bound(math.log(energy) if energy else -math.inf, correlation)


def homogeneous_energy_testing_certificate(*, count, energy, correlation):
    """Two-direction minimax testing bound with exact fixed total energy.

    Every sign vector has r_i=+/-sqrt(energy/count), hence sum r_i^2=energy.
    Product KL<=N*log(1+chi2_one). Affinity and Cauchy--Schwarz give
    TV<=sqrt(1-exp(-KL_upper)); equal-weight binary error is bounded below.
    A finite log error bound is retained even if its exponential underflows.
    """
    count, energy, correlation = _count(count), _energy(energy), _correlation(correlation)
    log_count = math.log(count)
    log_per_reading = math.log(energy)-log_count if energy else -math.inf
    one = _one_reading_bound(log_per_reading, correlation)
    # Retain the logarithm before multiplication by N; a tiny per-reading
    # moment can be rescued by a large count or by the final square root.
    log_kl = log_count+_log_softplus(one["log_chi_square_upper"])
    result = _testing_bound(log_kl)
    result.update(per_reading_energy=_exp_or_infinity(log_per_reading),
                  log_per_reading_energy=log_per_reading,
                  one_reading_chi_square_upper=one["chi_square_upper"],
                  kl_upper=_exp_or_infinity(log_kl), log_kl_upper=log_kl,
                  small_signal_kl_upper=(_exp_or_infinity(
                      log_count+one["log_small_signal_chi_square_upper"])
                      if log_per_reading <= 0 else None))
    return result


def scatter_minimax_upper_bound(*, count, energy, angle_tolerance):
    """Uniform TLS line-angle probability bound over an energy-only class.

    For the principal line of sum Y_iY_i',
    P(projective angle error >=epsilon)<=min(1,24*(N+Delta)/(Delta^2*sin^2 epsilon)).
    E||scatter-E scatter||_F^2=6*(N+Delta), uniformly in signs/allocation.
    The estimator needs observations only; Delta is a theorem parameter.
    """
    count, energy = _count(count), _energy(energy)
    angle_tolerance = float(angle_tolerance)
    if not math.isfinite(angle_tolerance) or not 0 < angle_tolerance <= math.pi/2:
        raise ValueError("angle_tolerance must lie in (0,pi/2]")
    information = energy/(1+count/energy) if energy else 0.
    sine = math.sin(angle_tolerance)
    denominator = (information*sine)*sine
    probability = 24/denominator if denominator > 24 else 1.
    return dict(information=information,
                scatter_fluctuation_second_moment=6*(count+energy),
                line_error_probability_upper=probability)

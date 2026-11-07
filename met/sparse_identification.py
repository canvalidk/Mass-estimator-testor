"""Research-only sparse-allocation testing certificates for T52.

Each proof mixture averages allowed fixed worlds with mean either zero or
amplitude*u. Activations are an auxiliary device for a lower bound, not a
signal population declared to an estimator. No estimator, prior or selected
identification criterion is implemented here. TV outputs are analytic upper
bounds evaluated in floating point, not exact TV or certified enclosures.
"""
import math

import numpy as np


def _amplitude(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError("amplitude must be finite and nonnegative")
    squared = value*value
    if not math.isfinite(squared) or squared > math.log(np.finfo(float).max):
        raise ValueError("Gaussian second moment exceeds the floating-point range")
    return value, squared


def _probabilities(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or np.any(~np.isfinite(values))
            or np.any(values < 0) or np.any(values > 1)):
        raise ValueError("probabilities must be a finite one-dimensional sequence in [0,1]")
    return values


def _scale(value):
    value = float(value)
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError("activation scale must lie in (0,1]")
    return value


def _count(value):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < 1):
        raise ValueError("count must be a positive integer")
    return int(value)


def _log_expm1(value):
    if value == 0:
        return -math.inf
    if value < math.log(2):
        return math.log(math.expm1(value))
    return value+math.log1p(-math.exp(-value))


def _log_softplus(log_value):
    # log(log(1+exp(x))). Keep x itself when exp(x) would lose a small
    # moment that can become representable after taking a square root.
    if log_value < -36:
        return log_value
    return math.log(float(np.logaddexp(0., log_value)))


def _log_expm1_from_log(log_value):
    """log(expm1(x)) when log(x), rather than x, is representable."""
    if log_value < -36:
        return log_value
    if log_value > math.log(np.finfo(float).max):
        return math.inf
    return _log_expm1(math.exp(log_value))


def _log_sum(values):
    return float(np.logaddexp.reduce(np.asarray(values, dtype=float), initial=-math.inf))


def _tv_from_log_squared_moment(log_value):
    log_bound = (log_value-math.log(2))/2
    return 1. if log_bound >= 0 else math.exp(log_bound)


def positive_sparse_log_ratio(points, *, amplitude, angle, probability):
    """Log ratio of (1-p)N(0,I2)+p N(b*u_angle,I2) to N(0,I2).

    Both coordinates of the active mean are positive for an acute angle.
    There is no random-sign mixture. Any finite angle is accepted by this
    density helper; the product certificate below uses nonnegative angle
    correlation for its direction-independent envelope.
    """
    amplitude, squared = _amplitude(amplitude)
    angle, probability = float(angle), float(probability)
    if not math.isfinite(angle) or not math.isfinite(probability) or not 0 <= probability <= 1:
        raise ValueError("angle must be finite and probability must lie in [0,1]")
    points = np.asarray(points, dtype=float)
    if points.ndim == 0 or points.shape[-1] != 2 or np.any(~np.isfinite(points)):
        raise ValueError("points must be finite vectors with final dimension two")
    if probability == 0 or amplitude == 0:
        result = np.zeros(points.shape[:-1])
    else:
        direction = np.array([math.cos(angle), math.sin(angle)])
        with np.errstate(over="ignore", invalid="ignore"):
            shifted = amplitude*(points@direction)-squared/2
        if np.any(~np.isfinite(shifted)):
            raise ValueError("density ratio exceeds the floating-point range")
        result = shifted if probability == 1 else np.logaddexp(
            math.log1p(-probability), math.log(probability)+shifted)
    return float(result) if np.ndim(result) == 0 else result


def positive_sparse_product_certificate(probabilities, *, amplitude, correlation):
    """Exact product second moments and Cauchy--Schwarz testing bounds.

    For directions with rho=u'v in [0,1], the one-reading cross moment is
    1+p_i^2*(exp(rho*b^2)-1). Let A and B be the products of the self and
    cross moments. TV<=sqrt((A-B)/2), and also
    TV<=sqrt(expm1((exp(b^2)-1)*sum p_i^2)/2). Bounds are capped at one.
    Binary error uses equal testing weights on two hypotheses only.
    """
    probabilities = _probabilities(probabilities)
    amplitude, squared = _amplitude(amplitude)
    correlation = float(correlation)
    if not math.isfinite(correlation) or not 0 <= correlation <= 1:
        raise ValueError("correlation must lie in [0,1]")
    log_energy = 2*math.log(amplitude) if amplitude else -math.inf
    same_coefficient = _log_expm1_from_log(log_energy)
    cross_coefficient = (_log_expm1_from_log(log_energy+math.log(correlation))
                         if correlation else -math.inf)
    # log(exp(b^2)-exp(rho*b^2)), retaining tiny b and near-unit rho.
    difference = (correlation*squared+_log_expm1_from_log(log_energy+math.log1p(-correlation))
                  if correlation < 1 else -math.inf)
    with np.errstate(divide="ignore"):
        log_squares = 2*np.log(probabilities)
    cross_terms = np.logaddexp(0., log_squares+cross_coefficient)
    log_cross = float(np.sum(cross_terms))
    # Each item is log(log(1+p_i^2*(e^b²-e^rho*b²)/(1+p_i²*(e^rho*b²-1)))).
    # Summing those logarithmically retains even an underflowing log(A/B).
    log_log_gap = _log_sum([_log_softplus(value+difference-cross)
                           for value, cross in zip(log_squares, cross_terms)])
    gap = math.exp(log_log_gap)
    log_same = log_cross+gap
    direct = _tv_from_log_squared_moment(log_cross+_log_expm1_from_log(log_log_gap))
    log_square_sum = _log_sum(log_squares)
    square_sum = math.exp(log_square_sum)
    log_chi_upper = _log_expm1_from_log(same_coefficient+log_square_sum)
    envelope = _tv_from_log_squared_moment(log_chi_upper)
    tv = min(direct, envelope)
    return dict(log_second_moment=log_same, log_cross_moment=log_cross,
                probability_square_sum=square_sum,
                log_probability_square_sum=log_square_sum,
                expected_active_count=float(np.sum(probabilities)),
                tv_direct_upper=direct, tv_envelope_upper=envelope,
                tv_upper=tv, binary_error_lower=(1-tv)/2)


def harmonic_activation_probabilities(*, scale, count):
    """Proof schedule p_i=scale/(i+1), i=1,...,count; no sampling occurs."""
    scale, count = _scale(scale), _count(count)
    return scale/np.arange(2, count+2, dtype=float)


def harmonic_sparse_certificate(*, amplitude, activation_scale):
    """One uniform-in-N positive-angle TV envelope for the fixed schedule.

    Sum p_i diverges, while sum p_i^2=c^2*(pi^2/6-1). Almost every fixed
    activation stream therefore has infinitely many active readings, but
    the direction mixtures can retain strictly positive binary testing error.
    """
    amplitude, _ = _amplitude(amplitude)
    scale = _scale(activation_scale)
    log_square_sum = 2*math.log(scale)+math.log(math.pi**2/6-1)
    log_energy = 2*math.log(amplitude) if amplitude else -math.inf
    log_upper = _log_expm1_from_log(log_energy)+log_square_sum
    tv = _tv_from_log_squared_moment(_log_expm1_from_log(log_upper))
    return dict(probability_square_sum=math.exp(log_square_sum),
                log_probability_square_sum=log_square_sum,
                log_second_moment_upper=math.exp(log_upper),
                tv_upper=tv, binary_error_lower=(1-tv)/2)


def harmonic_expected_active_bounds(*, scale, count):
    """Integral bounds on sum p_i, without allocating a large schedule."""
    scale, count = _scale(scale), _count(count)
    return dict(lower=scale*(math.log(count+2)-math.log(2)),
                upper=scale*math.log(count+1))


def active_support_scatter_moments(*, amplitude, active_count):
    """Oracle known-support comparison, not a support estimator.

    Summing only K known active readings gives mean scatter K(I+b^2 uu'),
    eigenvalue gap K*b^2 and squared Frobenius fluctuation 6K(1+b^2).
    An infinite known support therefore admits a consistent principal line
    at every b>0, with independent sign invariance. The support is extra
    oracle information, absent from the counterexample's declared inputs.
    """
    _, squared = _amplitude(amplitude)
    active_count = _count(active_count)
    return dict(mean_eigenvalue_gap=active_count*squared,
                fluctuation_second_moment=6*active_count*(1+squared),
                oracle_information=active_count*squared*squared/(1+squared))

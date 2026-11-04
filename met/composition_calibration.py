"""Finite-law diagnostics for T51's conditional composition obstruction.

This compares a fixed T47 normal-score rank warp with *independent*
composition. It neither supplies a joint law from marginals nor treats that
commutation requirement as a VD ruling. A supplied dependent joint can also
be pushed forward explicitly. No estimator or production adapter is changed.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy import special

from met.cdf_calibration import normal_log_cdf_warp


@dataclass(frozen=True)
class FinitePositiveLaw:
    """A sorted proper finite positive law from support and nonnegative weights.

    The supplied ``probabilities`` are normalized as weights; repeated support
    values are combined and zero-weight values removed. Arrays are read-only.
    This normalization introduces no distributional or dependence assumption.
    """
    support: np.ndarray
    probabilities: np.ndarray

    def __post_init__(self):
        support = np.asarray(self.support, dtype=float)
        weights = np.asarray(self.probabilities, dtype=float)
        if (support.ndim != 1 or weights.shape != support.shape or not support.size
                or np.any(~np.isfinite(support)) or np.any(support <= 0)
                or np.any(~np.isfinite(weights)) or np.any(weights < 0)):
            raise ValueError("support must be finite and positive with matching nonnegative weights")
        maximum = float(np.max(weights))
        if maximum <= 0:
            raise ValueError("weights must have finite positive total")
        weights = weights/maximum
        total = math.fsum(float(x) for x in weights)
        unique, inverse = np.unique(support, return_inverse=True)
        combined = np.array([math.fsum(float(x) for x in weights[inverse == j])/total
                             for j in range(len(unique))])
        positive = combined > 0
        support, probabilities = unique[positive], combined[positive]
        support.setflags(write=False)
        probabilities.setflags(write=False)
        object.__setattr__(self, "support", support)
        object.__setattr__(self, "probabilities", probabilities)

    def cdf(self, value):
        """CDF at explicit evaluation points, retaining the atomic convention."""
        value = np.asarray(value, dtype=float)
        if np.any(np.isnan(value)):
            raise ValueError("CDF evaluation points must not be NaN")
        result = np.array([math.fsum(float(p) for x, p in
                                    zip(self.support, self.probabilities) if x <= v)
                           for v in value.flat]).reshape(value.shape)
        return float(result) if result.ndim == 0 else result


def _law(value):
    if not isinstance(value, FinitePositiveLaw):
        raise ValueError("expected a FinitePositiveLaw")
    return value


def _log_probability(value):
    return math.log(value) if value else -math.inf


def _difference(log_larger, log_smaller):
    if log_larger == -math.inf:
        return 0.
    return math.exp(log_larger)*(-math.expm1(log_smaller-log_larger))


def warp_finite_law(law, *, eta):
    """Apply T47's H_eta to an already normalized finite law's CDF.

    Eta is explicit and fixed in a composition comparison. Both cumulative
    and survival ranks are formed from the original weights, avoiding loss
    of a small outer atom when its complementary CDF rounds to one. This is
    floating-point diagnostic code, not an arbitrary-precision rank oracle.
    """
    law = _law(law)
    # Validate eta even for point masses and the identity case.
    normal_log_cdf_warp(0., eta)
    if float(eta) == 1. or len(law.support) == 1:
        return law
    probabilities = []
    for j in range(len(law.support)):
        lower = math.fsum(float(x) for x in law.probabilities[:j])
        upper = math.fsum(float(x) for x in law.probabilities[j+1:])
        if lower <= upper:
            end = math.fsum((lower, float(law.probabilities[j])))
            larger = normal_log_cdf_warp(_log_probability(end), eta)
            smaller = normal_log_cdf_warp(_log_probability(lower), eta)
        else:
            start = math.fsum((upper, float(law.probabilities[j])))
            larger = normal_log_cdf_warp(_log_probability(start), eta)
            smaller = normal_log_cdf_warp(_log_probability(upper), eta)
        probabilities.append(_difference(larger, smaller))
    return FinitePositiveLaw(law.support, probabilities)


def _operation(left, right, operation):
    if operation == "addition":
        with np.errstate(over="ignore"):
            result = left+right
    elif operation == "parallel":
        small, large = np.minimum(left, right), np.maximum(left, right)
        result = small/(1+small/large)
    else:
        raise ValueError("operation must be explicitly 'addition' or 'parallel'")
    if np.any(~np.isfinite(result)) or np.any(result <= 0):
        raise ValueError("composed positive support exceeds floating-point range")
    return result


def pushforward_joint(left_support, right_support, joint_weights, *, operation):
    """Push a *supplied* finite joint law through the selected operation.

    Joint weights are normalized, not replaced by a product. This function
    makes no claim that a joint distribution can be inferred from marginals.
    """
    left, right = np.asarray(left_support, dtype=float), np.asarray(right_support, dtype=float)
    weights = np.asarray(joint_weights, dtype=float)
    if (left.ndim != 1 or right.ndim != 1 or not left.size or not right.size
            or weights.shape != (len(left), len(right))
            or np.any(~np.isfinite(left)) or np.any(left <= 0)
            or np.any(~np.isfinite(right)) or np.any(right <= 0)):
        raise ValueError("positive finite support and a matching joint weight matrix are required")
    support = _operation(left[:, None], right[None, :], operation)
    return FinitePositiveLaw(support.ravel(), weights.ravel())


def compose_independent(first, second, *, operation):
    """Compose positive laws under an explicit *independence* assumption."""
    first, second = _law(first), _law(second)
    return pushforward_joint(first.support, second.support,
                             np.outer(first.probabilities, second.probabilities),
                             operation=operation)


def gaussian_linear_coverage(actual_covariance, quoted_covariance, gradient, *, probability):
    """Gaussian coverage of a quoted central interval for a linearized quantity.

    Both positive-definite covariances and the nonzero gradient are explicit.
    Correct marginal intervals do not assert a correct joint covariance. This
    diagnostic neither estimates dependence nor constructs a joint mass law.
    """
    actual = np.asarray(actual_covariance, dtype=float)
    quoted = np.asarray(quoted_covariance, dtype=float)
    gradient = np.asarray(gradient, dtype=float)
    probability = float(probability)
    if (gradient.ndim != 1 or not gradient.size or np.any(~np.isfinite(gradient))
            or not np.any(gradient != 0.)
            or actual.shape != (len(gradient), len(gradient)) or quoted.shape != actual.shape
            or np.any(~np.isfinite(actual)) or np.any(~np.isfinite(quoted))
            or not np.array_equal(actual, actual.T) or not np.array_equal(quoted, quoted.T)
            or not math.isfinite(probability) or not 0 < probability < 1):
        raise ValueError("finite symmetric covariances, a nonzero gradient and probability in (0,1) are required")
    try:
        np.linalg.cholesky(actual)
        np.linalg.cholesky(quoted)
    except np.linalg.LinAlgError as error:
        raise ValueError("actual and quoted covariances must be positive definite") from error
    actual_variance = float(gradient@actual@gradient)
    quoted_variance = float(gradient@quoted@gradient)
    if (not math.isfinite(actual_variance) or not math.isfinite(quoted_variance)
            or actual_variance <= 0 or quoted_variance <= 0):
        raise ValueError("propagated variances must be finite and positive")
    z = -float(special.ndtri((1-probability)/2))
    standardized = z*math.sqrt(quoted_variance/actual_variance)
    coverage = float(special.erf(standardized/math.sqrt(2)))
    return dict(coverage=coverage, actual_variance=actual_variance,
                quoted_variance=quoted_variance)

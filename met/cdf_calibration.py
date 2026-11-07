"""Research normal-score CDF recalibration for T47, not a mass estimator.

For an already NORMALIZED raw CDF F, the corrected CDF is H_eta(F), where
H_eta(u)=Phi(sqrt(eta)*Phi^-1(u)). The scalar eta is an explicit calibration
input to these diagnostics, not a new physical input or a selected estimate.
This differs from powering an unnormalized likelihood and preserves medians.
Log-probability routines retain extreme ranks needed when eta is small.
"""
import math

import numpy as np
from scipy import special


def _eta(value):
    value = float(value)
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError("eta must be finite and lie in (0, 1]")
    return value


def _probability(value, *, interior=False):
    value = np.asarray(value, dtype=float)
    invalid = (value <= 0) | (value >= 1) if interior else (value < 0) | (value > 1)
    if np.any(~np.isfinite(value)) or np.any(invalid):
        raise ValueError("probabilities must lie in (0, 1)" if interior
                         else "probabilities must lie in [0, 1]")
    return value


def _scalar_or_array(value):
    return float(value) if value.ndim == 0 else value


def _log_normal_cdf(z):
    """Keep representable subnormal survival probabilities at positive z."""
    z = np.asarray(z, dtype=float)
    result = np.array(special.log_ndtr(z), copy=True)
    positive = z > 0
    with np.errstate(under="ignore"):
        result[positive] = np.log1p(-np.exp(special.log_ndtr(-z[positive])))
    return result


def normal_cdf_warp(probability, eta):
    """H_eta(u), including its exact endpoint values 0 and 1.

    Compose with a raw CDF only after normalizing that raw law. The input
    eta must be positive; no proper universal eta=0 extension exists for
    raw laws with full support on the positive mass line.
    """
    eta = _eta(eta)
    probability = _probability(probability)
    with np.errstate(under="ignore"):
        result = probability.copy() if eta == 1 else np.exp(_log_normal_cdf(math.sqrt(eta)*special.ndtri(probability)))
    return _scalar_or_array(result)


def normal_log_cdf_warp(log_probability, eta):
    """log H_eta(exp(log_probability)), without exponentiating the input.

    The same routine warps a raw upper-tail probability because H_eta is
    symmetric about probability 1/2. Minus infinity and zero are valid
    logarithmic endpoints. This avoids losing e.g. an exp(-1000) tail.
    """
    eta = _eta(eta)
    values = np.asarray(log_probability, dtype=float)
    if np.any(np.isnan(values)) or np.any(values > 0):
        raise ValueError("log probabilities must lie in [-infinity, 0]")
    result = values.copy() if eta == 1 else _log_normal_cdf(math.sqrt(eta)*special.ndtri_exp(values))
    return _scalar_or_array(result)


def raw_quantile_log_levels(probability, eta):
    """Raw CDF and survival log-ranks for a corrected probability quantile.

    Q_eta^-1(p)=F_raw^-1(Phi(Phi^-1(p)/sqrt(eta))). Return both log ranks
    rather than a rounded probability: at eta=1/90 the upper corrected
    97.5% quantile calls a raw survival probability around 10^-77, although
    its raw CDF rounds to exactly 1 in binary64. The caller must supply a
    raw inverse able to resolve the requested tail; no grid is certified.
    """
    eta = _eta(eta)
    probability = _probability(probability, interior=True)
    z = special.ndtri(probability)/math.sqrt(eta)
    return dict(log_cdf=_scalar_or_array(_log_normal_cdf(z)),
                log_sf=_scalar_or_array(_log_normal_cdf(-z)))

"""Independent moment integrals for O14's individually profiled Gram residual.

This is a research reference, not a scale correction in the estimator. Right-column
rotation applies to Q on the whole projective line, not to the positivity cost Q+.
Quantile quadrature convergence differences are estimates, not rigorous error bounds.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.special import roots_legendre
from scipy.stats import chi2, ncx2


def _arguments(k, noncentrality, nodes):
    if type(k) is not int or k not in (2, 3):
        raise ValueError("pooled Q moment reference requires dimension 2 or 3")
    if (isinstance(noncentrality, (bool, np.bool_))
            or not isinstance(noncentrality, (int, float, np.number))
            or not math.isfinite(noncentrality) or noncentrality < 0):
        raise ValueError("noncentrality must be finite and nonnegative")
    if type(nodes) is not int or not 8 <= nodes <= 64:
        raise ValueError("reference nodes must be an integer from 8 through 64")


@lru_cache(maxsize=256)
def _integrated_moments(k, noncentrality, nodes):
    grid, weights = roots_legendre(nodes)
    probabilities, weights = (grid + 1) / 2, weights / 2
    a = ncx2.ppf(probabilities, k, noncentrality)[:, None, None]
    w = chi2.ppf(probabilities, 1)[None, :, None]
    v = chi2.ppf(probabilities, k - 1)[None, None, :]
    total = a + w + v
    gap = np.hypot(a - w - v, 2 * np.sqrt(a * w))
    q = 2 * a * v / (total + gap)
    first = np.einsum("i,j,l,ijl->", weights, weights, weights, q, optimize=True)
    second = np.einsum("i,j,l,ijl->", weights, weights, weights, q * q, optimize=True)
    if not (np.isfinite(first) and np.isfinite(second)):
        raise ValueError("nonfinite pooled Q quadrature; reference range exceeded")
    return float(first), float(second)


def moments(k, noncentrality, nodes=32):
    """Moments of Q/v given total true squared SNR lambda (O14's exact integral).

    A~noncentral-chi-square(k,lambda), W~chi-square(1), V~chi-square(k-1)
    are independent, and Q/v=2AV/[A+W+V+hypot(A-W-V,2sqrt(AW))].
    The caller's investigator supplies lambda; the estimator is never given it.
    Null moments are exact. Other results use nodes and twice nodes in each variable.
    """
    _arguments(k, noncentrality, nodes)
    noncentrality = float(noncentrality)
    if noncentrality == 0:
        mean, second = ((2 - math.pi / 2, 10 - 3 * math.pi) if k == 2 else (1.0, 2.0))
        return dict(mean=mean, variance=second - mean * mean,
                    mean_error_estimate=0.0, variance_error_estimate=0.0,
                    method="exact central moments", fine_nodes=0)
    coarse, coarse_second = _integrated_moments(k, noncentrality, nodes)
    mean, second = _integrated_moments(k, noncentrality, 2 * nodes)
    variance, coarse_variance = max(0.0, second - mean * mean), max(0.0, coarse_second - coarse * coarse)
    return dict(mean=mean, variance=variance,
                mean_error_estimate=abs(mean - coarse),
                variance_error_estimate=abs(variance - coarse_variance),
                method="independent chi-square quantile quadrature; convergence differences are not bounds",
                fine_nodes=2 * nodes)

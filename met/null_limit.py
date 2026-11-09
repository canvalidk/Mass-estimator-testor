"""Deterministic research quadrature for T42's pure-null limiting law.

Use the doubled angle x=2*theta. Its random evidence is exp(t*cos(x-phi)),
where phi is uniform on the circle and t is chi-square with two degrees of
freedom. The optional beta multiplies this evidence strength; beta=1 is T42.
The prior is supplied as quadrature nodes and probability weights. This is
not a production estimator or a proof of convergence for finite data.
"""
import math

import numpy as np
from scipy.special import logsumexp, roots_laguerre


def average_null_update(support, weights, *, beta=1., queries=None,
                        radial_order=64, phase_order=512, fixed_strength=None):
    """Average *normalized* random updates, with independent response queries.

    Return averaged support probabilities, g(x)=E[L_x/Z], and dg/dx at the
    supplied queries (default: support). x may cover either the positive-mass
    half-circle or the whole projective circle. A fixed_strength replaces
    random t by that nonnegative value; beta still multiplies it. Independent
    query points far outside the prior support can have infinite response;
    finite quadrature must not be taken as a certificate of integrability.
    """
    x, p = np.asarray(support, dtype=float), np.asarray(weights, dtype=float)
    q = x if queries is None else np.asarray(queries, dtype=float)
    beta = float(beta)
    if (x.ndim != 1 or x.size == 0 or p.shape != x.shape or q.ndim != 1
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(q))
            or not np.all(np.isfinite(p)) or np.any(p < 0) or p.sum() <= 0
            or not math.isfinite(beta) or beta < 0):
        raise ValueError("finite 1D angles, nonnegative weights and beta required")
    if radial_order < 1 or phase_order < 4:
        raise ValueError("quadrature orders must be positive (phase_order >= 4)")
    p = p/p.sum()
    if fixed_strength is None:
        t, radial_weights = roots_laguerre(radial_order)
        t = 2*t
    else:
        fixed_strength = float(fixed_strength)
        if not math.isfinite(fixed_strength) or fixed_strength < 0:
            raise ValueError("fixed_strength must be finite and nonnegative")
        t, radial_weights = np.array([fixed_strength]), np.array([1.])
    phi = 2*math.pi*np.arange(phase_order)/phase_order
    support_cos = np.cos(x[:, None]-phi)
    query_cos, query_sin = np.cos(q[:, None]-phi), np.sin(q[:, None]-phi)
    log_p = np.full(p.shape, -np.inf)
    np.log(p, out=log_p, where=p > 0)
    averaged = np.zeros_like(p)
    response, derivative = np.zeros_like(q), np.zeros_like(q)
    for strength, weight in zip(beta*t, radial_weights):
        if weight == 0:
            continue
        log_joint = log_p[:, None]+strength*support_cos
        log_z = logsumexp(log_joint, axis=0)
        averaged += weight*np.exp(log_joint-log_z).mean(axis=1)
        weighted_response = np.exp(math.log(weight)+strength*query_cos-log_z)
        response += weighted_response.mean(axis=1)
        derivative -= strength*(weighted_response*query_sin).mean(axis=1)
    return dict(weights=averaged, response=response, derivative=derivative)

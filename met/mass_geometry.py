"""E9 research geometry; these helpers choose no production prior or size law.

See explore/e9-fisher-geometry.md. Structural Fisher information below is
conditional on an explicitly added Gaussian signal population.
"""

import math

import numpy as np
from scipy.linalg import eigh


def mass_curve(mass, covariance, subspace=None):
    """Scalar graph geometry in a fixed ambient covariance, including N.

    Returns P, dP/dm, the latent Gram G, transverse Gram E, and canonical
    squared speed g. The ambient data order is (force, acceleration).
    subspace is any finite full-column-rank basis; only its span affects P,g.
    No observation or latent signal is used to define this geometry.
    """
    m = float(mass)
    if not math.isfinite(m) or m <= 0:
        raise ValueError("a positive finite mass is required")
    cov = _covariance(covariance)
    if len(cov) % 2:
        raise ValueError("covariance must have equal force and acceleration dimensions")
    dimension = len(cov)//2
    n = _subspace(subspace, dimension)
    chol = np.linalg.cholesky(cov)
    basis = np.linalg.solve(chol, np.vstack([m*n, n]))
    derivative = np.linalg.solve(chol, np.vstack([n, np.zeros_like(n)]))
    u, r = np.linalg.qr(basis, mode="reduced")
    transverse = derivative-u @ (u.T @ derivative)
    half = transverse @ np.linalg.solve(r, u.T)
    dprojector = half+half.T
    gram = basis.T @ basis
    e = transverse.T @ transverse
    return {"projector": u @ u.T, "derivative": dprojector,
            "gram": gram, "transverse_gram": e,
            "metric": float(np.sum(dprojector*dprojector)/2)}


def functional_information(mass, signal, covariance, subspace=None):
    """Efficient mass Fisher information at a fixed physical latent signal.

    This Schur complement is not a normalized mass density or a completed
    multiparameter reference-prior construction.
    """
    state = mass_curve(mass, covariance, subspace)
    alpha = np.asarray(signal, dtype=float)
    if alpha.shape != (len(state["gram"]),) or not np.isfinite(alpha).all():
        raise ValueError("one finite latent coordinate per subspace basis column required")
    return float(alpha @ state["transverse_gram"] @ alpha)


def structural_fisher(mass, size, covariance, subspace=None):
    """Fisher matrix for whitened Z~N(0,I+size**2 P(m)), order (m,size).

    This expressly specified structural comparison model makes the latent
    acceleration covariance size**2 * G(m)**-1. It is not the functional
    model's Fisher information and is not licensed by L10 as a default.
    """
    omega = float(size)
    if not math.isfinite(omega) or omega < 0:
        raise ValueError("a finite nonnegative comparison size is required")
    state = mass_curve(mass, covariance, subspace)
    size2 = omega*omega
    weight = size2*(size2/(1+size2))
    nuisance = 2*len(state["gram"])*(omega/(1+size2))**2
    return np.diag([weight*state["metric"], nuisance])


def scalar_channel_metric(mass, ratios, correlations=None, weights=None):
    """Weighted product-curve metric for scalar instrument declarations.

    weights are explicit positive metric weights, not inferred reliability
    or a production population. Input arrays for mass are supported.
    """
    m = np.asarray(mass, dtype=float)
    s = np.asarray(ratios, dtype=float)
    if (not np.isfinite(m).all() or np.any(m <= 0) or s.ndim != 1 or len(s) == 0
            or not np.isfinite(s).all() or np.any(s <= 0)):
        raise ValueError("positive finite masses and a nonempty positive ratio vector required")
    rho = np.zeros_like(s) if correlations is None else np.asarray(correlations, dtype=float)
    w = np.ones_like(s) if weights is None else np.asarray(weights, dtype=float)
    if (rho.shape != s.shape or w.shape != s.shape or not np.isfinite(rho).all()
            or not np.isfinite(w).all() or np.any(np.abs(rho) >= 1) or np.any(w <= 0)):
        raise ValueError("matching correlations in (-1,1) and positive finite weights required")
    # Sum of squares avoids subtracting nearly equal quadratic terms.
    width = s*np.sqrt((1-rho)*(1+rho))
    denominator = (m[..., None]-rho*s)**2+width**2
    return np.sum(w*(width/denominator)**2, axis=-1)


def independent_channel_ratios(force_covariance, acceleration_covariance, subspace=None):
    """Effective scalar ratios from block-independent anisotropic channels.

    Uses restricted precision N.T Sigma_channel^-1 N, not the inverse of a
    projected marginal covariance. Valid for a fixed known ambient model.
    """
    force, acceleration = _covariance(force_covariance), _covariance(acceleration_covariance)
    if force.shape != acceleration.shape:
        raise ValueError("channel covariances must have the same dimension")
    n = _subspace(subspace, len(force))
    a = n.T @ np.linalg.solve(force, n)
    c = n.T @ np.linalg.solve(acceleration, n)
    return 1/np.sqrt(eigh(a, c, eigvals_only=True))


def _covariance(value):
    value = np.asarray(value, dtype=float)
    if (value.ndim != 2 or value.shape[0] == 0 or value.shape[0] != value.shape[1]
            or not np.isfinite(value).all()
            or not np.allclose(value, value.T, rtol=1e-12, atol=0)):
        raise ValueError("a finite symmetric square covariance is required")
    try:
        np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError("a positive-definite covariance is required") from error
    return value


def _subspace(value, dimension):
    n = np.eye(dimension) if value is None else np.asarray(value, dtype=float)
    if (n.ndim != 2 or n.shape[0] != dimension or n.shape[1] == 0
            or not np.isfinite(n).all() or np.linalg.matrix_rank(n) != n.shape[1]):
        raise ValueError("a finite nonempty full-column-rank subspace basis is required")
    return n

"""Research candidate law from the exact conditional Gaussian vector score.

The law is a supplied joint prior times exp(-T/2), T=S'C^+S. Exact
candidate tests do not force this conversion to a law or prove its sampling
calibration. In particular false stationary score roots remain present.
No production estimator, default rank gate or signal population is selected.
Finite atoms can charge the explicit exceptional rank-zero convention.
"""
import math

import numpy as np
from scipy import special, stats

from met.joint_mass import (FiniteJointMassLaw, gaussian_joint_geometry,
                            scalar_composite_law)


def _masses(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not values.size or np.any(~np.isfinite(values))
            or np.any(values <= 0)):
        raise ValueError("masses must be a nonempty finite positive vector")
    return values


def _rank_tolerances(relative, absolute):
    relative, absolute = float(relative), float(absolute)
    if (not math.isfinite(relative) or not math.isfinite(absolute)
            or relative < 0 or absolute < 0):
        raise ValueError("rank tolerances must be finite and nonnegative")
    return relative, absolute


def gaussian_joint_score_statistic(readings, masses, *, rank_rtol, rank_atol):
    """Conditional vector score and stable SVD statistic at one candidate.

    Stack G_i=Qperp_i' W_i(dC_i/dm) xhat_i and residual coordinates v_i.
    Then S=G'v, C=G'G and T is the residual's squared projection onto G's
    retained left singular space. Retain singular values strictly above
    max(rank_atol,rank_rtol*s_max). Rank zero gives T=0 and p-value=1.

    Mathematical exactness uses exact rank. Explicit tolerances resolve
    floating-point rank and can affect invariance at unresolved singular
    values; they are not a selected statistical signal threshold. The
    candidate fit uses only each event's own observations, never true signals.
    """
    masses = _masses(masses)
    relative, absolute = _rank_tolerances(rank_rtol, rank_atol)
    maps, residuals = [], []
    for reading in readings:
        geometry = gaussian_joint_geometry(reading, masses)
        if not np.any(geometry["whitened_data"]):
            # An exactly zero event contributes no score-map or residual
            # coordinate. Omitting its zero rows also preserves numerical Z.
            continue
        # Project the declared directions before multiplying fitted sizes.
        # This retains their common small tangent factor across readings,
        # instead of repeating a cancellation at each observed amplitude.
        directions = np.array([geometry["perpendicular"].T@item
                               for item in geometry["whitened_derivatives"]])
        maps.append(np.einsum("arq,q->ra", directions, geometry["fitted_latent"]))
        residuals.append(geometry["residual"])
    design = np.vstack(maps) if maps else np.empty((0, len(masses)))
    residual = np.concatenate(residuals) if residuals else np.empty(0)
    if np.any(~np.isfinite(design)) or np.any(~np.isfinite(residual)):
        raise ValueError("conditional score map exceeds the floating-point range")
    left, singular, _ = np.linalg.svd(design, full_matrices=False)
    cutoff = max(absolute, relative*singular[0]) if singular.size else absolute
    keep = singular > cutoff
    rank = int(np.sum(keep))
    coordinates = left[:, keep].T@residual
    with np.errstate(over="ignore", invalid="ignore"):
        statistic = float(coordinates@coordinates)
        score, covariance = design.T@residual, design.T@design
    if (not math.isfinite(statistic) or np.any(~np.isfinite(score))
            or np.any(~np.isfinite(covariance))):
        raise ValueError("score statistic exceeds the floating-point range")
    return dict(score=score, covariance=covariance, statistic=statistic,
                rank=rank, p_value=float(stats.chi2.sf(statistic, rank)) if rank else 1.,
                degenerate=rank == 0, rank_cutoff=float(cutoff), singular_values=singular)


def gaussian_joint_score_posterior(readings, prior, *, rank_rtol, rank_atol):
    """Normalize a supplied finite joint prior times exp(-T/2).

    The supplied joint dependence is retained. This proper finite research
    law is not a certified approximation to a continuous prior, a confidence
    distribution, or an adopted estimator. Log evidence is retained, and a
    common large statistic is removed before adding log prior weights.
    """
    if not isinstance(prior, FiniteJointMassLaw):
        raise ValueError("prior must be a FiniteJointMassLaw")
    relative, absolute = _rank_tolerances(rank_rtol, rank_atol)
    readings = tuple(readings)
    statistics = np.array([gaussian_joint_score_statistic(readings, masses,
                            rank_rtol=relative, rank_atol=absolute)["statistic"]
                           for masses in prior.atoms])
    minimum = float(np.min(statistics))
    log_weights = np.log(prior.probabilities)-(statistics-minimum)/2
    shifted_evidence = float(special.logsumexp(log_weights))
    return FiniteJointMassLaw(prior.atoms, np.exp(log_weights-shifted_evidence),
                             log_evidence=shifted_evidence-minimum/2)


def composite_quantile(law, probability, *, operation):
    """Finite generalized-inverse quantiles of the same joint pushforward.

    Probability is explicit; no median or interval level is selected. This
    is an atomic-law quantile, not a continuous-grid accuracy claim.
    """
    probability = np.asarray(probability, dtype=float)
    if (np.any(~np.isfinite(probability)) or np.any(probability <= 0)
            or np.any(probability >= 1)):
        raise ValueError("probability must lie strictly between zero and one")
    scalar = scalar_composite_law(law, operation=operation)
    cumulative = np.cumsum(scalar.probabilities)
    cumulative[-1] = 1.
    result = scalar.support[np.searchsorted(cumulative, probability, side="left")]
    return float(result) if result.ndim == 0 else result


def scalar_common_score_from_scatter(scatter, mass, *, covariance, rank_rtol, rank_atol):
    """Fast scalar/common-covariance score, vectorized over positive masses.

    Scatter is sum Y_iY_i' in physical (F,a) coordinates and covariance is
    the common declared 2x2 source covariance. No theoretical signal or count
    is supplied. In whitened coordinates let u be the candidate unit line,
    v its perpendicular, d its signed speed, A=u'Su and B=u'Sv. Then score
    is dB, variance d²A and T=B²/A on retained rank one. Exact rank zero is
    accepted. Scatter-only floating arithmetic does not certify every exact
    cancellation in singular input scatter; the generic own-event helper
    has the additional physical normal-equation safeguard.
    """
    relative, absolute = _rank_tolerances(rank_rtol, rank_atol)
    mass = np.asarray(mass, dtype=float)
    scatter, covariance = np.asarray(scatter, dtype=float), np.asarray(covariance, dtype=float)
    if (np.any(~np.isfinite(mass)) or np.any(mass <= 0)
            or scatter.shape != (2, 2) or covariance.shape != (2, 2)
            or np.any(~np.isfinite(scatter)) or np.any(~np.isfinite(covariance))
            or not np.allclose(scatter, scatter.T, rtol=1e-12, atol=0.)
            or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
        raise ValueError("positive finite mass, symmetric scatter and 2x2 covariance are required")
    scatter = scatter/2+scatter.T/2
    covariance = covariance/2+covariance.T/2
    try:
        chol = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError("covariance must be positive definite") from error
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    if np.any(eigenvalues < 0):
        raise ValueError("scatter must be positive semidefinite")
    # A square-root factor keeps A nonnegative without subtracting large
    # nearly equal quadratic-form terms near a low-rank scatter direction.
    factor = np.linalg.solve(chol, eigenvectors*np.sqrt(eigenvalues))
    basis = np.linalg.solve(chol, np.vstack([mass.ravel(), np.ones(mass.size)]))
    length = np.hypot(basis[0], basis[1])
    along_unit = basis/length
    across_unit = np.vstack([-along_unit[1], along_unit[0]])
    projected = along_unit.T@factor
    transverse = across_unit.T@factor
    along = np.einsum("ij,ij->i", projected, projected)
    cross = np.einsum("ij,ij->i", projected, transverse)
    # det(W[m,1],W[1,0])=-det(W). This avoids subtracting two O(1)
    # terms when a correlated covariance and a large mass make d tiny.
    speed = -np.exp(-np.log(np.diag(chol)).sum()-2*np.log(length))
    singular = np.abs(speed)*np.sqrt(along)
    cutoff = np.maximum(absolute, relative*singular)
    keep = singular > cutoff
    normalized = np.divide(cross, np.sqrt(along), out=np.zeros_like(cross), where=along > 0)
    with np.errstate(over="ignore", invalid="ignore"):
        statistic = np.where(keep, normalized*normalized, 0.)
        score, variance = speed*cross, singular*singular
    if (np.any(~np.isfinite(statistic)) or np.any(~np.isfinite(score))
            or np.any(~np.isfinite(variance))):
        raise ValueError("scatter score exceeds the floating-point range")
    p_value = np.where(keep, stats.chi2.sf(statistic, 1), 1.)

    def result(values):
        values = np.asarray(values).reshape(mass.shape)
        return values.item() if values.ndim == 0 else values

    return dict(score=result(score), variance=result(variance), statistic=result(statistic),
                rank=result(keep.astype(int)), p_value=result(p_value),
                degenerate=result(~keep), rank_cutoff=result(cutoff))

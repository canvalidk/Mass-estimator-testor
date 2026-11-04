"""Conditional Gaussian tests of fixed linear operators, a T50 reference.

The model is F_i = M a_i plus declared Gaussian observation error, with an
arbitrary fixed free vector a_i for each reading. No signal population is
estimated. This is a matrix-candidate test, not an extension of the production
mass output type. Only full reading subspaces are supported.

Tangent directions and numerical singular-value rank tolerances are explicit.
The latter resolve floating-point rank; they are not statistical signal gates.
The exact theorem concerns mathematical rank, including singular covariance.
"""
import math

import numpy as np
from scipy import stats

from met.full import Reading


def _matrix(value):
    value = np.asarray(value, dtype=float)
    if (value.ndim != 2 or value.shape[0] != value.shape[1]
            or value.shape[0] == 0 or np.any(~np.isfinite(value))):
        raise ValueError("candidate must be a finite nonempty square matrix")
    return value


def operator_tangent_basis(dimension, *, symmetric):
    """An explicit Frobenius-orthonormal full or symmetric matrix basis.

    Selecting the symmetric family tests its joint score at a symmetric
    candidate. No symmetry or positive-definiteness premise is imposed on an
    unrestricted candidate by this helper.
    """
    if (not isinstance(dimension, (int, np.integer))
            or isinstance(dimension, (bool, np.bool_)) or dimension < 1):
        raise ValueError("dimension must be a positive integer")
    if not isinstance(symmetric, (bool, np.bool_)):
        raise ValueError("symmetric must be an explicit boolean")
    basis = []
    for row in range(dimension):
        for col in range(row if symmetric else 0, dimension):
            direction = np.zeros((dimension, dimension))
            if symmetric and row != col:
                direction[row, col] = direction[col, row] = 1/math.sqrt(2.)
            else:
                direction[row, col] = 1.
            basis.append(direction)
    return np.asarray(basis)


def _geometry(reading, matrix):
    if not isinstance(reading, Reading):
        raise ValueError("each reading must be a met.full.Reading")
    d = matrix.shape[0]
    if reading.dimension != d or reading.subspace.shape[1] != d:
        raise ValueError("operator readings must have the candidate dimension and full subspace")
    data = np.r_[reading.force, reading.acceleration]
    if np.any(~np.isfinite(data)) or np.any(~np.isfinite(reading.cov)):
        raise ValueError("reading data and covariance must be finite")
    chol = np.linalg.cholesky(reading.cov)
    whitened = np.linalg.solve(chol, data)
    basis = np.linalg.solve(chol, np.vstack([matrix, np.eye(d)]))
    rotation, triangular = np.linalg.qr(basis, mode="complete")
    parallel, perpendicular = rotation[:, :d], rotation[:, d:]
    coefficients = parallel.T @ whitened
    # The same optional exact-floating-zero safeguard as T45. This is not a
    # general exact-arithmetic cancellation certificate for ill-conditioned
    # declarations; the stable whitened QR remains the ordinary path.
    try:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            precision_data = np.linalg.solve(reading.cov, data)
            normal_rhs = matrix.T @ precision_data[:d]+precision_data[d:]
        if np.all(np.isfinite(normal_rhs)) and np.all(normal_rhs == 0.):
            coefficients = np.zeros_like(coefficients)
    except np.linalg.LinAlgError:
        pass
    fitted_latent = np.linalg.solve(triangular[:d], coefficients)
    residual = perpendicular.T @ whitened
    return chol, perpendicular, fitted_latent, residual


def gaussian_operator_profile(readings, matrix):
    """Gaussian profile log factor, omitting matrix-independent constants.

    This is -sum ||(I-P_i) Sigma_i^-1/2 Y_i||**2 / 2 after maximizing each
    reading's free amplitude. It integrates neither amplitudes nor a matrix
    prior. In common isotropic noise it is a Grassmann Bingham kernel up to
    a matrix-independent constant, with the usual TLS principal-plane fit.
    """
    matrix = _matrix(matrix)
    costs = []
    for reading in readings:
        _, _, _, residual = _geometry(reading, matrix)
        costs.append(float(residual @ residual))
    result = -.5*math.fsum(costs)
    if not math.isfinite(result):
        raise ValueError("profile exceeds the floating-point range")
    return result


def gaussian_operator_score_test(readings, matrix, *, tangents,
                                 rank_rtol, rank_atol):
    """Exact conditional joint score test of a fixed graph operator.

    Given each reading's own candidate projection U_i, the score vector is
    S_a=sum V_i' G_ia and C_ab=sum <G_ia,G_ib>, where
    G_ia=(I-P_i) Sigma_i^-1/2 [E_a;0] B_i^+ U_i.
    At the true fixed candidate, S|U is N(0,C). The statistic S'C^+S has a
    chi-square law with the observed rank of C, not the number of tangents.

    Numerically, retain singular values of the stacked conditional score map
    strictly above max(rank_atol, rank_rtol*s_max). The statistic is computed
    from its left singular vectors, avoiding a second conditioning loss in
    forming/inverting C. Rank zero has statistic 0 and p-value 1. This null
    convention is total but is not claimed continuous at singular data.
    """
    matrix = _matrix(matrix)
    d = matrix.shape[0]
    tangents = np.asarray(tangents, dtype=float)
    if (tangents.ndim != 3 or tangents.shape[1:] != (d, d)
            or tangents.shape[0] == 0 or np.any(~np.isfinite(tangents))):
        raise ValueError("tangents must be a nonempty finite array of d by d matrices")
    rank_rtol, rank_atol = float(rank_rtol), float(rank_atol)
    if (not math.isfinite(rank_rtol) or not math.isfinite(rank_atol)
            or rank_rtol < 0 or rank_atol < 0):
        raise ValueError("rank tolerances must be finite and nonnegative")
    maps, residuals = [], []
    for reading in readings:
        chol, perpendicular, latent, residual = _geometry(reading, matrix)
        derivatives = np.einsum("aij,j->ia", tangents, latent)
        derivatives = np.vstack([derivatives, np.zeros_like(derivatives)])
        maps.append(perpendicular.T @ np.linalg.solve(chol, derivatives))
        residuals.append(residual)
    p = tangents.shape[0]
    linear_map = np.vstack(maps) if maps else np.empty((0, p))
    residual = np.concatenate(residuals) if residuals else np.empty(0)
    if np.any(~np.isfinite(linear_map)) or np.any(~np.isfinite(residual)):
        raise ValueError("conditional score exceeds the floating-point range")
    left, singular, _ = np.linalg.svd(linear_map, full_matrices=False)
    cutoff = max(rank_atol, rank_rtol*singular[0]) if singular.size else rank_atol
    keep = singular > cutoff
    rank = int(np.sum(keep))
    coordinates = left[:, keep].T @ residual
    with np.errstate(over="ignore", invalid="ignore"):
        statistic = float(coordinates @ coordinates)
        score = linear_map.T @ residual
        covariance = linear_map.T @ linear_map
    if (not math.isfinite(statistic) or np.any(~np.isfinite(score))
            or np.any(~np.isfinite(covariance))):
        raise ValueError("score, covariance or statistic exceeds the floating-point range")
    return dict(score=score, covariance=covariance, statistic=statistic,
                rank=rank, p_value=float(stats.chi2.sf(statistic, rank)) if rank else 1.,
                degenerate=rank == 0, rank_cutoff=float(cutoff),
                singular_values=singular)

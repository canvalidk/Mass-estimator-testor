"""Geometry of compatible corrections, independently of an estimator's mass law.

An observation is any object exposing force, acceleration, covariance and subspace
arrays. Corrections are observed minus fitted compatible values. They are minimum
cost corrections, not recovered instrument errors.

The Newton residual alone imposes F*=m a*. The declared fit additionally imposes
a* in span(N), which can cost more. Its mass-independent outside-subspace cost is
kept separately from the cost of codirectionality inside the declared span.
"""
import numpy as np


def profile(observation, masses):
    """Return geometric diagnostics at finite positive candidate masses.

    Fit the declared candidates by Cholesky whitening and QR least squares, without
    using the estimator's normal equations. A second projection onto diag(N,N)
    separates the cost of the subspace declaration from its mass-dependent cost.
    Costs are computed from residual vectors, avoiding subtraction of large energies.

    Scalar energies and costs use the declared covariance. Per-mass arrays have
    first axis equal to len(masses); conditional means/covariances are coordinates
    of the compatible acceleration in the supplied orthonormal basis N.
    """
    masses = np.asarray(masses, dtype=float)
    if (masses.ndim != 1 or not masses.size
            or not np.all(np.isfinite(masses) & (masses > 0))):
        raise ValueError("masses must be a nonempty 1D array of finite positive values")

    force = np.asarray(observation.force, dtype=float)
    acceleration = np.asarray(observation.acceleration, dtype=float)
    if (force.ndim != 1 or not force.size or acceleration.shape != force.shape
            or not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration))):
        raise ValueError("force and acceleration must be finite vectors of the same length")
    dimension = force.size
    covariance = np.asarray(observation.covariance, dtype=float)
    if (covariance.shape != (2 * dimension, 2 * dimension)
            or not np.all(np.isfinite(covariance))
            or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0)):
        raise ValueError("covariance must be a finite symmetric 2d by 2d matrix")
    try:
        cholesky = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError("covariance must be positive definite") from error
    subspace = (np.eye(dimension) if observation.subspace is None
                else np.asarray(observation.subspace, dtype=float))
    if (subspace.ndim != 2 or subspace.shape[0] != dimension
            or not 1 <= subspace.shape[1] <= dimension
            or not np.all(np.isfinite(subspace))
            or not np.allclose(subspace.T @ subspace, np.eye(subspace.shape[1]),
                               rtol=1e-10, atol=1e-12)):
        raise ValueError("subspace must have 1 to d orthonormal columns")

    y = np.concatenate((force, acceleration))
    whitened_y = np.linalg.solve(cholesky, y)
    zero = np.zeros_like(subspace)
    span = np.block([[subspace, zero], [zero, subspace]])
    span_q, _ = np.linalg.qr(np.linalg.solve(cholesky, span), mode="reduced")
    span_y = span_q @ (span_q.T @ whitened_y)
    outside_residual = whitened_y - span_y
    count, k = masses.size, subspace.shape[1]

    result = dict(
        residual=np.empty((count, dimension)),
        residual_covariance=np.empty((count, dimension, dimension)),
        residual_cost=np.empty(count),
        declared_cost=np.empty(count),
        within_subspace_cost=np.empty(count),
        candidate_force=np.empty((count, dimension)),
        candidate_acceleration=np.empty((count, dimension)),
        force_correction=np.empty((count, dimension)),
        acceleration_correction=np.empty((count, dimension)),
        conditional_mean=np.empty((count, k)),
        conditional_covariance=np.empty((count, k, k)),
        explained_energy=np.empty(count),
        data_energy=float(whitened_y @ whitened_y),
        span_energy=float(span_y @ span_y),
        outside_subspace_cost=float(outside_residual @ outside_residual),
    )
    eye = np.eye(dimension)
    for index, mass in enumerate(masses):
        constraint = np.concatenate((eye, -mass * eye), axis=1)
        residual = force - mass * acceleration
        residual_covariance = constraint @ covariance @ constraint.T
        result["residual"][index] = residual
        result["residual_covariance"][index] = residual_covariance
        result["residual_cost"][index] = residual @ np.linalg.solve(residual_covariance, residual)

        # Scaling all candidate columns leaves their span unchanged and avoids a
        # large mass needlessly worsening the linear solve's overall scale.
        scale = max(1.0, mass)
        candidate = np.concatenate(((mass / scale) * subspace, subspace / scale), axis=0)
        whitened_candidate = np.linalg.solve(cholesky, candidate)
        candidate_q, candidate_r = np.linalg.qr(whitened_candidate, mode="reduced")
        coordinates = np.linalg.solve(candidate_r, candidate_q.T @ whitened_y)
        fitted = candidate @ coordinates
        whitened_fitted = candidate_q @ (candidate_q.T @ whitened_y)
        fit_residual = whitened_y - whitened_fitted
        within_residual = span_y - whitened_fitted
        inverse_r = np.linalg.solve(candidate_r, np.eye(k)) / scale

        result["declared_cost"][index] = fit_residual @ fit_residual
        result["within_subspace_cost"][index] = within_residual @ within_residual
        result["candidate_force"][index] = fitted[:dimension]
        result["candidate_acceleration"][index] = fitted[dimension:]
        result["force_correction"][index] = force - fitted[:dimension]
        result["acceleration_correction"][index] = acceleration - fitted[dimension:]
        result["conditional_mean"][index] = coordinates / scale
        result["conditional_covariance"][index] = inverse_r @ inverse_r.T
        result["explained_energy"][index] = whitened_fitted @ whitened_fitted
    return result

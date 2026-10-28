"""O16 geometry of common sensor tilt; diagnostics do not alter the mass estimator.

Vectors remain in one physical frame. The local tilt fit has noisy regressors and
is not assigned a chi-square p-value or used as an instrument correction.
"""
import math

import numpy as np


def _vector(value, name):
    value = np.asarray(value, dtype=float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be a finite length-three vector")
    return value


def _covariance(value, dimension, name):
    value = np.asarray(value, dtype=float)
    if (value.shape != (dimension, dimension) or not np.all(np.isfinite(value))
            or not np.allclose(value, value.T, rtol=1e-12, atol=0)):
        raise ValueError(f"{name} must be a finite symmetric covariance")
    try:
        np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be positive definite") from error
    return value


def skew(vector):
    """The matrix satisfying skew(a) b = a cross b."""
    x, y, z = _vector(vector, "vector")
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotation(axis, angle_radians):
    """Right-handed Rodrigues rotation about a nonzero axis."""
    axis = _vector(axis, "axis")
    norm = np.linalg.norm(axis)
    if (norm == 0 or isinstance(angle_radians, bool) or not np.isscalar(angle_radians)
            or not np.isfinite(angle_radians)):
        raise ValueError("rotation needs a nonzero axis and a finite angle")
    matrix = skew(axis / norm)
    return np.eye(3) + math.sin(angle_radians) * matrix + (1 - math.cos(angle_radians)) * matrix @ matrix


def _singular_values(matrix):
    values = np.linalg.svd(matrix, compute_uv=False)
    return np.pad(values, (0, 3 - len(values)))


def excitation_spectrum(acceleration, covariance, alpha=0.05):
    """Physical excitation spectra and conservative conditional rank lower bounds.

    E||error||_F^2=sum trace(C_i), even with cross-reading dependence. Markov
    and Weyl bound every singular-value error by sqrt(sum trace(C_i)/alpha)
    with probability at least 1-alpha under correct zero-mean declarations.
    Centering is a contraction, so the same radius is valid for the centered
    spectrum without assuming independent errors. The bound is not a rank estimate.
    """
    acceleration = np.asarray(acceleration, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if (acceleration.ndim != 2 or acceleration.shape[1] != 3 or not len(acceleration)
            or not np.all(np.isfinite(acceleration))):
        raise ValueError("acceleration must be a nonempty finite N by 3 array")
    if covariance.shape != (len(acceleration), 3, 3):
        raise ValueError("covariance must contain one 3 by 3 matrix per reading")
    for value in covariance:
        _covariance(value, 3, "acceleration covariance")
    if isinstance(alpha, bool) or not np.isscalar(alpha) or not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be a finite number strictly between zero and one")
    raw = _singular_values(acceleration)
    centered = _singular_values(acceleration - acceleration.mean(axis=0))
    gram = acceleration.T @ acceleration
    noise = covariance.sum(axis=0)
    expected_squared = float(np.trace(noise))
    radius = math.sqrt(expected_squared / alpha)
    return dict(
        singular_values=raw, centered_singular_values=centered,
        gram_eigenvalues=np.linalg.eigvalsh(gram)[::-1],
        noise_subtracted_gram_eigenvalues=np.linalg.eigvalsh(gram - noise)[::-1],
        noise_frobenius_expected_squared=expected_squared, alpha=float(alpha),
        noise_radius=radius, singular_value_lower=np.maximum(0, raw - radius),
        singular_value_upper=raw + radius,
        centered_singular_value_lower=np.maximum(0, centered - radius),
        centered_singular_value_upper=centered + radius,
        certified_uncentered_rank_lower_bound=int(np.count_nonzero(raw > radius)),
        certified_centered_rank_lower_bound=int(np.count_nonzero(centered > radius)),
        assumptions=["Conditional on correct zero-mean covariance declarations; shared biases can violate them.",
                     "Physical common-frame spectra, not separately whitened row directions.",
                     "Noise-subtracted eigenvalues are signed moment diagnostics, not recovered latent energies."])


def shared_tilt_diagnostic(readings, reference_mass):
    """Fit r_i=delta_mass*a_i - m*skew(a_i)*omega in declared residual units.

    Each covariance is propagated through [I,-mI], including channel correlation.
    Rotation observability is measured after projecting out the scalar mass column.
    A relative gain is exactly confounded with mass; this fit does not identify it.
    """
    readings = list(readings)
    if (not readings or isinstance(reference_mass, bool) or not np.isscalar(reference_mass)
            or not np.isfinite(reference_mass) or reference_mass <= 0):
        raise ValueError("tilt diagnostic needs readings and a finite positive reference mass")
    m = float(reference_mass)
    residuals, covariances, designs, whitened_designs, whitened_residuals = [], [], [], [], []
    constraint = np.concatenate((np.eye(3), -m * np.eye(3)), axis=1)
    for reading in readings:
        a, f = _vector(reading.acceleration, "acceleration"), _vector(reading.force, "force")
        covariance = _covariance(reading.cov, 6, "reading covariance")
        residual, design = f - m * a, np.column_stack((a, -m * skew(a)))
        propagated = constraint @ covariance @ constraint.T
        chol = np.linalg.cholesky(propagated)
        residuals.append(residual)
        covariances.append(propagated)
        designs.append(design)
        whitened_designs.append(np.linalg.solve(chol, design))
        whitened_residuals.append(np.linalg.solve(chol, residual))
    design = np.vstack(whitened_designs)
    target = np.concatenate(whitened_residuals)
    beta, _, fit_rank, _ = np.linalg.lstsq(design, target, rcond=None)
    fitted = design @ beta
    mass_column, rotation_columns = design[:, 0], design[:, 1:]
    mass_norm = float(mass_column @ mass_column)
    mass_fit = mass_column * (mass_column @ target / mass_norm) if mass_norm > 0 else 0 * target
    projected = rotation_columns.copy()
    if mass_norm > 0:
        projected -= np.outer(mass_column, mass_column @ rotation_columns / mass_norm)
    _, values, vh = np.linalg.svd(projected, full_matrices=False)
    tolerance = float(np.finfo(float).eps * max(projected.shape) * (values[0] if len(values) else 0))
    rotation_rank = int(np.count_nonzero(values > tolerance))
    raw_cost, fitted_cost = float(target @ target), float((target - fitted) @ (target - fitted))
    return dict(delta_mass=float(beta[0]), rotation_vector=beta[1:], fit_rank=int(fit_rank),
                rotation_singular_values=values, rotation_rank=rotation_rank,
                rotation_rank_tolerance=tolerance,
                rotation_projector=vh[:rotation_rank].T @ vh[:rotation_rank],
                rotation_nullspace=vh[rotation_rank:].T,
                residual=np.asarray(residuals), residual_covariance=np.asarray(covariances),
                predicted_residual=np.einsum("nij,j->ni", np.asarray(designs), beta),
                per_reading_raw_cost=np.sum(target.reshape(-1, 3) ** 2, axis=1),
                per_reading_fitted_cost=np.sum((target - fitted).reshape(-1, 3) ** 2, axis=1),
                raw_cost=raw_cost, fitted_cost=fitted_cost,
                mass_only_cost=float((target - mass_fit) @ (target - mass_fit)),
                rotation_cost_reduction=float((target - mass_fit) @ (target - mass_fit) - fitted_cost),
                assumptions=["Local small-angle diagnostic about the supplied reference mass, with noisy acceleration regressors.",
                             "Costs and fitted parameters have no assigned p-values or causal interpretation.",
                             "Numerical rotation rank describes this observed design; noisy rows can lift a true null direction.",
                             "No instrument correction or estimator input is changed; relative gain and mass remain confounded."])


def score_tilt(readings, reference_mass, delta_mass, rotation_vector):
    """Score supplied tilt parameters without fitting them on these readings.

    This supports fit/score splits. The held-out covariance remains the declared
    [I,-mI] covariance, and fitted-parameter uncertainty is not added; improvements
    are predictive descriptive costs, not p-values. A bad prediction can increase cost.
    """
    readings = list(readings)
    if (not readings or isinstance(reference_mass, bool) or not np.isscalar(reference_mass)
            or not np.isfinite(reference_mass) or reference_mass <= 0):
        raise ValueError("tilt scoring needs readings and a finite positive reference mass")
    if isinstance(delta_mass, bool) or not np.isscalar(delta_mass) or not np.isfinite(delta_mass):
        raise ValueError("delta_mass must be finite")
    omega = _vector(rotation_vector, "rotation_vector")
    m = float(reference_mass)
    constraint = np.concatenate((np.eye(3), -m * np.eye(3)), axis=1)
    residuals, predicted, before, after = [], [], [], []
    for reading in readings:
        a, f = _vector(reading.acceleration, "acceleration"), _vector(reading.force, "force")
        covariance = _covariance(reading.cov, 6, "reading covariance")
        residual = f - m * a
        prediction = delta_mass * a + m * np.cross(omega, a)
        chol = np.linalg.cholesky(constraint @ covariance @ constraint.T)
        whitened = np.linalg.solve(chol, residual)
        corrected = np.linalg.solve(chol, residual - prediction)
        residuals.append(residual)
        predicted.append(prediction)
        before.append(float(whitened @ whitened))
        after.append(float(corrected @ corrected))
    return dict(cost_before=float(sum(before)), cost_after=float(sum(after)),
                improvement=float(sum(before) - sum(after)),
                residual=np.asarray(residuals), predicted_residual=np.asarray(predicted),
                per_reading_cost_before=np.asarray(before), per_reading_cost_after=np.asarray(after),
                assumptions=["Parameters supplied from a separate fit; this function fits nothing.",
                             "Prediction uses noisy acceleration and a local tilt approximation.",
                             "Declared residual covariance excludes fitted-parameter uncertainty; no p-value is assigned."])

"""Independent O17 diagnostics for a declared tensor-mass world.

Tensor truth belongs to the investigator. These functions neither choose a mass
readout nor alter ours29. All measured-data corrections use full 3D candidate space;
declared restricted subspaces require different constrained geometry.
"""
import math

import numpy as np
from scipy.stats import chi2


def validate_tensor(tensor):
    """Return a finite symmetric positive-definite 3-by-3 tensor."""
    tensor = np.asarray(tensor, dtype=float)
    if (tensor.shape != (3, 3) or not np.all(np.isfinite(tensor))
            or not np.allclose(tensor, tensor.T, rtol=1e-12, atol=0)):
        raise ValueError("tensor must be a finite symmetric 3 by 3 matrix")
    tensor = (tensor + tensor.T) / 2
    try:
        np.linalg.cholesky(tensor)
    except np.linalg.LinAlgError as error:
        raise ValueError("tensor must be positive definite") from error
    return tensor


def _readings(readings):
    readings = list(readings)
    if not readings:
        raise ValueError("geometry requires at least one reading")
    output = []
    for reading in readings:
        f = np.asarray(reading.force, dtype=float)
        a = np.asarray(reading.acceleration, dtype=float)
        covariance = np.asarray(reading.cov, dtype=float)
        if (f.shape != (3,) or a.shape != (3,) or not np.all(np.isfinite(f))
                or not np.all(np.isfinite(a))):
            raise ValueError("readings must have finite three-dimensional force and acceleration")
        if (covariance.shape != (6, 6) or not np.all(np.isfinite(covariance))
                or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0)):
            raise ValueError("reading covariance must be finite symmetric 6 by 6")
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError("reading covariance must be positive definite") from error
        subspace = getattr(reading, "subspace", None)
        if subspace is not None:
            subspace = np.asarray(subspace, dtype=float)
            if (subspace.shape != (3, 3) or not np.all(np.isfinite(subspace))
                    or not np.allclose(subspace.T @ subspace, np.eye(3), atol=1e-10)):
                raise ValueError("tensor geometry requires the full 3D candidate subspace")
        output.append((f, a, covariance))
    return output


def scalar_representability(tensor, latent_acceleration, tolerance=1e-10):
    """Test whether one scalar acts like G on the supplied latent excitation span.

    F*=G a* admits one scalar iff (G-mI) annihilates the span of all a*. For SPD
    G, a nonzero such span lies in one eigenspace. Decisions here use the stated
    relative Frobenius tolerance, not an inference from noisy measured vectors.
    The physical least-squares scalar is a diagnostic, not a true scalar mass when
    the result is incompatible. All-zero excitation admits every positive scalar.
    """
    tensor = validate_tensor(tensor)
    acceleration = np.asarray(latent_acceleration, dtype=float)
    if (acceleration.ndim != 2 or acceleration.shape[1] != 3 or not len(acceleration)
            or not np.all(np.isfinite(acceleration))):
        raise ValueError("latent_acceleration must be a nonempty finite N by 3 array")
    if (isinstance(tolerance, bool) or not np.isscalar(tolerance)
            or not np.isfinite(tolerance) or tolerance <= 0):
        raise ValueError("tolerance must be finite and positive")
    force = acceleration @ tensor.T
    eigenvalues = np.linalg.eigvalsh(tensor)
    singular = np.linalg.svd(acceleration, compute_uv=False)
    singular = np.pad(singular, (0, max(0, 3 - len(singular))))
    energy = float(np.sum(acceleration * acceleration))
    common = dict(tolerance=float(tolerance), tensor_eigenvalues=eigenvalues,
                  excitation_singular_values=singular, latent_force=force)
    if energy == 0:
        return dict(common, status="every", scalar_exact=True, unique_scalar=None,
                    least_squares_scalar=None, residual_norm=0.0, relative_residual=0.0,
                    residual_tolerance=0.0, excitation_rank=0)
    mass = float(np.sum(force * acceleration) / energy)
    residual = float(np.linalg.norm(force - mass * acceleration))
    force_norm = float(np.linalg.norm(force))
    threshold = float(tolerance * force_norm)
    compatible = residual <= threshold
    return dict(common, status="one" if compatible else "none", scalar_exact=bool(compatible),
                unique_scalar=mass if compatible else None, least_squares_scalar=mass,
                residual_norm=residual, relative_residual=residual / force_norm,
                residual_tolerance=threshold,
                excitation_rank=int(np.count_nonzero(singular > tolerance * singular[0])))


def scalar_profile(readings, masses):
    """Closest scalar-compatible pairs from independent Cholesky/QR latent fits.

    Each covariance may be fully correlated, but the candidate subspace is R^3.
    No mass is fitted by this function: it evaluates the supplied positive grid.
    Costs are squared corrections under declarations, without p-values.
    """
    observed = _readings(readings)
    masses = np.asarray(masses, dtype=float)
    if (masses.ndim != 1 or not masses.size
            or not np.all(np.isfinite(masses) & (masses > 0))):
        raise ValueError("masses must be a nonempty vector of finite positive values")
    eye = np.broadcast_to(np.eye(3), (len(masses), 3, 3))
    design = np.concatenate((masses[:, None, None] * eye, eye), axis=1)
    constraint = np.concatenate((eye, -masses[:, None, None] * eye), axis=2)
    rows = []
    for force, acceleration, covariance in observed:
        y = np.r_[force, acceleration]
        chol = np.linalg.cholesky(covariance)
        white_y = np.linalg.solve(chol, y)
        white_design = np.linalg.solve(chol, design)
        q, r = np.linalg.qr(white_design, mode="reduced")
        rhs = np.einsum("nij,i->nj", q, white_y)
        latent = np.linalg.solve(r, rhs[..., None])[..., 0]
        compatible = np.einsum("nij,nj->ni", design, latent)
        correction = y[None, :] - compatible
        white_error = white_y[None, :] - np.einsum("nij,nj->ni", white_design, latent)
        rows.append(dict(
            correction_cost=np.sum(white_error * white_error, axis=1),
            residual=force[None, :] - masses[:, None] * acceleration[None, :],
            residual_covariance=constraint @ covariance @ constraint.transpose(0, 2, 1),
            force_correction=correction[:, :3], acceleration_correction=correction[:, 3:],
            candidate_force=compatible[:, :3], candidate_acceleration=compatible[:, 3:]))
    return dict(mass=masses.copy(), readings=rows,
                total_cost=np.sum([row["correction_cost"] for row in rows], axis=0))


def common_scalar_minimum(readings):
    """Exact common-mass cost infimum for one independent isotropic noise pair.

    Every reading has the same sigma_F, sigma_a and full R^3 subspace. Stacking
    all whitened vectors gives the common Gram problem. Its positive domain has
    limiting endpoints; a boundary infimum never becomes a finite fitted mass.
    """
    observed = _readings(readings)
    covariance = observed[0][2]
    sf2, sa2 = float(covariance[0, 0]), float(covariance[3, 3])
    expected = np.diag(np.r_[np.full(3, sf2), np.full(3, sa2)])
    for _, _, current in observed:
        if not np.allclose(current, expected, rtol=1e-12, atol=0):
            raise ValueError("common minimum requires one shared independent isotropic noise pair")
    x = np.array([f for f, _, _ in observed]) / math.sqrt(sf2)
    y = np.array([a for _, a, _ in observed]) / math.sqrt(sa2)
    a, b, c = float(np.sum(x * x)), float(np.sum(y * y)), float(np.sum(x * y))
    singular = np.linalg.svd(np.column_stack((x.ravel(), y.ravel())), compute_uv=False)
    projective, maximum = float(singular[-1] ** 2), float(singular[0] ** 2)
    individual = float(sum(np.linalg.svd(np.column_stack((f, acc)), compute_uv=False)[-1] ** 2
                           for f, acc in zip(x, y)))
    mass = None
    gap = math.hypot(a - b, 2 * c)
    if c > 0:
        slope = ((a - b + gap) / (2 * c) if a >= b else 2 * c / (b - a + gap))
        mass = math.sqrt(sf2 / sa2) * slope
        minimum, location = projective, "interior"
    elif a < b:
        minimum, location = a, "mass_to_zero"
    elif b < a:
        minimum, location = b, "mass_to_infinity"
    else:
        minimum, location = a, "all_masses" if c == 0 else "both_boundaries"
    supremum = maximum if c < 0 else max(a, b)
    return dict(projective_min_cost=projective, positive_min_cost=float(minimum),
                positive_min_mass=mass, positive_min_location=location,
                positive_cost_supremum=float(supremum), positive_cost_span=(gap + abs(a - b)) / 2,
                sum_individual_projective_min_cost=individual,
                extra_common_mass_cost=max(0.0, projective - individual),
                common_gram=np.array([[a, c], [c, b]]),
                assumptions=["One common positive mass, independently profiled latent vectors, and one isotropic noise pair.",
                             "Infimum costs are diagnostics; no p-value or mass readout is assigned."])


def tensor_known_cost(readings, tensor):
    """Profile the pairs against a supplied known tensor without fitting that tensor.

    For independent Gaussian measured readings generated by this exact known G and
    correct covariance, the total cost is chi-square with 3N degrees of freedom.
    Fitting G from these readings would invalidate that unfitted reference.
    """
    tensor = validate_tensor(tensor)
    constraint = np.concatenate((np.eye(3), -tensor), axis=1)
    residuals, covariances, costs, corrections, compatible = [], [], [], [], []
    observed = _readings(readings)
    for force, acceleration, covariance in observed:
        y = np.r_[force, acceleration]
        residual = constraint @ y
        propagated = constraint @ covariance @ constraint.T
        multiplier = np.linalg.solve(propagated, residual)
        correction = covariance @ constraint.T @ multiplier
        residuals.append(residual)
        covariances.append(propagated)
        costs.append(float(residual @ multiplier))
        corrections.append(correction)
        compatible.append(y - correction)
    corrections, compatible = np.array(corrections), np.array(compatible)
    return dict(total_cost=float(sum(costs)), per_reading_cost=np.array(costs),
                residual=np.array(residuals), residual_covariance=np.array(covariances),
                force_correction=corrections[:, :3], acceleration_correction=corrections[:, 3:],
                candidate_force=compatible[:, :3], candidate_acceleration=compatible[:, 3:],
                residual_degrees_of_freedom=3 * len(observed),
                assumptions=["The tensor is supplied, not fitted; latent compatible pairs are profiled independently.",
                             "The chi-square reference requires correct independent zero-mean Gaussian channel errors."])


def scalar_gaussian_compatibility_bound(readings, alpha=0.01):
    """A conservative investigator diagnostic, never an estimator rejection gate.

    If one true positive scalar m0 exists, C_min <= K(m0), and K(m0) is
    chi-square(3N) under the correct independent Gaussian declarations. Therefore
    sf(C_min,3N) is a conservative p-value, without subtracting fitted degrees of
    freedom. The exact minimum here requires one shared isotropic noise pair.
    """
    readings = list(readings)
    if (isinstance(alpha, bool) or not np.isscalar(alpha)
            or not np.isfinite(alpha) or not 0 < alpha < 1):
        raise ValueError("alpha must be finite and strictly between zero and one")
    minimum = common_scalar_minimum(readings)["positive_min_cost"]
    degrees = 3 * len(readings)
    cutoff = float(chi2.isf(alpha, degrees))
    return dict(positive_common_cost=minimum, reference_degrees_of_freedom=degrees,
                conservative_p=float(chi2.sf(minimum, degrees)), alpha=float(alpha),
                cutoff=cutoff, exceeds_cutoff=bool(minimum > cutoff),
                assumptions=["One true positive scalar, correct independent zero-mean Gaussian channel errors.",
                             "The fitted minimum is bounded by a fixed-true-mass chi-square cost; this is not its exact distribution.",
                             "Per-case conditional diagnostic only; no change to the estimator or physical acceptance premise."])

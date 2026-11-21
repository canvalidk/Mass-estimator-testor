"""Same-record mass laws for a separable Gaussian channel covariance.

pairs[i, :, 0] is force and pairs[i, :, 1] is acceleration, after
whitening the spatial covariance only. The mass endpoints describe the
family of objects being measured, e.g. (b, infinity) for a known load b.
The returned log factor multiplies the given physical mass prior once.
"""
import math

import numpy as np

from explore.e12_residual_join import residual_log_evidence
from met.full import Reading


def standard_sensor_record(data, covariance, standard_mass, *, lump=False, subspace=None):
    """Reduce (unknown force, acceleration, standard force) to a T36 Reading.

    The return value also includes the observed calibration residual, which
    is independent of the returned reading in this declared Gaussian model.
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 1 or data.size not in (3, 6, 9) or not np.isfinite(data).all():
        raise ValueError("expected three finite vector channels of dimension 1, 2 or 3")
    if not math.isfinite(standard_mass) or standard_mass <= 0:
        raise ValueError("the known mass must be finite and positive")
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (data.size, data.size) or not np.isfinite(covariance).all():
        raise ValueError("invalid covariance for the declared channels")
    if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.0):
        raise ValueError("the three-channel covariance must be symmetric")
    try:
        root = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("the three-channel covariance must be positive definite") from exc

    d = data.size // 3
    eye, zero = np.eye(d), np.zeros((d, d))
    defect = np.column_stack((zero, -standard_mass * eye, eye))
    selected = np.block([[eye, zero, eye if lump else zero], [zero, eye, zero]])
    # Remove the components of the whitened noise observed by the calibration defect.
    basis, _ = np.linalg.qr((defect @ root).T, mode="complete")
    complement = basis[:, d:]
    effective_root = (selected @ root) @ complement
    effective = effective_root @ (complement.T @ np.linalg.solve(root, data))
    reading = Reading(effective[:d], effective[d:], effective_root @ effective_root.T, subspace)
    return reading, defect @ data


def _rays(masses):
    # Scaled homogeneous coordinates also represent the infinite endpoint.
    scale = np.maximum(1.0, masses)
    first = np.divide(masses, scale, out=np.ones_like(scale), where=np.isfinite(scale))
    return np.stack((first, 1.0 / scale), axis=-1)


def known_mass_log_evidence(pairs, covariance, masses, *, bounds=(0.0, np.inf), tau=0.0,
                            atol=1e-9, fit_domain="physical"):
    """Fit over the physical arc or all projective lines, then evaluate at physical masses."""
    if fit_domain not in ("physical", "projective"):
        raise ValueError("fit_domain must be physical or projective")
    pairs = np.asarray(pairs, dtype=float)
    if pairs.ndim != 3 or pairs.shape[-1] != 2 or pairs.shape[1] == 0:
        raise ValueError("pairs must have shape (number of readings, spatial dimension, 2)")
    if not np.isfinite(pairs).all():
        raise ValueError("readings must be finite")
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (2, 2) or not np.isfinite(covariance).all():
        raise ValueError("expected a finite 2-by-2 channel covariance")
    if not np.allclose(covariance, covariance.T, rtol=1e-14, atol=0):
        raise ValueError("channel covariance must be symmetric")
    try:
        root = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("channel covariance must be positive definite") from exc

    lower, upper = bounds
    if not 0 <= lower < upper:
        raise ValueError("mass bounds must be ordered and nonnegative")
    masses = np.asarray(masses, dtype=float)
    if np.any(np.isnan(masses)) or np.any((masses < lower) | (masses > upper)):
        raise ValueError("candidate lies outside the physical mass interval")

    endpoints = np.linalg.solve(root, _rays(np.array([lower, upper], dtype=float)).T)
    start = endpoints[:, 0] / math.hypot(*endpoints[:, 0])
    inward = np.array([start[1], -start[0]])
    frame = np.stack((inward, start))
    aperture = math.atan2(inward @ endpoints[:, 1], start @ endpoints[:, 1])
    if not 0 < aperture < math.pi:
        raise ValueError("the declared aperture cannot be resolved in floating point")

    chart = np.linalg.solve(root.T, frame.T).T
    chart_data = pairs @ chart.T
    grams = chart_data.transpose(0, 2, 1) @ chart_data
    directions = _rays(masses) @ chart.T
    angles = np.clip(np.arctan2(directions[..., 0], directions[..., 1]), 0.0, aperture)
    angles = np.where(masses == lower, 0.0, np.where(masses == upper, aperture, angles))
    fit_aperture = aperture if fit_domain == "physical" else math.pi
    return residual_log_evidence(grams, angles, dimension=pairs.shape[1], tau=tau,
                                 atol=atol, aperture=fit_aperture)

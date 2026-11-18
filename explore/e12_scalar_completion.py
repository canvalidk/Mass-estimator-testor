"""Numerical log-evidence for the boundary and equal-eigenvalue construction in E12.

Rows are scalar (force, acceleration) observations whitened by the declared
standard deviations. Angles are atan(m / instrument_ratio). The caller's
angular prior multiplies this factor once. This module is exploration code.
"""
import math

import numpy as np

from explore.e12_scalar_witness import pooled_witness


def completed_log_evidence(rows, angles, *, alpha, count_rate=1.0):
    """Evaluate the joint log factor, up to an angle-independent constant."""
    rows = np.asarray(rows, dtype=float)
    fit = pooled_witness(rows, alpha=alpha, count_rate=count_rate)
    angles = np.asarray(angles, dtype=float)
    if np.any(~np.isfinite(angles)) or np.any((angles < 0) | (angles > math.pi / 2)):
        raise ValueError("angles must be in the closed positive quadrant")
    if fit.largest == 0:
        return np.zeros_like(angles)
    rows = rows[np.any(rows != 0, axis=1)]

    # This cap only matters when one of the two observed channels has very little energy.
    channel_energy = np.einsum("ij,ij->j", rows, rows)
    smaller_channel = float(np.min(channel_energy))
    precision = fit.precision
    if smaller_channel < math.sqrt(precision):
        precision = min(precision, 1.0 + smaller_channel ** 2)
    log_factor = -0.5 * precision * np.sin(angles - fit.angle) ** 2

    spectral_precision = (fit.gap / fit.largest) * fit.gap
    gate = math.exp(-spectral_precision) * (fit.positive_residual / (1.0 + fit.positive_residual))
    if gate == 0:
        return log_factor

    # Normalization bounds the entire extra term by one even for very long series.
    energies = np.einsum("ij,ij->i", rows, rows)
    normalizer = 1.0 + np.sum(energies / (1.0 + energies))
    cosine, sine = np.cos(angles), np.sin(angles)
    shape = np.zeros_like(angles)
    for force, acceleration in rows / np.sqrt(1.0 + energies[:, None]):
        shape += np.log1p((force * cosine - acceleration * sine) ** 2)
    return log_factor - (gate / normalizer) * shape

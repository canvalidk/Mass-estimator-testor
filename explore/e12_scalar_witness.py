"""Pooled scalar comparison in e12-pooled-witness.md (whitened F,a rows).

Every choice of alpha and count_rate specifies a separate mathematical witness.
Only the uniform-angle intervals use the existing T42 continuous quadrature.
"""

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class PooledWitness:
    angle: float
    gap: float
    largest: float
    positive_residual: float
    count: float
    discount: float
    precision: float


def effective_count(row_energies, *, rate=1.0):
    """Sum h_b(t) with t the squared whitened length of one observed row."""
    if not math.isfinite(rate) or rate <= 0:
        raise ValueError("count rate must be finite and positive")
    energies = np.asarray(row_energies, dtype=float)
    if np.any(~np.isfinite(energies)) or np.any(energies < 0):
        raise ValueError("row energies must be finite and nonnegative")
    with np.errstate(over="ignore"):
        exponent = -(0.5 * energies) * rate
    u = np.exp(exponent)
    complement = -np.expm1(exponent)
    # Evaluate (1-u)(1+(2+1/b)u) without forming 1/b at a tiny b.
    return float(np.sum(complement * (1.0 + 2.0 * u) + (complement / rate) * u))


def pooled_witness(rows, *, alpha, count_rate=1.0):
    """Return the fitted angular line, effective count, discount and precision.

    Columns are force / sigma_F and acceleration / sigma_a. The evidence is
    proportional to exp(-precision * sin(theta-angle)**2/2) on positive
    angles. An exact spectral tie has constant evidence for any choice of angle.
    """
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be finite and positive")
    rows = np.asarray(rows, dtype=float)
    if rows.shape == (0,):
        rows = np.empty((0, 2))
    if rows.ndim != 2 or rows.shape[1] != 2 or np.any(~np.isfinite(rows)):
        raise ValueError("rows must be a finite N by 2 array")
    rows = rows[np.any(rows != 0, axis=1)]
    with np.errstate(over="ignore", invalid="ignore"):
        scatter = rows.T @ rows
        energies = np.einsum("ij,ij->i", rows, rows)
    if np.any(~np.isfinite(scatter)) or np.any(~np.isfinite(energies)):
        raise ValueError("whitened scatter exceeds floating-point range")
    count = effective_count(energies, rate=count_rate)
    values, vectors = np.linalg.eigh(scatter)
    largest = float(values[1])
    if not math.isfinite(largest):
        raise ValueError("largest eigenvalue exceeds floating-point range")
    gap = largest - max(0.0, float(values[0]))
    force, accel = vectors[:, 1]
    if accel < 0:
        force, accel = -force, -accel
    angle = math.atan2(force, accel) if gap else 0.0

    # Minimize the residual over the CLOSED positive arc. A negatively
    # sloped principal line has its restricted optimum at one of the axes.
    if scatter[0, 1] <= 0:
        positive_residual = float(min(scatter[0, 0], scatter[1, 1]))
    else:
        errors = rows[:, 0] * accel - rows[:, 1] * force
        positive_residual = min(largest, float(errors @ errors))
    discount = 1.0
    if positive_residual > count and count > 0:
        log_residual = math.log(positive_residual)
        log_weight = (math.log(count) - log_residual
                      + alpha * (math.log(largest) - log_residual))
        discount = math.exp(min(0.0, log_weight))
    precision = discount * gap * (gap / largest) if largest else 0.0
    return PooledWitness(angle, gap, largest, positive_residual, count, discount, precision)

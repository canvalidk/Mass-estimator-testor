"""Four additive real coordinates for the common-instrument E12 construction.

One declared reading contains k spatial components in each whitened channel.
Angles follow the project convention theta = atan(m / instrument_ratio).
"""
from dataclasses import dataclass
import math

import numpy as np


def _residual_count(k, cost):
    scale = k / (1.0 - 0.5 ** (0.5 * k))
    if cost < 1e-16:
        return (0.5 * scale) * cost
    return scale * (-math.expm1(-0.5 * cost))


@dataclass(frozen=True)
class JoinState:
    dimension: int
    force_energy: float
    cross: float
    acceleration_energy: float
    count: float

    @classmethod
    def from_reading(cls, force, acceleration):
        """Retain one reading, after whitening each channel by its standard error."""
        force = np.atleast_1d(np.asarray(force, dtype=float))
        acceleration = np.atleast_1d(np.asarray(acceleration, dtype=float))
        if force.ndim != 1 or force.size == 0 or force.shape != acceleration.shape:
            raise ValueError("both channels must have the same positive spatial dimension")
        ff = float(force @ force)
        aa = float(acceleration @ acceleration)
        if not math.isfinite(ff + aa):
            raise ValueError("whitened channel energy must be finite")
        # At m = sigma_F / sigma_a the unit transverse direction is (1, -1)/sqrt(2).
        defect = (force - acceleration) / math.sqrt(2.0)
        return cls(force.size, ff, float(force @ acceleration), aa,
                   _residual_count(force.size, float(defect @ defect)))

    def merge(self, other):
        """Combine two independent blocks on the same declared instrument."""
        if self.dimension != other.dimension:
            raise ValueError("spatial dimension must agree")
        ff = self.force_energy + other.force_energy
        aa = self.acceleration_energy + other.acceleration_energy
        count = self.count + other.count
        if not math.isfinite(ff + aa) or not math.isfinite(count):
            raise ValueError("combined energy and count must be finite")
        return JoinState(self.dimension, ff, self.cross + other.cross, aa, count)

    def log_evidence(self, angles, *, tau=0.0):
        """Return the angle-dependent log factor for a fixed family member."""
        reference_cost = max(0.0, 0.5 * (self.force_energy + self.acceleration_energy) - self.cross)
        return self._log_evidence(angles, tau, reference_cost)

    def _log_evidence(self, angles, tau, reference_cost, aperture=math.pi / 2):
        """A missing reference cost selects the minimum total residual on the mass arc."""
        angles = np.asarray(angles, dtype=float)
        if np.any(~np.isfinite(angles)) or np.any((angles < 0) | (angles > aperture)):
            raise ValueError("angles must lie in the declared mass arc")
        if not math.isfinite(tau) or tau < 0:
            raise ValueError("tau must be finite and nonnegative")

        ff, fa, aa = self.force_energy, self.cross, self.acceleration_energy
        trace = ff + aa
        if trace == 0:
            return np.zeros_like(angles)
        eigenvalues, eigenvectors = np.linalg.eigh([[ff, fa], [fa, aa]])
        small, large = max(0.0, eigenvalues[0]), eigenvalues[1]
        gap = large - small
        centre = math.atan2(eigenvectors[0, 1], eigenvectors[1, 1]) % math.pi
        # The first coordinate points along increasing angle at the lower endpoint.
        if aperture == math.pi:
            # A fit over all unoriented channel lines has no marked endpoints.
            sine, cosine, end_cost = 0.0, -1.0, ff
        elif aperture == math.pi / 2:
            sine, cosine, end_cost = 1.0, 0.0, aa
        else:
            sine, cosine = math.sin(aperture), math.cos(aperture)
            end_cost = small + gap * math.sin(aperture - centre) ** 2
        endpoint_cost = min(ff, end_cost)
        residual = small if centre <= aperture else endpoint_cost
        discount = 1.0
        if residual > 0 and self.count > 0:
            log_discount = math.log(self.count) - math.log(residual)
            log_discount += (math.log(large) - math.log(residual)) / 8.0
            discount = math.exp(min(0.0, log_discount))
        spectral = (gap / large) * gap
        precision = discount * spectral
        if aperture < math.pi and endpoint_cost < math.sqrt(precision):
            precision = min(precision, 1.0 + endpoint_cost ** 2)

        # The excess vanishes on every vector singleton, including full-rank ones.
        baseline = residual if reference_cost is None else reference_cost
        excess = max(0.0, self.count - _residual_count(self.dimension, baseline))
        excess /= 1.0 + self.count
        area = 4.0 * (large / (1.0 + trace)) * (small / (1.0 + trace))
        penalty = area * excess / (1.0 + trace)
        precision /= 1.0 + tau * penalty

        log_factor = -0.5 * precision * np.sin(angles - centre) ** 2
        if aperture == math.pi:
            return log_factor
        gate = math.exp(-spectral) * (residual / (1.0 + residual))
        # Weight each oriented endpoint by the residual at the opposite endpoint.
        first = math.sqrt(ff) * sine
        second = math.sqrt(ff) * cosine + math.sqrt(end_cost)
        norm = math.hypot(1.0, first, second)
        projection = (first / norm) * np.cos(angles) - (second / norm) * np.sin(angles)
        return log_factor - gate * projection ** 2

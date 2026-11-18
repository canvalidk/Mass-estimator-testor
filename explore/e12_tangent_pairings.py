"""Orthogonal components of T36's residual at a specified mass (research only).

T45 differentiates with respect to phi = atan(m). Scores, weights and local
corrections here use that chart. Energies and the normalized score are
independent of a positive reparameterization. See e12-tangent-pairings.md.
"""
from dataclasses import dataclass
import math

import numpy as np

from met.conditional_score import _geometry, _mass
from met.full import Reading


def _declared_span_residual(reading, mass):
    """Drop orthogonal noise, expressing the result in T45's original residual frame."""
    basis = reading.subspace
    zeros = np.zeros_like(basis)
    span = np.block([[basis, zeros], [zeros, basis]])
    data = np.r_[reading.force, reading.acceleration]
    chol = np.linalg.cholesky(reading.cov)
    orthogonal, _ = np.linalg.qr(np.linalg.solve(chol, span))
    projection = chol @ (orthogonal @ (orthogonal.T @ np.linalg.solve(chol, data)))
    d = reading.dimension
    projected = Reading(projection[:d], projection[d:], reading.cov, basis)
    _, residual, _ = _geometry(projected, mass)
    return residual


@dataclass(frozen=True)
class TangentPairings:
    candidate_mass: float
    score: float = 0.
    weight: float = 0.
    within: float = 0.
    between: float = 0.
    residual_dimension: int = 0
    active: int = 0

    @property
    def correction(self):
        """Signed displacement in phi given by the pooled reciprocal tangent."""
        return self.score / self.weight if self.weight else None

    @property
    def common(self):
        if self.weight == 0.:
            return 0.
        signed = self.score / math.sqrt(self.weight)
        return signed * signed

    @property
    def degrees(self):
        """Conditional degrees of freedom in the order within, between, common."""
        rank = int(self.active > 0)
        return self.residual_dimension - self.active, self.active - rank, rank

    def merge(self, other):
        """Combine blocks at the same candidate, adding their between-block cost."""
        if self.candidate_mass != other.candidate_mass:
            raise ValueError("blocks must use the same candidate mass")
        weight = self.weight + other.weight
        difference = 0.
        if self.weight and other.weight:
            # The difference of the two normalized scores in their orthogonal
            # contrast direction. This avoids division by a tiny squared norm.
            difference = (math.sqrt(other.weight / weight) * (self.score / math.sqrt(self.weight))
                          - math.sqrt(self.weight / weight) * (other.score / math.sqrt(other.weight)))
        state = TangentPairings(
            self.candidate_mass, self.score + other.score, weight, self.within + other.within,
            self.between + other.between + difference * difference,
            self.residual_dimension + other.residual_dimension, self.active + other.active)
        if not all(math.isfinite(x) for x in (state.score, state.weight, state.within, state.between)):
            raise ValueError("pairing energies exceed the floating-point range")
        return state


def tangent_pairings(readings, mass):
    """Compute the candidate-dependent three-way decomposition for independent readings.

    Measured-zero rows are omitted, including from the realized rank convention.
    Empty, all-zero, and nonzero zero-tangent inputs have a zero common statistic.
    """
    mass = _mass(mass)
    result = TangentPairings(mass)
    for reading in readings:
        # The original observation sets the tangent, including the exact-zero
        # normal-equation convention. Reprojecting first can change its rank.
        _, residual, tangent = _geometry(reading, mass)
        if not np.any(reading.force) and not np.any(reading.acceleration):
            continue
        weight = float(tangent @ tangent)
        score = float(residual @ tangent)
        if reading.k < reading.dimension:
            residual = _declared_span_residual(reading, mass)
        if weight:
            unit = tangent / math.sqrt(weight)
            residual = residual - (residual @ unit) * unit
        within = float(residual @ residual)
        result = result.merge(TangentPairings(mass, score, weight, within, 0.,
                                              reading.k, int(weight > 0)))
    return result

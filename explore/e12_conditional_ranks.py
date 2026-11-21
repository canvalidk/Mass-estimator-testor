"""Increasing variation of the scalar conditional-score rank on the positive mass arc."""
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr


HALF_PI = np.pi / 2


@dataclass(frozen=True)
class RankMeasurement:
    major_angle: float
    largest: float
    smallest: float
    discount: float

    def increasing_mass(self, angle):
        """Mass of increasing rank on (0,angle], in noise-angle radians."""
        angle = np.asarray(angle, dtype=float)
        if np.any(~np.isfinite(angle) | (angle < 0) | (angle > HALF_PI)):
            raise ValueError("angle must lie on the closed positive mass arc")
        return self._branch_mass(0., angle)

    def _branch_mass(self, start, angle):
        if self.largest == self.smallest:
            return np.zeros_like(angle)
        width = np.arctan2(self.largest**0.25, self.smallest**0.25)
        result = np.zeros_like(angle)
        for center in (self.major_angle - np.pi, self.major_angle, self.major_angle + np.pi):
            lo = np.clip(start-center, -width, width)
            hi = np.clip(angle-center, -width, width)
            a = self._normal_argument(lo)
            b = self._normal_argument(hi)
            result += np.where(a > 0., ndtr(-a)-ndtr(-b), ndtr(b)-ndtr(a))
        return np.clip(result, 0., 1.)

    def _normal_argument(self, offset):
        sine, cosine = np.sin(offset), np.cos(offset)
        if self.smallest == 0:
            signed = np.sqrt(self.largest)*sine
        else:
            size = np.hypot(np.sqrt(self.largest)*cosine, np.sqrt(self.smallest)*sine)
            signed = (self.largest-self.smallest)*sine*cosine/size
        return self.discount*signed

    @property
    def prior_weight(self):
        """Unassigned tail and out-of-domain rank, without subtracting from one."""
        separation = np.sqrt(self.largest)-np.sqrt(self.smallest)
        tails = 2*ndtr(-self.discount*separation)
        return np.clip(tails+self._branch_mass(HALF_PI, np.pi), 0., 1.)

    def cdf(self, angle, prior_cdf):
        """Return the CDF for a supplied prior CDF in the same angle coordinate."""
        value = self.increasing_mass(angle) + self.prior_weight*prior_cdf(angle)
        return np.clip(value, 0., 1.)


def rank_measurement(rows, *, discount=True):
    """Use (force,acceleration) rows already standardized to independent unit noise."""
    data = np.asarray(rows, dtype=float)
    if data.size == 0:
        data = np.empty((0, 2))
    if data.ndim != 2 or data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError("expected finite (force, acceleration) pairs")
    scatter = data.T @ data
    vals, vectors = np.linalg.eigh(scatter)
    smallest, largest = np.maximum(vals, 0.)
    axis = vectors[:, -1]
    center = np.arctan2(axis[0], axis[1]) % np.pi
    factor = 1.
    if discount and largest > 0:
        residual = smallest if center <= HALF_PI else min(scatter.diagonal())
        count = -4*np.expm1(-np.einsum("ij,ij->i", data, data)/2).sum()
        if residual > 0:
            factor = min(1., np.sqrt(count)*np.sqrt(largest)/residual)
        axis_cap = np.hypot(1., min(scatter.diagonal())) / np.sqrt(largest)
        factor = min(factor, axis_cap)
    return RankMeasurement(center, largest, smallest, factor)

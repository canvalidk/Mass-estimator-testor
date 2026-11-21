"""E12 readout using each reading's residual curve and a fitted positive mass."""
import heapq
import math

import numpy as np

from explore.e12_compact_join import JoinState


def _checked_grams(grams, dimension):
    if not isinstance(dimension, (int, np.integer)) or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    grams = np.asarray(grams, dtype=float)
    if grams.shape == (0,):
        grams = grams.reshape(0, 2, 2)
    if grams.ndim != 3 or grams.shape[1:] != (2, 2) or not np.isfinite(grams).all():
        raise ValueError("expected a finite list of 2-by-2 Gram matrices")
    if not np.array_equal(grams[:, 0, 1], grams[:, 1, 0]):
        raise ValueError("channel Grams must be symmetric")
    return grams


def _minimum(grams, dimension, atol, arc_end=math.pi / 2):
    if not math.isfinite(atol) or atol <= 0:
        raise ValueError("atol must be positive and finite")
    if not len(grams):
        return 0.0
    spectrum, vectors = np.linalg.eigh(grams)
    tolerance = 16 * np.finfo(float).eps * np.max(np.abs(grams), axis=(1, 2))
    if np.any(spectrum[:, 0] < -tolerance):
        raise ValueError("channel Grams must be positive semidefinite")
    if dimension == 1 and np.any(spectrum[:, 0] > tolerance):
        raise ValueError("a scalar reading must have a Gram of rank at most one")
    # The scalar declaration gives the exact rank even when a large rank-one
    # product rounds to a matrix with a small positive eigenvalue.
    small = np.zeros(len(grams)) if dimension == 1 else np.maximum(spectrum[:, 0], 0.0)
    gaps = np.maximum(spectrum[:, 1] - small, 0.0)
    centres = np.mod(np.arctan2(vectors[:, 0, 1], vectors[:, 1, 1]), np.pi)
    scale = dimension / (1.0 - 0.5 ** (0.5 * dimension))

    def costs(theta):
        return small + gaps * np.sin(theta - centres) ** 2

    def objective(cost):
        return float(scale * (-np.expm1(-0.5 * cost)).sum())

    # Individual minima seed sharp features. Bound the initial grid on long lists.
    seeds = np.unique(centres[(gaps > 0) & (centres < arc_end)])
    if len(seeds) > 64:
        seeds = seeds[np.linspace(0, len(seeds) - 1, 64, dtype=int)]
    knots = np.unique(np.r_[0.0, seeds, arc_end])
    values = [objective(costs(angle)) for angle in knots]
    best = min(values)
    queue = []

    def enqueue(left, right, left_cost, right_cost, left_value, right_value):
        # Each reading has at most one minimum inside the closed positive arc.
        row_min = np.where((left <= centres) & (centres <= right),
                           small, np.minimum(left_cost, right_cost))
        # h'' <= 0, so the full second derivative is at most this local curvature.
        curvature = scale * np.dot(gaps, np.exp(-0.5 * row_min))
        chord_bound = min(left_value, right_value) - curvature * (right - left) ** 2 / 8
        bound = max(objective(row_min), chord_bound)
        if bound + atol < best:
            heapq.heappush(queue, (bound, left, right, left_value, right_value))

    left_cost = costs(knots[0])
    for i in range(len(knots) - 1):
        right_cost = costs(knots[i + 1])
        enqueue(knots[i], knots[i + 1], left_cost, right_cost, values[i], values[i + 1])
        left_cost = right_cost
    while queue and queue[0][0] + atol < best:
        _, left, right, left_value, right_value = heapq.heappop(queue)
        middle = (left + right) / 2
        if middle in (left, right):
            raise ArithmeticError("minimum tolerance exceeds angular floating-point resolution")
        middle_cost = costs(middle)
        middle_value = objective(middle_cost)
        best = min(best, middle_value)
        enqueue(left, middle, costs(left), middle_cost, left_value, middle_value)
        enqueue(middle, right, middle_cost, costs(right), middle_value, right_value)
    return best


def residual_minimum(grams, *, dimension=1, atol=1e-9):
    """Minimize the summed saturated costs over the closed arc to absolute tolerance."""
    return _minimum(_checked_grams(grams, dimension), dimension, atol)


def residual_log_evidence(grams, angles, *, tau=0.0, dimension=1, atol=1e-9,
                          aperture=math.pi / 2):
    """Evaluate the joined residual-curve law for a list of whitened channel Grams."""
    grams = _checked_grams(grams, dimension)
    if not 0 < aperture <= math.pi:
        raise ValueError("aperture must be positive and at most pi")
    count = _minimum(grams, dimension, atol, aperture)
    pooled = grams.sum(axis=0)
    state = JoinState(dimension, pooled[0, 0], pooled[0, 1], pooled[1, 1], count)
    return state._log_evidence(angles, tau, reference_cost=None, aperture=aperture)

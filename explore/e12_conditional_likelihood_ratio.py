"""Conditional profile-drop ranks on a scalar instrument with whitened channel data."""

from numbers import Integral

import numpy as np
from scipy.integrate import quad
from scipy.special import erfc, gammaincc, gammaincinv, ndtri


def conditional_lr_pvalue(scatter, angle, n):
    """Return the exact-test tail (numerical quadrature) at one or more mass angles.

    scatter is Y.T @ Y for n independent whitened (force, acceleration) pairs.
    The mass line has unit vector (sin(angle), cos(angle)). Angles must already
    reflect the whitening. The best comparison line ranges over the full
    projective circle. Very small probabilities can underflow in float64.
    """
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer")
    gram = np.asarray(scatter, dtype=float)
    angles = np.asarray(angle, dtype=float)
    if gram.shape != (2, 2) or not np.all(np.isfinite(gram)):
        raise ValueError("expected a finite two-channel scatter")
    if not np.isfinite(angles).all():
        raise ValueError("angles must be finite")
    tolerance = 32 * np.finfo(float).eps * np.max(np.abs(gram))
    if abs(gram[0, 1] - gram[1, 0]) > tolerance:
        raise ValueError("scatter must be symmetric")
    eigenvalues, directions = np.linalg.eigh(gram)
    small, largest = eigenvalues
    if small < -tolerance or (n == 1 and small > tolerance):
        raise ValueError("scatter must be positive semidefinite with rank at most n")
    gap = largest - max(0., small)
    # Resolves very small drops near the optimum without subtracting two
    # nearly equal large explained energies.
    transverse = directions[0, 0] * np.sin(angles) + directions[1, 0] * np.cos(angles)
    drop = gap * np.minimum(1., np.abs(transverse))**2
    if largest == 0 or gap == 0:
        result = np.ones_like(angles)
    elif n == 1:
        result = erfc(np.sqrt(drop / 2))
    else:
        result = np.array([_conditional_tail(value, largest, n)
                           for value in drop.flat]).reshape(angles.shape)
    return float(result) if angles.ndim == 0 else result


def _conditional_tail(drop, largest, n):
    if drop == 0:
        return 1.
    ratio = drop / largest
    if ratio >= 1:
        return gammaincc(n / 2, drop / 2)
    # Condition on |N(0,1)|, retaining its upper tail exactly. In extreme
    # data the contribution below the endpoint can itself be a tiny tail.
    baseline = erfc(np.sqrt(drop / 2))
    if gammaincc(n / 2, drop / 2) == 0.:
        return 0.
    endpoint = np.sqrt(drop)
    shape = (n - 1) / 2

    def contribution(z):
        remaining = largest * ((endpoint-z)/endpoint) * (1+z/endpoint)
        return np.sqrt(2/np.pi) * np.exp(-.5*z*z) * gammaincc(shape, remaining/2)

    # The integrand can be confined to a thin boundary layer near sqrt(drop).
    # Locate its transition using reference chi-square quantiles.
    quantiles = 2 * gammaincinv(shape, [.001, .5, .999])
    points = [endpoint*np.sqrt(1-j/largest) for j in quantiles if j < largest]
    correction, _ = quad(contribution, 0., endpoint, points=points,
                         epsabs=1e-305, epsrel=1e-10, limit=150)
    return np.clip(baseline + correction, 0., 1.)


def gaussian_rank_weight(scatter, angle, n):
    """Relative mass weight under the note's two-sided Gaussian rank convention."""
    tail = conditional_lr_pvalue(scatter, angle, n)
    z = ndtri(np.asarray(tail) / 2)
    return np.exp(-.5 * z**2)

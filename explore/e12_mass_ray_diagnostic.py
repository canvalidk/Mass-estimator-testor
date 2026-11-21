"""Test a prespecified mass allocation with independent Gaussian force pairs."""
from dataclasses import dataclass
import math

import numpy as np
from scipy.special import gammaincc


@dataclass(frozen=True)
class MassRayDiagnostic:
    standardized_contrasts: np.ndarray
    scores: np.ndarray
    mean_score: float
    upper_discrepancy: float
    chi_square: float
    null_p_value: float


def mass_ray_diagnostic(forces, covariances, allocation, *, confidence=0.95):
    """Measure departure from a fixed ray using separately measured forces.

    forces[i] contains (F1, F2) at push i. Supply a common 2-by-2 force
    covariance or one per push. The specified allocation r must precede
    observation. The targets are conditional on the true masses and drives.

    Let c=(1-r,-r), v_i=c.T @ covariance_i @ c and d_i=c.T @ true_force_i.
    A shared acceleration a_i gives d_i=a_i*((1-r)*m1_i-r*m2_i).
    mean_score is unbiased for mean_i[1-exp(-d_i**2/(4*v_i))]. Its
    one-sided Hoeffding upper confidence bound is clipped to [0, 1].
    The chi-square upper-tail p-value tests d_i=0 for every push.
    A zero-drive push has d_i=0 for every mass pair.
    """
    if (np.ndim(allocation) != 0 or not 0 < allocation < 1
            or np.ndim(confidence) != 0 or not 0 < confidence < 1):
        raise ValueError("allocation and confidence must be scalars strictly between 0 and 1")
    r = float(allocation)
    forces = np.asarray(forces, dtype=float)
    if forces.shape == (2,):
        forces = forces[None, :]
    if (forces.ndim != 2 or forces.shape[1] != 2 or not len(forces)
            or not np.isfinite(forces).all()):
        raise ValueError("expected one or more finite force pairs")
    n = len(forces)
    cov = np.asarray(covariances, dtype=float)
    if cov.shape == (2, 2):
        cov = cov[None, :, :]
    if (cov.shape not in ((1, 2, 2), (n, 2, 2)) or not np.isfinite(cov).all()
            or not np.allclose(cov, cov.swapaxes(-1, -2), rtol=1e-12, atol=0)):
        raise ValueError("expected symmetric finite 2-by-2 force covariances")
    diagonal = np.diagonal(cov, axis1=1, axis2=2)
    if np.any(diagonal < 0):
        raise ValueError("force covariance must be positive semidefinite")
    sd1, sd2 = np.sqrt(diagonal).T
    # Divide by the larger channel standard deviation first. The product of
    # the two variances could overflow even for a valid covariance.
    larger, smaller = np.maximum(sd1, sd2), np.minimum(sd1, sd2)
    cross = cov[:, 0, 1]
    if np.any((smaller == 0) & (cross != 0)):
        raise ValueError("zero-variance force channel must have zero cross covariance")
    nonzero = smaller > 0
    rho = np.zeros_like(smaller)
    rho[nonzero] = (cross[nonzero] / larger[nonzero]) / smaller[nonzero]
    if np.any(np.abs(rho) > 1 + 4*np.finfo(float).eps):
        raise ValueError("force covariance must be positive semidefinite")
    # Var(c.F) = (u-v)^2 + 2*u*v*(1-rho), including rank-one instruments.
    u, v = (1-r)*sd1, r*sd2
    spread = np.sqrt(u) * np.sqrt(v) * np.sqrt(2*(1-np.clip(rho, -1., 1.)))
    contrast_sd = np.hypot(u-v, spread)
    if np.any(contrast_sd == 0):
        raise ValueError("the specified force contrast must have positive variance")
    # Beyond floating range the squared statistic has a well-defined limiting
    # score of one and a chi-square upper-tail probability of zero.
    with np.errstate(over="ignore"):
        z = ((1-r)*forces[:, 0] - r*forces[:, 1]) / contrast_sd
        square = np.square(z)
        chi_square = float(square.sum())
    scores = -np.expm1(0.5 * (math.log(2) - square))
    mean_score = float(scores.mean())
    radius = math.sqrt(-math.log1p(-confidence) / n)
    upper = min(1., max(0., mean_score + radius))
    z.setflags(write=False)
    scores.setflags(write=False)
    return MassRayDiagnostic(z, scores, mean_score, upper, chi_square,
                             float(gammaincc(0.5*n, 0.5*chi_square)))

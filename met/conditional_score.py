"""Exact conditional Gaussian score tests, a research confidence-set reference.

This uses the declared full-covariance Gaussian model and arbitrary fixed
signals in each reading's subspace. It returns tests of candidate masses,
not a probability law of mass, a point estimate, or a production adapter.
No signal population, true amplitudes, or theoretical sample count enters.
"""
import math

import numpy as np
from scipy import stats

from met.full import Reading


def _mass(value):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("mass must be finite and positive")
    return value


def _geometry(reading, mass):
    if not isinstance(reading, Reading):
        raise ValueError("each reading must be a met.full.Reading")
    data = np.concatenate([reading.force, reading.acceleration])
    if (np.any(~np.isfinite(data)) or np.any(~np.isfinite(reading.cov))
            or np.any(~np.isfinite(reading.subspace))):
        raise ValueError("reading data and declarations must be finite")
    # Use phi=atan(m), so the mean subspace is unchanged but the basis
    # and its derivative stay bounded for both large and small masses.
    scale = math.hypot(1., mass)
    sine, cosine = mass/scale, 1./scale
    n = reading.subspace
    basis = np.vstack([sine*n, cosine*n])
    derivative = np.vstack([cosine*n, -sine*n])
    chol = np.linalg.cholesky(reading.cov)
    whitened = np.linalg.solve(chol, data)
    basis = np.linalg.solve(chol, basis)
    derivative = np.linalg.solve(chol, derivative)
    q, r = np.linalg.qr(basis, mode="complete")
    k = n.shape[1]
    parallel, perpendicular = q[:, :k], q[:, k:]
    coefficients = parallel.T @ whitened
    # Preserve exactly vanishing physical normal equations when recognized
    # in the supplied floating data. Trigonometric basis normalization and
    # QR can otherwise turn e.g. (F,a)=(1,-2), m=2 into a tiny nonzero U,
    # whose self-normalized statistic does not become small with U.
    try:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            precision_data = np.linalg.solve(reading.cov, data)
            d = reading.dimension
            normal_rhs = n.T @ (precision_data[:d]+precision_data[d:]/mass
                                 if mass >= 1 else mass*precision_data[:d]+precision_data[d:])
        if np.all(np.isfinite(normal_rhs)) and np.all(normal_rhs == 0.):
            coefficients = np.zeros_like(coefficients)
    except np.linalg.LinAlgError:
        # The optional normal-equation check must not defeat the stable
        # whitened QR path when solving with Sigma is less well conditioned.
        pass
    residual = perpendicular.T @ whitened
    # Tangent to the fitted mean, in residual coordinates. Only the same
    # reading's own candidate projection is used, never another's signal.
    tangent = perpendicular.T @ derivative @ np.linalg.solve(r[:k], coefficients)
    return coefficients, residual, tangent


def gaussian_score_test(readings, mass):
    """Test one positive mass using the full declared Gaussian input type.

    Let P_i be the whitened candidate mean projector. The profile score is
    S=sum V_i' P_i' U_i, and its exact conditional variance given all U_i
    is C=sum ||P_i' U_i||**2 (derivative in phi=atan(m)). At the true mass,
    S/sqrt(C) is standard normal conditional on U whenever C>0.

    The deterministic C=0 convention accepts the candidate (p_value=1).
    This exceptional event has probability zero at any fixed candidate in
    every nonempty, nonsingular declared Gaussian experiment of this type.
    At all-zero observed data or an empty series it makes the test total.
    It is not claimed to be a continuous extension at such exceptional data.
    Exact zero normal equations are recognized separately from QR. This is
    floating-point research code, not a certified detector of every exact
    cancellation for arbitrary ill-conditioned declarations.
    """
    mass = _mass(mass)
    scores, variances = [], []
    for reading in readings:
        _, residual, tangent = _geometry(reading, mass)
        scores.append(float(residual @ tangent))
        variances.append(float(tangent @ tangent))
    score, variance = math.fsum(scores), math.fsum(variances)
    if not math.isfinite(score) or not math.isfinite(variance):
        raise ValueError("score exceeds the floating-point range")
    degenerate = variance == 0.
    normal_score = 0. if degenerate else score/math.sqrt(variance)
    # Multiplication saturates to +inf for an extreme test statistic;
    # Python's exponentiation would instead raise OverflowError.
    statistic = normal_score*normal_score
    return dict(score=score, variance=variance, normal_score=normal_score,
                statistic=statistic, p_value=float(stats.chi2.sf(statistic, 1)),
                degenerate=degenerate)


def gaussian_score_accepts(readings, masses, *, probability):
    """Membership of the exact confidence set, at explicit candidate masses.

    These are pointwise tests, not a CDF or a probability assigned to mass.
    The set can be unbounded or disconnected. The probability is required;
    no confidence level is selected by this reference.
    """
    probability = float(probability)
    if not math.isfinite(probability) or not 0 < probability < 1:
        raise ValueError("probability must lie strictly between zero and one")
    readings = tuple(readings)
    masses = np.asarray(masses, dtype=float)
    if np.any(~np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError("candidate masses must be finite and positive")
    threshold = float(stats.chi2.ppf(probability, 1))
    result = np.array([gaussian_score_test(readings, value)["statistic"] <= threshold
                       for value in masses.flat], dtype=bool).reshape(masses.shape)
    return bool(result) if result.ndim == 0 else result

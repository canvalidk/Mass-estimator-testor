"""Research graph-subspace geometry and measure diagnostics for T49.

Matrices here describe mathematical operator graphs. They do not change
the production scalar-mass hypothesis or choose a prior/integration measure.
The canonical metric lives on hypothesis subspaces; T48's latent-vector
volume lives inside a fixed subspace and is a different object.
"""
import math

import numpy as np
from scipy import special


def _matrix(value):
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or min(value.shape) < 1 or np.any(~np.isfinite(value)):
        raise ValueError("operator must be a finite nonempty matrix")
    return value


def _basis(matrix, covariance):
    rows, columns = matrix.shape
    basis = np.vstack([matrix, np.eye(columns)])
    if covariance is None:
        return basis, None
    covariance = np.asarray(covariance, dtype=float)
    if (covariance.shape != (rows+columns, rows+columns)
            or np.any(~np.isfinite(covariance))
            or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
        raise ValueError("covariance must be finite, symmetric and match the graph dimension")
    try:
        chol = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError("covariance must be positive definite") from error
    return np.linalg.solve(chol, basis), chol


def graph_projector(matrix, *, covariance=None):
    """Orthogonal projector onto the graph in declared whitened coordinates."""
    matrix = _matrix(matrix)
    basis, _ = _basis(matrix, covariance)
    q, _ = np.linalg.qr(basis, mode="reduced")
    return q@q.T


def graph_projector_derivative(matrix, direction, *, covariance=None):
    """Derivative for the path M+tA with fixed declared covariance."""
    matrix, direction = _matrix(matrix), np.asarray(direction, dtype=float)
    if direction.shape != matrix.shape or np.any(~np.isfinite(direction)):
        raise ValueError("direction must be finite and have the operator shape")
    basis, chol = _basis(matrix, covariance)
    derivative = np.vstack([direction, np.zeros((matrix.shape[1], matrix.shape[1]))])
    if chol is not None:
        derivative = np.linalg.solve(chol, derivative)
    q, r = np.linalg.qr(basis, mode="reduced")
    transverse = derivative-q@(q.T@derivative)
    half = transverse@np.linalg.solve(r, q.T)
    return half+half.T


def graph_metric(matrix, direction, other=None, *, covariance=None):
    """Canonical pullback metric tr(P'_A P'_B)/2 on operator graphs.

    It is signal independent and is NOT the efficient Fisher information
    of a measured experiment. Whitening covariance is held fixed.
    """
    first = graph_projector_derivative(matrix, direction, covariance=covariance)
    second = first if other is None else graph_projector_derivative(matrix, other, covariance=covariance)
    return float(np.sum(first*second)/2)


def graph_chart_log_volume(matrix, *, covariance=None):
    """Canonical volume in a declared-noise graph chart, up to a constant.

    For M of shape p by q: det(I+M.T M)^(-(p+q)/2). SVD avoids squaring
    extreme singular values. For nonidentity covariance Omega=Sigma^-1,
    the result is det(Omega)^(q/2) det(B.T Omega B)^(-(p+q)/2).
    No normalization constant or prior is selected.
    """
    matrix = _matrix(matrix)
    if covariance is not None:
        basis, chol = _basis(matrix, covariance)
        _, r = np.linalg.qr(basis, mode="reduced")
        log_det_omega = -2*float(np.sum(np.log(np.diag(chol))))
        log_det_g = 2*float(np.sum(np.log(np.abs(np.diag(r)))))
        return (matrix.shape[1]/2)*log_det_omega-(sum(matrix.shape)/2)*log_det_g
    singular = np.linalg.svd(matrix, compute_uv=False)
    log_determinant = 2*float(np.sum(np.log(np.hypot(1., singular))))
    return -(sum(matrix.shape)/2)*log_determinant


def _scalar_inputs(mass, dimension, noise_ratio):
    mass, ratio = float(mass), float(noise_ratio)
    if (not math.isfinite(mass) or mass <= 0 or not math.isfinite(ratio) or ratio <= 0
            or isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 1):
        raise ValueError("mass and noise_ratio must be finite positive; dimension a positive integer")
    return mass, ratio


def scalar_curve_volume(mass, *, dimension, noise_ratio):
    """Intrinsic scalar-graph speed sqrt(d)*s/(s^2+m^2), not normalized."""
    mass, ratio = _scalar_inputs(mass, dimension, noise_ratio)
    return math.exp(.5*math.log(dimension)+math.log(ratio)-2*math.log(math.hypot(mass, ratio)))


def scalar_chart_substitution_density(mass, *, dimension, noise_ratio):
    """Diagnostic NAIVE substitution of M=(m/s)I into ambient chart density.

    This normalized 1D function is NOT conditioning Haar measure on scalar
    matrices or inducing canonical curve volume. For d>1 it gives a concrete
    different answer and violates the scalar reciprocal-prior requirement.
    """
    mass, ratio = _scalar_inputs(mass, dimension, noise_ratio)
    exponent = dimension**2
    log_constant = math.log(2)-.5*math.log(math.pi)+special.gammaln(exponent)-special.gammaln(exponent-.5)
    log_kernel = -2*exponent*(math.log(math.hypot(mass, ratio))-math.log(ratio))
    return math.exp(log_constant-math.log(ratio)+log_kernel)

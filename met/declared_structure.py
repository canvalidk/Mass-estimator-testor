"""E6 research profiles for explicitly declared relationships between signals.

These are Gaussian least-squares diagnostics, not an estimator or a new input
to ours30/32/33. No distribution on the nuisance signals is introduced.
See explore/e6-declared-structure.md for the declarations and their scope.
"""

import math

import numpy as np
from scipy.special import logsumexp


def gaussian_profile(observation, covariance, design):
    """Minimum Mahalanobis cost over the columns of a fixed mean design.

    covariance must be SPD. A copied observation with singular covariance
    must first be represented by its unique acquisition channels; treating
    those copies as independent is a different experiment.
    """
    y = np.asarray(observation, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    x = np.asarray(design, dtype=float)
    if y.ndim != 1 or x.ndim != 2 or x.shape[0] != len(y):
        raise ValueError("one observation vector and a matching design matrix required")
    if cov.shape != (len(y), len(y)) or not np.allclose(cov, cov.T, rtol=1e-12, atol=0.0):
        raise ValueError("a matching symmetric covariance required")
    if not all(np.isfinite(a).all() for a in (y, cov, x)):
        raise ValueError("finite inputs required")
    chol = np.linalg.cholesky(cov)
    z, w = np.linalg.solve(chol, y), np.linalg.solve(chol, x)
    coefficients = np.linalg.lstsq(w, z, rcond=None)[0]
    residual = z - w @ coefficients
    return float(residual @ residual)


def linear_mass_profile(force, acceleration, mass, basis, covariance=None):
    """Profile with declared stacked accelerations a* = basis @ t.

    Stacking order is (all force components; all acceleration components).
    basis is fixed before the data. Its column space is the declaration;
    changing its basis leaves the experiment unchanged.
    """
    f, a = _pairs(force, acceleration)
    m = _mass(mass)
    b = np.asarray(basis, dtype=float)
    if b.ndim != 2 or b.shape[0] != f.size:
        raise ValueError("basis must have one row per stacked acceleration component")
    y = np.r_[f.ravel(), a.ravel()]
    cov = np.eye(len(y)) if covariance is None else covariance
    return gaussian_profile(y, cov, np.vstack((m * b, b)))


def compress_gaussian_mean(observation, covariance, mean_basis):
    """Sufficient GLS coordinates for a fixed known Gaussian mean subspace.

    mean_basis has full column rank and does not depend on mass or data.
    Returns estimate, covariance, and ancillary residual cost. This is an
    exact statistical reduction conditional on the supplied mean constraint;
    it does not establish that constraint's physical or VD provenance.
    """
    y = np.asarray(observation, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    design = np.asarray(mean_basis, dtype=float)
    # Reuse the profile's finite/SPD/shape checks, without imposing mass.
    gaussian_profile(y, cov, design)
    if design.shape[1] == 0 or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("a nonempty full-column-rank mean basis required")
    chol = np.linalg.cholesky(cov)
    z, w = np.linalg.solve(chol, y), np.linalg.solve(chol, design)
    q, r = np.linalg.qr(w, mode="reduced")
    estimate = np.linalg.solve(r, q.T @ z)
    inverse_r = np.linalg.solve(r, np.eye(len(r)))
    reduced_cov = inverse_r @ inverse_r.T
    residual = z - q @ (q.T @ z)
    return estimate, reduced_cov, float(residual @ residual)


def compress_linear_experiment(force, acceleration, basis, covariance=None):
    """GLS pair assuming BOTH channel means lie in col(basis).

    No candidate mass is used. If the force constraint was inferred using
    Newton II from acceleration's constraint, that provenance still belongs
    to this reduction; mass-independent algebra does not remove it.
    Output coefficient dimension need not be a physical Reading dimension.
    """
    f, a = _pairs(force, acceleration)
    b = np.asarray(basis, dtype=float)
    if b.ndim != 2 or b.shape[0] != f.size:
        raise ValueError("basis must have one row per stacked channel component")
    zero = np.zeros_like(b)
    mean_basis = np.block([[b, zero], [zero, b]])
    y = np.r_[f.ravel(), a.ravel()]
    cov = np.eye(len(y)) if covariance is None else covariance
    estimate, reduced_cov, cost = compress_gaussian_mean(y, cov, mean_basis)
    count = b.shape[1]
    return estimate[:count], estimate[count:], reduced_cov, cost


def direction_log_factors(force, acceleration, masses, directions, covariances=None):
    """Own-reading Gaussian profile factors on explicit finite hypotheses.

    Result axes: reading, mass, direction. These finite references are for
    exact E6 witnesses, not a quadrature/default direction prior. Each push
    amplitude remains a separate unrestricted nuisance. Full ambient costs
    retain evidence about the shared direction outside its candidate line.
    """
    f, a = _pairs(force, acceleration)
    ms = np.asarray(masses, dtype=float)
    ns = np.asarray(directions, dtype=float)
    if ms.ndim != 1 or len(ms) == 0:
        raise ValueError("a nonempty vector of candidate masses required")
    ms = np.array([_mass(m) for m in ms])
    if ns.ndim != 2 or ns.shape[1] != f.shape[1] or not len(ns) or not np.isfinite(ns).all():
        raise ValueError("directions must be a finite nonempty J by d array")
    lengths = np.linalg.norm(ns, axis=1)
    if np.any(lengths == 0):
        raise ValueError("directions must be nonzero")
    ns = ns / lengths[:, None]
    size = 2*f.shape[1]
    covs = np.tile(np.eye(size), (len(f), 1, 1)) if covariances is None else np.asarray(covariances)
    if covs.shape != (len(f), size, size):
        raise ValueError("one covariance per reading required")
    result = np.empty((len(f), len(ms), len(ns)))
    for i, (fi, ai, cov) in enumerate(zip(f, a, covs)):
        y = np.r_[fi, ai]
        for j, m in enumerate(ms):
            c, s = 1/math.hypot(1, m), m/math.hypot(1, m)
            for k, n in enumerate(ns):
                result[i, j, k] = -0.5*gaussian_profile(y, cov, np.r_[s*n, c*n][:, None])
    return result


def update_discrete_joint(joint, log_factor):
    """Normalize one joint-state update; the shared reference enters once."""
    base, factor = np.asarray(joint, dtype=float), np.asarray(log_factor, dtype=float)
    if (base.ndim != 2 or factor.shape != base.shape or
            not np.isfinite(base).all() or np.any(base < 0) or not np.any(base > 0) or
            not np.isfinite(factor).all()):
        raise ValueError("a nonnegative nonzero finite joint and matching finite log factor required")
    with np.errstate(divide="ignore"):
        log_joint = np.log(base) + factor
    return np.exp(log_joint - logsumexp(log_joint))


def calibrated_discrete_joint(mass_weights, direction_weights, log_factors):
    """Joint (m,n) working law on explicit finite references.

    direction_weights already contain independent calibration evidence.
    They describe ONE shared fixture parameter, not a population of pushes.
    log_factors has axes reading, mass, direction; zero readings can be empty.
    """
    mw, nw = np.asarray(mass_weights, dtype=float), np.asarray(direction_weights, dtype=float)
    factors = np.asarray(log_factors, dtype=float)
    if mw.ndim != 1 or nw.ndim != 1 or factors.ndim != 3 or factors.shape[1:] != (len(mw), len(nw)):
        raise ValueError("weight vectors and reading by mass by direction factors required")
    if (not np.isfinite(mw).all() or not np.isfinite(nw).all() or
            np.any(mw < 0) or np.any(nw < 0) or not np.any(mw > 0) or not np.any(nw > 0)):
        raise ValueError("nonnegative nonzero finite reference weights required")
    # Normalize first so the outer product cannot overflow for large weights.
    mw, nw = mw / mw.max(), nw / nw.max()
    return update_discrete_joint(np.outer(mw, nw), factors.sum(axis=0))


def isotropic_profile(force, acceleration, mass, structure="free", direction=None):
    """Closed profiles on the unit, independent, isotropic Gaussian slice.

    free: every acceleration vector arbitrary.
    known_line: each acceleration is an arbitrary signed multiple of supplied n.
    unknown_line: one unknown common line, arbitrary signed multiples.
    equal: one unknown vector is the acceleration at every acquisition.

    The coordinates must be transported to a common physical frame before
    calling this helper. The relation must travel with any frame recoding.
    """
    f, a = _pairs(force, acceleration)
    m = _mass(mass)
    c, s = 1.0 / math.hypot(1.0, m), m / math.hypot(1.0, m)
    parallel = s * f + c * a
    perpendicular = c * f - s * a
    cost = float(np.sum(perpendicular**2))
    if structure == "free":
        return cost
    if structure == "known_line":
        n = np.asarray(direction, dtype=float)
        if n.shape != (f.shape[1],) or not np.isfinite(n).all() or np.linalg.norm(n) == 0:
            raise ValueError("a finite nonzero direction of the spatial dimension required")
        n = n / np.linalg.norm(n)
        remainder = parallel - np.outer(parallel @ n, n)
        return cost + float(np.sum(remainder**2))
    if structure == "unknown_line":
        # Eckart--Young: distance of the parallel data matrix to rank <= 1.
        singular = np.linalg.svd(parallel, compute_uv=False)
        return cost + float(singular[1:] @ singular[1:])
    if structure == "equal":
        remainder = parallel - parallel.mean(axis=0)
        return cost + float(np.sum(remainder**2))
    raise ValueError("structure must be free, known_line, unknown_line, or equal")


def _pairs(force, acceleration):
    f, a = np.asarray(force, dtype=float), np.asarray(acceleration, dtype=float)
    if f.ndim != 2 or a.shape != f.shape or not all(f.shape):
        raise ValueError("matching nonempty N by d arrays required")
    if not np.isfinite(f).all() or not np.isfinite(a).all():
        raise ValueError("finite inputs required")
    return f, a


def _mass(mass):
    m = float(mass)
    if not math.isfinite(m) or m <= 0:
        raise ValueError("mass must be positive and finite")
    return m

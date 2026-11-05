"""B15: a plug-in model score variance, for studies only.

``flat`` is T45's known-unit-Gaussian benchmark. ``ours33`` uses the
Student residual kernel obtained by T48's specified flat latent integration.
The latter is not the fixed-signal sampling conditional under a shared
inverse-gamma variance: that conditional also depends on the unknown signal.
Neither function selects a calibration rule. The working equation imports
nothing here; ours33 can opt into this study through an explicit switch.
"""
import math

import numpy as np

from met.full import Reading, bartlett_eta as flat_bartlett_eta
from met.ours31 import _per_reading
from met.ours32 import bartlett_eta
from met.ours33 import _l_terms


def _inputs(readings, c, nu_eff, factors):
    c = float(c)
    if not math.isfinite(c) or c <= 0:
        raise ValueError("c must be finite and positive")
    if factors not in ("flat", "ours33"):
        raise ValueError("factors must be 'flat' or 'ours33'")
    readings = list(readings)
    if np.ndim(nu_eff) == 0:
        nu_eff = [float(nu_eff)] * len(readings)
    else:
        nu_eff = [float(value) for value in nu_eff]
    if len(nu_eff) != len(readings) or not all(math.isfinite(x) and x > 0 for x in nu_eff):
        raise ValueError("nu_eff must have one positive finite value per reading")
    active, degrees = [], []
    for reading, nu in zip(readings, nu_eff):
        if not isinstance(reading, Reading):
            raise ValueError("readings must be met.full.Reading objects")
        if any(np.any(~np.isfinite(value)) for value in
               (reading.force, reading.acceleration, reading.cov, reading.subspace)):
            raise ValueError("reading data and declarations must be finite")
        # Remove zero summands before both the fit and its reductions. Appending
        # exact-zero readings must not perturb summation/search rounding.
        if np.any(reading.force != 0) or np.any(reading.acceleration != 0):
            active.append(reading)
            degrees.append(nu)
    return active, c, degrees


def residual_score_multiplier(k, nu):
    """Student residual score variance / T45's Gaussian variance (T45, B15)."""
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k < 1:
        raise ValueError("k must be a positive integer")
    nu = float(nu)
    if not math.isfinite(nu) or nu <= 0:
        raise ValueError("nu must be finite and positive")
    return 1.0 - 2.0 / (nu + k + 2.0)


def _tangent_squared(reading, c, phi):
    """||R P' U||² in phi=atan(m/c), allowing the whole projective line."""
    n = reading.subspace
    sine, cosine = math.sin(phi), math.cos(phi)
    chol = np.linalg.cholesky(reading.cov)
    data = np.linalg.solve(chol, np.r_[reading.force, reading.acceleration])
    basis = np.linalg.solve(chol, np.vstack([c * sine * n, cosine * n]))
    derivative = np.linalg.solve(chol, np.vstack([c * cosine * n, -sine * n]))
    parallel, upper = np.linalg.qr(basis, mode="reduced")
    coefficients = np.linalg.solve(upper, parallel.T @ data)
    tangent = derivative @ coefficients
    tangent -= parallel @ (parallel.T @ tangent)
    return float(tangent @ tangent)


def _variance(readings, c, phi, degrees, factors):
    return math.fsum(
        _tangent_squared(reading, c, phi)
        * (1.0 if factors == "flat" else residual_score_multiplier(reading.k, nu))
        for reading, nu in zip(readings, degrees))


def score_variance(readings, c, phi, nu_eff, *, factors="ours33"):
    """Fixed-candidate conditional variance under the specified working model.

    ``flat`` assumes fixed Gaussian V=1; its nu_eff argument is validated but
    does not affect its variance. ``ours33`` uses T48's residual Student law,
    giving sum (nu+k)/(nu+k+2) ||R P' U||². Data outside the declared span
    do not contribute. This expression can be evaluated at a fitted phi,
    but that does not make it an exact conditional post-fit sampling pivot.
    """
    readings, c, degrees = _inputs(readings, c, nu_eff, factors)
    phi = float(phi)
    if not math.isfinite(phi):
        raise ValueError("phi must be finite")
    return _variance(readings, c, phi, degrees, factors)


def eta_c(readings, c, nu_eff, *, factors="ours33"):
    """Return (eta_C, observed curvature, C, phi_hat), with the cap unchanged.

    Use the same whole-line fit as the chosen factor family, then evaluate
    its fixed-candidate variance at that fit. eta_C=min(1,H/C), with the
    existing eta=1 convention when C<=1e-300 or H<=0. Exact-zero observations
    are removed before fitting, so padding leaves all four outputs bitwise
    unchanged. An empty/all-zero series returns (1,0,0,0).
    """
    readings, c, degrees = _inputs(readings, c, nu_eff, factors)
    if not readings:
        return 1.0, 0.0, 0.0, 0.0
    if factors == "flat":
        _, curvature, _, phi = flat_bartlett_eta(readings, c, [1.0] * len(readings))
    else:
        states = _per_reading(readings, c, degrees)
        _, curvature, _, phi = bartlett_eta(states, c, terms=_l_terms)
    variance = _variance(readings, c, phi, degrees, factors)
    if not math.isfinite(curvature) or not math.isfinite(variance):
        raise ValueError("score variance or curvature exceeds floating-point range")
    eta = 1.0 if variance <= 1e-300 or curvature <= 0 else min(1.0, curvature / variance)
    return float(eta), float(curvature), float(variance), float(phi)

"""Single-pair baselines for known, independent diagonal measurement errors.

Inputs and returned masses retain physical units: F_true = m * a_true.
Each row is one 3D estimation task. Axes are never whitened independently in
the two channels and then treated as if they still had a common scalar mass.

Profiling the unrestricted latent acceleration coordinates from
    sum_i [(F_i - m*x_i)^2 / sf_i^2 + (a_i - x_i)^2 / sa_i^2]
gives Q(m) = sum_i (F_i - m*a_i)^2 / (sf_i^2 + m^2*sa_i^2).
This is a through-origin, independent-error specialization of weighted
errors-in-variables regression, not a reproduction of York's intercept fit.
See Mikkonen et al. (2019), Appendix A4, equations A13-A14.
"""

from copy import deepcopy
import math

import numpy as np

from comparison_estimators import METHOD_METADATA as ISOTROPIC_METADATA
from comparison_estimators import _inputs, _positive_profile
from comparison_estimators import baseline_points as isotropic_baseline_points


SOURCE_URL = "https://acp.copernicus.org/articles/19/12531/2019/"
YORK_DOI = "https://doi.org/10.1119/1.1632486"
_DOMAIN_RISK = (
    "Nondegenerate Gaussian measurements give positive-probability boundary, "
    "nonpositive, or unavailable masses. Under the declared infinite penalty "
    "for invalid positive masses, every unbounded criterion has infinite "
    "population risk, even if a finite sample happens to contain no failures."
)
METHOD_METADATA = deepcopy(ISOTROPIC_METADATA)
METHOD_METADATA["norm_floor"].update(
    formula="max(||F||, sqrt(2/pi)*rms(sigma_F)) / max(||a||, sqrt(2/pi)*rms(sigma_a))",
    assumptions="Existing explicit-floor heuristic, with one RMS noise scale per channel; it does not model unequal coordinate errors.")
METHOD_METADATA["positive_profile"].update(
    label="RMS-isotropic profile (misspecified for unequal axes)",
    formula="argmin_(m>0) ||F-m*a||^2 / (rms(sigma_F)^2 + m^2*rms(sigma_a)^2)",
    assumptions="Replaces each known diagonal covariance by its mean variance times the identity. Correct when each channel is isotropic; deliberately misspecified otherwise.")
METHOD_METADATA["noise_corrected_ols"].update(
    formula="(F dot a) / (||a||^2 - sum_i sigma_ai^2), only if the denominator is positive",
    assumptions="Corrects total acceleration energy using its known noise trace. The estimating equation is unbiased under independent channel errors; the one-pair ratio is not unbiased.")
METHOD_METADATA.update({
    "covariance_profile": {
        "label": "Diagonal-covariance Gaussian profile",
        "formula": "argmin_(m>0) sum_i (F_i-m*a_i)^2 / (sigma_Fi^2 + m^2*sigma_ai^2)",
        "assumptions": "Profiles unrestricted latent acceleration under known independent diagonal Gaussian errors. All positive stationary points and both boundaries are compared; zero/infinity are retained and a flat objective returns NaN. No log-determinant term belongs to this profile likelihood.",
        "source": SOURCE_URL, "source_locator": "Appendix A4, Eqs. A13-A14; independent errors and fixed zero intercept, followed here by exact elimination of adjusted coordinates.",
        "related_source": YORK_DOI, "uses_coordinate_covariance": True,
        "risk_note": _DOMAIN_RISK,
    },
    "weighted_forward_ols": {
        "label": "Force-weighted forward OLS",
        "formula": "sum_i F_i*a_i/sigma_Fi^2 / sum_i a_i^2/sigma_Fi^2",
        "assumptions": "Uses force-coordinate precision weights while treating observed acceleration as exact; acceleration noise still makes that assumption misspecified.",
        "source": SOURCE_URL, "source_locator": "Weighted through-origin specialization of the vertical residual criterion; Appendix A1.",
        "uses_coordinate_covariance": True, "risk_note": _DOMAIN_RISK,
    },
    "weighted_reverse_ols": {
        "label": "Acceleration-weighted reverse OLS",
        "formula": "sum_i F_i^2/sigma_ai^2 / sum_i F_i*a_i/sigma_ai^2",
        "assumptions": "Fits inverse mass with acceleration-coordinate precision weights while treating force as exact, then inverts. Force noise makes that assumption misspecified; a zero weighted cross-product is a pole.",
        "source": SOURCE_URL, "source_locator": "Channel-swapped weighted through-origin vertical residual criterion; Appendix A1.",
        "uses_coordinate_covariance": True, "risk_note": _DOMAIN_RISK,
    },
    "weighted_corrected_ols": {
        "label": "Force-weighted corrected moment",
        "formula": "sum_i F_i*a_i/sigma_Fi^2 / sum_i (a_i^2-sigma_ai^2)/sigma_Fi^2, only if the denominator is positive",
        "assumptions": "Known-variance correction to the forward weighted moment equation. Its expectation vanishes at true mass, but the ratio itself is biased and has a pole at zero corrected energy. This is a one-pair reduction, not pooled BCES.",
        "source": ISOTROPIC_METADATA["noise_corrected_ols"]["source"],
        "source_locator": "Local known-diagonal-variance extension of the corrected moment equation; Kelly Sec. 7.1.",
        "uses_coordinate_covariance": True, "risk_note": _DOMAIN_RISK,
    },
})


def _sigmas(sigma_force, sigma_acceleration, rows=None):
    result = []
    for name, value in (("sigma_force", sigma_force), ("sigma_acceleration", sigma_acceleration)):
        array = np.asarray(value, dtype=float)
        shapes = ((3,),) if rows is None else ((3,), (rows, 3))
        if array.shape not in shapes or not np.all(np.isfinite(array)) or np.any(array <= 0):
            raise ValueError(f"{name} must contain finite positive SDs of shape (3,) or (N,3)")
        result.append(array)
    return tuple(result)


def _rms(values):
    scale = np.max(values, axis=-1, keepdims=True)
    return np.squeeze(scale, axis=-1) * np.sqrt(np.mean((values / scale)**2, axis=-1))


def _row_scales(values):
    """Keep scalar scale broadcasting and give per-trial scales an axis."""
    return values if np.ndim(values) == 0 else values[:, None]


def _objective(z, force, acceleration, u, v):
    """Scaled Q; the reciprocal evaluation avoids squaring large masses."""
    if z > 1:
        inverse = 1 / z
        return np.sum((inverse*force - acceleration)**2 / (inverse*inverse*u + v))
    return np.sum((force - z*acceleration)**2 / (u + z*z*v))


def _polish(z, force, acceleration, u, v):
    """Newton refinement of the rational derivative near an enumerated root."""
    reciprocal = z > 1
    if reciprocal:
        z = 1 / z
        force, acceleration, u, v = acceleration, force, v, u
    product = force * acceleration
    linear = u*acceleration**2 - v*force**2
    for _ in range(12):
        denominator = u + v*z*z
        numerator = -u*product + linear*z + v*product*z*z
        gradient = np.sum(numerator / denominator**2)
        curvature = np.sum((linear + 2*v*product*z) / denominator**2
                           - 4*v*z*numerator / denominator**3)
        if curvature == 0 or not math.isfinite(curvature):
            break
        step = gradient / curvature
        if not math.isfinite(step):
            break
        candidate = z - step
        if candidate <= 0 or not math.isfinite(candidate):
            break
        z = candidate
        if abs(step) <= 8*np.finfo(float).eps*z:
            break
    return 1/z if reciprocal else z


def _profile_single(force, acceleration, u, v, products):
    scale = max(np.max(np.abs(force)), np.max(np.abs(acceleration)))
    if scale == 0:
        return np.nan
    force, acceleration = force/scale, acceleration/scale
    fa = force * acceleration
    terms = []
    for axis in range(3):
        # Q'/2 = sum_i N_i(m)/(u_i+v_i*m^2)^2, where N_i is quadratic.
        numerator = np.array([-u[axis]*fa[axis],
                              u[axis]*acceleration[axis]**2 - v[axis]*force[axis]**2,
                              v[axis]*fa[axis]])
        terms.append(np.convolve(numerator, products[axis]))
    coefficients = np.sum(terms, axis=0)
    magnitude = np.max(np.abs(coefficients))
    cancellation_scale = np.max(np.sum(np.abs(terms), axis=0))
    if magnitude <= 32*np.finfo(float).eps*cancellation_scale:
        # Exact flatness or cancellation beyond resolvable precision must not
        # be turned into an arbitrary finite positive mass.
        return np.nan
    coefficients /= magnitude
    candidates = []
    # Solve in m and 1/m: reciprocal roots improve conditioning near either
    # boundary. Positive denominators introduce no spurious positive roots.
    for reciprocal, polynomial in ((False, coefficients), (True, coefficients[::-1])):
        for root in np.polynomial.polynomial.polyroots(polynomial):
            if root.real > 0 and abs(root.imag) <= 1e-7*abs(root.real):
                candidate = 1/root.real if reciprocal else root.real
                if math.isfinite(candidate) and candidate > 0:
                    candidates.append(_polish(candidate, force, acceleration, u, v))
    # Both global boundary optima remain explicit failures of a positive-mass
    # point estimate. Infinity is the legacy convention for tied boundaries.
    candidates = sorted(set(candidates)) + [np.inf, 0.]
    values = [_objective(z, force, acceleration, u, v) for z in candidates]
    return candidates[int(np.argmin(values))]


def covariance_profile(force, acceleration, sigma_force, sigma_acceleration):
    """Globally compare all stationary points of the diagonal Gaussian profile.

    Eliminating positive denominators from Q' gives a degree-at-most-ten
    polynomial in three coordinates. Its roots are enumerated in both mass
    and reciprocal mass, then polished against Q' and compared with 0/infinity.
    This does not assume Q is unimodal. It is numerical polynomial solving,
    not a certified symbolic root isolation algorithm.
    """
    force, acceleration = _inputs(force, acceleration)
    sf, sa = _sigmas(sigma_force, sigma_acceleration, len(force))
    # Covariance can be supplied independently for every Monte Carlo trial.
    # The existing global stationary-point comparison is itself per row;
    # retaining it here also preserves its exact scalar specialization.
    if sf.ndim == 2 or sa.ndim == 2:
        sf = np.broadcast_to(sf, force.shape)
        sa = np.broadcast_to(sa, acceleration.shape)
        return np.array([covariance_profile(f, a, row_sf, row_sa)[0]
                         for f, a, row_sf, row_sa in zip(force, acceleration, sf, sa)])
    rms_f, rms_a = _rms(sf), _rms(sa)
    mass_scale = rms_f / rms_a
    if not math.isfinite(mass_scale) or mass_scale <= 0:
        raise ValueError("Noise scales must have a representable positive ratio")
    normalized_f, normalized_a = force/rms_f, acceleration/rms_a
    u, v = (sf/rms_f)**2, (sa/rms_a)**2
    if np.any(u <= 0) or np.any(v <= 0):
        raise ValueError("Noise contrast is too large for double-precision profile calculation")
    if np.all(sf == sf[0]) and np.all(sa == sa[0]):
        with np.errstate(invalid="ignore", divide="ignore"):
            return mass_scale * _positive_profile(np.sum(normalized_f**2, axis=1),
                                                   np.sum(normalized_a**2, axis=1),
                                                   np.sum(normalized_f*normalized_a, axis=1))
    # Normalize each denominator by a constant, which changes neither Q nor
    # physical mass. This is not a separate whitening of the observed vectors.
    scale = u + v
    normalized_f, normalized_a = normalized_f/np.sqrt(scale), normalized_a/np.sqrt(scale)
    u, v = u/scale, v/scale
    squares = [np.array([u[i]**2, 0., 2*u[i]*v[i], 0., v[i]**2]) for i in range(3)]
    products = [np.convolve(squares[(i+1) % 3], squares[(i+2) % 3]) for i in range(3)]
    result = [_profile_single(f, a, u, v, products) for f, a in zip(normalized_f, normalized_a)]
    return mass_scale * np.asarray(result)


def baseline_points(force, acceleration, sigma_force, sigma_acceleration):
    """Return ten (N,) point arrays, using SDs of shape (3,) or (N,3)."""
    force, acceleration = _inputs(force, acceleration)
    sf, sa = _sigmas(sigma_force, sigma_acceleration, len(force))
    result = isotropic_baseline_points(force, acceleration)
    rms_f, rms_a = _rms(sf), _rms(sa)
    floor = math.sqrt(2/math.pi)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        p, q = np.sum(force**2, axis=1), np.sum(acceleration**2, axis=1)
        dot = np.sum(force*acceleration, axis=1)
        result["norm_floor"] = np.maximum(np.sqrt(p), floor*rms_f) / np.maximum(np.sqrt(q), floor*rms_a)
        normalized_f, normalized_a = force/_row_scales(rms_f), acceleration/_row_scales(rms_a)
        result["positive_profile"] = rms_f/rms_a * _positive_profile(
            np.sum(normalized_f**2, axis=1), np.sum(normalized_a**2, axis=1),
            np.sum(normalized_f*normalized_a, axis=1))
        correction = q - np.sum(sa**2, axis=-1)
        result["noise_corrected_ols"] = np.full(q.shape, np.nan)
        np.divide(dot, correction, out=result["noise_corrected_ols"], where=correction > 0)
        # A common scaling of all weights cancels out of each weighted ratio.
        wf, wa = (np.min(sf, axis=-1, keepdims=True)/sf)**2, (np.min(sa, axis=-1, keepdims=True)/sa)**2
        cross_f = np.sum(force*acceleration*wf, axis=1)
        energy_a = np.sum(acceleration**2*wf, axis=1)
        result["weighted_forward_ols"] = cross_f / energy_a
        result["weighted_reverse_ols"] = np.sum(force**2*wa, axis=1) / np.sum(force*acceleration*wa, axis=1)
        corrected = energy_a - np.sum(sa**2*wf, axis=-1)
        result["weighted_corrected_ols"] = np.full(q.shape, np.nan)
        np.divide(cross_f, corrected, out=result["weighted_corrected_ols"], where=corrected > 0)
    result["covariance_profile"] = covariance_profile(force, acceleration, sf, sa)
    return result

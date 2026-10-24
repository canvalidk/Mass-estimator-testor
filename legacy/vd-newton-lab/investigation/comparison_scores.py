"""External-truth scores for the mass-estimator comparison study (NumPy only).

Every array element is a trial, and every estimator is scored against the same
positive true mass. Nonfinite and nonpositive estimates are failures: they get
infinite unbounded losses, the maximum capped loss, and failed tolerances.
Unbounded means are never calculated only on the surviving estimates.

Reported means, quantiles, variances, and Monte Carlo standard errors describe
the supplied finite sample; their existence does not establish finite
population moments. MCSEs assume independent trials (paired trials for
``paired_summary``). No observational agreement or compatibility gate is used.

``absolute_mass`` and ``squared_mass`` use the supplied mass units, or squared
units. They are kg MAE/MSE only when inputs are in kilograms. All other losses
are dimensionless. ``reciprocal_root`` is the mass-only normalized loss, not
the magnitude-weighted pair loss that elicits E[F]/E[A].
"""

import numpy as np


TOLERANCE_FACTORS = (1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0)
CAP_FACTORS = (2.0, 5.0, 10.0, 20.0)
ABS_LOG_QUANTILES = (0.5, 0.9, 0.95, 0.99)
UNBOUNDED_LOSSES = (
    "abs_log", "squared_log", "reciprocal_root", "abs_relative",
    "squared_relative", "absolute_mass", "squared_mass",
)


def _factor_label(value):
    return format(value, "g")


def _inputs(estimate, truth):
    estimates = np.asarray(estimate, dtype=float)
    if estimates.size == 0:
        raise ValueError("At least one estimate is required.")
    if np.ndim(truth) != 0:
        raise ValueError("truth must be one finite positive scalar mass.")
    true_mass = float(truth)
    if not np.isfinite(true_mass) or true_mass <= 0:
        raise ValueError("truth must be one finite positive scalar mass.")
    valid = np.isfinite(estimates) & (estimates > 0)
    return estimates, true_mass, valid


def _log_errors(estimates, truth, valid):
    """Avoid overflowing a ratio, while preserving small relative errors."""
    result = np.full(estimates.shape, np.inf)
    values = estimates[valid]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        relative = (values - truth) / truth
        log_error = np.log(values) - np.log(truth)
        safe_relative = np.isfinite(relative) & (relative > -1)
        log_error[safe_relative] = np.log1p(relative[safe_relative])
    result[valid] = log_error
    return result


def point_losses(estimate, truth):
    """Return name -> per-trial loss arrays, preserving the input shape.

    All returned quantities are losses: smaller is better. Capped squared-log
    loss is min(log(estimate/truth)**2, log(cap_factor)**2), without rescaling.
    ``outside_factor_*`` is zero on success and one on failure, including every
    invalid estimate. Tolerance endpoints are inclusive.
    """
    estimates, true_mass, valid = _inputs(estimate, truth)
    errors = _log_errors(estimates, true_mass, valid)
    absolute_log = np.abs(errors)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        squared_log = absolute_log ** 2
        # This identity is stable close to zero, unlike cosh(z) - 1.
        reciprocal_root = 4 * np.sinh(absolute_log / 4) ** 2
        absolute_mass = np.full(estimates.shape, np.inf)
        absolute_mass[valid] = np.abs(estimates[valid] - true_mass)
        absolute_relative = absolute_mass / true_mass
        losses = {
            "abs_log": absolute_log,
            "squared_log": squared_log,
            "reciprocal_root": reciprocal_root,
            "abs_relative": absolute_relative,
            "squared_relative": absolute_relative ** 2,
            "absolute_mass": absolute_mass,
            "squared_mass": absolute_mass ** 2,
        }
        for factor in CAP_FACTORS:
            losses["capped_squared_log_factor_" + _factor_label(factor)] = (
                np.minimum(squared_log, np.log(factor) ** 2)
            )
        for factor in TOLERANCE_FACTORS:
            success = valid & (estimates >= true_mass / factor) & (
                estimates <= true_mass * factor
            )
            losses["outside_factor_" + _factor_label(factor)] = (~success).astype(float)
    return losses


def _finite_mean_mcse(values):
    """Finite-input sample statistics without overflowing sums or squares."""
    values = np.asarray(values, dtype=float).ravel()
    scale = float(np.max(np.abs(values)))
    if scale == 0:
        return 0.0, 0.0 if values.size > 1 else None
    normalized = values / scale
    mean = float(np.mean(normalized) * scale)
    mcse = None
    if values.size > 1:
        mcse = float((np.std(normalized, ddof=1) / np.sqrt(values.size)) * scale)
    return mean, mcse


def _extended_quantile(values, probability):
    """Linear empirical quantile with +inf failures, avoiding inf-inf NaNs."""
    ordered = np.sort(np.asarray(values, dtype=float).ravel())
    position = (ordered.size - 1) * probability
    low, high = int(np.floor(position)), int(np.ceil(position))
    if low == high:
        return float(ordered[low])
    if not np.isfinite(ordered[high]):
        return np.inf
    weight = position - low
    return float((1 - weight) * ordered[low] + weight * ordered[high])


def score_point(estimate, truth):
    """Summarize all trials with JSON-safe numbers and explicit statuses.

    ``means`` contains the empirical mean of each loss. For an unbounded loss,
    any invalid estimate makes its mean and MCSE null; numerical overflow does
    likewise. Bounded scores always retain all trials. An MCSE is null when
    only one trial is supplied.

    Log variance uses ddof=0, so squared-log mean equals signed-log-bias squared
    plus log variance. Quantiles include invalid estimates as +infinite errors;
    a quantile reaching these errors is null with its own explicit status.
    """
    estimates, true_mass, valid = _inputs(estimate, truth)
    losses = point_losses(estimates, true_mass)
    count = int(estimates.size)
    valid_count = int(np.count_nonzero(valid))
    invalid_count = count - valid_count
    result = {
        "n": count,
        "valid_count": valid_count,
        "invalid_count": invalid_count,
        "finite_positive_fraction": valid_count / count,
        "status": "all_estimates_valid" if invalid_count == 0 else "invalid_estimates_present",
        "mcse_basis": "sample_sd_over_sqrt_n_population_moments_not_assessed",
        "means": {},
        "mcse": {},
        "mean_status": {},
    }
    for name, values in losses.items():
        if name in UNBOUNDED_LOSSES and invalid_count:
            mean, mcse, status = None, None, "invalid_estimates"
        elif not np.isfinite(values).all():
            mean, mcse, status = None, None, "nonfinite_loss"
        else:
            mean, mcse = _finite_mean_mcse(values)
            status = "finite_sample_mean"
        result["means"][name] = mean
        result["mcse"][name] = mcse
        result["mean_status"][name] = status

    diagnostics = {
        "signed_log_bias": None,
        "log_variance": None,
        "variance_ddof": 0,
        "status": "invalid_estimates" if invalid_count else "finite_sample_diagnostics",
    }
    if invalid_count == 0:
        errors = _log_errors(estimates, true_mass, valid)
        diagnostics["signed_log_bias"] = _finite_mean_mcse(errors)[0]
        diagnostics["log_variance"] = float(np.var(errors, ddof=0))
    result["log_diagnostics"] = diagnostics

    result["absolute_log_quantiles"] = {}
    result["absolute_log_quantile_status"] = {}
    for probability in ABS_LOG_QUANTILES:
        label = _factor_label(probability)
        quantile = _extended_quantile(losses["abs_log"], probability)
        result["absolute_log_quantiles"][label] = quantile if np.isfinite(quantile) else None
        result["absolute_log_quantile_status"][label] = (
            "finite_sample_quantile" if np.isfinite(quantile) else "infinite_failure_error"
        )

    result["tolerance_success"] = {}
    for factor in TOLERANCE_FACTORS:
        label = _factor_label(factor)
        key = "outside_factor_" + label
        result["tolerance_success"][label] = {
            "fraction": 1 - result["means"][key],
            "mcse": result["mcse"][key],
        }
    return result


def paired_summary(loss_a, loss_b):
    """Mean paired loss A-B and its sample MCSE; negative favors A.

    Inputs must contain the same trials in the same order and have the same
    shape. Any nonfinite input or difference makes the overall comparison
    unavailable; no successful-trial subset is silently substituted.
    """
    first = np.asarray(loss_a, dtype=float)
    second = np.asarray(loss_b, dtype=float)
    if first.shape != second.shape or first.size == 0:
        raise ValueError("Paired losses must have the same nonempty shape.")
    finite_pairs = np.isfinite(first) & np.isfinite(second)
    finite_count = int(np.count_nonzero(finite_pairs))
    result = {
        "n": int(first.size),
        "finite_pair_count": finite_count,
        "nonfinite_pair_count": int(first.size) - finite_count,
        "difference_definition": "A_minus_B",
        "mcse_basis": "sample_sd_over_sqrt_n_population_moments_not_assessed",
        "mean_difference": None,
        "mcse": None,
        "status": "nonfinite_pairs",
    }
    if not finite_pairs.all():
        return result
    with np.errstate(over="ignore", invalid="ignore"):
        difference = first - second
    if not np.isfinite(difference).all():
        result["status"] = "nonfinite_differences"
        return result
    result["mean_difference"], result["mcse"] = _finite_mean_mcse(difference)
    result["status"] = "finite_sample_comparison" if first.size > 1 else "single_pair_no_mcse"
    return result


def interval_score(log_lower, log_upper, logtruth, alpha):
    """Central (1-alpha) interval score on natural-log mass, lower is better.

    Returns width + (2/alpha) times the distance outside the interval. Arguments
    broadcast as NumPy arrays. The log truth must be finite, and 0 < alpha < 1.
    Unbounded, reversed, or NaN predicted intervals receive infinite loss.
    This scores reported uncertainty, independently of any point estimate.
    """
    if np.ndim(alpha) != 0 or not np.isfinite(alpha) or not 0 < float(alpha) < 1:
        raise ValueError("alpha must be one finite scalar strictly between zero and one.")
    lower, upper, truth = np.broadcast_arrays(
        np.asarray(log_lower, dtype=float),
        np.asarray(log_upper, dtype=float),
        np.asarray(logtruth, dtype=float),
    )
    if not np.isfinite(truth).all():
        raise ValueError("logtruth must be finite.")
    good = np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
    result = np.full(lower.shape, np.inf)
    with np.errstate(over="ignore", invalid="ignore"):
        outside_distance = (
            np.maximum(lower[good] - truth[good], 0)
            + np.maximum(truth[good] - upper[good], 0)
        )
        result[good] = upper[good] - lower[good] + 2 * (outside_distance / float(alpha))
    return result

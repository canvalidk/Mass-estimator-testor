"""Single-pair point comparators for the standardized 3D mass benchmark.

Inputs are F/sigma_F and a/sigma_a, with independent unit-variance errors in
all six coordinates. Returned points estimate the dimensionless mass
m*sigma_a/sigma_F. These are through-origin scalar reductions, not reproductions
of complete regression or robot-identification algorithms in the cited papers.

Negative values, zero, infinity, and NaN remain visible to the evaluator. Only
the explicitly named ``norm_floor`` comparator imposes a positive floor.
"""

import math

import numpy as np


KELLY_URL = "https://arxiv.org/pdf/0705.2774"
KUBUS_URL = (
    "https://www.researchgate.net/publication/224339569_"
    "On-Line_Estimation_of_Inertial_Parameters_Using_a_"
    "Recursive_Total_Least-Squares_Approach"
)

METHOD_METADATA = {
    "norm_ratio": {
        "label": "Vector-length ratio",
        "formula": "||F|| / ||a||",
        "assumptions": "Plug-in magnitudes; does not correct vector-length noise bias.",
        "source": "Existing local crossover comparator.",
        "risk_note": (
            "Under independent nondegenerate 3D Gaussian errors, squared mass loss "
            "has finite mean but infinite variance: inverse fourth powers of ||a|| "
            "are not integrable near zero. An empirical standard error for that "
            "mean does not have its usual finite-variance justification."
        ),
    },
    "norm_floor": {
        "label": "Vector-length ratio with explicit noise floor",
        "formula": "max(||F||, sqrt(2/pi)) / max(||a||, sqrt(2/pi))",
        "assumptions": (
            "Existing heuristic baseline; floor fixed in unit-noise coordinates. "
            "The floor is not the expected length of a 3D Gaussian vector."
        ),
        "source": "Existing local crossover comparator.",
        "risk_note": "The explicit denominator floor bounds all inverse powers.",
    },
    "positive_profile": {
        "label": "Positive-slope Gaussian profile / TLS",
        "formula": "argmin_(m>0) ||F-m*a||^2/(1+m^2)",
        "assumptions": (
            "Profiles unrestricted latent acceleration under independent isotropic "
            "unit Gaussian errors in both channels. Boundaries return 0 or infinity; "
            "a completely flat objective returns NaN. This is not recursive robot TLS."
        ),
        "source": KUBUS_URL,
        "source_locator": "Section IV-D, Eq.23 (errors in both matrix and wrench).",
        "risk_note": (
            "Nondegenerate Gaussian data give positive-probability 0/infinity "
            "boundary outputs. Unbounded loss with infinite penalty for invalid "
            "positive masses consequently has infinite population risk."
        ),
    },
    "forward_ols": {
        "label": "Forward through-origin OLS",
        "formula": "(F dot a) / ||a||^2",
        "assumptions": (
            "Minimizes ||F-m*a||^2, treating observed acceleration as exact. "
            "This assumption is deliberately misspecified in the noisy-both-channels benchmark."
        ),
        "source": KELLY_URL,
        "source_locator": "Section 3, pp.3-4; Section 7.1, pp.19-20 (OLS attenuation).",
        "risk_note": (
            "Can be nonpositive. Even before domain penalties, the second moment "
            "is finite in 3D but the fourth is not, so usual squared-loss MCSE "
            "does not have a finite-variance justification."
        ),
    },
    "reverse_ols": {
        "label": "Inverted reverse through-origin OLS",
        "formula": "||F||^2 / (F dot a)",
        "assumptions": (
            "Fits a=b*F treating force as exact, then reports 1/b. "
            "This assumption is deliberately misspecified; it can give negative or infinite mass."
        ),
        "source": KELLY_URL,
        "source_locator": "Through-origin channel-swapped OLS reduction; not Kelly's model.",
        "risk_note": (
            "The pole at F dot a=0 causes divergent absolute and squared mass "
            "risk under nondegenerate independent Gaussian errors, in addition "
            "to nonpositive outputs. Finite simulation means can conceal this."
        ),
    },
    "noise_corrected_ols": {
        "label": "Noise-corrected moment slope",
        "formula": "(F dot a) / (||a||^2 - 3), if ||a||^2 > 3; otherwise NaN",
        "assumptions": (
            "Known independent unit errors imply E[F dot a]=m*||a_true||^2 "
            "and E[||a||^2-3]=||a_true||^2. Solves the resulting corrected "
            "moment equation for one pair. The estimating equation is unbiased; "
            "the ratio is not. This fragile one-pair reduction is not pooled BCES. "
            "A nonpositive corrected energy is flagged as unavailable, not floored."
        ),
        "source": KELLY_URL,
        "source_locator": "Section 7.1, p.20: corrected-moment/BCES slope and instability.",
        "risk_note": (
            "Besides unavailable/nonpositive outcomes, approaching ||a||^2=3 "
            "from above gives a pole with divergent absolute and squared mass "
            "loss. Finite simulation averages are only sample diagnostics."
        ),
    },
}


def _inputs(force, acceleration):
    force = np.asarray(force, dtype=float)
    acceleration = np.asarray(acceleration, dtype=float)
    if force.ndim == 1:
        force = force[None, :]
    if acceleration.ndim == 1:
        acceleration = acceleration[None, :]
    if (force.ndim != 2 or acceleration.shape != force.shape
            or force.shape[1] != 3):
        raise ValueError("force and acceleration must have matching shape (N, 3) or (3,)")
    if not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration)):
        raise ValueError("input measurements must be finite")
    return force, acceleration


def _positive_profile(p, q, dot):
    """Same positive-profile rule and boundary conventions as the crossover study."""
    difference = p - q
    disc = np.hypot(difference, 2 * dot)
    result = np.where(p >= q, np.inf, 0.)
    positive = dot > 0
    high = positive & (p >= q)
    low = positive & (p < q)
    result[high] = (difference[high] + disc[high]) / (2 * dot[high])
    result[low] = 2 * dot[low] / (-difference[low] + disc[low])
    # For negative dot and p=q both unattained boundaries tie; infinity is the
    # existing reporting convention. For dot=0 and p=q every m is optimal.
    result[(p == q) & (dot == 0)] = np.nan
    return result


def baseline_points(force, acceleration):
    """Return six named (N,) arrays, one estimate per paired 3D observation.

    Rows are independent estimation tasks, not samples pooled into one mass.
    Input sigmas are fixed to one; standardize each channel before calling and
    multiply outputs by sigma_F/sigma_a to return to physical mass units.
    Undefined ratios retain NumPy's NaN/inf conventions. Nonpositive finite
    slopes are retained for the scorer to classify as invalid positive masses.
    """
    force, acceleration = _inputs(force, acceleration)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        p = np.sum(force * force, axis=1)
        q = np.sum(acceleration * acceleration, axis=1)
        dot = np.sum(force * acceleration, axis=1)
        nf, na = np.sqrt(p), np.sqrt(q)
        floor = math.sqrt(2 / math.pi)
        corrected = np.full(q.shape, np.nan)
        np.divide(dot, q - 3., out=corrected, where=q > 3.)
        return {
            "norm_ratio": nf / na,
            "norm_floor": np.maximum(nf, floor) / np.maximum(na, floor),
            "positive_profile": _positive_profile(p, q, dot),
            "forward_ols": dot / q,
            "reverse_ols": p / dot,
            "noise_corrected_ols": corrected,
        }

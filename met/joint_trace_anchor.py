"""T60/T61 research finite joint law with a frame-invariant trace anchor.

The law is a supplied joint prior times exp[-T/2-lambda*K/2], with exact
conditional-score statistic T, Gaussian profile cost K, and explicit
lambda=kappa/sqrt(1+sum ||W_i Y_i||^2). Kappa and numerical rank tolerances
are required inputs to this research comparison, not selected defaults.

The finite construction accepts per-event SPD covariance declarations and
is invariant under independent event frames. T61's asymptotic extension
allows arbitrary overall precisions and uniformly comparable covariance
shapes under its information-growth condition. T63 shows why unrestricted
degenerating shapes need a different argument. No atomic-prior TV limit is claimed.
No production estimator, signal population or prior coupling is selected.
"""
import math

import numpy as np
from scipy import linalg, special

from met.joint_mass import FiniteJointMassLaw, JointReading, gaussian_profile_cost
from met.joint_score_law import gaussian_joint_score_statistic


def _nonnegative(value, name):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _state(readings, mass_count, anchor_prefactor, rank_rtol, rank_atol):
    prefactor = float(anchor_prefactor)
    if not math.isfinite(prefactor) or prefactor <= 0:
        raise ValueError("anchor_prefactor must be finite and positive")
    relative = _nonnegative(rank_rtol, "rank_rtol")
    absolute = _nonnegative(rank_atol, "rank_atol")
    dimension, squared_norms = None, []
    all_zero = True
    for reading in readings:
        if (not isinstance(reading, JointReading)
                or reading.protocol != "common_acceleration"
                or reading.body_indices != tuple(range(mass_count))):
            raise ValueError("canonical full-body common-acceleration events are required")
        if dimension is None:
            dimension = reading.dimension
        elif reading.dimension != dimension:
            raise ValueError("all events must have the same spatial dimension")
        whitened = np.linalg.solve(np.linalg.cholesky(reading.covariance), reading.data)
        with np.errstate(over="ignore", invalid="ignore"):
            squared = float(whitened@whitened)
        if not math.isfinite(squared):
            raise ValueError("observed whitened trace exceeds the floating-point range")
        squared_norms.append(squared)
        all_zero = all_zero and not np.any(reading.data)
    try:
        trace = math.fsum(squared_norms)
    except OverflowError as error:
        raise ValueError("observed whitened trace exceeds the floating-point range") from error
    weight = prefactor/math.sqrt(1+trace)
    if not math.isfinite(weight) or weight <= 0:
        raise ValueError("anchor weight exceeds the supported floating-point range")
    return dict(anchor_prefactor=prefactor, anchor_weight=weight,
                observed_whitened_trace=trace, body_count=mass_count,
                spatial_dimension=dimension, rank_rtol=relative, rank_atol=absolute,
                all_zero=bool(all_zero),
                finite_scope="declared per-event SPD common-acceleration geometry",
                asymptotic_scope="T60/T61 comparable covariance shapes and stated information growth; continuous regular prior")


def joint_trace_anchor(readings, prior, *, anchor_prefactor, rank_rtol, rank_atol):
    """Return (finite joint law, observed state) under the explicit trace rule.

    Each source event is used once. Covariances need not be numerically
    identical: independently rotated representations are valid, and the finite
    law does not infer an asymptotic covariance-shape bound from matrix equality.
    Body labels must be canonical and all bodies must share one acceleration
    in every event. Other declared protocols are outside this bounded helper.

    Minima of score statistics and profile costs are removed separately,
    then the combined contrast is shifted again before adding log prior
    probabilities. Its minimum can still be huge when the two individual
    minima occur at different atoms. Finite support is not a continuous-grid
    accuracy or Gaussian-limit assertion. Empty/all-zero events return the same law.
    """
    if not isinstance(prior, FiniteJointMassLaw):
        raise ValueError("prior must be a FiniteJointMassLaw")
    readings = tuple(readings)
    state = _state(readings, prior.mass_count, anchor_prefactor, rank_rtol, rank_atol)
    if state["all_zero"]:
        return FiniteJointMassLaw(prior.atoms, prior.probabilities, log_evidence=0.), state
    statistics, costs = [], []
    for masses in prior.atoms:
        statistics.append(gaussian_joint_score_statistic(readings, masses,
            rank_rtol=state["rank_rtol"], rank_atol=state["rank_atol"])["statistic"])
        try:
            costs.append(math.fsum(gaussian_profile_cost(reading, masses) for reading in readings))
        except OverflowError as error:
            raise ValueError("total profile cost exceeds the floating-point range") from error
    statistics, costs = np.asarray(statistics), np.asarray(costs)
    minimum_statistic, minimum_cost = float(np.min(statistics)), float(np.min(costs))
    with np.errstate(over="ignore", invalid="ignore"):
        shifted = .5*(statistics-minimum_statistic)+.5*state["anchor_weight"]*(costs-minimum_cost)
    if np.any(~np.isfinite(shifted)):
        raise ValueError("anchored contrast exceeds the floating-point range")
    minimum_combined = float(np.min(shifted))
    log_weights = np.log(prior.probabilities)-(shifted-minimum_combined)
    shifted_evidence = float(special.logsumexp(log_weights))
    log_evidence = (shifted_evidence-minimum_combined-minimum_statistic/2
                    -state["anchor_weight"]*minimum_cost/2)
    if not math.isfinite(log_evidence):
        raise ValueError("log evidence exceeds the floating-point range")
    state.update(minimum_score_statistic=minimum_statistic, minimum_profile_cost=minimum_cost,
                 minimum_combined_contrast=minimum_combined)
    law = FiniteJointMassLaw(prior.atoms, np.exp(log_weights-shifted_evidence),
                             log_evidence=log_evidence)
    return law, state


def _positive_definite(value, name):
    value = np.asarray(value, dtype=float)
    if (value.ndim != 2 or value.shape[0] != value.shape[1] or not value.size
            or np.any(~np.isfinite(value))
            or not np.allclose(value, value.T, rtol=1e-12, atol=0.)):
        raise ValueError(f"{name} must be finite, symmetric and square")
    value = value/2+value.T/2
    try:
        np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be positive definite") from error
    return value


def anchored_local_precision(sensitivity, variability, *, anchor_weight):
    """Theoretical Gaussian comparison A B^-1 A + lambda*A.

    A and B are explicit regular-experiment benchmarks, not inputs supplied
    to joint_trace_anchor or learned latent populations. This function does
    not fit masses, regularize weak identification or calibrate finite atoms.
    """
    sensitivity = _positive_definite(sensitivity, "sensitivity")
    variability = _positive_definite(variability, "variability")
    weight = _nonnegative(anchor_weight, "anchor_weight")
    if sensitivity.shape != variability.shape:
        raise ValueError("sensitivity and variability must have matching shapes")
    with np.errstate(over="ignore", invalid="ignore"):
        score = sensitivity@np.linalg.solve(variability, sensitivity)
        anchor = weight*sensitivity
        total = score+anchor
    if np.any(~np.isfinite(total)):
        raise ValueError("local precision exceeds the floating-point range")
    score, anchor, total = (value/2+value.T/2 for value in (score, anchor, total))
    covariance = np.linalg.solve(total, np.eye(len(total)))
    relative = linalg.eigvalsh(anchor, score)
    if np.any(~np.isfinite(covariance)) or np.any(~np.isfinite(relative)):
        raise ValueError("local Gaussian benchmark exceeds the floating-point range")
    return dict(score_precision=score, anchor_precision=anchor, total_precision=total,
                covariance=covariance/2+covariance.T/2,
                relative_anchor_eigenvalues=relative)

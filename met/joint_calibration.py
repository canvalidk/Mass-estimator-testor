"""Gaussian joint-profile score moments and local covariance for T55.

True latent signals are explicit theorem parameters, never fitting inputs.
These references distinguish sensitivity from score variability when a free
signal is profiled out of each acquisition. They select no production law,
global transport, dependence rule or replacement for L19/property N.
"""
import numpy as np
from scipy import linalg

from met.joint_mass import JointReading


def _masses(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not values.size or np.any(~np.isfinite(values))
            or np.any(values <= 0)):
        raise ValueError("masses must be a nonempty finite positive vector")
    return values


def _geometry(reading, masses):
    if not isinstance(reading, JointReading):
        raise ValueError("each acquisition must be a JointReading")
    design = reading.design(masses)
    derivatives = reading.design_derivatives(masses)
    chol = np.linalg.cholesky(reading.covariance)
    basis = np.linalg.solve(chol, design)
    directions = np.array([np.linalg.solve(chol, item) for item in derivatives])
    rotation, full_triangular = np.linalg.qr(basis, mode="complete")
    latent_dimension = basis.shape[1]
    parallel, perpendicular = rotation[:, :latent_dimension], rotation[:, latent_dimension:]
    triangular = full_triangular[:latent_dimension]
    residual_directions = np.array([perpendicular.T@item for item in directions])
    transverse = np.array([perpendicular@item for item in residual_directions])
    right_inverse = np.linalg.solve(triangular, parallel.T)
    halves = transverse@right_inverse
    projector_derivatives = halves+halves.transpose(0, 2, 1)
    if (np.any(~np.isfinite(transverse))
            or np.any(~np.isfinite(projector_derivatives))):
        raise ValueError("joint graph derivatives exceed the floating-point range")
    return dict(chol=chol, basis=basis, directions=directions,
                parallel=parallel, perpendicular=perpendicular, triangular=triangular,
                residual_directions=residual_directions, transverse=transverse,
                projector_derivatives=projector_derivatives)


def _observed_event(reading, masses):
    geometry = _geometry(reading, masses)
    parallel, perpendicular = geometry["parallel"], geometry["perpendicular"]
    triangular = geometry["triangular"]
    whitened = np.linalg.solve(geometry["chol"], reading.data)
    latent = np.linalg.solve(triangular, parallel.T@whitened)
    residual = perpendicular.T@whitened
    with np.errstate(over="ignore", invalid="ignore"):
        fitted_directions = np.einsum("anq,q->na", geometry["directions"], latent)
        longitudinal = parallel.T@fitted_directions
        transverse = np.einsum("arq,q->ra", geometry["residual_directions"], latent)
        noise_rhs = np.einsum("arq,r->qa", geometry["residual_directions"], residual)
        normalized_noise = np.linalg.solve(triangular.T, noise_rhs)
        score = transverse.T@residual
        # Expand D'D-(a-p)'(a-p) in orthogonal coordinates. The p'p
        # terms cancel analytically, retaining small transverse curvature.
        sensitivity = (transverse.T@transverse-normalized_noise.T@normalized_noise
                       +normalized_noise.T@longitudinal+longitudinal.T@normalized_noise)
    return score, sensitivity


def gaussian_joint_profile_score(readings, masses):
    """Observed gradient of the joint Gaussian profile log factor.

    At each fixed candidate its component is z'P'_j z/2. Its sampling
    covariance is not generally the negative expected profile Hessian.
    """
    masses = _masses(masses)
    score = np.zeros(len(masses))
    for reading in readings:
        contribution, _ = _observed_event(reading, masses)
        score += contribution
    if np.any(~np.isfinite(score)):
        raise ValueError("profile score exceeds the floating-point range")
    return score


def gaussian_joint_observed_moments(readings, masses):
    """Observed profile sensitivity and score outer products, using sums.

    Evaluated at a consistent fitted mass in T55's repeated-block regime,
    these provide a consistent local sandwich covariance. This helper does
    not fit a mass or claim finite-sample calibration. The observed Hessian
    may be indefinite; no eigenvalue floor or regularization is inserted.
    An independent exactly zero acquisition contributes zero to every sum.
    """
    masses = _masses(masses)
    score_sum = np.zeros(len(masses))
    sensitivity = np.zeros((len(masses), len(masses)))
    variability = np.zeros_like(sensitivity)
    for reading in readings:
        score, observed_sensitivity = _observed_event(reading, masses)
        with np.errstate(over="ignore", invalid="ignore"):
            sensitivity += observed_sensitivity
            variability += np.outer(score, score)
            score_sum += score
    if (np.any(~np.isfinite(sensitivity)) or np.any(~np.isfinite(variability))
            or np.any(~np.isfinite(score_sum))):
        raise ValueError("observed profile moments exceed the floating-point range")
    return dict(score=score_sum, sensitivity=sensitivity, variability=variability)


def gaussian_joint_score_moments(readings, masses, latent_signals):
    """Exact score sensitivity A and variability B=A+G at the true mass.

    Each signal has that acquisition's declared latent dimension. Acquisition
    errors are independent; covariance within an acquisition is unrestricted.
    The signal is fixed separately for each event. No signal population or
    shared size across events is inferred. Data values do not enter A or B.
    """
    readings, signals = list(readings), list(latent_signals)
    masses = _masses(masses)
    if len(readings) != len(signals):
        raise ValueError("one explicit true latent signal is required per acquisition")
    sensitivity = np.zeros((len(masses), len(masses)))
    geometry = np.zeros_like(sensitivity)
    for reading, signal in zip(readings, signals):
        model = _geometry(reading, masses)
        transverse, derivatives = model["transverse"], model["projector_derivatives"]
        signal = np.asarray(signal, dtype=float)
        if signal.shape != (transverse.shape[-1],) or np.any(~np.isfinite(signal)):
            raise ValueError("each true signal must match its acquisition's latent dimension")
        with np.errstate(over="ignore", invalid="ignore"):
            signal_map = np.einsum("anq,q->na", transverse, signal)
            sensitivity += signal_map.T@signal_map
            geometry += .5*np.einsum("aij,bij->ab", derivatives, derivatives)
    variability = sensitivity+geometry
    if np.any(~np.isfinite(variability)):
        raise ValueError("score moments exceed the floating-point range")
    return dict(sensitivity=sensitivity, noise_geometry=geometry,
                variability=variability)


def _positive_definite(value, name):
    value = np.asarray(value, dtype=float)
    if (value.ndim != 2 or value.shape[0] != value.shape[1] or not value.size
            or np.any(~np.isfinite(value))
            or not np.allclose(value, value.T, rtol=1e-12, atol=0.)):
        raise ValueError(f"{name} must be finite, symmetric and square")
    value = .5*value+.5*value.T
    try:
        np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be positive definite") from error
    return value


def joint_covariance_benchmark(sensitivity, variability):
    """Raw and sandwich local Gaussian covariance per protocol block.

    For n repeats both covariance matrices are divided by n. The generalized
    eigenvalues of (B,A) are the Gaussian variance-inflation spectrum. They
    coincide exactly iff one scalar inflation supplies the full covariance.
    Their maximum is the least conservative inflation valid for all linear
    contrasts in this Gaussian benchmark, not a selected physical rule.
    These are covariances of Gaussian limits, not necessarily moments of a
    finite-data mass law; heavy-tailed input priors can keep those infinite.
    Positive sensitivity is required; singular weak-identification cases
    need a separate theorem and are not regularized here.
    """
    sensitivity = _positive_definite(sensitivity, "sensitivity")
    variability = _positive_definite(variability, "variability")
    if sensitivity.shape != variability.shape:
        raise ValueError("sensitivity and variability must have matching shape")
    raw = np.linalg.solve(sensitivity, np.eye(len(sensitivity)))
    with np.errstate(over="ignore", invalid="ignore"):
        sampling = raw@variability@raw
    if np.any(~np.isfinite(raw)) or np.any(~np.isfinite(sampling)):
        raise ValueError("covariance benchmark exceeds the floating-point range")
    spectrum = linalg.eigvalsh(variability, sensitivity)
    if np.any(~np.isfinite(spectrum)) or np.any(spectrum <= 0):
        raise ValueError("inflation spectrum exceeds the floating-point range")
    return dict(raw_covariance=.5*raw+.5*raw.T,
                sampling_covariance=.5*sampling+.5*sampling.T,
                inflation_eigenvalues=spectrum,
                minimal_conservative_inflation=float(spectrum[-1]))

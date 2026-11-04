"""Research score/profile law anchored by prior projector variance.

For each own-event whitened observation z, average the candidate graph
projector under the supplied joint prior. The state is the prior variance
of P(m)z, equivalently z'(Pbar-Pbar**2)z. It measures variation between
mass hypotheses, with no fitted mass, prior density or parameter metric.

This finite-atom reference does not certify continuous-prior quadrature,
finite-sample coverage or an adopted production rule. The prior coupling,
positive anchor prefactor and numerical rank controls are explicit.
"""
import math

import numpy as np
from scipy import integrate, special

from met.joint_mass import (FiniteJointMassLaw, JointReading,
                            gaussian_joint_geometry, gaussian_profile_cost)
from met.joint_score_law import gaussian_joint_score_statistic


def _controls(prefactor, relative, absolute):
    prefactor, relative, absolute = map(float, (prefactor, relative, absolute))
    if not math.isfinite(prefactor) or prefactor <= 0:
        raise ValueError("anchor_prefactor must be finite and positive")
    if (not math.isfinite(relative) or not math.isfinite(absolute)
            or relative < 0 or absolute < 0):
        raise ValueError("rank tolerances must be finite and nonnegative")
    return prefactor, relative, absolute


def _events(readings, prior):
    if not isinstance(prior, FiniteJointMassLaw):
        raise ValueError("prior must be a FiniteJointMassLaw")
    readings = tuple(readings)
    dimension = None
    for reading in readings:
        if (not isinstance(reading, JointReading)
                or reading.protocol != "common_acceleration"
                or reading.body_indices != tuple(range(prior.mass_count))):
            raise ValueError("canonical full-body common-acceleration events are required")
        if dimension is None:
            dimension = reading.dimension
        elif reading.dimension != dimension:
            raise ValueError("all events must have the same spatial dimension")
    return readings, dimension


def observed_mass_contrast_energy(readings, prior):
    """Observed sum of own-event prior variances of projected data.

    Projected-vector differences are centered on the first atom before
    averaging. This retains exact zero when all represented projections
    agree, including huge common projections whose weighted mean might
    otherwise introduce a spurious difference. The reference atom only
    changes this numerical centering; it changes no mathematical state.

    The supplied prior is averaged jointly, not replaced by a product of
    its marginals. Independent events are kept separate. In exact
    arithmetic each variance is between zero and ||z||**2/4.
    """
    readings, _ = _events(readings, prior)
    contributions = []
    for reading in readings:
        if not np.any(reading.data):
            continue
        projections = []
        for masses in prior.atoms:
            geometry = gaussian_joint_geometry(reading, masses)
            # Honor the geometry helper's recognized exact zero projection.
            projections.append(geometry["whitened_design"]@geometry["fitted_latent"])
        projections = np.asarray(projections)
        with np.errstate(over="ignore", invalid="ignore"):
            offsets = projections-projections[0]
            center = prior.probabilities@offsets
            differences = offsets-center
            squared = np.einsum("ai,ai->a", differences, differences)
        if np.any(~np.isfinite(squared)):
            raise ValueError("mass-contrast energy exceeds the floating-point range")
        contributions.append(math.fsum(float(p*v) for p, v in zip(prior.probabilities, squared)))
    try:
        result = math.fsum(contributions)
    except OverflowError as error:
        raise ValueError("total mass-contrast energy exceeds the floating-point range") from error
    if not math.isfinite(result):
        raise ValueError("total mass-contrast energy exceeds the floating-point range")
    return result


def joint_mass_contrast_anchor(readings, prior, *, anchor_prefactor, rank_rtol, rank_atol):
    """Return (finite joint law, observed prior-contrast anchor state).

    The kernel is exp[-S'C^+S/2-lambda*profile_cost/2], with
    lambda=kappa/sqrt(1+U). U uses the original supplied prior once, not
    the posterior or learned latent-signal distribution. Per-event SPD
    covariances are accepted; asymptotic results have separate conditions.
    """
    prefactor, relative, absolute = _controls(anchor_prefactor, rank_rtol, rank_atol)
    readings, dimension = _events(readings, prior)
    contrast = observed_mass_contrast_energy(readings, prior)
    weight = prefactor/math.sqrt(1+contrast)
    if not math.isfinite(weight) or weight <= 0:
        raise ValueError("anchor weight exceeds the supported floating-point range")
    all_zero = all(not np.any(reading.data) for reading in readings)
    state = dict(observed_mass_contrast_energy=contrast, anchor_weight=weight,
                 anchor_prefactor=prefactor, body_count=prior.mass_count,
                 spatial_dimension=dimension, rank_rtol=relative, rank_atol=absolute,
                 all_zero=all_zero, state_measure="supplied joint prior projector variance",
                 finite_scope="declared per-event SPD common-acceleration geometry",
                 asymptotic_scope="separate theorem conditions; finite atoms do not assert a Gaussian TV limit")
    if all_zero:
        return FiniteJointMassLaw(prior.atoms, prior.probabilities, log_evidence=0.), state
    statistics, costs = [], []
    for masses in prior.atoms:
        statistics.append(gaussian_joint_score_statistic(readings, masses,
            rank_rtol=relative, rank_atol=absolute)["statistic"])
        try:
            costs.append(math.fsum(gaussian_profile_cost(reading, masses) for reading in readings))
        except OverflowError as error:
            raise ValueError("total profile cost exceeds the floating-point range") from error
    statistics, costs = np.asarray(statistics), np.asarray(costs)
    if np.any(~np.isfinite(costs)):
        raise ValueError("total profile cost exceeds the floating-point range")
    minimum_statistic, minimum_cost = float(np.min(statistics)), float(np.min(costs))
    with np.errstate(over="ignore", invalid="ignore"):
        shifted = .5*(statistics-minimum_statistic)+.5*weight*(costs-minimum_cost)
    if np.any(~np.isfinite(shifted)):
        raise ValueError("anchored contrast exceeds the floating-point range")
    minimum_combined = float(np.min(shifted))
    log_weights = np.log(prior.probabilities)-(shifted-minimum_combined)
    shifted_evidence = float(special.logsumexp(log_weights))
    log_evidence = shifted_evidence-minimum_combined-minimum_statistic/2-weight*minimum_cost/2
    if not math.isfinite(log_evidence):
        raise ValueError("log evidence exceeds the floating-point range")
    state.update(minimum_score_statistic=minimum_statistic, minimum_profile_cost=minimum_cost,
                 minimum_combined_contrast=minimum_combined)
    return FiniteJointMassLaw(prior.atoms, np.exp(log_weights-shifted_evidence),
                             log_evidence=log_evidence), state


def scalar_prior_projector_contrast(covariance, *, prior_center, mass_bounds,
                                    epsabs, epsrel, limit):
    """Continuous scalar prior's projector-variance coefficient and diagnostics.

    The prior is the explicitly supplied half-Cauchy conditioned on fixed
    positive-domain bounds, permitting (0, infinity). For two-dimensional
    rank-one projectors, Pbar-Pbar**2=det(Pbar)*I, so a grouped whitened
    scatter F contributes coefficient*trace(F) to U. No data, count or
    latent signal enters this covariance/prior-only coefficient.

    Compute det(Pbar) as Var(sin(delta)**2)+Var(sin(delta)*cos(delta))
    about a fixed reference line. Normalize these coordinates before
    adaptive integration, so epsabs does not conceal a tiny contrast.
    Direct subtraction Pbar-Pbar**2 would lose a small positive variance
    when an instrument shape degenerates. Full-prior endpoint layers for
    arbitrarily extreme shapes are not certified by floating quadrature.
    The returned error diagnostics are numerical estimates, not certificates.
    """
    covariance = np.asarray(covariance, dtype=float)
    if (covariance.shape != (2, 2) or np.any(~np.isfinite(covariance))
            or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
        raise ValueError("a finite symmetric scalar two-channel covariance is required")
    center, epsabs, epsrel = map(float, (prior_center, epsabs, epsrel))
    bounds = np.asarray(mass_bounds, dtype=float)
    if (bounds.shape != (2,) or np.any(np.isnan(bounds)) or bounds[0] < 0
            or not bounds[0] < bounds[1] or not math.isfinite(bounds[0])
            or not math.isfinite(center) or center <= 0
            or not math.isfinite(epsabs) or epsabs <= 0
            or not math.isfinite(epsrel) or epsrel <= 0
            or type(limit) is not int or limit < 1):
        raise ValueError("positive prior/integration controls and increasing nonnegative mass bounds required")
    chol = np.linalg.cholesky(covariance)
    whitening = np.linalg.solve(chol, np.eye(2))
    # The graph projector ignores any common whitening scale. Remove it
    # before determinant/norm calculations to retain extreme precisions.
    whitening /= np.max(np.abs(whitening))
    start, stop = map(lambda mass: math.atan(mass/center), bounds)
    width, reference = stop-start, (start+stop)/2
    if width <= 0:
        raise ValueError("the supplied prior angle interval is numerically unresolved")
    reference_vector = whitening@np.array([center*math.sin(reference), math.cos(reference)])
    reference_norm = math.hypot(*reference_vector)
    determinant = float(np.linalg.det(whitening))

    def line_coordinates(fraction):
        theta = start+width*fraction
        vector = whitening@np.array([center*math.sin(theta), math.cos(theta)])
        if theta == math.pi/2:
            vector = whitening@np.array([center, 0.])
        length = math.hypot(*vector)
        sine = determinant*center*math.sin(width*(fraction-.5))/(length*reference_norm)
        cosine = float((vector/length)@(reference_vector/reference_norm))
        return sine, cosine

    def coordinates(fraction):
        sine, cosine = line_coordinates(fraction)
        return sine*sine, sine*cosine

    # Bound the coordinate magnitudes over the monotone candidate-line arc.
    # The extrema of sin²(delta) and sin(delta)cos(delta) occur at the
    # endpoints or their standard pi/2 and pi/4 stationary angles.
    endpoints = [line_coordinates(fraction) for fraction in (0., 1.)]
    deltas = [math.atan2(sine, cosine) for sine, cosine in endpoints]
    scales = [max(sine*sine for sine, _ in endpoints),
              max(abs(sine*cosine) for sine, cosine in endpoints)]
    if any(deltas[0] <= angle <= deltas[1] for angle in (-math.pi/2, math.pi/2)):
        scales[0] = 1.
    if any(deltas[0] <= angle <= deltas[1]
           for angle in (-3*math.pi/4, -math.pi/4, math.pi/4, 3*math.pi/4)):
        scales[1] = .5
    controls = dict(epsabs=epsabs, epsrel=epsrel, limit=limit)
    means, errors, variance_terms, variance_errors = [], [], [], []
    for column in (0, 1):
        scale = scales[column]
        if scale == 0:
            means.append(0.)
            errors.append(0.)
            variance_terms.append(0.)
            variance_errors.append(0.)
            continue
        normalized = lambda fraction: coordinates(fraction)[column]/scale
        value, error = integrate.quad(normalized, 0., 1., **controls)
        means.append(value*scale)
        errors.append(error*scale)
        variance, variance_error = integrate.quad(
            lambda fraction: (normalized(fraction)-value)**2, 0., 1., **controls)
        variance_terms.append(variance*scale*scale)
        variance_errors.append(variance_error*scale*scale)
    coefficient = math.fsum(variance_terms)
    if not math.isfinite(coefficient) or coefficient < 0:
        raise ValueError("projector contrast coefficient exceeds the supported numerical range")
    return dict(coefficient=coefficient, variance_integral_error=math.fsum(variance_errors),
                mean_integral_errors=errors, prior_center=center,
                normalized_coordinate_scales=scales,
                mass_bounds=bounds.tolist(), prior_angle_bounds=[start, stop],
                prior="half-Cauchy conditioned on the fixed supplied mass bounds")

"""T61 checks from physical residuals, not projector differentiation.

These are Gaussian moment, conditional-distribution and finite-law checks;
they do not assert finite-sample coverage or simulate a signal population.
"""
import math

import numpy as np
import pytest

from met.joint_calibration import (gaussian_joint_observed_moments,
                                   gaussian_joint_score_moments)
from met.joint_mass import FiniteJointMassLaw, common_acceleration
from met.joint_score_law import gaussian_joint_score_statistic
from met.joint_trace_anchor import anchored_local_precision, joint_trace_anchor


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _physical_geometry(masses, covariance, dimension):
    """Constraint F_j-m_j*A and its covariance, with no QR/projector API."""
    masses = np.asarray(masses)
    count = len(masses)
    constraint = np.column_stack([np.eye(count*dimension),
                                  -np.kron(masses[:, None], np.eye(dimension))])
    design = np.vstack([np.kron(masses[:, None], np.eye(dimension)),
                        np.eye(dimension)])
    residual_covariance = constraint@covariance@constraint.T
    precision = np.linalg.inv(residual_covariance)
    latent_covariance = np.linalg.inv(design.T@np.linalg.solve(covariance, design))
    blocks = precision.reshape(count, dimension, count, dimension).transpose(0, 2, 1, 3)
    return design, constraint, residual_covariance, latent_covariance, blocks


def _physical_moments(masses, covariance, signal):
    design, constraint, residual_covariance, gamma, blocks = _physical_geometry(
        masses, covariance, len(signal))
    signal = np.asarray(signal)
    sensitivity = np.einsum("r,abrs,s->ab", signal, blocks, signal)
    geometry = np.einsum("abrs,sr->ab", blocks, gamma)
    # C_jk is a quadratic form in the fitted acceleration. Its variance
    # contains only constant and signal-quadratic terms, never ||a||^4.
    variance_frobenius = 0.
    for block in blocks.reshape(-1, len(signal), len(signal)):
        symmetric = (block+block.T)/2
        variance_frobenius += (2*np.trace(symmetric@gamma@symmetric@gamma)
            +4*signal@symmetric@gamma@symmetric@signal)
    return dict(sensitivity=sensitivity, noise_geometry=geometry,
                variability=sensitivity+geometry,
                conditional_covariance_variance=variance_frobenius,
                signal_energy=signal@design.T@np.linalg.solve(covariance, design)@signal)


def _dense_covariance(size):
    factor = np.eye(size)+np.arange(size*size).reshape(size, size)/(3*size*size)
    return factor@factor.T


def _event(point, covariance, dimension):
    return common_acceleration(point[:-dimension].reshape(-1, dimension),
                               point[-dimension:], covariance)


@pytest.mark.parametrize("dimension,signal", [(1, [.7]), (2, [.7, -.4])])
def test_physical_residual_moments_include_the_exact_conditional_covariance(dimension, signal):
    masses = np.array([1.2, .8])
    covariance = _dense_covariance(3*dimension)
    signal = np.asarray(signal)
    reference = _physical_moments(masses, covariance, signal)
    declaration = _event(np.zeros(3*dimension), covariance, dimension)
    theoretical = gaussian_joint_score_moments([declaration], masses, [signal])
    for key in ("sensitivity", "noise_geometry", "variability"):
        np.testing.assert_allclose(theoretical[key], reference[key], rtol=3e-13, atol=3e-14)

    design = _physical_geometry(masses, covariance, dimension)[0]
    chol = np.linalg.cholesky(covariance)
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    size = 3*dimension
    mean_score = np.zeros(2)
    mean_hessian = np.zeros((2, 2))
    mean_covariance = np.zeros((2, 2))
    second_score = np.zeros((2, 2))
    covariance_variance = 0.
    # Three-node product Gaussian cubature integrates every polynomial
    # needed here exactly (degree at most four in each noise coordinate).
    for index in np.ndindex(*(3,)*size):
        probability = np.prod(weights[list(index)])/math.pi**(size/2)
        point = design@signal+chol@(math.sqrt(2)*nodes[list(index)])
        event = _event(point, covariance, dimension)
        observed = gaussian_joint_observed_moments([event], masses)
        conditional = gaussian_joint_score_statistic([event], masses, **RANK)
        mean_score += probability*conditional["score"]
        mean_hessian += probability*observed["sensitivity"]
        mean_covariance += probability*conditional["covariance"]
        second_score += probability*np.outer(conditional["score"], conditional["score"])
        covariance_variance += probability*np.sum(
            (conditional["covariance"]-reference["variability"])**2)
    np.testing.assert_allclose(mean_score, 0., atol=5e-15)
    np.testing.assert_allclose(mean_hessian, reference["sensitivity"], rtol=4e-13, atol=4e-14)
    np.testing.assert_allclose(mean_covariance, reference["variability"], rtol=4e-13, atol=4e-14)
    np.testing.assert_allclose(second_score, reference["variability"], rtol=4e-13, atol=4e-14)
    assert covariance_variance == pytest.approx(
        reference["conditional_covariance_variance"], rel=5e-13, abs=4e-14)


def test_different_event_geometries_have_one_exact_conditional_gaussian_score():
    masses = np.array([1.2, .8])
    covariances = [_dense_covariance(3),
                   np.array([[1.3, -.2, .1], [-.2, 2., .4], [.1, .4, .9]])]
    fitted_signals = [-.7, 1.1]
    geometries = [_physical_geometry(masses, covariance, 1) for covariance in covariances]
    expected_covariance = sum(a*a*g[4][:, :, 0, 0]
                              for a, g in zip(fitted_signals, geometries))
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    mean, second, statistic_moments = np.zeros(2), np.zeros((2, 2)), np.zeros(2)
    for index in np.ndindex(3, 3, 3, 3):
        probability = np.prod(weights[list(index)])/math.pi**2
        residual_noise = math.sqrt(2)*nodes[list(index)].reshape(2, 2)
        events, expected_score = [], np.zeros(2)
        for covariance, fitted, geometry, noise in zip(
                covariances, fitted_signals, geometries, residual_noise):
            design, constraint, residual_covariance, _, blocks = geometry
            residual = np.linalg.cholesky(residual_covariance)@noise
            # This reconstruction fixes the GLS fitted acceleration while
            # varying independent physical residual coordinates.
            point = (design[:, 0]*fitted
                     +covariance@constraint.T@np.linalg.solve(residual_covariance, residual))
            events.append(_event(point, covariance, 1))
            expected_score += fitted*blocks[:, :, 0, 0]@residual
        conditional = gaussian_joint_score_statistic(events, masses, **RANK)
        np.testing.assert_allclose(conditional["score"], expected_score, atol=8e-15)
        np.testing.assert_allclose(conditional["covariance"], expected_covariance,
                                   rtol=5e-14, atol=5e-15)
        assert conditional["rank"] == 2
        mean += probability*conditional["score"]
        second += probability*np.outer(conditional["score"], conditional["score"])
        statistic_moments += probability*np.array([
            conditional["statistic"], conditional["statistic"]**2])
    np.testing.assert_allclose(mean, 0., atol=4e-15)
    np.testing.assert_allclose(second, expected_covariance, rtol=3e-13, atol=4e-14)
    np.testing.assert_allclose(statistic_moments, [2., 8.], rtol=3e-13, atol=4e-14)


def _physical_kernel(events, masses, prefactor):
    score, conditional, cost, energy = np.zeros(2), np.zeros((2, 2)), 0., 0.
    for event in events:
        design, constraint, residual_covariance, gamma, blocks = _physical_geometry(
            masses, event.covariance, 1)
        fitted = (gamma@design.T@np.linalg.solve(event.covariance, event.data))[0]
        residual = constraint@event.data
        precision = blocks[:, :, 0, 0]
        score += fitted*precision@residual
        conditional += fitted*fitted*precision
        cost += residual@precision@residual
        energy += event.data@np.linalg.solve(event.covariance, event.data)
    statistic = score@np.linalg.solve(conditional, score)
    return -.5*(statistic+prefactor/math.sqrt(1+energy)*cost)


def test_finite_heterogeneous_law_matches_physical_residual_likelihood_ratios():
    events = [_event(np.array([2.1, .4, .7]), _dense_covariance(3), 1),
              _event(np.array([-.3, 1.4, -.2]), np.diag([.7, 2.3, 1.1]), 1)]
    prior = FiniteJointMassLaw([[.4, .8], [1.2, .8], [1.3, 2.], [3., 1.]],
                               [.1, .2, .3, .4])
    law, state = joint_trace_anchor(events, prior, anchor_prefactor=.7, **RANK)
    log_factors = np.array([_physical_kernel(events, masses, .7) for masses in prior.atoms])
    expected = prior.probabilities*np.exp(log_factors-log_factors.max())
    expected /= expected.sum()
    np.testing.assert_allclose(law.probabilities, expected, rtol=3e-14, atol=3e-15)
    permuted, changed = joint_trace_anchor(events[::-1], prior, anchor_prefactor=.7, **RANK)
    np.testing.assert_allclose(permuted.probabilities, law.probabilities, atol=3e-15)
    assert changed["observed_whitened_trace"] == state["observed_whitened_trace"]


def test_arbitrary_event_noise_scales_cancel_with_their_actual_signal_and_data_scales():
    masses = np.array([1.2, .8])
    covariance = _dense_covariance(6)
    signal = np.array([.7, -.4])
    reference = _physical_moments(masses, covariance, signal)
    point = np.array([2.1, .4, .7, -.2, .5, .9])
    original = _event(point, covariance, 2)
    prior = FiniteJointMassLaw([[.4, .8], [1.2, .8], [1.3, 2.], [3., 1.]],
                               [.1, .2, .3, .4])
    law, state = joint_trace_anchor([original, original], prior, anchor_prefactor=.7, **RANK)
    events, signals = [], []
    for multiplier in (1e-80, 1e80):
        scaled_covariance = multiplier*covariance
        scaled_signal = math.sqrt(multiplier)*signal
        moments = _physical_moments(masses, scaled_covariance, scaled_signal)
        for key in ("sensitivity", "noise_geometry", "variability",
                    "conditional_covariance_variance", "signal_energy"):
            np.testing.assert_allclose(moments[key], reference[key], rtol=8e-14, atol=2e-14)
        events.append(_event(math.sqrt(multiplier)*point, scaled_covariance, 2))
        signals.append(scaled_signal)
    observed = gaussian_joint_score_moments(events, masses, signals)
    for key in ("sensitivity", "noise_geometry", "variability"):
        np.testing.assert_allclose(observed[key], 2*reference[key], rtol=8e-14, atol=3e-14)
    rescaled, changed = joint_trace_anchor(events, prior, anchor_prefactor=.7, **RANK)
    np.testing.assert_allclose(rescaled.probabilities, law.probabilities, rtol=8e-14, atol=3e-15)
    assert changed["observed_whitened_trace"] == pytest.approx(
        state["observed_whitened_trace"], rel=5e-14)
    assert changed["anchor_weight"] == pytest.approx(state["anchor_weight"], rel=5e-14)


def test_changing_geometry_and_concentrated_signal_preserve_the_anchor_rate_window():
    masses = np.array([1.2, .8])
    distortions, relative_errors, separation = [], [], []
    for count in (16, 256, 4096):
        # A bounded, nonrepeating covariance history; signals here form a
        # deterministic triangular-array check, concentrated in one event.
        target_energy = count**.9
        total_a, total_g = np.zeros((2, 2)), np.zeros((2, 2))
        covariance_variance, energy = 0., 0.
        for i in range(count):
            angle = (i+1)**.7
            vector = np.array([math.cos(angle), math.sin(angle),
                               math.cos(2*angle), math.sin(2*angle), .3, -.2])
            vector /= np.linalg.norm(vector)
            covariance = np.eye(6)+np.outer(vector, vector)
            signal = np.array([math.sqrt(target_energy), 0.]) if i == 0 else np.zeros(2)
            moments = _physical_moments(masses, covariance, signal)
            total_a += moments["sensitivity"]
            total_g += moments["noise_geometry"]
            covariance_variance += moments["conditional_covariance_variance"]
            energy += moments["signal_energy"]
        total_b = total_a+total_g
        weight = 1/math.sqrt(1+6*count+energy)
        comparison = anchored_local_precision(total_a, total_b, anchor_weight=weight)
        distortions.append(comparison["relative_anchor_eigenvalues"][-1])
        relative_errors.append(math.sqrt(covariance_variance)/np.linalg.eigvalsh(total_b)[0])
        separation.append(weight*energy)
        # Uniformly bounded geometry gives information comparable to total
        # signal energy, even though a single event has all the signal.
        assert np.linalg.eigvalsh(total_a/energy)[0] > .02
        assert np.linalg.eigvalsh(total_a/energy)[-1] < 1.
    assert separation[0] < separation[1] < separation[2]
    assert distortions[0] > distortions[1] > distortions[2]
    assert relative_errors[0] > relative_errors[1] > relative_errors[2]
    assert distortions[-1] < .2
    assert relative_errors[-1] < .05


def test_conditional_covariance_concentration_needs_no_fourth_signal_moment_bound():
    masses, covariance = np.array([1.2, .8]), _dense_covariance(6)
    relative_errors = []
    for amplitude in (10., 100., 1000.):
        moments = _physical_moments(masses, covariance, [amplitude, -.4*amplitude])
        relative_errors.append(math.sqrt(moments["conditional_covariance_variance"])
                               /np.linalg.eigvalsh(moments["variability"])[0])
    assert relative_errors[1] < relative_errors[0]/9
    assert relative_errors[2] < relative_errors[1]/9
    assert relative_errors[2] < .01

"""E6: normalized future events, retained state and a common-null obstruction.

Independent checks of explore/e6-predictive-state.md. No estimator is added.
"""

import itertools

import numpy as np
import pytest
from scipy.special import ndtr
from scipy.stats import norm

from met.declared_structure import (calibrated_discrete_joint,
                                    compress_gaussian_mean, direction_log_factors)


def _outside_probability(mean, threshold=1.0, sigma=1.0):
    mean = np.asarray(mean)
    return ndtr((mean-threshold)/sigma)+ndtr((-mean-threshold)/sigma)


def _mass_neutral_states():
    masses, directions = np.array([0.5, 1.0, 2.0]), np.eye(2)
    prior, calibration = np.array([1, 2, 1])/4, np.array([1, 3])/4
    before = np.outer(prior, calibration)
    history = direction_log_factors([[1, 0], [1, 0]], [[1, 0], [-1, 0]],
                                    masses, directions)
    after = calibrated_discrete_joint(prior, calibration, history)
    return before, after


@pytest.mark.parametrize("amplitude", [-3.0, -1.0, -0.1, 0.0, 0.1, 1.0, 3.0])
def test_mass_neutral_history_changes_a_normalized_future_raw_event(amplitude):
    before, after = _mass_neutral_states()
    # Y=(F,A)|m,n,t ~ N((m*t*n,t*n), I4), with t fixed, not random.
    # The event |A_x|>1 and its complement are exhaustive for every t.
    event = np.tile(_outside_probability([amplitude, 0.0]), (3, 1))
    assert np.all((0 <= event) & (event <= 1))
    np.testing.assert_allclose(event+(1-event), 1)
    np.testing.assert_allclose(before.sum(axis=1), after.sum(axis=1), atol=1e-15)
    change = np.sum((after-before)*event)
    expected = (np.e/(np.e+3)-0.25)*(_outside_probability(amplitude)-_outside_probability(0))
    assert change == pytest.approx(expected, abs=1e-15)
    if amplitude == 0:
        assert change == pytest.approx(0, abs=1e-15)
    else:
        assert change > 0


def test_raw_event_witness_values_and_amplitude_dependence():
    before, after = _mass_neutral_states()
    event = np.tile(_outside_probability([1.0, 0.0]), (3, 1))
    assert np.sum(before*event) == pytest.approx(0.3686704138842304)
    assert np.sum(after*event) == pytest.approx(0.41496970231134894)
    # Even a full Q does not supply one forecast when the next t is unspecified.
    conditional = [np.sum(after*np.tile(_outside_probability([t, 0]), (3, 1)))
                   for t in [0.0, 0.5, 1.0, 2.0]]
    assert np.all(np.diff(conditional) > 0)


@pytest.mark.parametrize("mass", [0.3, 1.0, 5.0])
def test_common_null_event_limit_and_gaussian_tv_contraction(mass):
    covariance = np.array([[2.0, 0.4], [0.4, 1.0]])
    probe = np.array([0.8, -0.3])
    threshold = 0.4
    sigma = np.sqrt(probe @ covariance @ probe)
    p0 = norm.sf(threshold/sigma)
    differences, tv_bounds = [], []
    for amplitude in [1.0, 0.1, 0.01, 0.001, 0.0]:
        mean = amplitude*np.array([mass, 1.0])
        # One fixed event, not a residual map chosen using the candidate.
        probability = norm.sf((threshold-probe @ mean)/sigma)
        distance = np.sqrt(mean @ np.linalg.solve(covariance, mean))
        tv = 2*ndtr(distance/2)-1  # exact TV of equal-covariance Gaussian shifts
        assert abs(probability-p0) <= tv+1e-15
        differences.append(abs(probability-p0))
        tv_bounds.append(tv)
    assert differences[-1] == 0 and tv_bounds[-1] == 0
    assert np.all(np.diff(differences) <= 0)
    assert np.all(np.diff(tv_bounds) <= 0)


def test_candidate_residual_is_a_pivot_only_at_the_matching_candidate():
    true_mass = 2.0
    candidate_masses = [0.5, 1.0, 2.0, 4.0]
    at_zero, at_one = [], []
    for candidate in candidate_masses:
        normal = np.array([1.0, -candidate])/np.hypot(1, candidate)
        # The candidate-dependent map is a different event for each candidate.
        mean_at_one = normal @ np.array([true_mass, 1.0])
        at_zero.append(_outside_probability(0))
        at_one.append(_outside_probability(mean_at_one))
    assert at_one[2] == pytest.approx(at_zero[2])
    assert all(at_one[i] > at_zero[i] for i in [0, 1, 3])


def test_one_step_prediction_is_not_closed_under_conditioning():
    p = np.array([0.2, 0.5, 0.8])
    first, second = np.array([0, 1, 0]), np.array([0.5, 0, 0.5])
    assert first @ p == pytest.approx(second @ p)
    assert first @ (p*p) == pytest.approx(0.25)
    assert second @ (p*p) == pytest.approx(0.34)
    assert (first*p) @ p/(first @ p) == pytest.approx(0.5)
    assert (second*p) @ p/(second @ p) == pytest.approx(0.68)
    assert np.linalg.matrix_rank(np.array([np.ones(3), p, p*p])) == 3
    # Fresh noise repeats; redrawing the hidden bias would give a different law.
    assert (second @ p)**2 == pytest.approx(0.25)


def test_same_both_marginals_and_next_calibration_law_can_hide_mass_update():
    first = np.array([[0.4, 0.1], [0.1, 0.4]])
    second = np.array([[0.1, 0.4], [0.4, 0.1]])
    # Independent future fixture readout Z|[n] ~ N(n_x**2,1).
    # The event Z>1/2 has no direct m-dependence.
    effect = ndtr(np.array([0.5, -0.5]))
    np.testing.assert_array_equal(first.sum(axis=0), second.sum(axis=0))
    np.testing.assert_array_equal(first.sum(axis=1), second.sum(axis=1))
    p1, p2 = np.sum(first*effect), np.sum(second*effect)
    assert p1 == pytest.approx(0.5) and p2 == pytest.approx(0.5)
    post1, post2 = (first*effect).sum(axis=1)/p1, (second*effect).sum(axis=1)/p2
    np.testing.assert_allclose(post1, post2[::-1], atol=1e-15)
    assert abs(post1[0]-post2[0]) > 0.2


@pytest.mark.parametrize("mass_count", [1, 2, 4])
def test_two_axis_calibration_and_mass_readout_close_to_the_full_joint(mass_count):
    a, b = ndtr(0.5), ndtr(-0.5)
    mass_indicators = np.repeat(np.eye(mass_count), 2, axis=1)
    effect = np.tile([a, b], mass_count)
    closed = np.vstack([mass_indicators, mass_indicators*effect])
    assert np.linalg.matrix_rank(closed) == 2*mass_count
    np.testing.assert_allclose(effect**2, effect-a*b, atol=1e-15)
    rng = np.random.default_rng(6510+mass_count)
    joint = rng.uniform(0.1, 1.0, (mass_count, 2))
    joint /= joint.sum()
    s = joint.sum(axis=1)
    t = (joint*np.array([a, b])).sum(axis=1)
    recovered = np.column_stack([(t-b*s)/(a-b), (a*s-t)/(a-b)])
    np.testing.assert_allclose(recovered, joint, atol=2e-16)


def test_separate_one_shot_sensors_do_not_preserve_hidden_association():
    # Hidden (bit1,bit2), with separately available one-shot noiseless sensors.
    bits = np.array(list(itertools.product([0, 1], repeat=2)))
    first, second = np.array([0.5, 0, 0, 0.5]), np.array([0, 0.5, 0.5, 0])
    for bit in bits.T:
        assert first @ bit == pytest.approx(second @ bit)
    joint_event = bits[:, 0]*bits[:, 1]
    assert first @ joint_event == pytest.approx(0.5)
    assert second @ joint_event == 0
    assert np.linalg.matrix_rank(np.vstack([np.ones(4), bits.T])) == 3
    assert np.linalg.matrix_rank(np.vstack([np.ones(4), bits.T, joint_event])) == 4


def test_a_true_duplicate_axis_can_be_quotiented_for_every_calibration_history():
    # The first two axis labels have identical experiments and mass readouts.
    effect = np.tile(ndtr([0.5, 0.5, -0.5]), (2, 1))
    first = np.array([[0.1, 0.2, 0.1], [0.15, 0.25, 0.2]])
    second = np.array([[0.2, 0.1, 0.1], [0.35, 0.05, 0.2]])
    for yes in range(5):
        for no in range(5):
            word = effect**yes*(1-effect)**no
            numerator1, numerator2 = (first*word).sum(axis=1), (second*word).sum(axis=1)
            np.testing.assert_allclose(numerator1, numerator2, atol=1e-16)
            np.testing.assert_allclose(numerator1/numerator1.sum(),
                                       numerator2/numerator2.sum(), atol=1e-15)


def test_mass_only_future_operations_admit_the_mass_marginal_as_exact_state():
    first = np.array([[0.05, 0.2], [0.3, 0.45]])
    second = np.array([[0.15, 0.1], [0.05, 0.7]])
    effect = np.repeat(np.array([[0.3], [0.8]]), 2, axis=1)
    np.testing.assert_allclose(first.sum(axis=1), second.sum(axis=1), atol=1e-15)
    assert not np.allclose(first.sum(axis=0), second.sum(axis=0))
    for yes in range(5):
        for no in range(5):
            word = effect**yes*(1-effect)**no
            numerator1, numerator2 = (first*word).sum(axis=1), (second*word).sum(axis=1)
            np.testing.assert_allclose(numerator1, numerator2, atol=1e-15)
            np.testing.assert_allclose(numerator1/numerator1.sum(),
                                       numerator2/numerator2.sum(), atol=1e-15)


@pytest.mark.parametrize("amplitude", [1.0, 2.0])
def test_normalized_symmetric_gaussian_events_separate_the_finite_joint(amplitude):
    masses, axes = [0.5, 1.0, 2.0], np.eye(2)
    means = np.array([np.r_[mass*n, n] for mass in masses for n in axes])*amplitude
    rng = np.random.default_rng(6518)
    probes = rng.normal(size=(40, 4))
    probes /= np.linalg.norm(probes, axis=1)[:, None]
    projected = probes @ means.T
    effects = np.vstack([_outside_probability(projected, threshold=c)
                         for c in [0.5, 1.0, 2.0]])
    assert np.all((0 <= effects) & (effects <= 1))
    assert np.linalg.matrix_rank(np.vstack([np.ones(6), effects])) == 6
    # These effects refer to |probe.Y|>c, so reversing axis orientation is invisible.
    np.testing.assert_allclose(effects, np.vstack([
        _outside_probability(-projected, threshold=c) for c in [0.5, 1.0, 2.0]]))
    zero_effects = np.vstack([_outside_probability(np.zeros_like(projected), threshold=c)
                             for c in [0.5, 1.0, 2.0]])
    assert np.linalg.matrix_rank(np.vstack([np.ones(6), zero_effects])) == 1


def test_gaussian_symmetric_density_gram_detects_only_genuine_axis_duplicates():
    means = np.array([np.r_[m*n, n] for m in [0.5, 1, 2] for n in np.eye(2)])*2
    minus = means[:, None]-means[None, :]
    plus = means[:, None]+means[None, :]
    gram = 0.5*(np.exp(-np.sum(minus**2, axis=2)/4)
                + np.exp(-np.sum(plus**2, axis=2)/4))
    assert np.linalg.eigvalsh(gram).min() > 1e-3
    # Negating one representative leaves its entire symmetric density unchanged.
    antipodal = means.copy()
    antipodal[1] *= -1
    cross_minus = means[:, None]-antipodal[None, :]
    cross_plus = means[:, None]+antipodal[None, :]
    cross_gram = 0.5*(np.exp(-np.sum(cross_minus**2, axis=2)/4)
                      + np.exp(-np.sum(cross_plus**2, axis=2)/4))
    np.testing.assert_allclose(cross_gram, gram, atol=1e-15)


def test_fixed_linear_gls_has_a_parameter_independent_raw_data_regenerator():
    rng = np.random.default_rng(6520)
    design = rng.normal(size=(6, 2))
    raw = rng.normal(size=(6, 6))
    covariance = raw @ raw.T+np.eye(6)
    reduced_cov = np.linalg.inv(design.T @ np.linalg.solve(covariance, design))
    regression = reduced_cov @ np.linalg.solve(covariance, design).T
    residual_cov = covariance-design @ reduced_cov @ design.T
    assert np.linalg.eigvalsh(residual_cov).min() >= -1e-13
    np.testing.assert_allclose(regression @ residual_cov, 0, atol=2e-14)
    for theta in [np.zeros(2), np.array([1, 2]), np.array([-3, 0.5])]:
        mean = design @ theta
        estimate, actual_cov, _ = compress_gaussian_mean(mean, covariance, design)
        np.testing.assert_allclose(estimate, theta, atol=2e-15)
        np.testing.assert_allclose(actual_cov, reduced_cov, atol=2e-14)
        # S~N(theta,V); Y*=design*S+independent N(0,residual_cov).
        # Same mean and covariance as the raw normalized Gaussian experiment.
        np.testing.assert_allclose(design @ estimate, mean, atol=5e-15)
        np.testing.assert_allclose(design @ actual_cov @ design.T+residual_cov,
                                   covariance, atol=3e-14)

"""O2: physical composition, shared evidence and Gaussian coarsening.

These deterministic proof checks use an explicitly supplied joint model, not a
new ours29 readout or an adopted signal-size law. All instrument covariances are
positive definite. The finite priors are counterexamples to information-loss
claims; their reachability through the current ours29 adapter is not asserted.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit, logsumexp
from scipy.stats import multivariate_normal, norm

from met.full import Reading, isotropic


def _sum_and_contrast(dimension, weight):
    eye, zero = np.eye(dimension), np.zeros((dimension, dimension))
    # (F1,F2,A) -> (S,A,C), with S=F1+F2 and C=F1-weight*S.
    return np.block([[eye, eye, zero], [zero, zero, eye],
                     [(1 - weight) * eye, -weight * eye, zero]])


def _means(states):
    """Each state is (m1,m2,a), sharing the same positive scalar acceleration."""
    states = np.asarray(states, float)
    return np.column_stack((states[:, 0] * states[:, 2],
                            states[:, 1] * states[:, 2], states[:, 2]))


def _posterior(states, prior, observation, covariance, coarse=False):
    means = _means(states)
    if coarse:
        linear = _sum_and_contrast(1, .5)[:2]
        means = means @ linear.T
        observation = linear @ observation
        covariance = linear @ covariance @ linear.T
    log_weights = np.log(prior) + multivariate_normal.logpdf(
        observation - means, mean=np.zeros(means.shape[1]), cov=covariance)
    return np.exp(log_weights - logsumexp(log_weights))


def _readouts(states, probability):
    m1, m2, acceleration = np.asarray(states, float).T
    denominator = probability @ np.abs(acceleration)
    return np.array([probability @ (m * np.abs(acceleration)) / denominator
                     for m in (m1, m2, m1 + m2)])


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_shared_observation_copies_have_forced_cross_covariance(dimension):
    eye, zero = np.eye(dimension), np.zeros((dimension, dimension))
    factor = np.eye(3 * dimension) + .05 * np.tril(np.ones((3 * dimension, 3 * dimension)), -1)
    source_covariance = factor @ factor.T
    copy = np.block([[eye, zero, zero], [zero, zero, eye],
                     [zero, eye, zero], [zero, zero, eye]])
    copied_covariance = copy @ source_covariance @ copy.T
    difference = np.block([zero, eye, zero, -eye])
    assert np.linalg.matrix_rank(copied_covariance) == 3 * dimension
    np.testing.assert_allclose(difference @ copied_covariance @ difference.T, 0, atol=1e-14)
    np.testing.assert_array_equal(copied_covariance[dimension:2 * dimension, 3 * dimension:],
                                  source_covariance[2 * dimension:, 2 * dimension:])

    # Each standalone pair is a legitimate Reading, but their independent product
    # cannot represent a single acceleration observation copied twice.
    pair_covariances = []
    for force in (0, 1):
        indices = list(range(force * dimension, (force + 1) * dimension)) + list(range(2 * dimension, 3 * dimension))
        covariance = source_covariance[np.ix_(indices, indices)]
        Reading(np.zeros(dimension), np.zeros(dimension), covariance)
        pair_covariances.append(covariance)
    independent = np.block([[pair_covariances[0], np.zeros((2 * dimension, 2 * dimension))],
                            [np.zeros((2 * dimension, 2 * dimension)), pair_covariances[1]]])
    np.testing.assert_allclose(difference @ independent @ difference.T,
                               2 * source_covariance[2 * dimension:, 2 * dimension:])


def test_current_separate_clouds_cannot_preserve_one_shared_acceleration():
    variances, means, precisions = [], [], []
    for mass, force in ((1., 1.), (2., 2.)):
        reading = isotropic([force], [1.], 1., 1.)
        pre = reading.pre()
        precision = mass * mass * pre["k_ff"] + mass * pre["k_s"] + pre["k_aa"]
        means.append(float(np.linalg.solve(precision, mass * pre["u"] + pre["v"])[0]))
        variances.append(float(np.linalg.inv(precision)[0, 0]))
        precisions.append(float(precision[0, 0]))
    np.testing.assert_allclose(means, [1, 1])
    np.testing.assert_allclose(variances, [1 / 2, 1 / 5])
    # A shared random variable cannot retain two different marginal variances.
    # The joint likelihood counts the acceleration instrument once.
    design, data = np.array([1., 2., 1.]), np.array([1., 2., 1.])
    joint_precision = design @ design
    assert (design @ data) / joint_precision == pytest.approx(1)
    assert 1 / joint_precision == pytest.approx(1 / 6)
    assert 1 / sum(precisions) == pytest.approx(1 / 7)  # duplicated A likelihood


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_correlated_gaussian_factors_into_coarse_data_and_innovation(dimension):
    factor = np.eye(3 * dimension) + .07 * np.tril(np.ones((3 * dimension, 3 * dimension)), -1)
    covariance = factor @ factor.T
    acceleration = np.linspace(.4, 1.1, dimension)
    mean = np.concatenate((.7 * acceleration, 1.9 * acceleration, acceleration))
    observation = mean + np.linspace(-.3, .5, 3 * dimension)
    transform = _sum_and_contrast(dimension, 1 / 3)
    assert np.linalg.slogdet(transform)[1] == pytest.approx(0, abs=1e-14)
    gamma = transform @ covariance @ transform.T
    transformed, predicted = transform @ observation, transform @ mean
    size = 2 * dimension
    coarse_covariance = gamma[:size, :size]
    regression = np.linalg.solve(coarse_covariance, gamma[:size, size:]).T
    conditional_covariance = gamma[size:, size:] - regression @ gamma[:size, size:]
    innovation = transformed[size:] - regression @ transformed[:size]
    innovation_mean = predicted[size:] - regression @ predicted[:size]
    direct = multivariate_normal.logpdf(observation, mean=mean, cov=covariance)
    factored = multivariate_normal.logpdf(transformed[:size], mean=predicted[:size], cov=coarse_covariance)
    factored += multivariate_normal.logpdf(innovation, mean=innovation_mean, cov=conditional_covariance)
    assert factored == pytest.approx(direct, abs=1e-12)
    np.testing.assert_allclose(gamma[size:, :size] - regression @ coarse_covariance, 0, atol=1e-14)


@pytest.mark.parametrize("contrast", [-2., 0., .5, 1., 3.])
def test_positive_noise_counterexample_preserves_composition_but_not_refitting(contrast):
    states, prior = np.array([[1., 1., 1.], [3., 1., 1.]]), np.array([.5, .5])
    covariance = np.eye(3)  # acceleration SD is 1, not zero
    observation = np.array([1.5 + contrast, 1.5 - contrast, 1.])
    fine = _posterior(states, prior, observation, covariance)
    coarse = _posterior(states, prior, observation, covariance, coarse=True)
    np.testing.assert_allclose(fine, [1 - expit(2 * contrast - 1), expit(2 * contrast - 1)], atol=1e-14)
    np.testing.assert_allclose(coarse, [.5, .5], atol=1e-14)
    first, second, total = _readouts(states, fine)
    assert total == pytest.approx(first + second, abs=1e-14)
    assert total == pytest.approx(2 + 2 * expit(2 * contrast - 1), abs=1e-14)
    assert _readouts(states, coarse)[2] == pytest.approx(3)
    if contrast != .5:
        assert abs(total - 3) > .4


def test_same_coarse_measurement_and_prior_can_have_different_fine_points():
    states, prior, covariance = np.array([[1., 1., 1.], [3., 1., 1.]]), [.5, .5], np.eye(3)
    data = [np.array([1.5, 1.5, 1.]), np.array([2.5, .5, 1.])]
    linear = _sum_and_contrast(1, .5)[:2]
    np.testing.assert_array_equal(linear @ data[0], linear @ data[1])
    np.testing.assert_array_equal(linear @ covariance @ linear.T, np.diag([2, 1]))
    points = [_readouts(states, _posterior(states, prior, y, covariance))[2] for y in data]
    np.testing.assert_allclose(points, [2.5378828427399904, 3.4621171572600096], rtol=1e-14)


@pytest.mark.parametrize("contrast", [0., 1.])
def test_shared_force_dual_composes_inverse_mass_on_the_same_joint_law(contrast):
    # Now the coefficients are inverse masses and the shared channel is force:
    # states=(mu1,mu2,F*), source data=(a1,a2,F). All three SDs are still 1.
    states, prior = np.array([[1., 1., 1.], [3., 1., 1.]]), [.5, .5]
    observation, covariance = np.array([1.5 + contrast, 1.5 - contrast, 1.]), np.eye(3)
    probability = _posterior(states, prior, observation, covariance)
    mu1, mu2, shared_force = states.T
    mean_force = probability @ shared_force
    masses = np.array([mean_force / (probability @ (mu * shared_force))
                       for mu in (mu1, mu2, mu1 + mu2)])
    assert 1 / masses[2] == pytest.approx(1 / masses[0] + 1 / masses[1], abs=1e-14)
    assert 1 / masses[2] == pytest.approx(2 + 2 * expit(2 * contrast - 1), abs=1e-14)
    coarse = _posterior(states, prior, observation, covariance, coarse=True)
    coarse_mass = (coarse @ shared_force) / (coarse @ ((mu1 + mu2) * shared_force))
    assert coarse_mass == pytest.approx(1 / 3)
    assert abs(coarse_mass - masses[2]) > .04


def test_omitted_contrast_factor_is_exact_with_propagated_joint_prior():
    states = np.array([[1., 1., 1.], [.5, 1.5, 1.], [3., 1., 1.5], [1., 3., 1.5]])
    prior, observation, covariance = np.array([.1, .3, .2, .4]), np.array([2., 1., .8]), np.eye(3)
    fine = _posterior(states, prior, observation, covariance)
    coarse = _posterior(states, prior, observation, covariance, coarse=True)
    fine_q = np.array([fine[:2].sum(), fine[2:].sum()])
    coarse_q = np.array([coarse[:2].sum(), coarse[2:].sum()])
    contrast = observation[0] - .5 * (observation[0] + observation[1])
    contrast_mean = (states[:, 0] - .5 * states[:, :2].sum(axis=1)) * states[:, 2]
    split_likelihood = norm.pdf(contrast, loc=contrast_mean, scale=math.sqrt(.5))
    omitted = np.array([prior[:2] @ split_likelihood[:2] / prior[:2].sum(),
                        prior[2:] @ split_likelihood[2:] / prior[2:].sum()])
    predicted = coarse_q * omitted
    predicted /= predicted.sum()
    np.testing.assert_allclose(fine_q, predicted, rtol=1e-14)
    assert np.max(np.abs(fine_q - coarse_q)) > .01


@pytest.mark.parametrize("force_variances", [(1., 1.), (1., 4.), (4., 1.)])
def test_variance_matched_split_makes_discarded_contrast_ancillary(force_variances):
    v1, v2 = force_variances
    weight = v1 / (v1 + v2)
    total_mass, acceleration = np.array([.5, 1., 2., 4.]), np.array([.4, 1., 1.5, 2.])
    prior, covariance = np.array([.1, .2, .3, .4]), np.diag([v1, v2, .7])
    states = np.column_stack((weight * total_mass, (1 - weight) * total_mass, acceleration))
    for contrast in (-1., 0., 2.):
        observation = np.array([weight * 1.3 + contrast, (1 - weight) * 1.3 - contrast, .9])
        fine = _posterior(states, prior, observation, covariance)
        coarse = _posterior(states, prior, observation, covariance, coarse=True)
        np.testing.assert_allclose(fine, coarse, rtol=1e-12, atol=1e-14)
    # A different fixed physical split is not generally ancillary to this instrument.
    other_weight = weight + .1
    states[:, 0], states[:, 1] = other_weight * total_mass, (1 - other_weight) * total_mass
    observation = np.array([weight * 1.3, (1 - weight) * 1.3, .9])
    assert np.max(np.abs(_posterior(states, prior, observation, covariance)
                         - _posterior(states, prior, observation, covariance, coarse=True))) > 1e-4


def test_tower_uses_size_weighted_fine_readouts_and_own_predictive_law():
    # Coarse observation S=3.5,A=1 lies midway between both candidate means.
    # Its posterior weights are 1/2, with a=1/2 or 3/2 and total mass 2 or 4.
    states = np.array([[1., 1., .5], [3., 1., 1.5]])
    observation, covariance = np.array([1.75, 1.75, 1.]), np.eye(3)
    coarse = _posterior(states, [.5, .5], observation, covariance, coarse=True)
    np.testing.assert_allclose(coarse, [.5, .5])
    assert _readouts(states, coarse)[2] == pytest.approx(3.5)

    def integrands(contrast):
        second = expit(3 * contrast - 2.25)
        acceleration_size = .5 + second
        readout = (1 + 5 * second) / acceleration_size
        predictive = .5 * (norm.pdf(contrast, scale=math.sqrt(.5))
                            + norm.pdf(contrast, loc=1.5, scale=math.sqrt(.5)))
        return acceleration_size * readout * predictive, acceleration_size * predictive, readout * predictive

    integrals = [quad(lambda c: integrands(c)[i], -np.inf, np.inf, epsabs=1e-11)[0] for i in range(3)]
    np.testing.assert_allclose(integrals[:2], [3.5, 1], rtol=1e-12)
    assert integrals[0] / integrals[1] == pytest.approx(_readouts(states, coarse)[2], abs=1e-12)
    assert abs(integrals[2] - 3.5) > .05  # ordinary averaging of points is a different operation


def test_part_prior_marginals_do_not_determine_the_composite_prior():
    aligned = np.array([[1., 1., 1.], [2., 2., 1.]])
    opposed = np.array([[1., 2., 1.], [2., 1., 1.]])
    for part in (0, 1):
        np.testing.assert_array_equal(np.sort(aligned[:, part]), np.sort(opposed[:, part]))
    np.testing.assert_array_equal(aligned[:, :2].sum(axis=1), [2, 4])
    np.testing.assert_array_equal(opposed[:, :2].sum(axis=1), [3, 3])
    # Both use equal probabilities and the same shared a and instrument model.
    # Their pushed-forward priors differ despite identical part prior marginals.

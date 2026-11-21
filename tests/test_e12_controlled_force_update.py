"""Finite Gaussian conditioning on known signed total loads."""
import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from explore.e12_mass_measurement import controlled_force_posterior


def test_heterogeneous_correlated_instruments_match_gaussian_pdf_and_mass_moments():
    pairs = np.array([[1., 3.], [3., 1.], [.9, .6], [4., 2.]])
    prior = np.array([2., .7, .4, 0.])
    loads = np.array([2.5, -.7])
    readings = np.array([[1.1, 1.6, 1.3], [-.1, -.9, -.6]])
    cov = np.array([[[1.2, .3, .2], [.3, .7, -.1], [.2, -.1, 1.]],
                    [[.3, -.1, 0.], [-.1, 1.4, .5], [0., .5, .8]]])
    posterior = controlled_force_posterior(pairs, prior, loads, readings, cov)
    expected = prior.copy()
    for j, (x, y) in enumerate(pairs):
        for s, record, covariance in zip(loads, readings, cov):
            mean = s * np.array([x, y, 1]) / (x + y)
            expected[j] *= multivariate_normal.pdf(record, mean, covariance)
    expected /= expected.sum()
    np.testing.assert_allclose(posterior.probabilities, expected)
    np.testing.assert_allclose(posterior.mean, expected @ pairs)
    np.testing.assert_allclose(posterior.covariance,
                               np.cov(pairs.T, aweights=expected, ddof=0))
    assert np.isneginf(posterior.log_probabilities[-1])
    np.testing.assert_array_equal(posterior.mass_pairs, pairs)


def test_remote_common_acceleration_does_not_erase_allocation_odds():
    # Tiny common absolute likelihood: subtracting two full log densities
    # would round their difference to zero. Force log odds are exactly 1.
    posterior = controlled_force_posterior([(1, 3), (2, 2)], [1., 2.],
                                           4., [1., 3., 1e150], np.eye(3))
    assert np.diff(posterior.log_probabilities).item() == pytest.approx(np.log(2) - 1)


def test_joint_and_sequential_updates_agree_after_probability_underflow():
    pairs = [[1., 3.], [2., 2.]]
    cov = np.eye(3) * .001
    prior = [2., 3.]
    rows = [[1., 3., 1.], [2., 2., 1.]]
    first = controlled_force_posterior(pairs, prior, 4., rows[0], cov)
    assert first.probabilities[1] == 0
    assert np.isfinite(first.log_probabilities[1])
    last = controlled_force_posterior(pairs, first.log_probabilities, 4., rows[1],
                                      cov, log_prior=True)
    batch = controlled_force_posterior(pairs, prior, [4., 4.], rows, cov)
    np.testing.assert_allclose(last.probabilities, [.4, .6], atol=1e-12)
    np.testing.assert_allclose(batch.probabilities, last.probabilities)


def test_small_odds_between_leading_candidates_survive_a_remote_reference():
    posterior = controlled_force_posterior([(1, 3), (2, 2), (4, 4)], [3., 2., 1.],
                                           4., [2., 2., 1.], np.diag([1e-20, 1e-20, 1.]))
    assert posterior.probabilities[0] == 0.
    assert posterior.log_probabilities[1] - posterior.log_probabilities[2] == pytest.approx(
        np.log(2) + .125
    )


@pytest.mark.parametrize("loads,readings,cov", [
    ([], [], []), ([0.], [[2e200, -1., 0.]], np.eye(3)),
])
def test_no_drive_preserves_prior_even_for_large_weights(loads, readings, cov):
    pairs = [[2., 2.], [2., 2.], [1., 3.]]
    posterior = controlled_force_posterior(pairs, [1e308, 0., 1e308], loads, readings, cov)
    np.testing.assert_allclose(posterior.probabilities, [.5, 0., .5])
    np.testing.assert_allclose(posterior.covariance, [[.25, -.25], [-.25, .25]])


def test_large_total_and_covariance_contribution_from_tiny_posterior_atom():
    # At the supplied load both accelerations are O(1) although the first
    # exact total overflows double precision.
    loads = 1e308
    pairs = [[1e308, 1e308], [4e307, 4e307]]
    posterior = controlled_force_posterior(pairs, [1., 1.], loads,
                                           [5e307, 5e307, .5], np.eye(3))
    expected_logs = np.array([0., -.5 * (.5 - 1.25)**2])
    np.testing.assert_allclose(posterior.log_probabilities, expected_logs - logsumexp(expected_logs))
    tiny = controlled_force_posterior([[1., 1.], [1e160, 1.]], [0., -750.],
                                      [], [], [], log_prior=True)
    assert tiny.probabilities[1] == 0.
    assert tiny.covariance[0, 0] == pytest.approx(np.exp(320*np.log(10) - 750))


@pytest.mark.parametrize("field,value", [
    ("mass_pairs", [[1., 0.]]), ("prior_weights", [0.]),
    ("loads", [1., 2.]), ("readings", [1., 2.]),
    ("covariances", [[1., 0., 0.], [.1, 1., 0.], [0., 0., 1.]]),
    ("covariances", np.diag([1., 1., 0.])),
])
def test_ill_defined_experiment_is_rejected(field, value):
    arguments = dict(mass_pairs=[[1., 2.]], prior_weights=[1.], loads=0.,
                     readings=[0., 0., 0.], covariances=np.eye(3))
    arguments[field] = value
    with pytest.raises(ValueError):
        controlled_force_posterior(**arguments)

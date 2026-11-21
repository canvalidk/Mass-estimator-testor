"""Exact and numerical checks of route agreement for deliberate Gaussian noise."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erf, expit, roots_hermitenorm
from scipy.stats import chi2, ncx2, norm


def direction(mass):
    vector = np.array([mass, 1.])
    return vector / np.linalg.norm(vector)


def projected_cost(rows, mass):
    unit = direction(mass)
    return np.sum((rows - np.outer(rows @ unit, unit))**2)


def log_reference_factor(rows, mass, variance, drive_variance):
    coefficient = drive_variance / (2 * variance * (variance + drive_variance))
    return -coefficient * projected_cost(rows, mass)


@pytest.mark.parametrize("drive_variance,log_odds", [(1., .45), (4., .72)])
def test_distinct_odds_and_measured_zero_neutrality(drive_variance, log_odds):
    rows = np.array([[3., 1.]])
    factors = [log_reference_factor(rows, m, 1., drive_variance) for m in (1, 2)]
    padded = np.vstack((rows, np.zeros((1, 2))))
    assert factors[1] - factors[0] == pytest.approx(log_odds)
    for i, mass in enumerate((1, 2)):
        assert log_reference_factor(padded, mass, 1., drive_variance) == factors[i]


@pytest.mark.parametrize("drive_variance", [1., 4.])
def test_proper_predictive_tower_through_full_gaussian_bridge(drive_variance):
    masses = (1., 2.)
    prior = .4
    s, t = 1., 2.5
    observation = np.array([3., 1.])
    nodes, weights = roots_hermitenorm(28)
    standard = np.array(np.meshgrid(nodes, nodes)).reshape(2, -1).T
    weights = np.outer(weights, weights).ravel() / (2*np.pi)

    def prob_at_mass2(rows, variance):
        k0 = np.sum(rows**2, axis=-1) - (rows @ direction(masses[0]))**2
        k1 = np.sum(rows**2, axis=-1) - (rows @ direction(masses[1]))**2
        gain = drive_variance * (k0-k1) / (2*variance*(variance+drive_variance))
        return expit(np.log(prior/(1-prior)) + gain)

    coarse_mass2 = prob_at_mass2(observation, t)
    fine_average = 0.
    for mass, posterior_weight in zip(masses, (1-coarse_mass2, coarse_mass2)):
        u = direction(mass)
        fine = s*np.eye(2) + drive_variance*np.outer(u, u)
        bridge_gain = np.linalg.solve(fine + (t-s)*np.eye(2), fine)
        bridge_cov = (t-s)*bridge_gain
        draws = bridge_gain @ observation + standard @ np.linalg.cholesky(bridge_cov).T
        fine_average += posterior_weight * (weights @ prob_at_mass2(draws, s))
    assert fine_average == pytest.approx(coarse_mass2, abs=2e-9)


def test_positive_preparations_give_same_bridge_from_different_coarse_records():
    fine_fraction = .5
    y, other_y = np.array([3., 1.]), np.array([1., 2.])
    mu, other_mu = np.array([4., 4.]), np.array([6., 3.])
    center = fine_fraction * y + (1-fine_fraction) * mu
    alternative = fine_fraction * other_y + (1-fine_fraction) * other_mu
    np.testing.assert_array_equal(center, [3.5, 2.5])
    np.testing.assert_array_equal(center, alternative)


def test_two_preparations_at_one_record_realize_any_bridge_center_difference():
    centers = np.array([[7., 2.], [3., -2.]])
    physical_lines = np.array([[1., 1.], [2., 1.]]).T
    fraction = .4
    signed_drives = np.linalg.solve(physical_lines, centers[0]-centers[1]) / (1-fraction)
    mu1, mu2 = signed_drives[0]*physical_lines[:, 0], -signed_drives[1]*physical_lines[:, 1]
    common_record = (centers[0] - (1-fraction)*mu1) / fraction
    np.testing.assert_allclose(fraction*common_record + (1-fraction)*mu2, centers[1])


def test_rank_two_conflict_nearest_to_a_signal_even_under_a_singular_preparation():
    graph_signal = np.array([[2.5, 2.5], [2.5, 2.5]])
    rows = graph_signal + np.array([[1.5, -1.5], [-1.5, 1.5]])
    assert np.linalg.matrix_rank(rows) == 2
    np.testing.assert_allclose(rows @ np.outer(direction(1), direction(1)), graph_signal)
    costs = [projected_cost(rows, m) for m in (1, 2)]
    assert costs == pytest.approx([9., 10.6])
    # The second atom may use any drive. Even its best possible position loses.
    gain = []
    for variance in (1., .1, .01):
        posterior_mass1 = expit(np.log(.1/.9) + (costs[1]-costs[0])/(2*variance))
        gain.append(posterior_mass1)
    assert gain[0] < .5 < gain[1] < gain[2]
    assert gain[-1] == pytest.approx(1.)


@pytest.mark.parametrize("degrees", [1, 3, 6])
def test_transported_covector_rejection_matches_fine_level(degrees):
    fine_fraction = .4
    level = .1
    cutoff = chi2.isf(level, degrees)

    def probability(coarse_cost):
        noncentrality = fine_fraction*coarse_cost/(1-fine_fraction)
        return ncx2.sf(cutoff/(1-fine_fraction), degrees, noncentrality)

    assert probability(3.) > probability(1.)
    result, _ = quad(lambda value: probability(value)*chi2.pdf(value, degrees),
                     0., np.inf, epsabs=2e-10)
    assert result == pytest.approx(level, abs=2e-9)


def test_anisotropic_added_noise_exposes_the_hidden_drive_in_full_reverse_average():
    fine = np.eye(2)
    added = np.diag([3., 1.])
    bridge_gain = np.linalg.solve(fine + added, fine)
    bridge_cov = fine - fine @ bridge_gain
    contrast = np.array([1., -1.])/np.sqrt(2)
    variance = contrast @ bridge_cov @ contrast
    assert variance == pytest.approx(5/8)
    for drive in (0., 1., 5.):
        mean = contrast @ (np.eye(2)-bridge_gain) @ [drive, drive]
        expected = 1-4/np.sqrt(13)*np.exp(-drive**2/104)
        result, _ = quad(lambda x: (1-np.sqrt(2)*np.exp(-x*x/2)) *
                         norm.pdf(x, mean, np.sqrt(variance)), -np.inf, np.inf)
        assert result == pytest.approx(expected, abs=1e-10)


def test_raw_gaussian_mixture_convolves_but_padding_changes_its_odds():
    s, t, extra = .3, .8, 2.
    broad_weight = .4

    def spread(v, variance, dims):
        return (2*np.pi*variance)**(-dims/2) * np.exp(-v/(2*variance))

    def factor(v, variance, dims):
        return spread(v, variance, dims) + broad_weight*spread(v, variance+extra, dims)

    for residual in (0., .5, 3.):
        summed, _ = quad(lambda change: factor((residual+change)**2, s, 1) *
                         norm.pdf(change, scale=np.sqrt(t-s)), -np.inf, np.inf)
        assert summed == pytest.approx(factor(residual**2, t, 1))
    odds = factor(1., s, 1)/factor(.2, s, 1)
    padded_odds = factor(1., s, 2)/factor(.2, s, 2)
    assert abs(odds - padded_odds) > .02


def test_residual_likelihood_ratio_transports_without_a_push_law():
    s, t, shift, coarse = .8, 2.7, 1.3, -.9
    fraction = s/t
    scale = np.sqrt(s*(1-fraction))

    def e_value(residual, variance):
        return np.exp((shift*residual - shift**2/2)/variance)

    # Integrate over the centered candidate's bridge for the exact covector.
    nodes, weights = roots_hermitenorm(24)
    bridge = fraction*coarse + scale*nodes
    average = (weights @ e_value(bridge, s))/np.sqrt(2*np.pi)
    assert average == pytest.approx(e_value(coarse, t))
    under_null = (weights @ e_value(np.sqrt(s)*nodes, s))/np.sqrt(2*np.pi)
    assert under_null == pytest.approx(1.)
    assert e_value(5., t) > e_value(1., t)


def test_covariance_probes_need_the_off_diagonal_entries_for_sign_symmetry():
    # 1+c*sgn(F)*sgn(A) is invariant under joint reversal of the row.
    # Independent noise cannot see the second term in the coarse zero.
    def probe(c, corr):
        # Condition on the first channel of a correlated standard Gaussian.
        signal = corr/np.sqrt(2*(1-corr*corr))
        correction, _ = quad(lambda force: 2*norm.pdf(force)*erf(signal*force), 0, np.inf)
        return 1 + c*correction

    assert probe(.2, 0.) == probe(.7, 0.) == 1.
    assert probe(.2, .6) != pytest.approx(probe(.7, .6))

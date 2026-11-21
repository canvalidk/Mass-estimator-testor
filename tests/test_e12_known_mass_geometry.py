"""Known-mass addition through a single recorded acceleration channel."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from explore.e12_tangent_pairings import tangent_pairings
from met.conditional_score import gaussian_score_accepts, gaussian_score_test
from met.full import Reading, _h_terms, bartlett_eta, gram_eta, isotropic, reading_state


def add_known_mass(reading, added_mass):
    dimension = reading.dimension
    identity = np.eye(dimension)
    transform = np.block([[identity, added_mass * identity], [np.zeros_like(identity), identity]])
    observation = transform @ np.r_[reading.force, reading.acceleration]
    return Reading(observation[:dimension], observation[dimension:],
                   transform @ reading.cov @ transform.T, reading.subspace, reading.window)


def candidate_cost(reading, masses):
    pre = reading.pre()
    return pre["energy"] - _h_terms(pre, 1.0, np.arctan(masses))


def sample_readings():
    result = []
    rng = np.random.default_rng(129)
    for dimension, k in [(1, 1), (3, 2), (2, 2)]:
        matrix = rng.normal(size=(2 * dimension, 2 * dimension))
        cov = matrix @ matrix.T + np.eye(2 * dimension)
        basis, _ = np.linalg.qr(rng.normal(size=(dimension, k)))
        result.append(Reading(rng.normal(size=dimension), rng.normal(size=dimension), cov, basis))
    return result


def test_absolute_cost_and_fitted_latent_drive_on_correlated_restricted_graphs():
    masses = np.array([0.02, 0.55, 1.4, 4.0, 15.0])
    added_mass = 1.75
    for reading in sample_readings():
        total = add_known_mass(reading, added_mass)
        np.testing.assert_allclose(candidate_cost(reading, masses),
                                   candidate_cost(total, masses + added_mass), atol=3e-13)
        for m in masses:
            fits = []
            for record, mass in [(reading, m), (total, m + added_mass)]:
                root = np.linalg.cholesky(record.cov)
                design = np.linalg.solve(root, np.vstack([mass * record.subspace, record.subspace]))
                data = np.linalg.solve(root, np.r_[record.force, record.acceleration])
                fits.append(np.linalg.lstsq(design, data, rcond=None)[0])
            np.testing.assert_allclose(*fits, atol=5e-14)


@pytest.mark.parametrize("mass", [0.08, 0.6, 2.1, 11.0])
def test_reciprocal_tangent_energies_and_actual_conditional_score_tests_transport(mass):
    added_mass = 2.5
    original = sample_readings()
    total = [add_known_mass(reading, added_mass) for reading in original]
    old = tangent_pairings(original, mass)
    new = tangent_pairings(total, mass + added_mass)
    factor = (1 + (mass + added_mass)**2) / (1 + mass**2)
    np.testing.assert_allclose([new.score, new.weight, new.correction],
                               [old.score * factor, old.weight * factor**2, old.correction / factor])
    np.testing.assert_allclose([new.within, new.between, new.common], [old.within, old.between, old.common])
    assert old.degrees == new.degrees
    score = gaussian_score_test(original, mass)
    translated = gaussian_score_test(total, mass + added_mass)
    assert translated["normal_score"] == pytest.approx(score["normal_score"], abs=3e-13)
    assert translated["p_value"] == pytest.approx(score["p_value"], abs=3e-13)
    masses = np.array([0.2, mass, 10])
    np.testing.assert_array_equal(gaussian_score_accepts(total, masses + added_mass, probability=0.6),
                                  gaussian_score_accepts(original, masses, probability=0.6))


def test_spatial_reciprocal_and_exterior_pairing():
    reading = isotropic([4.0, -1.0, 3.0], [1.0, 0.0, 2.0], 1.5, 1.0)
    added_mass = 1.7
    total = add_known_mass(reading, added_mass)
    metric = np.array([[3., 1., 0.], [1., 2., 0.], [0., 0., 0.5]])
    covector = (reading.acceleration @ metric) / (reading.acceleration @ metric @ reading.acceleration)
    old, new = covector @ reading.force, covector @ total.force
    assert new == pytest.approx(old + added_mass)
    np.testing.assert_allclose(total.force - new * total.acceleration,
                               reading.force - old * reading.acceleration)
    np.testing.assert_allclose(np.cross(total.force, total.acceleration),
                               np.cross(reading.force, reading.acceleration))


def test_preserving_only_channel_variances_changes_the_experiment():
    reading = isotropic([2], [1], 1, 1)
    total = add_known_mass(reading, 1)
    diagonalized = Reading(total.force, total.acceleration, np.diag(np.diag(total.cov)))
    np.testing.assert_allclose(total.cov, [[2, 1], [1, 1]])
    assert candidate_cost(reading, [1])[0] == pytest.approx(1 / 2)
    assert candidate_cost(total, [2])[0] == pytest.approx(1 / 2)
    assert candidate_cost(diagonalized, [2])[0] == pytest.approx(1 / 6)


def test_translated_candidate_cone_changes_the_reachable_minimum():
    reading = isotropic([-1], [1], 1, 1)
    added_mass = 2
    total = add_known_mass(reading, added_mass)
    assert reading_state(reading, 1)["Q_plus"] == pytest.approx(1)
    assert reading_state(total, 1)["Q_plus"] == pytest.approx(0, abs=1e-14)
    old_domain = np.r_[0, np.geomspace(1e-3, 1e3, 101), np.inf]
    assert candidate_cost(total, old_domain + added_mass).min() == pytest.approx(1)


def test_zero_tangent_and_zero_reading_acceptance():
    reading = isotropic([1], [-1], 1, 1)
    translated = add_known_mass(reading, 1)
    for records, mass in [([reading], 1), ([translated], 2)]:
        assert gaussian_score_test(records, mass)["degenerate"]
        assert tangent_pairings(records, mass).within == pytest.approx(2)
    zero = isotropic([0], [0], 1, 1)
    assert gaussian_score_test([add_known_mass(zero, 1.5)], 2.5)["p_value"] == 1


def test_two_existing_series_discounts_both_pass_addition_yet_disagree():
    records = [isotropic([4], [1], 1, 1), isotropic([1], [4], 1, 1)]
    for b in [0, 1.5, 4]:
        total = [add_known_mass(row, b) for row in records]
        assert all(reading_state(row, 1)["Q_plus"] == pytest.approx(0, abs=1e-12) for row in total)
        eta, curvature, score_square, phi = bartlett_eta(total, 1, [1, 1])
        other_eta, residual, explained, other_phi = gram_eta(total, 1, [1, 1])
        chart_factor = (1 + (1 + b)**2) / 2
        assert math.tan(phi) == pytest.approx(1 + b)
        assert phi == pytest.approx(other_phi)
        np.testing.assert_allclose([curvature, score_square],
                                   np.array([16, 225 / 2]) * chart_factor**2)
        np.testing.assert_allclose([eta, other_eta, residual, explained], [32 / 225, 16 / 25, 9, 25])


def test_liseo_ratio_and_reference_both_translate_and_remain_distinct():
    original = isotropic([2], [0], 1, 1)
    added_mass = 1.75
    total = add_known_mass(original, added_mass)

    def reference_density(reading, mass):
        normal = np.array([1., -mass])
        variance = normal @ reading.cov @ normal
        defect = reading.force[0] - mass * reading.acceleration[0]
        return math.sqrt(np.linalg.det(reading.cov)) / variance * math.exp(-0.5 * defect**2 / variance)

    # The positive-event probabilities follow by direct integration of the posterior Gaussian
    # over original positive denominators and negative denominators separately.
    def ratio_probability(reading, low, high):
        ff, cross = reading.cov[0]
        aa = reading.cov[1, 1]
        conditional_sd = math.sqrt(ff - cross**2 / aa)

        def integrand(y):
            conditional_mean = reading.force[0] + cross / aa * (y - reading.acceleration[0])
            p_low = norm.cdf((low * y - conditional_mean) / conditional_sd)
            p_high = norm.cdf((high * y - conditional_mean) / conditional_sd)
            return abs(p_high - p_low) * norm.pdf(y, loc=reading.acceleration[0], scale=math.sqrt(aa))

        return sum(quad(integrand, left, right)[0] for left, right in [(-np.inf, 0), (0, np.inf)])

    probabilities = []
    for reading, lower in [(original, 0), (total, added_mass)]:
        reference = quad(lambda mass: reference_density(reading, mass), lower, lower + 1)[0]
        reference /= quad(lambda mass: reference_density(reading, mass), lower, np.inf)[0]
        ratio = ratio_probability(reading, lower, lower + 1) / ratio_probability(reading, lower, np.inf)
        probabilities.append((reference, ratio))
    np.testing.assert_allclose(probabilities[0], probabilities[1], atol=1e-12)
    reference, ratio = probabilities[0]
    assert ratio == pytest.approx(2 * norm.cdf(math.sqrt(2)) * norm.cdf(-math.sqrt(2)), abs=1e-12)
    assert reference > ratio + 0.05

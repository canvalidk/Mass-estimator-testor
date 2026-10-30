"""B06: exact flat-factor information moments and conditional median limits."""
import math

import numpy as np
import pytest
from scipy.special import roots_hermitenorm

from met.full import Reading, _h_terms, isotropic, ours30


def test_t17_median_concentrates_with_correlated_covariance():
    covariance = np.array([[2, .3, .2, -.1], [.3, 1, .1, .2],
                           [.2, .1, 1.5, .4], [-.1, .2, .4, 1.2]])
    pushes = (np.array([1, .3]), np.array([-.4, 1]))
    errors = (np.array([2, -1, .5, 1.2]), np.array([-1, .4, .8, -1.5]))
    results = []
    for scale in (100, 1000):
        readings = []
        for push, error in zip(pushes, errors):
            y = scale * np.r_[2.5 * push, push] + error
            readings.append(Reading(y[:2], y[2:], covariance))
        results.append(ours30(readings, 1, readout=False))
    assert abs(results[-1]["median"] - 2.5) < .004
    assert abs(results[-1]["median"] - 2.5) < .11 * abs(results[0]["median"] - 2.5)
    widths = [np.log(r["interval"][1] / r["interval"][0]) for r in results]
    assert widths[-1] < .11 * widths[0]


def test_t17_fixed_noise_need_not_send_g_or_eta_to_one():
    u, b = np.array([0., 40., -30.]), np.array([8., -6., 4.])
    limit_g = 7.5 / b ** 2
    mean_u = np.dot(limit_g, u) / limit_g.sum()
    limit_eta = 5 * limit_g.sum() / np.sum(limit_g ** 2 * (u - mean_u) ** 2)
    assert limit_eta == pytest.approx(0.018772296552679514)
    deviations = []
    for scale in (1000, 10000):
        readings = [isotropic([2 * scale + parallel, transverse], [scale, 0], 1, 1)
                    for parallel, transverse in zip(u, b)]
        result = ours30(readings, 1, readout=False)
        assert max(result["g"]) < .5
        assert 0 < result["eta"] < .03
        deviations.append(abs(result["median"] - 2))
    np.testing.assert_allclose(result["g"], limit_g, rtol=.004)
    assert result["eta"] == pytest.approx(limit_eta, rel=.006)
    assert deviations[-1] < .001
    assert deviations[-1] < .13 * deviations[0]


def _score_moments(mass, ratio, pair_length):
    """Three-node Gaussian quadrature is exact for the degree-four moments."""
    nodes, weights = roots_hermitenorm(3)
    weights /= math.sqrt(2 * math.pi)
    angle = math.atan(mass / ratio)
    derivative = ratio / (ratio ** 2 + mass ** 2)
    second = -2 * ratio * mass / (ratio ** 2 + mass ** 2) ** 2
    mean_score, curvature, score_square = 0., 0., 0.
    for i, force_noise in enumerate(nodes):
        for j, acceleration_noise in enumerate(nodes):
            force = ratio * (pair_length * math.sin(angle) + force_noise)
            acceleration = pair_length * math.cos(angle) + acceleration_noise
            reading = isotropic([force], [acceleration], ratio, 1)
            _, h1, h2 = _h_terms(reading.pre(), ratio, np.array([angle]), 2)
            score = .5 * h1[0] * derivative
            negative_curvature = -.5 * (h2[0] * derivative ** 2 + h1[0] * second)
            probability = weights[i] * weights[j]
            mean_score += probability * score
            curvature += probability * negative_curvature
            score_square += probability * score ** 2
    return mean_score, curvature, score_square, derivative


@pytest.mark.parametrize("mass,ratios,lengths", [
    (2., [1, 1, 1], [0., .5, 3.]),
    (.7, [.4, 1.5, 3.], [1., 2., .2]),
])
def test_t18_score_information_moments_are_exact_gaussian_integrals(mass, ratios, lengths):
    information, variance, derivatives = [], [], []
    for ratio, length in zip(ratios, lengths):
        score, a, b, derivative = _score_moments(mass, ratio, length)
        assert score == pytest.approx(0, abs=3e-15)
        assert a == pytest.approx(derivative ** 2 * length ** 2, abs=3e-15)
        assert b == pytest.approx(derivative ** 2 * (1 + length ** 2), abs=3e-15)
        information.append(a)
        variance.append(b)
        derivatives.append(derivative)
    a, b = sum(information), sum(variance)
    expected = 1 + np.dot(derivatives, derivatives) / a
    assert b / a == pytest.approx(expected, rel=2e-14)
    assert (b / a ** 2) / (1 / a) == pytest.approx(expected)
    eta_limit = a / b
    assert 0 < eta_limit <= 1
    assert 1 / (eta_limit * a) == pytest.approx(b / a ** 2)


def test_t18_constant_power_changes_the_unweighted_width_ratio():
    _, a, b, derivative = _score_moments(1.7, .8, 1.2)
    power = .4
    weighted_a, weighted_b = power * a, power ** 2 * b
    raw_ratio = 1 + 1 / 1.2 ** 2
    assert b / a == pytest.approx(raw_ratio)
    assert weighted_b / weighted_a == pytest.approx(power * raw_ratio)
    assert weighted_b / weighted_a != pytest.approx(raw_ratio)

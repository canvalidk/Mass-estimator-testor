"""Finite witnesses for forgotten observations and changes of mass prior."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.integrate import trapezoid

from explore.e12_residual_join import residual_log_evidence


def update(prior, factor):
    weights = [p * value for p, value in zip(prior, factor)]
    total = sum(weights)
    return tuple(w / total for w in weights)


def test_fair_sign_forgets_at_initial_prior_but_changes_a_reached_prior():
    prior = (F(1, 2), F(1, 2))
    factors = [(3, 1), (1, 3)]
    first = [update(prior, factor) for factor in factors]
    assert [post[0] for post in first] == [F(3, 4), F(1, 4)]
    assert sum(post[0] for post in first) / 2 == prior[0]

    second = [update(first[0], factor) for factor in factors]
    assert [post[0] for post in second] == [F(9, 10), F(1, 2)]
    second_mean = sum(post[0] for post in second) / 2
    assert first[0][0] - second_mean == F(1, 20)


def test_exact_reciprocal_moment_identity_detects_hidden_preparation():
    prior, other = (F(1, 2), F(1, 2)), (F(3, 4), F(1, 4))
    ratios = [b / a for a, b in zip(prior, other)]
    values = []
    for factor in [(3, 1), (1, 3)]:
        post = update(prior, factor)
        z = sum(h * p for h, p in zip(ratios, post))
        other_post = update(other, factor)
        assert tuple(h * p / z for h, p in zip(ratios, post)) == other_post
        assert sum(p / h for h, p in zip(ratios, other_post)) == 1 / z
        values.append(z)
    assert sum(values) / 2 == 1
    gap = sum(1 / z for z in values) / 2 - 1
    assert gap == sum((z - 1)**2 / z for z in values) / 2 == F(1, 15)


def test_two_half_cauchy_priors_allow_changes_on_an_inversion_invariant_band():
    # Inversion m -> 4/m interchanges scales 1 and 4, and preserves [1,4].
    def band_probability(scale):
        return 2 / np.pi * (np.arctan(4 / scale) - np.arctan(1 / scale))

    p = band_probability(1)
    assert 0 < p < 1
    assert band_probability(4) == pytest.approx(p)
    for scale in [1, 4]:
        prior = (band_probability(scale), 1 - band_probability(scale))
        posts = []
        for sign in [-1, 1]:
            # The factor is constant separately inside and outside the band.
            factor = (1 + sign * .7 * (1 - p), 1 - sign * .7 * p)
            assert np.dot(prior, factor) == pytest.approx(1.)
            posts.append(update(prior, factor))
        assert abs(posts[0][0] - posts[1][0]) > .1
        np.testing.assert_allclose(np.mean(posts, axis=0), prior, atol=2e-15)


def test_e12_factor_has_nonunit_posterior_ratio_moment_at_positive_signal():
    # Work at noise scale 1. Under the scale-1 half-Cauchy, theta=arctan(m)
    # is uniform. This form of dPi_4/dPi_1 is finite including at theta=pi/2.
    theta = np.linspace(0., np.pi / 2, 4097)
    h = 4 / (16 * np.cos(theta)**2 + np.sin(theta)**2)
    assert 2 / np.pi * trapezoid(h, theta) == pytest.approx(1.)
    assert 4 * (1 + 2**2) / (16 + 2**2) == 1.
    target_moment = 4 * (1 + 3**2) / (16 + 3**2)
    for strength in [20., 100.]:
        for delta in [0., .01]:
            y = strength * np.array([3., 1.]) + [delta, -delta]
            log_factor = residual_log_evidence([np.outer(y, y)], theta)
            weight = np.exp(log_factor - log_factor.max())
            z = trapezoid(h * weight, theta) / trapezoid(weight, theta)
            assert z == pytest.approx(target_moment, abs=.007)
            assert z > 1.5


def test_three_centres_force_positive_variance_even_at_the_compactified_extremizer():
    # U=1/(4+m²), so the extremizing U-values correspond to infinite mass
    # and m²=8/7. Both Stieltjes moments agree with the scale-1 prior.
    prob = [F(1, 7), F(6, 7)]
    values = [F(0), F(7, 36)]
    mean = sum(p * u for p, u in zip(prob, values))
    first = sum(p * u / (1 + 12*u) for p, u in zip(prob, values))
    assert mean == F(1, 6)
    assert first == F(1, 20)
    second = (mean - first) / 12
    variance_bound = second + 12 * second**2 / first - mean**2
    measured_variance = sum(p * (4*u - 4*mean)**2 for p, u in zip(prob, values))
    assert measured_variance == 16 * variance_bound == F(2, 27)


def test_reached_prior_mixture_identity_detects_positive_variance():
    # g_A(m) averages densities over a subset of outcomes. Its variance equals
    # the event-weighted excess of the pairwise density inner products over 1.
    prior = (F(1, 2), F(1, 2))
    probabilities = (F(1, 4), F(1, 4), F(1, 2))
    densities = ((F(3, 2), F(1, 2)), (F(5, 4), F(3, 4)),
                 (F(5, 8), F(11, 8)))
    for i in range(2):
        assert sum(q * f[i] for q, f in zip(probabilities, densities)) == 1
    event = (0, 1)
    p_event = sum(probabilities[y] for y in event)
    mixture = [sum(probabilities[y] * densities[y][i] for y in event)
               for i in range(2)]
    variance = sum(p * (g - p_event)**2 for p, g in zip(prior, mixture))
    assert variance == F(9, 256)
    pair_gap = sum(probabilities[y] * probabilities[z]
                   * (sum(p * a * b for p, a, b in
                          zip(prior, densities[y], densities[z])) - 1)
                   for y in event for z in event)
    assert pair_gap == variance

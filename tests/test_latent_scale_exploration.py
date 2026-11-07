"""E5: moment preservation and proper diffuse limits do not fix T48's measure.

Independent quadratures check the counterexamples documented in T48's
5 Nov exploration. These are conditional probability constructions.
"""
import math
from fractions import Fraction

import numpy as np
import pytest
from scipy import integrate, special, stats


def _polynomial(x):
    return x**3 - 18 * x**2 + 90 * x - 120


def _weight(x):
    """Here x=5/v, with the original x distribution Gamma(5, rate=1)."""
    return 1 + _polynomial(x) / 240


def _expectation(function):
    value, error = integrate.quad(
        lambda x: stats.gamma.pdf(x, 5) * function(x),
        0, np.inf, epsabs=2e-11, epsrel=2e-10, limit=200)
    assert error < 1e-8
    return value


def _density(cost, dimension):
    shape = 5 + dimension / 2
    rate = 1 + cost / 10
    posterior_polynomial = (
        shape * (shape + 1) * (shape + 2) / rate**3
        - 18 * shape * (shape + 1) / rate**2
        + 90 * shape / rate - 120)
    base = math.exp(
        special.gammaln(shape) - special.gammaln(5)
        - dimension / 2 * math.log(10 * math.pi)
        - shape * math.log(rate))
    return base * (1 + posterior_polynomial / 240)


@pytest.mark.parametrize("power, expected", [(-1, 1), (0, 1), (1, 1.25), (2, 65 / 32)])
def test_positive_tilt_preserves_total_precision_and_covariance(power, expected):
    actual = _expectation(lambda x: _weight(x) * (5 / x)**power)
    assert actual == pytest.approx(expected, rel=2e-10)


def test_preserved_moments_and_changed_fourth_moment_are_exact():
    def gamma_moment(power):
        return Fraction(math.factorial(4 + power), math.factorial(4))

    def polynomial_moment(power):
        return (gamma_moment(power + 3) - 18 * gamma_moment(power + 2)
                + 90 * gamma_moment(power + 1) - 120 * gamma_moment(power))

    assert [polynomial_moment(p) for p in [-1, 0, 1]] == [0, 0, 0]
    assert polynomial_moment(-2) == Fraction(-1, 2)
    assert 25 * (gamma_moment(-2) + polynomial_moment(-2) / 240) == Fraction(65, 32)


def test_weight_is_strictly_positive_at_all_possible_minima():
    # The cubic has positive leading coefficient. Its only nonnegative
    # critical points plus its left endpoint determine its global minimum.
    candidates = [0, 6 - math.sqrt(6), 6 + math.sqrt(6)]
    assert min(_weight(x) for x in candidates) == pytest.approx(0.5)


@pytest.mark.parametrize("dimension", [1, 3])
@pytest.mark.parametrize("cost", [0., 9., 100.])
def test_changed_residual_density_matches_independent_variance_integral(dimension, cost):
    actual = _expectation(
        lambda x: _weight(x) * (x / (10 * math.pi))**(dimension / 2)
        * math.exp(-cost * x / 10))
    assert actual == pytest.approx(_density(cost, dimension), rel=2e-10)


@pytest.mark.parametrize("dimension", [1, 3])
def test_changed_residual_law_is_normalized(dimension):
    area = 2 * math.pi**(dimension / 2) / special.gamma(dimension / 2)
    actual, error = integrate.quad(
        lambda radius: area * radius**(dimension - 1) * _density(radius**2, dimension),
        0, np.inf, epsabs=2e-10, epsrel=2e-10, limit=200)
    assert error < 2e-8
    assert actual == pytest.approx(1, abs=2e-9)


def test_same_noise_moments_give_different_relative_mass_evidence():
    # These are exact gamma-posterior expectations at K=0 and K=9.
    first_ratio, second_ratio = Fraction(639, 640), Fraction(110367, 109744)
    baseline = lambda cost: stats.t.pdf(math.sqrt(cost), df=10)
    assert _density(0, 1) / baseline(0) == pytest.approx(float(first_ratio))
    assert _density(9, 1) / baseline(9) == pytest.approx(float(second_ratio))
    assert abs(float(second_ratio / first_ratio) - 1) > 0.007


@pytest.mark.parametrize("nu,b", [(0.5, -1), (10, -1), (10, -0.5), (10, 0), (10, 1)])
def test_power_family_precision_selects_b_zero_within_that_family(nu, b):
    alpha = nu / 2
    # Normalizing v^b q gives IG(alpha-b, alpha).
    precision = stats.invgamma.expect(
        lambda v: 1 / v, args=(alpha - b,), scale=alpha)
    assert precision == pytest.approx((alpha - b) / alpha, rel=2e-8)
    assert (abs(precision - 1) < 1e-8) == (b == 0)


@pytest.mark.parametrize("weight", [lambda x: 1., lambda x: x / 5, _weight])
def test_proper_signal_prior_has_q_marginal_but_its_diffuse_likelihood_keeps_weight(weight):
    parallel, residual = 0.6, 1.1

    def scaled_likelihood(volume_radius):
        def integrand(x):
            variance = 5 / x
            h = weight(x)
            half_width = volume_radius / h
            probability = (special.ndtr((half_width - parallel) / math.sqrt(variance))
                           - special.ndtr((-half_width - parallel) / math.sqrt(variance)))
            return h * stats.norm.pdf(residual, scale=math.sqrt(variance)) * probability
        return _expectation(integrand)

    # Integral of p_A(t|v) is exactly one at every finite A,v, preserving q.
    for variance in [0.1, 1., 20.]:
        h = weight(5 / variance)
        half_width = 3 / h
        total = integrate.quad(lambda t: h / 6, -half_width, half_width)[0]
        assert total == pytest.approx(1)

    limit = _expectation(lambda x: weight(x) * stats.norm.pdf(residual, scale=math.sqrt(5 / x)))
    values = [scaled_likelihood(radius) for radius in [1., 10., 100.]]
    assert values[0] < values[1] <= values[2] + 1e-12
    assert values[-1] == pytest.approx(limit, rel=2e-7)


def test_nuisance_free_residual_uses_original_q_for_every_signal_prior():
    # Integrating a proper signal prior cancels its h dependent width,
    # leaving q in the residual pushforward at each finite A.
    residual = 1.1
    actual = _expectation(lambda x: stats.norm.pdf(residual, scale=math.sqrt(5 / x)))
    assert actual == pytest.approx(stats.t.pdf(residual, df=10), rel=2e-10)
    assert abs(actual - _density(residual**2, 1)) > 1e-5


@pytest.mark.parametrize("volume_radius", [1., 10., 100.])
def test_coupled_diffuse_model_changes_calibration_conditional_on_zero_signal(volume_radius):
    # h=1/v. At t=0 the proper uniform density is h/(2A), so
    # conditioning the original joint on that signal tilts q at every A.
    denominator = _expectation(lambda x: (x / 5) / (2 * volume_radius))
    numerator = _expectation(lambda x: (x / 5)**2 / (2 * volume_radius))
    assert numerator / denominator == pytest.approx(1.2, rel=2e-10)
    # Independent uniform signals give the declared precision at every A.
    independent = _expectation(lambda x: (x / 5) / (2 * volume_radius))
    assert independent / (1 / (2 * volume_radius)) == pytest.approx(1, rel=2e-10)

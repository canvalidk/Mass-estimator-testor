"""Deterministic checks of T51's conditional law-composition statements."""
import math

import numpy as np
import pytest
from scipy import integrate, special, stats

from met.cdf_calibration import normal_cdf_warp
from met.composition_calibration import (FinitePositiveLaw, compose_independent,
                                         gaussian_linear_coverage,
                                         pushforward_joint, warp_finite_law)


def _assert_law(actual, support, probabilities):
    assert actual.support == pytest.approx(support, abs=1e-14)
    assert actual.probabilities == pytest.approx(probabilities, abs=2e-14)


@pytest.mark.parametrize("operation,outputs", [("addition", [2., 3., 4.]),
                                               ("parallel", [.5, 2/3, 1.])])
def test_exact_t47_bernoulli_witness_and_witness_specific_dependent_coupling(operation, outputs):
    original = FinitePositiveLaw([1., 2.], [.5, .5])
    eta = .25
    marginal = warp_finite_law(original, eta=eta)
    _assert_law(marginal, [1., 2.], [.5, .5])
    independent = compose_independent(marginal, marginal, operation=operation)
    after = warp_finite_law(compose_independent(original, original, operation=operation), eta=eta)
    # Independent normal integration, rather than the production warp call.
    h, _ = integrate.quad(stats.norm.pdf, -np.inf, math.sqrt(eta)*stats.norm.ppf(.25), epsabs=1e-12)
    assert h == pytest.approx(.3679661556049961, abs=1e-14)
    _assert_law(independent, outputs, [.25, .5, .25])
    _assert_law(after, outputs, [h, 1-2*h, h])
    assert abs(after.cdf(outputs[0])-independent.cdf(outputs[0])) > .1

    joint = np.array([[h, .5-h], [.5-h, h]])
    assert joint.sum(axis=0) == pytest.approx([.5, .5], abs=0.)
    assert joint.sum(axis=1) == pytest.approx([.5, .5], abs=0.)
    supplied = pushforward_joint([1., 2.], [1., 2.], joint, operation=operation)
    _assert_law(supplied, after.support, after.probabilities)
    assert 4*h-1 == pytest.approx(.4718646224199845, abs=1e-14)


@pytest.mark.parametrize("operation", ["addition", "parallel"])
def test_independent_composition_matches_separate_joint_enumeration(operation):
    first = FinitePositiveLaw([.7, 2., 4.], [2., 3., 5.])
    second = FinitePositiveLaw([.4, 1.2], [1., 3.])
    enumerated = {}
    for x, px in zip(first.support, first.probabilities):
        for y, py in zip(second.support, second.probabilities):
            value = x+y if operation == "addition" else 1/(1/x+1/y)
            enumerated[value] = enumerated.get(value, 0.)+px*py
    expected = sorted(enumerated.items())
    actual = compose_independent(first, second, operation=operation)
    _assert_law(actual, [x for x, _ in expected], [p for _, p in expected])
    identity = warp_finite_law(actual, eta=1.)
    _assert_law(identity, actual.support, actual.probabilities)


def test_properness_point_masses_monotone_coordinate_changes_and_inverse():
    law = FinitePositiveLaw([3., 1., 3., 7.], [.2, .3, .1, .4])
    _assert_law(law, [1., 3., 7.], [.3, .3, .4])
    warped = warp_finite_law(law, eta=.1)
    assert np.all(warped.probabilities > 0)
    assert sum(warped.probabilities) == pytest.approx(1., abs=1e-15)
    composed_warp = warp_finite_law(warp_finite_law(law, eta=.5), eta=.2)
    _assert_law(composed_warp, warped.support, warped.probabilities)
    changed = warp_finite_law(FinitePositiveLaw(law.support**2+2, law.probabilities), eta=.1)
    _assert_law(changed, warped.support**2+2, warped.probabilities)
    inverted = warp_finite_law(FinitePositiveLaw(1/law.support, law.probabilities), eta=.1)
    expected_inverse = FinitePositiveLaw(1/warped.support, warped.probabilities)
    _assert_law(inverted, expected_inverse.support, expected_inverse.probabilities)
    point = warp_finite_law(FinitePositiveLaw([2.], [1.]), eta=.1)
    _assert_law(point, [2.], [1.])
    assert law.cdf([-np.inf, 0., 1., 3., np.inf]) == pytest.approx([0., 0., .3, .6, 1.])


def test_small_outer_atoms_are_kept_using_both_tail_ranks():
    law = FinitePositiveLaw([1., 2., 3.], [1e-100, 1., 1e-100])
    warped = warp_finite_law(law, eta=.25)
    assert warped.probabilities[0] > 1e-30
    assert warped.probabilities[-1] == pytest.approx(warped.probabilities[0], rel=1e-13)
    assert warped.probabilities.sum() == pytest.approx(1., abs=1e-15)


def test_smooth_strictly_positive_lognormal_density_is_also_a_witness():
    # H_eta maps a lognormal's normal log-coordinate variance to sigma^2/eta.
    # Both original and warped laws have smooth positive density on (0,infinity).
    eta, sigma, stop = .25, .2, 2.

    def sum_cdf(sd):
        value, error = integrate.quad(
            lambda x: stats.lognorm.pdf(x, sd)*stats.lognorm.cdf(stop-x, sd),
            0., stop, epsabs=1e-11, epsrel=1e-11)
        assert error < 1e-9
        return value

    raw = sum_cdf(sigma)
    independently_warped_sum = sum_cdf(sigma/math.sqrt(eta))
    warped_raw_sum = normal_cdf_warp(raw, eta)
    assert raw == pytest.approx(.4723820114869846, abs=2e-12)
    assert independently_warped_sum == pytest.approx(.44759993809258586, abs=2e-12)
    assert warped_raw_sum == pytest.approx(.4861827208446029, abs=2e-12)
    assert warped_raw_sum-independently_warped_sum > .038


def test_gaussian_addition_commutes_for_a_common_fixed_eta_by_normal_integrals():
    # These whole-real-line Gaussians are a local reference, not positive laws.
    eta, means, scales, stop = .4, [2., 3.], [.3, .5], 5.3

    def sum_cdf(temperature):
        sd0, sd1 = np.asarray(scales)/math.sqrt(temperature)
        return integrate.quad(
            lambda x: stats.norm.pdf(x, means[0], sd0)*stats.norm.cdf(stop-x, means[1], sd1),
            -np.inf, np.inf, epsabs=1e-11, epsrel=1e-11)[0]

    assert sum_cdf(eta) == pytest.approx(normal_cdf_warp(sum_cdf(1.), eta), abs=2e-12)


def test_correct_gaussian_marginals_do_not_supply_a_correct_sum_joint_covariance():
    probability, rho = .95, .9
    actual, quoted = np.array([[1., rho], [rho, 1.]]), np.eye(2)
    for gradient in ([1., 0.], [0., 1.]):
        assert gaussian_linear_coverage(actual, quoted, gradient, probability=probability)["coverage"] == pytest.approx(probability, abs=1e-15)
    result = gaussian_linear_coverage(actual, quoted, [1., 1.], probability=probability)
    z = stats.norm.ppf((1+probability)/2)
    reference, _ = integrate.quad(stats.norm.pdf, -z/math.sqrt(1+rho), z/math.sqrt(1+rho), epsabs=1e-12)
    assert result["coverage"] == pytest.approx(reference, abs=1e-14)
    assert .844 < result["coverage"] < .846
    assert gaussian_linear_coverage(actual, actual, [1., 1.], probability=probability)["coverage"] == pytest.approx(probability, abs=1e-15)
    # The same covariance propagation applies to the parallel-sum gradient.
    left, right = 2., 3.
    gradient = np.array([right**2, left**2])/(left+right)**2
    consistent = gaussian_linear_coverage(actual, actual, gradient, probability=probability)
    assert consistent["coverage"] == pytest.approx(probability, abs=1e-15)


def test_exact_centered_gaussian_part_and_sum_coverage_selects_the_cross_covariance():
    actual = np.array([[1., .35], [.35, 2.]])
    coverages = []
    for cross in (-.4, .35, .8):
        quoted = np.array([[1., cross], [cross, 2.]])
        for gradient in ([1., 0.], [0., 1.]):
            assert gaussian_linear_coverage(actual, quoted, gradient, probability=.9)["coverage"] == pytest.approx(.9, abs=1e-15)
        coverages.append(gaussian_linear_coverage(actual, quoted, [1., 1.], probability=.9)["coverage"])
    assert coverages[0] < .9
    assert coverages[1] == pytest.approx(.9, abs=1e-15)
    assert coverages[2] > .9


def test_validation_and_explicit_operation():
    for support, weights in (([0.], [1.]), ([1.], [-1.]), ([1.], [0.])):
        with pytest.raises(ValueError):
            FinitePositiveLaw(support, weights)
    law = FinitePositiveLaw([1.], [1.])
    with pytest.raises(ValueError):
        compose_independent(law, law, operation="unspecified")
    with pytest.raises(ValueError):
        pushforward_joint([1., 2.], [1.], [[1.]], operation="addition")
    with pytest.raises(ValueError):
        warp_finite_law(law, eta=0.)
    with pytest.raises(ValueError):
        gaussian_linear_coverage(np.eye(2), np.array([[1., 2.], [2., 1.]]),
                                  [1., 1.], probability=.95)

"""Independent Gaussian integration and fixed-stream certificates for T52."""
import math

import numpy as np
import pytest

from met.sparse_identification import (active_support_scatter_moments,
    harmonic_activation_probabilities, harmonic_expected_active_bounds,
    harmonic_sparse_certificate, positive_sparse_log_ratio,
    positive_sparse_product_certificate)


def _normal_quadrature(order, dimension):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    indices = np.asarray(list(np.ndindex(*(order,)*dimension)))
    return math.sqrt(2)*nodes[indices], np.prod(weights[indices], axis=1)/math.pi**(dimension/2)


@pytest.mark.parametrize("amplitude", [.25, .8, 1.6])
def test_positive_sparse_moments_against_independent_gaussian_quadrature(amplitude):
    points, weights = _normal_quadrature(28, 2)
    probability, first_angle, second_angle = .3, .2, 1.1

    def ratio(angle):
        direction = np.array([math.cos(angle), math.sin(angle)])
        return 1-probability+probability*np.exp(amplitude*(points@direction)-amplitude**2/2)

    first, second = ratio(first_angle), ratio(second_angle)
    assert weights@first == pytest.approx(1., abs=3e-14)
    np.testing.assert_allclose(np.exp(positive_sparse_log_ratio(points, amplitude=amplitude,
        angle=first_angle, probability=probability)), first, rtol=3e-14, atol=3e-14)
    certificate = positive_sparse_product_certificate([probability], amplitude=amplitude,
        correlation=math.cos(first_angle-second_angle))
    assert weights@(first*first) == pytest.approx(math.exp(certificate["log_second_moment"]), rel=3e-13)
    assert weights@(first*second) == pytest.approx(math.exp(certificate["log_cross_moment"]), rel=3e-13)
    chi_squared = weights@((first-1)**2)
    assert chi_squared == pytest.approx(probability**2*math.expm1(amplitude**2), rel=3e-13)


def test_two_reading_product_certificate_by_four_dimensional_gaussian_integration():
    points, weights = _normal_quadrature(12, 4)
    amplitude, probabilities, angles = .7, [.2, .35], [.3, 1.]
    ratios = []
    for angle in angles:
        direction = np.array([math.cos(angle), math.sin(angle)])
        product = np.ones(len(points))
        for i, probability in enumerate(probabilities):
            product *= 1-probability+probability*np.exp(
                amplitude*(points[:, 2*i:2*i+2]@direction)-amplitude**2/2)
        ratios.append(product)
    certificate = positive_sparse_product_certificate(probabilities, amplitude=amplitude,
        correlation=math.cos(angles[0]-angles[1]))
    integrated_bound = math.sqrt((weights@((ratios[0]-ratios[1])**2))/4)
    assert certificate["tv_direct_upper"] == pytest.approx(integrated_bound, rel=2e-12)
    integrated_tv = weights@np.abs(ratios[0]-ratios[1])/2
    assert integrated_tv < certificate["tv_upper"]


def test_one_fixed_harmonic_schedule_keeps_positive_testing_error_at_every_prefix():
    infinite = harmonic_sparse_certificate(amplitude=.5, activation_scale=.5)
    expected_exponent = math.expm1(.25)*.25*(math.pi**2/6-1)
    assert infinite["log_second_moment_upper"] == pytest.approx(expected_exponent)
    assert infinite["tv_upper"] == pytest.approx(math.sqrt(math.expm1(expected_exponent)/2))
    assert infinite["tv_upper"] < .16
    assert infinite["binary_error_lower"] > .42
    expected_active = []
    for count in (20, 200, 2000):
        probabilities = harmonic_activation_probabilities(scale=.5, count=count)
        assert probabilities[0] == .25
        finite = positive_sparse_product_certificate(probabilities, amplitude=.5, correlation=.8)
        assert finite["tv_upper"] <= finite["tv_envelope_upper"] <= infinite["tv_upper"]
        assert finite["binary_error_lower"] >= infinite["binary_error_lower"]
        bounds = harmonic_expected_active_bounds(scale=.5, count=count)
        assert bounds["lower"] <= sum(probabilities) <= bounds["upper"]
        expected_active.append(finite["expected_active_count"])
    assert expected_active[0] < expected_active[1] < expected_active[2]
    # Divergence of this lower bound uses no allocated large array.
    assert harmonic_expected_active_bounds(scale=.5, count=10**20)["lower"] > 22


@pytest.mark.parametrize("amplitude", [.3, 1.2])
def test_known_active_support_scatter_contrast_by_independent_fourth_moment_integral(amplitude):
    noise, weights = _normal_quadrature(5, 2)
    direction = np.array([math.cos(.6), math.sin(.6)])
    points = noise+amplitude*direction
    scatter = points[:, :, None]*points[:, None, :]
    mean = np.einsum("n,nij->ij", weights, scatter)
    eigenvalues = np.linalg.eigvalsh(mean)
    fluctuation = weights@np.sum((scatter-mean)**2, axis=(1, 2))
    result = active_support_scatter_moments(amplitude=amplitude, active_count=7)
    assert result["mean_eigenvalue_gap"] == pytest.approx(7*(eigenvalues[1]-eigenvalues[0]), abs=2e-14)
    assert result["fluctuation_second_moment"] == pytest.approx(7*fluctuation, rel=2e-14)
    assert result["oracle_information"] == pytest.approx(7*amplitude**4/(1+amplitude**2))
    # Relative scatter noise tends to zero with the oracle active count.
    larger = active_support_scatter_moments(amplitude=amplitude, active_count=700)
    assert (larger["fluctuation_second_moment"]/larger["mean_eigenvalue_gap"]**2
        == pytest.approx(result["fluctuation_second_moment"]/result["mean_eigenvalue_gap"]**2/100))


def test_tiny_moments_and_nearly_equal_directions_do_not_become_false_identical_laws():
    tiny = positive_sparse_product_certificate([.2, .3], amplitude=1e-8, correlation=.8)
    assert tiny["log_second_moment"] == pytest.approx(.13e-16, rel=1e-14, abs=0.)
    assert tiny["tv_upper"] > 0
    nearly_equal = positive_sparse_product_certificate([.4], amplitude=.5,
        correlation=np.nextafter(1., 0.))
    assert nearly_equal["tv_upper"] > 0
    assert positive_sparse_product_certificate([.4], amplitude=.5, correlation=1.)["tv_upper"] == 0
    for probabilities, amplitude in (([], .5), ([0., 0.], .5), ([.3], 0.)):
        result = positive_sparse_product_certificate(probabilities, amplitude=amplitude, correlation=.8)
        assert result["tv_upper"] == 0 and result["binary_error_lower"] == .5


@pytest.mark.parametrize("amplitude", [.5, 26.])
def test_tiny_activation_probability_retains_representable_testing_bounds(amplitude):
    probability, correlation = 1e-170, .8
    # Factor p out before the square root; squaring p would underflow.
    squared = amplitude*amplitude
    reference = probability*math.sqrt(math.exp(correlation*squared)
        *math.expm1((1-correlation)*squared)/2)
    result = positive_sparse_product_certificate([probability], amplitude=amplitude,
        correlation=correlation)
    assert result["tv_direct_upper"] == pytest.approx(reference, rel=2e-12, abs=0.)
    assert result["tv_upper"] > 0.
    if amplitude == 26.:
        assert result["tv_upper"] == pytest.approx(4.375421388703658e-24, rel=2e-12, abs=0.)


def test_underflowing_amplitude_square_and_harmonic_scale_keep_their_tv_bounds():
    result = positive_sparse_product_certificate([.3], amplitude=1e-200, correlation=.8)
    reference = .3e-200*math.sqrt((1-.8)/2)
    assert result["tv_upper"] == pytest.approx(reference, rel=2e-12, abs=0.)
    harmonic = harmonic_sparse_certificate(amplitude=.5, activation_scale=1e-170)
    reference = 1e-170*math.sqrt(math.expm1(.25)*(math.pi**2/6-1)/2)
    assert harmonic["tv_upper"] == pytest.approx(reference, rel=2e-12, abs=0.)


def test_density_endpoints_and_range_validation():
    assert positive_sparse_log_ratio([2., 3.], amplitude=.5, angle=.2, probability=0.) == 0.
    expected = .5*(2*math.cos(.2)+3*math.sin(.2))-.125
    assert positive_sparse_log_ratio([2., 3.], amplitude=.5, angle=.2, probability=1.) == pytest.approx(expected)
    for probabilities in ([-.1], [1.1], [math.nan], [[.2]]):
        with pytest.raises(ValueError):
            positive_sparse_product_certificate(probabilities, amplitude=.5, correlation=.8)
    for amplitude in (-.1, math.nan, 1e200):
        with pytest.raises(ValueError):
            harmonic_sparse_certificate(amplitude=amplitude, activation_scale=.5)
    with pytest.raises(ValueError):
        positive_sparse_product_certificate([.2], amplitude=.5, correlation=-.1)
    for scale, count in ((0., 10), (1.1, 10), (.5, 0), (.5, True)):
        with pytest.raises(ValueError):
            harmonic_activation_probabilities(scale=scale, count=count)
    with pytest.raises(ValueError):
        positive_sparse_log_ratio([1., 2., 3.], amplitude=.5, angle=.2, probability=.3)

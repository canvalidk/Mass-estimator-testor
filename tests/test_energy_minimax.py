"""Independent Gaussian integrals and uniform regime bounds for T53."""
import math

import numpy as np
import pytest

from met.energy_minimax import (homogeneous_energy_testing_certificate,
    one_reading_relative_chi_square_bound, scatter_minimax_upper_bound)


def _quadrature(order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    first, second = np.meshgrid(nodes, nodes, indexing="ij")
    return (math.sqrt(2)*np.column_stack([first.ravel(), second.ravel()]),
            np.outer(weights, weights).ravel()/math.pi)


@pytest.mark.parametrize("energy", [.03, .4, 1., 1.2])
@pytest.mark.parametrize("correlation", [.2, .8, .99])
def test_relative_chi_square_and_kl_bounds_against_independent_gaussian_quadrature(energy, correlation):
    points, weights = _quadrature(40)
    radius = math.sqrt(energy)
    first = math.exp(-energy/2)*np.cosh(radius*points[:, 0])
    direction = np.array([correlation, math.sqrt(1-correlation**2)])
    second = math.exp(-energy/2)*np.cosh(radius*(points@direction))
    chi_squared = weights@((second-first)**2/first)
    kl = weights@(second*np.log(second/first))
    raw_second_moment = weights@((second-first)**2)
    bound = one_reading_relative_chi_square_bound(energy=energy, correlation=correlation)
    assert bound["chi_square_upper"] == pytest.approx(math.exp(energy/2)*raw_second_moment, rel=2e-12)
    assert 0 <= chi_squared <= bound["chi_square_upper"]*(1+1e-12)
    assert 0 <= kl <= math.log1p(chi_squared)*(1+1e-10)
    if energy <= 1:
        assert bound["chi_square_upper"] <= bound["small_signal_chi_square_upper"]*(1+1e-14)
    else:
        assert bound["small_signal_chi_square_upper"] is None
    certificate = homogeneous_energy_testing_certificate(count=3, energy=3*energy, correlation=correlation)
    assert 3*kl <= certificate["kl_upper"]
    assert certificate["kl_upper"] == pytest.approx(3*math.log1p(bound["chi_square_upper"]))


def test_affinity_tv_and_binary_error_rationalization_remain_informative():
    certificate = homogeneous_energy_testing_certificate(count=1000, energy=20., correlation=.8)
    k = certificate["kl_upper"]
    expected_tv = math.sqrt(1-math.exp(-k))
    assert certificate["tv_upper"] == pytest.approx(expected_tv)
    assert certificate["binary_error_lower"] == pytest.approx((1-expected_tv)/2)
    # Finite KL may produce TV rounded to one; log error still states the bound.
    large = homogeneous_energy_testing_certificate(count=10**6, energy=10**5, correlation=.8)
    assert large["tv_upper"] == 1.
    assert math.isfinite(large["log_binary_error_lower"])
    assert large["log_binary_error_lower"] < -700


@pytest.mark.parametrize("scale", [.5, 2., 5.])
def test_every_bounded_energy_over_sqrt_count_subsequence_has_a_fixed_positive_error_bound(scale):
    correlation = .8
    uniform_kl = math.exp(.5)*math.cosh(1.)*(1-correlation**2)*scale**2
    fixed_error = math.exp(-uniform_kl)/(2*(1+math.sqrt(-math.expm1(-uniform_kl))))
    for count in (100, 1000, 10000):
        certificate = homogeneous_energy_testing_certificate(count=count,
            energy=scale*math.sqrt(count), correlation=correlation)
        assert certificate["kl_upper"] <= certificate["small_signal_kl_upper"] <= uniform_kl*(1+2e-15)
        assert certificate["binary_error_lower"] >= fixed_error*(1-2e-14)


def test_scatter_fluctuation_identity_and_allocation_independent_upper_bound():
    noise, weights = _quadrature(5)
    radii = [.3, 2., .8]
    mean_scatter, integrated_variance = np.zeros((2, 2)), 0.
    direction = np.array([math.cos(.7), math.sin(.7)])
    for radius in radii:
        points = noise+radius*direction
        products = points[:, :, None]*points[:, None, :]
        mean = np.einsum("n,nij->ij", weights, products)
        mean_scatter += mean
        integrated_variance += weights@np.sum((products-mean)**2, axis=(1, 2))
    energy = sum(r*r for r in radii)
    result = scatter_minimax_upper_bound(count=len(radii), energy=energy, angle_tolerance=.4)
    assert result["scatter_fluctuation_second_moment"] == pytest.approx(integrated_variance, rel=2e-14)
    assert np.diff(np.linalg.eigvalsh(mean_scatter))[0] == pytest.approx(energy, abs=1e-14)
    assert result["information"] == pytest.approx(energy**2/(len(radii)+energy))


def test_uniform_sufficient_bound_vanishes_in_the_identifying_energy_regime():
    upper = []
    for count in (10**4, 10**8, 10**12):
        energy = count**.75
        result = scatter_minimax_upper_bound(count=count, energy=energy, angle_tolerance=.5)
        independently = min(1., 24*(count+energy)/(energy**2*math.sin(.5)**2))
        assert result["line_error_probability_upper"] == pytest.approx(independently)
        upper.append(result["line_error_probability_upper"])
    assert upper[0] > upper[1] > upper[2]
    assert upper[2] < .001


def test_zero_signal_identical_directions_tiny_energy_and_validation():
    for energy, correlation in ((0., .8), (1., 1.)):
        result = homogeneous_energy_testing_certificate(count=10, energy=energy, correlation=correlation)
        assert result["tv_upper"] == 0. and result["binary_error_lower"] == .5
    tiny = one_reading_relative_chi_square_bound(energy=1e-10, correlation=.8)
    assert tiny["chi_square_upper"] == pytest.approx((1-.8**2)*1e-20, rel=1e-9, abs=0.)
    assert scatter_minimax_upper_bound(count=10, energy=0., angle_tolerance=.3)["line_error_probability_upper"] == 1.
    for count, energy, correlation in ((0, 1., .8), (True, 1., .8), (10, -1., .8),
                                       (10, math.nan, .8), (10, 1., -.1)):
        with pytest.raises(ValueError):
            homogeneous_energy_testing_certificate(count=count, energy=energy, correlation=correlation)
    with pytest.raises(ValueError):
        scatter_minimax_upper_bound(count=10, energy=1., angle_tolerance=0.)


def test_tiny_per_reading_information_survives_count_and_square_root_operations():
    for count, energy in ((10**300, 1e120), (1, 1e-170), (10**300, 1e-100)):
        result = homogeneous_energy_testing_certificate(count=count, energy=energy, correlation=.8)
        # Here t is tiny: KL_upper~(1-rho²)*Delta²/N, and TV_upper~sqrt(KL).
        log_reference = (math.log1p(-.8)+math.log1p(.8)+2*math.log(energy)-math.log(count))/2
        assert result["tv_upper"] == pytest.approx(math.exp(log_reference), rel=2e-12, abs=0.)
        assert result["tv_upper"] > 0
        assert result["log_kl_upper"] == pytest.approx(2*log_reference, abs=2e-12)
    smallest = one_reading_relative_chi_square_bound(energy=np.nextafter(0., 1.), correlation=.8)
    assert math.isfinite(smallest["log_chi_square_upper"])


def test_tiny_angle_does_not_divide_by_an_underflowed_squared_sine():
    result = scatter_minimax_upper_bound(count=1, energy=1., angle_tolerance=1e-200)
    assert result["line_error_probability_upper"] == 1.

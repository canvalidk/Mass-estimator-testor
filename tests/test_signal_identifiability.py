"""Independent Gaussian quadrature and deterministic T41 testing certificates."""
import math

import numpy as np
import pytest

from met.signal_identifiability import (harmonic_energy_bounds, harmonic_radii,
                                        harmonic_testing_certificate, log_cosh,
                                        log_cross_moment, product_testing_certificate,
                                        scatter_fluctuation_second_moment,
                                        sign_mixture_log_ratio)


def _gaussian_quadrature(order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    a, b = np.meshgrid(nodes, nodes, indexing="ij")
    points = math.sqrt(2)*np.stack([a.ravel(), b.ravel()], axis=-1)
    w = np.outer(weights, weights).ravel()/math.pi
    return points, w


@pytest.mark.parametrize("radius", [.25, .8, 1.6])
def test_exact_cross_moment_by_independent_normal_quadrature(radius):
    points, weights = _gaussian_quadrature(32)
    theta = .2

    def ratio(angle):
        direction = np.array([math.cos(angle), math.sin(angle)])
        z = radius*(points @ direction)
        return math.exp(-radius**2/2)*(np.exp(z)+np.exp(-z))/2

    first = ratio(theta)
    assert weights @ first == pytest.approx(1., abs=2e-13)
    for other in (theta, 1.1, theta+math.pi/2):
        actual = weights @ (first*ratio(other))
        expected = math.exp(log_cross_moment(radius, angle_difference=theta-other))
        assert actual == pytest.approx(expected, rel=3e-13, abs=2e-13)
    np.testing.assert_allclose(np.exp(sign_mixture_log_ratio(points, radius=radius, angle=theta)),
                               first, rtol=2e-13, atol=2e-14)


def test_log_operations_remain_accurate_in_small_and_large_regimes():
    assert log_cosh(1e-12) == pytest.approx(5e-25, rel=1e-14, abs=0)
    assert log_cosh(1000.) == pytest.approx(1000-math.log(2))
    assert log_cosh(-1000.) == log_cosh(1000.)
    certificate = product_testing_certificate([1e3]*4, angle_difference=.7)
    assert certificate["tv_upper"] == 1
    assert certificate["binary_error_lower"] == 0
    assert math.isfinite(certificate["log_second_moment"])


def test_finite_product_moments_and_both_testing_bounds():
    radii, difference = np.array([.1, .3, .6]), .7
    certificate = product_testing_certificate(radii, angle_difference=difference)
    same = math.prod(math.cosh(r*r) for r in radii)
    cross = math.prod(math.cosh(r*r*math.cos(difference)) for r in radii)
    assert certificate["chi_square_null"] == pytest.approx(same-1)
    assert certificate["tv_direct_upper"] == pytest.approx(math.sqrt((same-cross)/2))
    assert certificate["tv_direct_upper"] <= certificate["tv_triangle_upper"]
    assert certificate["binary_error_lower"] == pytest.approx(.5*(1-certificate["tv_upper"]))
    assert product_testing_certificate(radii, angle_difference=0)["tv_upper"] == 0
    assert product_testing_certificate([.5], angle_difference=1e-10)["tv_upper"] > 0
    assert product_testing_certificate([], angle_difference=.7)["binary_error_lower"] == .5


def test_harmonic_energy_grows_but_testing_error_stays_positive():
    infinite = harmonic_testing_certificate(scale=.5)
    assert infinite["binary_error_triangle_lower"] > .38
    assert infinite["binary_error_direct_lower"] > .41
    energies, moments = [], []
    for n in (100, 1000, 10000):
        radii = harmonic_radii(scale=.5, count=n)
        finite = product_testing_certificate(radii, angle_difference=math.atan(2)-math.atan(.5))
        bounds = harmonic_energy_bounds(scale=.5, count=n)
        assert bounds["lower"] <= finite["total_energy"] <= bounds["upper"]
        assert finite["log_second_moment"] <= infinite["log_second_moment_upper"]
        assert finite["tv_triangle_upper"] <= infinite["tv_triangle_upper"]
        assert finite["tv_direct_upper"] <= infinite["tv_direct_upper"]
        assert finite["binary_error_lower"] >= infinite["binary_error_direct_lower"]
        energies.append(finite["total_energy"])
        moments.append(finite["log_second_moment"])
    assert energies[0] < energies[1] < energies[2]
    assert moments[0] < moments[1] < moments[2]
    assert harmonic_energy_bounds(scale=.5, count=10**12)["lower"] > 6


@pytest.mark.parametrize("radius", [0., .4, 2.])
def test_scatter_fluctuation_by_independent_normal_quadrature(radius):
    noise, weights = _gaussian_quadrature(8)
    mean = radius*np.array([math.cos(.6), math.sin(.6)])
    points = noise+mean
    scatter = points[:, :, None]*points[:, None, :]
    centered = scatter-(np.eye(2)+np.outer(mean, mean))
    actual = weights @ np.sum(centered**2, axis=(1, 2))
    assert actual == pytest.approx(scatter_fluctuation_second_moment([radius]), rel=2e-14)
    assert scatter_fluctuation_second_moment([radius, radius]) == pytest.approx(2*actual)


def test_input_validation():
    for invalid in ([-1.], [math.nan], [[1.]]):
        with pytest.raises(ValueError):
            product_testing_certificate(invalid, angle_difference=.5)
    with pytest.raises(ValueError):
        harmonic_radii(scale=.5, count=0)
    with pytest.raises(ValueError):
        harmonic_testing_certificate(scale=0)
    with pytest.raises(ValueError):
        sign_mixture_log_ratio([1., 2., 3.], radius=.5, angle=.2)

"""Independent identities for changing declared priors on fixed measured inputs.

Expected relationships follow from direct scalar integrals, reciprocal symmetry
and the half-Cauchy CDF. They do not choose a desired physical response to fixed
incompatible measurements.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.inspection.case import Observation
from met.inspection.references import angle_reference, cloud_readout_reference, equal_perpendicular_point


def observation(force, acceleration, sigma=0.1):
    return Observation("r", force, acceleration, sigma ** 2 * np.eye(4))


@pytest.mark.parametrize("center", [0.1, 1.0, 3.0, 10.0])
def test_perpendicular_point_closed_form_matches_direct_mass_integrals(center):
    """Constant evidence gives pi_c; the independent cloud factor is (1+m²)^-1/2."""
    def denominator(mass):
        return center / (center ** 2 + mass ** 2) / math.sqrt(1 + mass ** 2)

    d = quad(denominator, 0, np.inf, epsabs=1e-12, epsrel=1e-11)[0]
    n = quad(lambda mass: mass * denominator(mass), 0, np.inf, epsabs=1e-12, epsrel=1e-11)[0]
    assert equal_perpendicular_point(center) == pytest.approx(n / d, rel=2e-11)


@pytest.mark.parametrize("center", [0.1, 1.0, 10.0])
def test_perpendicular_law_returns_each_prior_and_point_matches_closed_form(center):
    masses = center * np.array([0.1, 1.0, 10.0])
    r = observation([1, 0], [0, 1])
    cdf = angle_reference(r, center, 0.01, masses)
    np.testing.assert_allclose(cdf["cdf"], 2 / math.pi * np.arctan(masses / center), atol=2e-12)
    point = cloud_readout_reference(r, center, 0.01)
    assert point["point"] == pytest.approx(equal_perpendicular_point(center), rel=2e-10)
    assert point["quadrature_warnings"] == []
    assert point["numerator"] > 0 and point["denominator"] > 0


def test_perpendicular_cloud_point_is_independent_of_common_noise_scale():
    """Both conditional mean and SD factor by the same (1+m²)^-1/2."""
    points = [cloud_readout_reference(observation([1, 0], [0, 1], sigma), 0.1, 1)["point"]
              for sigma in (0.1, 0.4)]
    np.testing.assert_allclose(points, equal_perpendicular_point(0.1), rtol=2e-10)


@pytest.mark.parametrize("force,acceleration", [([1, 0], [-1, 0]), ([1, 0], [0, 1])])
def test_reciprocal_priors_exchange_quantiles_and_cloud_readouts(force, acceleration):
    """L(1/m)=L(m), pi_c(1/m)/m²=pi_(1/c)(m), A(1/m)=m*A(m)."""
    r = observation(force, acceleration)
    masses = np.array([0.001, 0.01, 0.1, 1.0, 10.0])
    small = angle_reference(r, 0.1, 1, masses)["cdf"]
    large = angle_reference(r, 10, 1, 1 / masses)["cdf"]
    np.testing.assert_allclose(small + large, 1.0, atol=2e-10, rtol=0)
    low = cloud_readout_reference(r, 0.1, 1)
    high = cloud_readout_reference(r, 10, 1)
    assert low["point"] * high["point"] == pytest.approx(1, rel=2e-10)


def test_opposite_reference_agrees_with_direct_unscaled_angle_integrals():
    """The opposite unit-pair excess cost is 2 sin(phi)cos(phi)/sigma²."""
    sigma, center = 0.1, 0.1
    r = observation([1, 0], [-1, 0], sigma)

    def kernel(phi):
        co, si = math.cos(phi), math.sin(phi)
        return center / (center ** 2 * co ** 2 + si ** 2) * math.exp(-si * co / sigma ** 2)

    z = quad(kernel, 0, math.pi / 4, epsabs=1e-12)[0] + quad(kernel, math.pi / 4, math.pi / 2, epsabs=1e-12)[0]
    expected = quad(kernel, 0, math.pi / 4, epsabs=1e-12)[0] / z
    got = angle_reference(r, center, 1, [1.0])
    assert got["cdf"][0] == pytest.approx(expected, abs=2e-10)
    assert 0.98 < expected < 1


def test_opposite_point_at_neutral_prior_is_exactly_one():
    """Reciprocal cloud-size symmetry fixes the ratio of integrals at c=1."""
    got = cloud_readout_reference(observation([1, 0], [-1, 0]), 1, 1)
    assert got["point"] == pytest.approx(1, rel=2e-10)
    assert got["quadrature_warnings"] == []


@pytest.mark.parametrize("center", [0.0, -1.0, float("nan"), float("inf")])
def test_closed_form_requires_a_finite_positive_prior_center(center):
    with pytest.raises(ValueError):
        equal_perpendicular_point(center)

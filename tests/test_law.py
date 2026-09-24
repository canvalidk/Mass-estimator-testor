"""The per-reading law: closed forms against direct integration, and its splits."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import i0e

from met.law import ReadingSet, log_radial_kernel, mean_radius


def _angular_scaled(z, d):
    """< exp(z u) >_u * exp(-z), stable for large z."""
    if d == 1:
        return 0.5 * (1 + math.exp(-2 * z))
    if d == 2:
        return i0e(z)
    return (1 - math.exp(-2 * z)) / (2 * z) if z > 1e-9 else 1.0


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("h", [0.0, 0.3, 1.0, 3.0, 8.0, 20.0])
def test_kernel_and_mean_radius_match_direct_integration(d, h):
    def integrand(r, k):
        return r**k * math.exp(-(r - h) ** 2 / 2) * _angular_scaled(r * h, d)

    k1 = quad(integrand, 0, np.inf, args=(1,), limit=400, epsabs=0, epsrel=1e-13)[0]
    k2 = quad(integrand, 0, np.inf, args=(2,), limit=400, epsabs=0, epsrel=1e-13)[0]
    assert math.log(k1) + h * h / 2 == pytest.approx(float(log_radial_kernel(h, d)), abs=1e-11)
    assert k2 / k1 == pytest.approx(float(mean_radius(h, d)), rel=1e-11)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_small_h_is_continuous(d):
    hs = np.array([0.0, 1e-9, 1e-7, 1e-6, 1.0001e-6, 1e-5])
    assert np.allclose(log_radial_kernel(hs, d), log_radial_kernel(hs, d)[0], atol=1e-9)
    assert np.allclose(mean_radius(hs, d), math.sqrt(math.pi / 2), rtol=1e-9)


def _random_readings(rng, b, n, d):
    return ReadingSet(rng.normal(size=(b, n, d)) * 2, rng.normal(size=(b, n, d)),
                      rng.uniform(0.5, 2, size=(b, n)), rng.uniform(0.5, 2, size=(b, n)))


@pytest.mark.parametrize("d", [1, 2, 3])
def test_nuisance_splits_give_the_same_single_reading_law(d):
    rng = np.random.default_rng(1)
    readings = _random_readings(rng, 4, 3, d)
    u = np.tile(np.linspace(-6, 6, 101), (4, 1))
    totals = []
    for nuisance in ("radius", "acceleration", "force"):
        like, ref, cond = readings.reading_terms(u, nuisance, reference="flat")
        totals.append(like + ref)
        cond_radius = readings.reading_terms(u, "radius", reference="flat")[2]
        assert np.allclose(cond, cond_radius)
    for other in totals[1:]:
        diff = other - totals[0]
        assert np.allclose(diff, diff[..., :1], atol=1e-10)


def test_pooled_is_the_inverse_variance_mean():
    rng = np.random.default_rng(2)
    readings = ReadingSet(rng.normal(size=(2, 5, 3)), rng.normal(size=(2, 5, 3)),
                          np.full((2, 5), 2.0), np.full((2, 5), 0.5))
    pair = readings.pooled()
    assert np.allclose(pair.force[:, 0], readings.force.mean(axis=1))
    assert np.allclose(pair.acceleration[:, 0], readings.acceleration.mean(axis=1))
    assert np.allclose(pair.force_sd, 2.0 / math.sqrt(5))
    assert np.allclose(pair.acceleration_sd, 0.5 / math.sqrt(5))
    unequal = ReadingSet(np.array([[[1.0], [4.0]]]), np.array([[[1.0], [1.0]]]),
                         np.array([[1.0, 2.0]]), np.array([[1.0, 1.0]]))
    pooled = unequal.pooled()
    assert pooled.force[0, 0, 0] == pytest.approx((1 * 1 + 4 * 0.25) / 1.25)
    assert pooled.force_sd[0, 0] == pytest.approx(1 / math.sqrt(1.25))

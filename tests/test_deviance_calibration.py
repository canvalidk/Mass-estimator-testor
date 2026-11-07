"""T40's research maps: algebra, continuous scales and deterministic limits."""
import math

import numpy as np
import pytest
from scipy import integrate, optimize, stats

from met.deviance_calibration import (deviance_loss, effective_scale,
                                      local_gaussian_reference, pair_angle_cdf,
                                      pair_angle_quantile)


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_bounds_endpoints_and_extensive_homogeneity(family):
    d = np.r_[0.0, np.logspace(-12, 6, 30)]
    eta, a = 0.13, 1.7
    values = deviance_loss(d, eta, family=family, scale=a)
    assert np.all(values >= eta*d - 1e-12)
    assert np.all(values <= d + 1e-12)
    assert np.all(values >= d + a*math.log(eta) - 1e-12)
    assert values[0] == 0
    np.testing.assert_array_equal(deviance_loss(d, 0, family=family, scale=a), np.zeros_like(d))
    np.testing.assert_array_equal(deviance_loss(d, 1, family=family, scale=a), d)
    np.testing.assert_allclose(deviance_loss(7*d, eta, family=family, scale=7*a), 7*values,
                               rtol=2e-13, atol=2e-12)


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_semigroup_small_deviance_slope_and_tail_offset(family):
    d, a = np.r_[0.0, 1e-8, 0.1, 1, 10, 1e4, 1e6], 1.0
    first = deviance_loss(d, 0.2, family=family, scale=a)
    np.testing.assert_allclose(deviance_loss(first, 0.3, family=family, scale=a),
                               deviance_loss(d, 0.06, family=family, scale=a),
                               rtol=2e-13, atol=2e-10)
    assert deviance_loss(1e-8, 0.1, family=family, scale=a) / 1e-8 == pytest.approx(.1, rel=2e-7)
    values = deviance_loss(np.array([1e4, 1e6]), 0.1, family=family, scale=a)
    assert np.all(np.isfinite(values))
    assert values[-1] - 1e6 == pytest.approx(math.log(.1), abs=3e-6)
    assert abs(values[-1] - 1e6 - math.log(.1)) <= abs(values[0] - 1e4 - math.log(.1)) + 1e-10
    for value in (1e-20, 1e-100):
        assert deviance_loss(value, .1, family=family, scale=1e300) == pytest.approx(.1*value, rel=2e-14, abs=0)


def test_maps_are_distinct_despite_shared_constraints():
    assert deviance_loss(1., .1, family="logistic", scale=1.) == pytest.approx(.15856507874042916)
    assert deviance_loss(1., .1, family="lambert", scale=1.) == pytest.approx(.21847923056999352)
    # An independent implicit equation checks the Lambert branch without W.
    expected = optimize.brentq(lambda f: f + math.log(f) - (1 + math.log(.1)), .01, 1.)
    assert deviance_loss(1., .1, family="lambert", scale=1.) == pytest.approx(expected)


@pytest.mark.parametrize("family", ["rational", "exponential"])
def test_effective_scale_is_bounded_continuous_and_exactly_zero_safe(family):
    values, unit = np.array([0., .1, 1., 100.]), .7
    base = effective_scale(values, unit=unit, family=family)
    assert 0 < base < len(values)
    assert effective_scale(np.r_[values, np.zeros(20)], unit=unit, family=family) == base
    assert effective_scale(np.zeros(5), unit=unit, family=family) == 0
    assert effective_scale([], unit=unit, family=family) == 0
    assert effective_scale(np.tile(values, 2), unit=unit, family=family) == pytest.approx(2*base)
    assert effective_scale(3*values, unit=3*unit, family=family) == pytest.approx(base)
    eps = 1e-9
    changed = effective_scale(np.r_[values, eps], unit=unit, family=family)
    assert 0 < changed-base < 2*eps/unit
    assert deviance_loss(2., .2, family="logistic", scale=changed) == pytest.approx(
        deviance_loss(2., .2, family="logistic", scale=base), abs=2e-9)
    tiny = effective_scale([eps], unit=unit, family=family)
    assert 0 <= deviance_loss(eps, .2, family="logistic", scale=tiny) <= eps
    assert deviance_loss(0., .2, family="logistic", scale=0.) == 0


def test_logistic_fixed_scale_local_limit_is_underdispersed():
    result = local_gaussian_reference(.1, family="logistic", scale=1., probability=.95)
    assert result["variance"] == pytest.approx(.20188767802530033, abs=2e-11)
    assert result["half_width"] == pytest.approx(.8145822491841591, abs=2e-11)
    assert result["variance"] < 1
    assert 2*stats.norm.cdf(result["half_width"])-1 < .95
    # Independent direct logistic-ratio quadrature, with no deviance operator.
    kernel = lambda z: 1/(.9 + .1*math.exp(z*z/.2)) if z < 10 else 0.
    normalizer = 2*integrate.quad(kernel, 0, 10, epsabs=2e-12)[0]
    variance = 2*integrate.quad(lambda z: z*z*kernel(z), 0, 10, epsabs=2e-12)[0]/normalizer
    assert result["variance"] == pytest.approx(variance, abs=2e-11)


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_growing_scale_recovers_the_eta_powered_local_reference(family):
    target = stats.norm.ppf(.975)
    small = local_gaussian_reference(.5, family=family, scale=1., probability=.95)
    large = local_gaussian_reference(.5, family=family, scale=1000., probability=.95)
    assert small["half_width"] < large["half_width"] < target
    assert abs(large["half_width"]-target) < .005
    assert abs(large["variance"]-1) < .004


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_extensive_scale_keeps_the_stated_null_rate_bounded(family):
    errors, ns = [], np.array([100., 1e4, 1e6])
    for n in ns:
        eta, d = 1/math.sqrt(n), math.sqrt(n)*np.linspace(0, 3, 101)
        transformed = deviance_loss(d, eta, family=family, scale=n)
        assert np.all(transformed >= eta*d - 1e-14)
        errors.append(float(np.max(np.abs(transformed - eta*d))))
    assert errors[0] > errors[1] > errors[2] >= 0
    assert errors[0] < 1 and errors[-1] < .01
    scaled = np.asarray(errors)*np.sqrt(ns)
    assert scaled[-1] == pytest.approx(scaled[-2], rel=.02)


@pytest.mark.parametrize("family", ["logistic", "lambert"])
@pytest.mark.parametrize("scale_source", ["fixed", "oscillation"])
def test_pair_fixed_scale_regains_concentration_at_ordinary_nu(family, scale_source):
    arguments = dict(truth_angle=math.atan(2), nu_eff=10., error_size=3., family=family)
    medians, probabilities = [], []
    for radius in (100., 1000.):
        oscillation = 5*math.log1p((radius**2+9)/10)
        scale = (2. if scale_source == "fixed" else
                 effective_scale([oscillation, oscillation], unit=1., family="rational"))
        params = dict(arguments, push=radius, scale=scale)
        medians.append(math.tan(pair_angle_quantile(.5, **params)))
        low, le = pair_angle_cdf(math.atan(1.9), **params)
        high, he = pair_angle_cdf(math.atan(2.1), **params)
        probabilities.append(high-low)
        assert max(le, he) < 1e-7
    # At these radii the asymmetric prior boundaries already contribute less
    # than quadrature precision to the median; concentration is checked below.
    np.testing.assert_allclose(medians, [2., 2.], atol=2e-9, rtol=0)
    assert probabilities[-1] > probabilities[0]
    assert probabilities[-1] > .999
    # eta=1 is unchanged, so the intrinsic nu<1 single-reading failure remains.
    np.testing.assert_array_equal(deviance_loss(np.arange(5.), 1., family=family, scale=2.), np.arange(5.))


def test_explicit_choices_and_valid_domains_are_required():
    with pytest.raises(TypeError):
        deviance_loss(1., .1)
    with pytest.raises(TypeError):
        effective_scale([1.])
    for kwargs in (dict(family="missing", scale=1), dict(family="logistic", scale=-1),
                   dict(family="logistic", scale=0)):
        with pytest.raises(ValueError):
            deviance_loss(1., .1, **kwargs)
    with pytest.raises(ValueError):
        deviance_loss(-1., .1, family="logistic", scale=1.)
    with pytest.raises(ValueError):
        deviance_loss(1., 1.1, family="logistic", scale=1.)
    with pytest.raises(ValueError):
        effective_scale([1.], unit=0., family="rational")

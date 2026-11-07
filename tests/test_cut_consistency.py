"""T39: independent population integrals and the fixed-N tail counterexample.

These are deterministic checks of a proof, not a coverage ranking. In
particular, a correct population optimum does not prove every L11 regime.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.cut_consistency import (aligned_angle_cdf, aligned_angle_quantile,
                                 heavy_tail_limit_cdf, heavy_tail_limit_quantile,
                                 opposed_pair_angle_cdf, opposed_pair_angle_quantile,
                                 opposed_pair_eta,
                                 population_loss_gap)
from met.full import isotropic
from met.ours32 import ours32


@pytest.mark.parametrize("k", [1, 2, 3])
@pytest.mark.parametrize("noise_scale", [0.4, 2.0])
def test_population_gap_matches_independent_noncentral_chi_square_density(k, noise_scale):
    nu, delta = 4.0, 1.7
    coefficient = (nu + k - 1) / 2
    loss = lambda q: coefficient * np.log1p(noise_scale * q / nu)
    # SciPy uses the Bessel density; the reference uses a Laplace transform.
    central = stats.chi2.expect(loss, args=(k,), epsabs=2e-10)
    shifted = stats.ncx2.expect(loss, args=(k, delta / noise_scale), epsabs=2e-10)
    gap, error = population_loss_gap(k, nu, delta, noise_scale=noise_scale)
    assert gap == pytest.approx(shifted - central, abs=2e-9)
    assert error < 1e-9


@pytest.mark.parametrize("k", [1, 2, 3])
@pytest.mark.parametrize("nu", [0.5, 10.0])
def test_population_optimum_is_at_truth_for_every_positive_scale(k, nu):
    # Integrating over these scales is a non-Gaussian finite scale mixture.
    # The theorem covers any positive scale law with a finite log expectation.
    theta0, radius = math.atan(2), 1.3
    for scale in (0.2, 1.0, 4.0):
        at_truth = population_loss_gap(k, nu, 0.0, noise_scale=scale)
        assert at_truth == (0.0, 0.0)
        for displacement in (-0.5, -0.2, 0.2):
            delta = radius**2 * math.sin(displacement)**2
            gap, error = population_loss_gap(k, nu, delta, noise_scale=scale)
            assert gap > 100 * error
        # A true zero-push reading has delta=0 for every hypothesis.
        for theta in (0.1, theta0, 1.4):
            assert population_loss_gap(k, nu, 0 * math.sin(theta - theta0)**2,
                                       noise_scale=scale) == (0.0, 0.0)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_population_curvature_has_positive_chi_square_reference(k):
    nu, radius, scale = 3.0, 1.7, 0.8
    coefficient = (nu + k - 1) / 2
    reciprocal = stats.chi2.expect(lambda q: 1 / (nu + scale * q), args=(k + 2,))
    expected_curvature = 2 * coefficient * radius**2 * reciprocal
    h = 1e-3
    gap, _ = population_loss_gap(k, nu, radius**2 * math.sin(h)**2, noise_scale=scale)
    assert 2 * gap / h**2 == pytest.approx(expected_curvature, rel=2e-6)
    assert expected_curvature > 0


def _singular_integral(length, nu):
    # x = length*u^(1/(1-nu)) removes the endpoint singularity. This
    # quadrature is independent of the incomplete-beta implementation.
    power = 1 / (1 - nu)
    return length**(1 - nu) * power * integrate.quad(
        lambda u: np.sinc(length * u**power / math.pi)**(-nu),
        0, 1, epsabs=2e-12, epsrel=2e-12)[0]


@pytest.mark.parametrize("nu", [0.2, 0.5, 0.9])
def test_tail_limit_beta_integral_matches_removed_singularity_quadrature(nu):
    theta0 = math.atan(2)
    left = _singular_integral(theta0, nu)
    right = _singular_integral(math.pi / 2 - theta0, nu)
    assert heavy_tail_limit_cdf(theta0, truth_angle=theta0, tail_index=nu) == pytest.approx(
        left / (left + right), abs=2e-11)


def test_allowed_heavy_tail_limit_has_wrong_point_at_infinite_signal():
    theta0 = math.atan(2)
    median = heavy_tail_limit_quantile(0.5, truth_angle=theta0, tail_index=0.5)
    assert math.tan(median) == pytest.approx(1.8219669037884, abs=2e-11)
    assert heavy_tail_limit_cdf(theta0, truth_angle=theta0, tail_index=0.5) == pytest.approx(
        0.6113338264869, abs=2e-11)
    # Positive probability remains outside a fixed neighbourhood of truth.
    band = heavy_tail_limit_cdf(theta0 + 0.05, truth_angle=theta0, tail_index=0.5) - \
        heavy_tail_limit_cdf(theta0 - 0.05, truth_angle=theta0, tail_index=0.5)
    assert band < 0.5


def test_finite_signal_independent_quadrature_approaches_the_wrong_limit():
    theta0 = math.atan(2)
    limit = math.tan(heavy_tail_limit_quantile(0.5, truth_angle=theta0, tail_index=0.5))
    medians = [math.tan(aligned_angle_quantile(0.5, truth_angle=theta0, push=push, nu_eff=0.5))
               for push in (10.0, 100.0, 1000.0, 10000.0)]
    np.testing.assert_allclose(medians, [1.5172050552803, 1.7301265688584,
                                       1.7942704379789, 1.8133707243065], atol=2e-11)
    errors = [abs(median - limit) for median in medians]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 0.01
    # At a resolved finite signal this is the production law too.
    push = 10.0
    reading = isotropic([push * math.sin(theta0)], [push * math.cos(theta0)], 1.0, 1.0)
    result = ours32([reading], c=1.0, nu_eff=0.5)
    assert result["eta"] == pytest.approx(1.0, abs=1e-8)
    assert result["median"] == pytest.approx(medians[0], abs=2e-6)


def test_subcritical_counterexample_respects_swap_and_is_not_a_null_world():
    theta0 = math.atan(2)
    direct = heavy_tail_limit_quantile(0.5, truth_angle=theta0, tail_index=0.5)
    reflected = heavy_tail_limit_quantile(0.5, truth_angle=math.pi / 2 - theta0, tail_index=0.5)
    assert math.tan(direct) * math.tan(reflected) == pytest.approx(1.0, abs=2e-12)
    assert direct < theta0
    assert reflected > math.pi / 2 - theta0
    assert heavy_tail_limit_quantile(0.5, truth_angle=math.pi / 4, tail_index=0.5) == pytest.approx(
        math.pi / 4, abs=2e-12)


def test_above_threshold_single_reading_law_concentrates():
    theta0, delta = math.atan(2), 0.05
    probabilities = []
    for radius in (10.0, 100.0, 1000.0):
        high = aligned_angle_cdf(theta0 + delta, truth_angle=theta0, push=radius, nu_eff=3.0)[0]
        low = aligned_angle_cdf(theta0 - delta, truth_angle=theta0, push=radius, nu_eff=3.0)[0]
        probabilities.append(high - low)
    assert probabilities[0] < probabilities[1] < probabilities[2]
    assert probabilities[-1] > 0.999


def test_bartlett_power_creates_a_subcritical_tail_at_ordinary_nu():
    theta0, nu, error = math.atan(2), 10.0, 3.0
    exponent = 2 * (nu - error**2) / error**2
    assert exponent == pytest.approx(2 / 9)
    limit = math.tan(heavy_tail_limit_quantile(0.5, truth_angle=theta0, tail_index=exponent))
    assert limit == pytest.approx(1.3107867539725, abs=2e-11)
    assert heavy_tail_limit_cdf(theta0, truth_angle=theta0, tail_index=exponent) == pytest.approx(
        0.6654942769223355, abs=2e-11)
    medians = []
    for radius in (100.0, 1000.0, 10000.0):
        eta = opposed_pair_eta(push=radius, nu_eff=nu, error_size=error)
        assert eta == pytest.approx(1 / 90 - 19 / (10 * radius**2), rel=2e-14)
        medians.append(math.tan(opposed_pair_angle_quantile(
            0.5, truth_angle=theta0, push=radius, nu_eff=nu,
            error_size=error, calibrate=True)))
    np.testing.assert_allclose(medians, [1.2553365823020, 1.3021637339378,
                                       1.3093548686541], atol=2e-10)
    assert abs(medians[-1] - limit) < 0.002


def test_same_readings_concentrate_without_eta_but_not_with_it():
    theta0, radius = math.atan(2), 1000.0
    arguments = dict(truth_angle=theta0, push=radius, nu_eff=10.0, error_size=3.0)
    low, high = math.atan(1.9), math.atan(2.1)
    probabilities = {}
    for calibrated in (False, True):
        probabilities[calibrated] = opposed_pair_angle_cdf(high, **arguments, calibrate=calibrated)[0] - \
            opposed_pair_angle_cdf(low, **arguments, calibrate=calibrated)[0]
    assert probabilities[False] > 0.999
    assert probabilities[True] < 0.06
    # A resolved finite-push production check uses the identical inputs.
    radius = 100.0
    u = np.array([math.sin(theta0), math.cos(theta0)])
    v = np.array([math.cos(theta0), -math.sin(theta0)])
    readings = [isotropic(*((radius * u + sign * 3 * v)[:, None]), 1.0, 1.0)
                for sign in (-1, 1)]
    raw = ours32(readings, c=1.0, nu_eff=10.0, calibrate=False)
    calibrated = ours32(readings, c=1.0, nu_eff=10.0)
    exact = math.tan(opposed_pair_angle_quantile(
        0.5, truth_angle=theta0, push=radius, nu_eff=10.0, error_size=3.0, calibrate=True))
    assert raw["median"] == pytest.approx(2.0, abs=0.002)
    assert calibrated["median"] == pytest.approx(exact, abs=2e-4)


def test_eta_obstruction_persists_on_an_open_error_interval():
    # For arbitrary errors e1,e2 their perpendicular half-difference is h.
    # Parallel components and the perpendicular mean do not affect this limit.
    nu = 10.0
    for h in (2.9, 3.0, 3.1):
        eta_limit = (nu - h*h) / (nu * h*h)
        tail_index = 2 * nu * eta_limit
        assert 0 < tail_index < 1
        limit = math.tan(heavy_tail_limit_quantile(
            0.5, truth_angle=math.atan(2), tail_index=tail_index))
        assert limit < 1.6


def test_opposed_pair_reference_requires_the_named_mode_regime():
    with pytest.raises(ValueError, match="central line"):
        opposed_pair_eta(push=1.0, nu_eff=10.0, error_size=3.0)


@pytest.mark.xfail(strict=True, reason="T39/O18: native uniform log grid misses the narrow integrable peak")
def test_native_ours32_grid_resolves_the_subcritical_high_signal_peak():
    theta0, radius = math.atan(2), 10000.0
    reading = isotropic([radius * math.sin(theta0)], [radius * math.cos(theta0)], 1.0, 1.0)
    exact = math.tan(aligned_angle_quantile(0.5, truth_angle=theta0, push=radius, nu_eff=0.5))
    assert ours32([reading], c=1.0, nu_eff=0.5)["median"] == pytest.approx(exact, abs=1e-3)


@pytest.mark.parametrize("kwargs", [dict(k=0), dict(nu_eff=0), dict(displacement_squared=-1),
                                    dict(noise_scale=math.inf)])
def test_population_reference_rejects_invalid_inputs(kwargs):
    arguments = dict(k=1, nu_eff=0.5, displacement_squared=1, noise_scale=1)
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        population_loss_gap(**arguments)


def test_nonconcentrating_limit_is_not_defined_at_or_above_critical_index():
    for nu in (1, 2):
        with pytest.raises(ValueError, match="tail_index < 1"):
            heavy_tail_limit_quantile(0.5, truth_angle=1.0, tail_index=nu)

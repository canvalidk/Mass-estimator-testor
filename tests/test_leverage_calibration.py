"""Independent references for score sampling and leverage-sensitive coverage.

These checks certify research identities and an unadopted spectral control;
they do not select a replacement mass equation or finite-N calibration.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.leverage_calibration import (equal_dominant_coverage,
                                      one_dominant_coverage, score_log_cf,
                                      score_moments, spectral_calibration)


@pytest.mark.parametrize("radii", [[0.], [.7], [.2, -1.3, 2.]])
def test_exact_score_cf_and_moments_against_independent_hermite_quadrature(radii):
    nodes, weights = np.polynomial.hermite.hermgauss(40)
    xi, w = np.meshgrid(math.sqrt(2)*nodes, math.sqrt(2)*nodes, indexing="ij")
    measure = np.outer(weights, weights)/math.pi
    variance, fourth, prior_variance = 0., 0., 0.
    ts = np.array([-.7, 0., .3, .6])
    characteristic = np.ones(ts.shape, dtype=complex)
    for radius in radii:
        score = (radius+xi)*w
        second = float(np.sum(measure*score**2))
        fourth += float(np.sum(measure*score**4)) + 6*prior_variance*second
        variance += second
        prior_variance += second
        assert np.sum(measure*score) == pytest.approx(0., abs=2e-14)
        characteristic *= np.array([np.sum(measure*np.exp(1j*t*score)) for t in ts])
    energy = sum(r*r for r in radii)
    moments = score_moments(len(radii), energy)
    assert moments["mean"] == 0
    assert moments["variance"] == pytest.approx(variance, rel=2e-14)
    assert moments["fourth"] == pytest.approx(fourth, rel=2e-14)
    np.testing.assert_allclose(np.exp(score_log_cf(len(radii), energy, ts)),
                               characteristic.real, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(characteristic.imag, 0., atol=2e-14)


def test_score_cf_is_even_and_stable_without_squaring_extreme_frequencies():
    assert score_log_cf(4, 3., 0.) == 0.
    assert score_log_cf(4, 3., .7) == score_log_cf(4, 3., -.7)
    actual = score_log_cf(4, 3., 1e200)
    assert actual == pytest.approx(-4*math.log(1e200)-1.5)
    # Small frequencies must not disappear through a rounded hypot(1,t).
    assert score_log_cf(1, 0., 1e-10) == pytest.approx(-.5e-20, rel=1e-14, abs=0.)


def _coverage_conditioning_on_fitted_error(tau, probability):
    """Separate h-integral; helper instead conditions on G and splits its cap."""
    z = stats.norm.ppf((1+probability)/2)
    h_sd = math.sqrt(1+tau)
    g_sd = 1/h_sd

    def integrand(h):
        threshold = math.sqrt(h*h/(z*z*tau)-1)
        mean = math.sqrt(tau)*h/(1+tau)
        tail = stats.norm.sf((threshold-mean)/g_sd) + stats.norm.cdf((-threshold-mean)/g_sd)
        return 2 * stats.norm.pdf(h/h_sd)/h_sd * tail

    value, error = integrate.quad(integrand, z, np.inf, epsabs=2e-12, epsrel=2e-11)
    return 2*stats.norm.cdf(z/h_sd)-1+value, error


@pytest.mark.parametrize("tau", [.1, .25, .5])
def test_mixed_regime_coverage_against_independent_conditional_on_h_integral(tau):
    actual, error = one_dominant_coverage(tau, .95)
    independent, independent_error = _coverage_conditioning_on_fitted_error(tau, .95)
    assert actual == pytest.approx(independent, abs=2e-11)
    assert error < 1e-9 and independent_error < 1e-9
    assert actual < .95


def test_one_dominant_coverage_pins_undercoverage_and_conservative_regimes():
    value, error = one_dominant_coverage(.25, .95)
    assert value == pytest.approx(.9272811852285586, abs=2e-12)
    assert error < 1e-9
    assert one_dominant_coverage(1., .95)[0] == pytest.approx(.9602930310231653, abs=2e-12)
    assert one_dominant_coverage(0., .95) == (.95, 0.)
    assert one_dominant_coverage(1e-8, .95)[0] == pytest.approx(.95, abs=2e-8)


@pytest.mark.parametrize("count, expected", [(2, .955011058313721),
                                              (3, .956673997626778),
                                              (5, .957807270511314)])
def test_equal_dominant_coverage_against_separate_chi_square_density_integral(count, expected):
    actual, error = equal_dominant_coverage(count, .95)
    z = stats.norm.ppf(.975)
    density_integral = integrate.quad(
        lambda x: (2*stats.norm.cdf(z*math.sqrt(x/count))-1)*stats.chi2.pdf(x, count-1),
        count, np.inf, epsabs=2e-12, epsrel=2e-11)[0]
    independent = .95*stats.chi2.cdf(count, count-1)+density_integral
    assert actual == pytest.approx(independent, abs=2e-12)
    assert actual == pytest.approx(expected, abs=2e-12)
    assert actual > .95 and error < 1e-9
    assert equal_dominant_coverage(1, .95) == (.95, 0.)


def _basis(theta):
    return np.array([math.sin(theta), math.cos(theta)]), np.array([math.cos(theta), -math.sin(theta)])


def test_spectral_control_uses_observed_spectrum_not_true_count_or_energy():
    radius, null_count = 3., 20
    u, v = _basis(math.atan(2))
    observations = np.array([radius*u]+[u+v, u-v]*(null_count//2))
    actual = spectral_calibration(observations)
    np.testing.assert_allclose(actual["eigenvalues"], [null_count, null_count+radius**2])
    assert actual["curv"] == pytest.approx(radius**2)
    assert actual["eta_sp"] == pytest.approx(radius**2/(radius**2+null_count))
    assert actual["precision_sp"] == pytest.approx(radius**4/(radius**2+null_count))
    assert actual["angle"] == pytest.approx(math.atan(2))
    assert not actual["prior_only"]


def test_spectral_control_zero_padding_sign_and_swap_invariance():
    data = np.array([[2., 1.], [-.4, .7], [1.2, -.1], [.3, .5]])
    base = spectral_calibration(data)
    variants = [np.vstack([data, np.zeros((7, 2))]),
                data*np.array([1., -1., -1., 1.])[:, None], data[:, ::-1]]
    for data_changed in variants:
        actual = spectral_calibration(data_changed)
        np.testing.assert_allclose(actual["eigenvalues"], base["eigenvalues"], atol=2e-14)
        assert actual["eta_sp"] == pytest.approx(base["eta_sp"], abs=2e-14)
        assert actual["precision_sp"] == pytest.approx(base["precision_sp"], abs=2e-14)
    # Changing physical units and their SDs leaves the standardized input fixed.
    units = np.array([7., .2])
    rescaled = spectral_calibration((data*units)/units)
    assert rescaled["eta_sp"] == pytest.approx(base["eta_sp"], abs=2e-14)


def test_spectral_control_zero_scatter_and_positive_isotropic_tie_are_prior_branches():
    zero = spectral_calibration(np.zeros((3, 2)))
    assert zero["prior_only"] and zero["eta_sp"] == 1.
    assert zero["precision_sp"] == zero["curv"] == 0.
    tied = spectral_calibration(np.eye(2))
    assert tied["prior_only"] and tied["eta_sp"] == 0.
    assert tied["precision_sp"] == tied["curv"] == 0.


@pytest.mark.parametrize("kwargs", [dict(count=0), dict(count=True), dict(energy=-1),
                                    dict(energy=math.inf), dict(t=math.nan)])
def test_score_reference_rejects_invalid_inputs(kwargs):
    arguments = dict(count=2, energy=1., t=.2)
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        score_log_cf(**arguments)


@pytest.mark.parametrize("probability", [0., 1., math.nan])
def test_coverage_requires_a_valid_explicit_probability(probability):
    with pytest.raises(ValueError):
        one_dominant_coverage(.25, probability)
    with pytest.raises(ValueError):
        equal_dominant_coverage(2, probability)


def test_spectral_and_coverage_input_validation():
    with pytest.raises(ValueError):
        spectral_calibration([1., 2.])
    with pytest.raises(ValueError):
        one_dominant_coverage(-1., .95)
    with pytest.raises(ValueError):
        equal_dominant_coverage(0, .95)

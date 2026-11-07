"""T42: independent moments, sparse geometry and continuous Gaussian laws."""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.deviance_calibration import deviance_loss, local_gaussian_reference
from met.full import isotropic, ours29
from met.gaussian_information import (flat_gaussian_fit, gaussian_angle_cdf,
                                       gaussian_angle_quantile, score_polynomial_moments)


@pytest.mark.parametrize("radius", [0., .4, -1.3, 2.])
def test_score_polynomial_moments_by_independent_hermite_quadrature(radius):
    nodes, weights = np.polynomial.hermite.hermgauss(8)
    a, b = np.meshgrid(radius+math.sqrt(2)*nodes, math.sqrt(2)*nodes, indexing="ij")
    w = np.outer(weights, weights)/math.pi
    u, v, z = (a*b)**2, (a*a-b*b)**2, a*b*(a*a-b*b)
    expected = score_polynomial_moments([radius])
    for key, values in (("U", u), ("V", v), ("W", z)):
        mean = float(np.sum(w*values))
        variance = float(np.sum(w*(values-mean)**2))
        assert mean == pytest.approx(expected["mean_"+key], abs=2e-12)
        assert variance == pytest.approx(expected["var_"+key], rel=2e-14)
    assert np.sum(w*(u-expected["mean_U"])*(v-expected["mean_V"])) == pytest.approx(expected["cov_UV"], rel=2e-14)
    assert np.sum(w*(u-expected["mean_U"])*z) == pytest.approx(0., abs=2e-11)
    assert np.sum(w*(v-expected["mean_V"])*z) == pytest.approx(0., abs=2e-11)


def _basis(theta):
    return np.array([math.sin(theta), math.cos(theta)]), np.array([math.cos(theta), -math.sin(theta)])


def test_sparse_exact_family_and_production_law_against_angle_quadrature():
    theta, radius, count = math.atan(2), 3., 20
    u, v = _basis(theta)
    observations = np.array([radius*u]+[u+v, u-v]*(count//2))
    fit = flat_gaussian_fit(observations)
    assert fit["angle"] == pytest.approx(theta)
    assert fit["curv"] == pytest.approx(radius**2)
    assert fit["J"] == pytest.approx(count)
    assert fit["eta"] == pytest.approx(.45)
    assert fit["precision"] == pytest.approx(4.05)
    readings = [isotropic([row[0]], [row[1]], 1., 1.) for row in observations]
    result = ours29(readings, 1., misfit=False, calibrate=True, readout=False)
    assert result["curv"] == pytest.approx(fit["curv"])
    assert result["J"] == pytest.approx(fit["J"])
    expected = [math.tan(gaussian_angle_quantile(p, fitted_angle=theta, precision=4.05))
                for p in (.025, .5, .975)]
    np.testing.assert_allclose([result["interval"][0], result["point"], result["interval"][1]],
                               expected, rtol=2e-4)
    # Separate direct integration checks the reference CDF, including normalization.
    kernel = lambda x: math.exp(-2.025*math.sin(x-theta)**2)
    actual, error = gaussian_angle_cdf(.8, fitted_angle=theta, precision=4.05)
    reference = integrate.quad(kernel, 0, .8)[0]/integrate.quad(kernel, 0, math.pi/2)[0]
    assert actual == pytest.approx(reference, abs=2e-12) and error < 1e-9


def test_sparse_information_growth_can_concentrate_while_eta_vanishes():
    theta = math.atan(2)
    etas, precisions, widths, bands = [], [], [], []
    for count in (100., 1e4, 1e8):
        delta = count**.75
        eta, precision = delta/count, delta*delta/count
        etas.append(eta)
        precisions.append(precision)
        low = gaussian_angle_quantile(.25, fitted_angle=theta, precision=precision)
        high = gaussian_angle_quantile(.75, fitted_angle=theta, precision=precision)
        widths.append(math.tan(high)-math.tan(low))
        bands.append(gaussian_angle_cdf(math.atan(2.1), fitted_angle=theta, precision=precision)[0]
                     - gaussian_angle_cdf(math.atan(1.9), fitted_angle=theta, precision=precision)[0])
    assert etas[0] > etas[1] > etas[2]
    assert precisions[0] < precisions[1] < precisions[2]
    assert widths[0] > widths[1] > widths[2]
    assert bands[0] < bands[1] < bands[2] and bands[-1] > .9


def test_opposed_pair_has_the_gaussian_precision_not_the_cut_precision():
    theta, radius, error = math.atan(2), 10., 3.
    u, v = _basis(theta)
    observations = np.array([radius*u+error*v, radius*u-error*v])
    fit = flat_gaussian_fit(observations)
    assert fit["angle"] == pytest.approx(theta)
    assert fit["curv"] == pytest.approx(2*(radius**2-error**2))
    assert fit["J"] == pytest.approx(2*radius**2*error**2)
    assert fit["eta"] == pytest.approx((radius**2-error**2)/(radius**2*error**2))


def test_single_row_self_fit_and_zero_data():
    fit = flat_gaussian_fit([[2., 1.]])
    assert fit["eta"] == 1 and fit["precision"] == pytest.approx(5.)
    assert flat_gaussian_fit(np.zeros((3, 2)))["precision"] == 0
    assert gaussian_angle_cdf(.3, fitted_angle=0., precision=0.) == (.3/(math.pi/2), 0.)


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_growing_scale_and_vanishing_eta_recover_the_local_reference(family):
    errors = []
    for n in (100., 1e4, 1e6):
        result = local_gaussian_reference(n**(-.25), family=family, scale=n, probability=.95)
        errors.append(abs(result["half_width"]-stats.norm.ppf(.975)))
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < .001


@pytest.mark.parametrize("family", ["logistic", "lambert"])
def test_diverging_scale_alone_does_not_give_the_gaussian_local_limit(family):
    # Scale diverges, but scale*eta=1 retains a nonlinear standard-score loss.
    eta, z = 1e-5, np.array([0., .5, 1., 1.5, 2.])
    x = z*z/2
    limit = np.expm1(x) if family == "logistic" else x*np.exp(x)
    actual = deviance_loss(x/eta, eta, family=family, scale=1/eta)
    # log(1+t) and W(t) give this explicit O(eta) pointwise error bound.
    assert np.all(actual <= limit+1e-12)
    assert np.all(limit-actual <= eta*limit**2+1e-12)
    assert actual[-1]-x[-1] > 1
    reference = local_gaussian_reference(eta, family=family, scale=1/eta, probability=.95)
    assert reference["half_width"] < stats.norm.ppf(.975)-.1


def test_invalid_inputs():
    with pytest.raises(ValueError):
        flat_gaussian_fit([1., 2.])
    with pytest.raises(ValueError):
        gaussian_angle_cdf(.5, fitted_angle=.5, precision=-1)
    with pytest.raises(ValueError):
        gaussian_angle_quantile(1., fitted_angle=.5, precision=1.)

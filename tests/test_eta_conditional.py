"""Independent checks for B15's fixed-candidate working score variance.

The Student calculation uses T48's explicitly integrated flat nuisance measure.
It is not the conditional law of a scale mixture with a fixed unknown signal.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.conditional_score import gaussian_score_test
from met.eta_conditional import eta_c, residual_score_multiplier, score_variance
from met.full import Reading, isotropic
from met.ours33 import ours33


def _scalar(data):
    return [Reading([f], [a], np.eye(2)) for f, a in data]


def _full_example():
    matrix = np.array([[1.4, .2, -.1, .3], [.1, 1.1, .4, 0.],
                       [.3, -.2, 1.7, .2], [0., .1, -.1, 1.2]])
    return Reading([2., -.7], [.5, 1.1], matrix @ matrix.T,
                   np.array([[.6], [.8]]))


def _pair(radius=1000., perpendicular=3.):
    theta = math.atan(2.)
    u = np.array([math.sin(theta), math.cos(theta)])
    v = np.array([math.cos(theta), -math.sin(theta)])
    return _scalar([radius*u + perpendicular*v, radius*u - perpendicular*v])


@pytest.mark.parametrize("k,nu", [(1, .5), (3, .5), (1, 10.), (3, 10.)])
def test_student_score_variance_from_independent_radial_density(k, nu):
    # For unit tangent, symmetry gives E[(r.t)^2 | ||r||^2=q] = q/k.
    # Integrate the F density of q/k, without using the beta-moment formula.
    def integrand(q):
        score_square = (nu+k)**2*q/(k*(nu+q)**2)
        return score_square*stats.f.pdf(q/k, k, nu)/k

    actual, error = integrate.quad(integrand, 0., np.inf,
                                   epsabs=2e-10, epsrel=2e-10)
    assert actual == pytest.approx(residual_score_multiplier(k, nu), abs=2e-9)
    assert error < 2e-8


@pytest.mark.parametrize("c,mass", [(.4, .3), (1., 2.), (3., 4.7)])
@pytest.mark.parametrize("factors", ["flat", "ours33"])
def test_full_covariance_variance_has_t45_geometry_and_correct_angle(c, mass, factors):
    reading = _full_example()
    phi = math.atan(mass/c)
    # T45 differentiates in theta=atan(m), while this module uses
    # phi=atan(m/c). Their variances transform by (dtheta/dphi)^2.
    derivative = (c*c+mass*mass)/(c*(1+mass*mass))
    expected = gaussian_score_test([reading], mass)["variance"]*derivative**2
    if factors == "ours33":
        expected *= (10.+reading.k)/(10.+reading.k+2.)
    assert score_variance([reading], c, phi, 10., factors=factors) == pytest.approx(expected, rel=2e-12)


@pytest.mark.parametrize("factors", ["flat", "ours33"])
def test_full_covariance_units_frame_and_reciprocal_variance(factors):
    reading, c, phi = _full_example(), 1.7, .6
    base = score_variance([reading], c, phi, 10., factors=factors)
    theta = .7
    rotation = np.array([[math.cos(theta), -math.sin(theta)],
                         [math.sin(theta), math.cos(theta)]])
    transform = np.zeros((4, 4))
    transform[:2, :2] = 7.*rotation
    transform[2:, 2:] = .2*rotation
    data = transform @ np.r_[reading.force, reading.acceleration]
    changed = Reading(data[:2], data[2:], transform @ reading.cov @ transform.T,
                      -rotation @ reading.subspace)
    assert score_variance([changed], c*7./.2, phi, 10., factors=factors) == pytest.approx(base, rel=2e-12)
    order = [2, 3, 0, 1]
    swapped = Reading(reading.acceleration, reading.force,
                      reading.cov[np.ix_(order, order)], reading.subspace)
    assert score_variance([swapped], 1./c, math.pi/2-phi, 10., factors=factors) == pytest.approx(base, rel=2e-12)


@pytest.mark.parametrize("data", [
    [[4., .2], [-.3, 1.4], [1., -.8]],
    [[.2, 3.], [.7, -.2], [-1., .3]],
    [[2., 1.]],
])
def test_flat_fit_matches_independent_scatter_spectrum(data):
    data = np.asarray(data)
    eigenvalues = np.linalg.eigvalsh(data.T @ data)
    gap, largest = eigenvalues[1]-eigenvalues[0], eigenvalues[1]
    eta, curvature, variance, phi = eta_c(_scalar(data), 1., 10., factors="flat")
    assert curvature == pytest.approx(gap, rel=2e-11)
    assert variance == pytest.approx(largest, rel=2e-11)
    assert eta == pytest.approx(gap/largest, rel=2e-11)
    u = np.array([math.sin(phi), math.cos(phi)])
    assert float(u @ data.T @ data @ u) == pytest.approx(largest, rel=2e-11)


def test_t39_pair_retains_small_conditional_eta():
    radius, d, nu = 1000., 3., 10.
    eta, curvature, variance, phi = eta_c(_pair(radius, d), 1., nu, factors="ours33")
    a, denominator = (nu+1.)/2., nu+d*d
    expected_h = 4*a*(radius*radius*(nu-d*d)-d*d*denominator)/denominator**2
    expected_c = 2.*radius*radius*(nu+1.)/(nu+3.)
    assert phi == pytest.approx(math.atan(2.), abs=2e-12)
    assert curvature == pytest.approx(expected_h, rel=2e-10)
    assert variance == pytest.approx(expected_c, rel=2e-12)
    assert eta == pytest.approx(expected_h/expected_c, rel=2e-10)
    assert .035 < eta < .037


@pytest.mark.parametrize("factors", ["flat", "ours33"])
def test_exact_zero_padding_is_bitwise_neutral_and_empty_has_convention(factors):
    readings = [_full_example(), isotropic([1.2], [.8], 1., 1.)]
    zero = isotropic([0., 0., 0.], [0., 0., 0.], 2., .3)
    base = eta_c(readings, 1.7, [10., .5], factors=factors)
    padded = eta_c([zero, *readings, zero], 1.7, [3., 10., .5, 100.], factors=factors)
    assert padded == base
    for values in ([], [zero, zero]):
        assert eta_c(values, 1.7, 10., factors=factors) == (1., 0., 0., 0.)
    assert score_variance([zero, *readings, zero], 1.7, .6,
                          [3., 10., .5, 100.], factors=factors) == score_variance(readings, 1.7, .6, [10., .5], factors=factors)
    if factors == "flat":
        assert eta_c(readings, 1.7, [.1, 100.], factors=factors) == base


def test_optional_ours33_switch_changes_only_calibration_and_default_stays_fitted():
    readings = _pair()
    base = ours33(readings, 1., 10.)
    fitted = ours33(readings, 1., 10., eta_method="fitted")
    conditional = ours33(readings, 1., 10., eta_method="conditional")
    assert base["eta"] == fitted["eta"]
    assert base["interval"] == fitted["interval"]
    assert conditional["eta"] == eta_c(readings, 1., 10., factors="ours33")[0]
    assert conditional["eta"] > base["eta"]
    assert conditional["point"] == base["point"]
    assert conditional["raw_interval"] == base["raw_interval"]
    np.testing.assert_array_equal(conditional["law"], base["law"])
    assert ours33(readings, 1., 10., calibrate=False, eta_method="conditional")["eta"] == 1.
    with pytest.raises(ValueError):
        ours33(readings, 1., 10., eta_method="unknown")


@pytest.mark.parametrize("kwargs", [
    {"c": 0.}, {"c": math.inf}, {"nu_eff": 0.},
    {"nu_eff": math.nan}, {"nu_eff": [10., 10.]}, {"factors": "unknown"},
])
def test_invalid_conditional_variance_inputs_raise(kwargs):
    arguments = dict(c=1., nu_eff=10., factors="ours33")
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        eta_c([_full_example()], **arguments)
    with pytest.raises(ValueError):
        score_variance([_full_example()], phi=.6, **arguments)


@pytest.mark.parametrize("k,nu", [(1, 10.), (3, .5)])
def test_fixed_unknown_signal_conditioning_is_a_different_residual_law(k, nu):
    # At an observed U equal to its true mean, V|U has shape (nu+k)/2.
    # The conditional radius q/k is nu/(nu+k) times F(k,nu+k), rather
    # than F(k,nu). Independently integrate the same nonlinear score.
    scale = nu/(nu+k)

    def integrand(q):
        density = stats.f.pdf(q/(k*scale), k, nu+k)/(k*scale)
        return (nu+k)**2*q/(k*(nu+q)**2)*density

    actual, _ = integrate.quad(integrand, 0., np.inf, epsabs=2e-10)
    expected = (nu+k)**3/(nu*(nu+2*k)*(nu+2*k+2))
    assert actual == pytest.approx(expected, abs=2e-9)
    assert abs(actual-residual_score_multiplier(k, nu)) > .01

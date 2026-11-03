"""T48: normalized Gaussian mixture, amplitude limits and surviving power defect."""
import math

import numpy as np
import pytest
from scipy import integrate, stats

from met.cut_consistency import (heavy_tail_limit_cdf, opposed_pair_angle_cdf,
                                 opposed_pair_eta)
from met.residual_mixture import (aligned_angle_cdf, aligned_angle_quantile,
                                  pair_angle_cdf, student_pair_parameters,
                                  student_residual_density, student_residual_log_density)


CONTROLS = dict(epsabs=2e-12, epsrel=2e-11, limit=200)


@pytest.mark.parametrize("dimension", [1, 2, 3])
@pytest.mark.parametrize("nu", [.2, .5, 10.])
def test_density_matches_independent_variance_integral_and_multivariate_t(dimension, nu):
    residual = np.arange(1, dimension+1, dtype=float)/(dimension+1)
    cost = float(residual@residual)

    def integrand(v):
        log_normal = -.5*dimension*math.log(2*math.pi*v)-cost/(2*v)
        return math.exp(log_normal+stats.invgamma.logpdf(v, nu/2, scale=nu/2))

    numerical, error = integrate.quad(integrand, 0, np.inf, epsabs=2e-11, epsrel=2e-10, limit=200)
    actual = student_residual_density(cost, dimension=dimension, nu_eff=nu)
    assert actual == pytest.approx(numerical, rel=3e-9, abs=2e-11)
    assert error < 1e-9
    expected = stats.multivariate_t.pdf(residual, loc=np.zeros(dimension), shape=np.eye(dimension), df=nu)
    assert actual == pytest.approx(expected, rel=2e-14)


@pytest.mark.parametrize("nu", [.2, .5, 2.])
def test_corrected_single_reading_band_mass_grows_for_every_tested_nu(nu):
    theta, probabilities, medians = math.atan(2), [], []
    for radius in (10., 100., 1000.):
        arguments = dict(observed_angle=theta, radius=radius, nu_eff=nu, **CONTROLS)
        low, le = aligned_angle_cdf(math.atan(1.9), **arguments)
        high, he = aligned_angle_cdf(math.atan(2.1), **arguments)
        probabilities.append(high-low)
        medians.append(math.tan(aligned_angle_quantile(.5, **arguments)))
        assert max(le, he) < 1e-7
    assert probabilities[0] < probabilities[1] < probabilities[2]
    assert abs(medians[-1]-2) < abs(medians[0]-2)
    if nu == .5:
        old_limit = heavy_tail_limit_cdf(math.atan(2.1), truth_angle=theta, tail_index=nu) - \
            heavy_tail_limit_cdf(math.atan(1.9), truth_angle=theta, tail_index=nu)
        assert probabilities[-1] > .8 and old_limit < .2


def test_global_bartlett_power_exactly_cancels_the_prefactor_change_on_the_pair():
    theta, nu, d = math.atan(2), 10., 3.
    for radius in (100., 1000.):
        parameters = student_pair_parameters(radius=radius, nu_eff=nu, error_size=d)
        old_eta = opposed_pair_eta(push=radius, nu_eff=nu, error_size=d)
        assert parameters["eta"] < 1
        assert parameters["eta"]*parameters["coefficient"] == pytest.approx(old_eta*nu/2, rel=2e-14)
        for angle in (math.atan(1.9), theta, math.atan(2.1)):
            new, error = pair_angle_cdf(angle, truth_angle=theta, radius=radius,
                                         nu_eff=nu, error_size=d, calibrate=True, **CONTROLS)
            old = opposed_pair_angle_cdf(angle, truth_angle=theta, push=radius,
                                          nu_eff=nu, error_size=d, calibrate=True)[0]
            assert new == pytest.approx(old, abs=2e-10) and error < 1e-7
    parameters = student_pair_parameters(radius=1000., nu_eff=nu, error_size=d)
    assert parameters["eta"] == pytest.approx(1/99-19/(11*1000.**2), rel=2e-13)
    assert parameters["quoted_tail_index"] < 1 and parameters["raw_tail_index"] == 22


def test_same_corrected_pair_concentrates_raw_but_not_under_the_global_power():
    arguments = dict(truth_angle=math.atan(2), radius=1000., nu_eff=10., error_size=3., **CONTROLS)
    probabilities = {}
    for calibrate in (False, True):
        low = pair_angle_cdf(math.atan(1.9), calibrate=calibrate, **arguments)[0]
        high = pair_angle_cdf(math.atan(2.1), calibrate=calibrate, **arguments)[0]
        probabilities[calibrate] = high-low
    assert probabilities[False] > .999 and probabilities[True] < .06


def test_fixed_antiparallel_reading_keeps_bounded_nonvanishing_relative_evidence():
    angles = np.unique(np.r_[np.linspace(0, math.pi/2, 41), math.atan(2)])
    observed = np.array([1., -2.])
    # Independent compatible-vector least squares supplies the physical cost.
    costs = []
    for angle in angles:
        design = np.array([[math.sin(angle)], [math.cos(angle)]])
        fitted = design @ np.linalg.lstsq(design, observed, rcond=None)[0]
        costs.append(float(np.sum((observed-fitted)**2)))
    costs = np.asarray(costs)
    np.testing.assert_allclose(costs, (np.cos(angles)+2*np.sin(angles))**2, atol=2e-14)
    assert costs.min() == pytest.approx(1.) and costs.max() == pytest.approx(5.)
    for nu in (.5, 10.):
        coefficient, errors, oscillations = (nu+1)/2, [], []
        limit = -coefficient*np.log(costs)
        for sigma in (1., .1, .01, 1e-4):
            log_factor = student_residual_log_density(costs/sigma**2, dimension=1, nu_eff=nu)
            best = student_residual_log_density(1/sigma**2, dimension=1, nu_eff=nu)
            relative = log_factor-best
            errors.append(float(np.max(np.abs(relative-limit))))
            oscillations.append(float(np.ptp(relative)))
        assert errors[0] > errors[1] > errors[2] > errors[3]
        assert errors[-1] < 1e-6
        assert oscillations[0] < oscillations[1] < oscillations[2] < oscillations[3]
        assert oscillations[-1] == pytest.approx(coefficient*math.log(5), abs=1e-6)


def test_model_and_numerical_choices_are_explicit():
    with pytest.raises(ValueError):
        student_residual_density(-1., dimension=1, nu_eff=.5)
    with pytest.raises(ValueError):
        student_residual_density(1., dimension=0, nu_eff=.5)
    with pytest.raises(ValueError):
        student_pair_parameters(radius=1., nu_eff=10., error_size=3.)
    with pytest.raises(TypeError):
        aligned_angle_cdf(.5, observed_angle=.7, radius=10., nu_eff=.5)

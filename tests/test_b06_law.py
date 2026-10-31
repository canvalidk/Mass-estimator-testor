"""B06: current-law qualifications for a tilt and the exact no-mean tail."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.full import Reading, _h_terms, _log_law, _log_prior, isotropic, ours30


@pytest.mark.parametrize("tilt_centre,lam,prior_centre", [
    (1.7, 0.8, 0.4), (0.3, -0.1, 2.0), (0.3, -1.5, 2.0),
])
def test_equal_energy_reading_has_the_actual_ours30_raw_tilt(tilt_centre, lam, prior_centre):
    base = [isotropic([1.2, -0.4], [0.3, 0.8], 0.7, 0.2)]
    x = math.sqrt(2 * abs(lam))
    extra = isotropic([tilt_centre * x], [math.copysign(x, lam)], tilt_centre, 1.0)
    before = ours30(base, prior_centre, calibrate=False, readout=False)
    after = ours30(base + [extra], prior_centre, calibrate=False, readout=False)
    np.testing.assert_allclose(before["g"], after["g"][:-1], atol=1e-13)
    extra_g = after["g"][-1]
    expected_g = 1.0 if lam > 0 or abs(lam) <= 0.25 else 1 / (4 * abs(lam))
    assert extra_g == pytest.approx(expected_g, abs=1e-13)
    lm = np.linspace(-6, 6, 301)
    difference = (_log_law(base + [extra], prior_centre, after["g"], 1.0, lm)
                  - _log_law(base, prior_centre, before["g"], 1.0, lm))
    tilt = np.sin(2 * np.arctan(np.exp(lm) / tilt_centre))
    assert np.ptp(difference - extra_g * lam * tilt) < 2e-12
    if lam < -0.25:
        assert extra_g * lam == pytest.approx(-0.25)
        assert np.ptp(difference - lam * tilt) > 1


def test_recomputing_eta_changes_the_unseen_reading_identity():
    base = [isotropic([3], [1], 1, 1), isotropic([1], [3], 1, 1)]
    lam = 0.8
    extra = isotropic([math.sqrt(2 * lam)], [math.sqrt(2 * lam)], 1, 1)
    before = ours30(base, 1.0, readout=False)
    after = ours30(base + [extra], 1.0, readout=False)
    assert before["eta"] == pytest.approx(0.375)
    assert after["eta"] == pytest.approx(0.475)
    lm = np.linspace(-5, 5, 241)
    tilt = np.sin(2 * np.arctan(np.exp(lm)))
    difference = (_log_law(base + [extra], 1, after["g"], after["eta"], lm)
                  - _log_law(base, 1, before["g"], before["eta"], lm))
    assert np.ptp(difference - 0.68 * tilt) < 2e-12
    assert np.ptp(difference - lam * tilt) > 0.1


@pytest.mark.parametrize("k", [1, 2, 3])
def test_full_covariance_ours30_tail_has_positive_inverse_square_coefficient(k):
    rng = np.random.default_rng(1906 + k)
    a = rng.normal(size=(6, 6))
    covariance = a @ a.T + np.eye(6)
    subspace = np.linalg.qr(rng.normal(size=(3, 3)))[0][:, :k]
    readings = [Reading(rng.normal(size=3), rng.normal(size=3), covariance, subspace),
                isotropic([0.7], [-1.1], 0.5, 0.8)]
    centre = 1.6
    result = ours30(readings, centre, readout=False)
    g, eta = result["g"], result["eta"]
    energies = [np.r_[r.force, r.acceleration] @ np.linalg.solve(
        r.cov, np.r_[r.force, r.acceleration]) for r in readings]
    upper = 0.5 * sum(energies)
    grid = np.linspace(-7, 7, 201)
    log_factor = _log_law(readings, centre, g, eta, grid) - _log_prior(grid, centre)
    assert log_factor.min() >= -2e-14
    assert log_factor.max() <= upper + 2e-14

    endpoint = sum(gi * r.pre()["u"] @ np.linalg.solve(r.pre()["k_ff"], r.pre()["u"])
                   for r, gi in zip(readings, g)) * eta / 2
    def angular(phi):
        return math.exp(0.5 * eta * sum(gi * _h_terms(r.pre(), centre, np.array([phi]))[0]
                                       for r, gi in zip(readings, g)))
    normalizer = quad(angular, 0, math.pi / 2, epsabs=1e-11, epsrel=1e-11)[0]
    assert math.pi / 2 <= normalizer <= math.exp(upper) * math.pi / 2
    coefficient = centre * math.exp(endpoint) / normalizer
    assert 0 < coefficient < np.inf
    x = np.linspace(20, 28, 81)
    # e^x p_x(x) is the first-moment integrand with respect to log mass.
    integrand = np.exp(_log_law(readings, centre, g, eta, x) + x) / normalizer
    np.testing.assert_allclose(integrand, coefficient, rtol=2e-8)
    left = np.trapezoid(integrand[:41], x[:41])
    right = np.trapezoid(integrand[40:], x[40:])
    assert left == pytest.approx(4 * coefficient, rel=2e-8)
    assert right == pytest.approx(left, rel=2e-8)

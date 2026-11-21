"""Known-standard Gaussian experiment and loss from retaining only a force sum."""
import math

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

from explore.e12_residual_join import residual_log_evidence
from explore.e12_known_mass import known_mass_log_evidence, standard_sensor_record
from met.conditional_score import gaussian_score_test
from met.full import Reading


def effective_record(data, covariance, standard, *, lump=False):
    d = len(data) // 3
    eye, zero = np.eye(d), np.zeros((d, d))
    selection = np.block([[eye, zero, eye if lump else zero], [zero, eye, zero]])
    defect = np.column_stack((zero, -standard * eye, eye))
    residual_noise = defect @ covariance @ defect.T
    beta = np.linalg.solve(residual_noise, defect @ covariance @ selection.T).T
    correction = selection - beta @ defect
    effective = correction @ data
    return effective, correction @ covariance @ correction.T, correction, defect


def profile_cost(data, covariance, design):
    chol = np.linalg.cholesky(covariance)
    y = np.linalg.solve(chol, data)
    a = np.linalg.solve(chol, design)
    residual = y - a @ np.linalg.lstsq(a, y, rcond=None)[0]
    return residual @ residual


def test_gaussian_reduction_and_shear_with_full_correlated_spatial_sensors():
    n, b = 2, 0.7
    rng = np.random.default_rng(19362)
    root = rng.normal(size=(3*n, 3*n))
    covariance = root @ root.T + np.eye(3*n)
    data = rng.normal(size=3*n)
    part, s, correction, defect = effective_record(data, covariance, b)
    total, t, _, _ = effective_record(data, covariance, b, lump=True)
    for lump, expected, noise in [(False, part, s), (True, total, t)]:
        reading, ancillary = standard_sensor_record(data, covariance, b, lump=lump)
        np.testing.assert_allclose(np.r_[reading.force, reading.acceleration], expected)
        np.testing.assert_allclose(reading.cov, noise)
        np.testing.assert_allclose(ancillary, defect @ data)
    eye, zero = np.eye(n), np.zeros((n, n))
    shear = np.block([[eye, b*eye], [zero, eye]])
    np.testing.assert_allclose(correction @ covariance @ defect.T, 0, atol=1e-14)
    np.testing.assert_allclose(total, shear @ part)
    np.testing.assert_allclose(t, shear @ s @ shear.T)

    nuisance = defect @ data
    ancillary_cost = nuisance @ np.linalg.solve(defect @ covariance @ defect.T, nuisance)
    for mass in (0.2, 1.0, 5.0):
        model = np.vstack((mass*eye, eye, b*eye))
        full_cost = profile_cost(data, covariance, model)
        part_cost = profile_cost(part, s, np.vstack((mass*eye, eye)))
        total_cost = profile_cost(total, t, np.vstack(((mass+b)*eye, eye)))
        assert full_cost == pytest.approx(part_cost + ancillary_cost)
        assert total_cost == pytest.approx(part_cost)
        fine = gaussian_score_test([Reading(part[:n], part[n:], s)], mass)
        summed = gaussian_score_test([Reading(total[:n], total[n:], t)], mass+b)
        assert summed["statistic"] == pytest.approx(fine["statistic"])


def test_force_sum_loses_observed_mass_information():
    # The sensor channels are the unknown force, common acceleration and standard force.
    pairs = []
    for y in (np.array([1., 1., 1.]), np.array([2., 1., 0.])):
        raw_total = [y[0] + y[2], y[1]]
        np.testing.assert_allclose(raw_total, [2, 1])
        part, s, _, _ = effective_record(y, np.eye(3), 1)
        np.testing.assert_allclose(s, np.diag([1, .5]))
        pairs.append(part)
    assert [f/a for f, a in pairs] == pytest.approx([1, 4])
    assert np.array([1, -1]) @ s @ [1, -1] == pytest.approx(1.5)
    raw_total_cov = np.diag([2., 1.])
    assert np.array([1, -2]) @ raw_total_cov @ [1, -2] == pytest.approx(6.)


@pytest.mark.parametrize("tau,expected", [(0., 2.528611), (64., 2.510203)])
def test_distinct_readouts_extend_to_measured_standard(tau, expected):
    covariance = np.diag([1., 2., 2.])
    effective = []
    combined = []
    for y in (np.array([3., 0., 2.]), np.array([1., 1.5, 2.5])):
        part, s, _, _ = effective_record(y, covariance, 1)
        total, t, _, _ = effective_record(y, covariance, 1, lump=True)
        np.testing.assert_allclose(s, np.eye(2))
        np.testing.assert_allclose(t, [[2, 1], [1, 1]])
        np.testing.assert_allclose(part, [total[0] - total[1], total[1]])
        effective.append(part)
        reading, _ = standard_sensor_record(y, covariance, 1, lump=True)
        combined.append(np.column_stack((reading.force, reading.acceleration)))
    np.testing.assert_allclose(effective, [[3, 1], [1, 2]])
    grams = [np.outer(row, row) for row in effective]
    angle = np.linspace(0, math.pi/2, 14001)
    log_density = residual_log_evidence(grams, angle, tau=tau)
    direct = known_mass_log_evidence(combined, reading.cov, 1 + np.tan(angle),
                                     bounds=(1, np.inf), tau=tau)
    np.testing.assert_allclose(direct, log_density, atol=1e-12)
    cdf = cumulative_trapezoid(np.exp(log_density - np.max(log_density)), angle, initial=0)
    median = 1 + math.tan(np.interp(0.5*cdf[-1], cdf, angle))
    assert median == pytest.approx(expected, abs=2e-6)


def test_irrelevant_calibration_error_changes_unprojected_nonlinear_evidence():
    effective = np.array([2., 1.])
    covariance = np.diag([1., 2., 2.])
    by_ancillary = []
    for r in (0., 4.):
        data = np.array([effective[0], effective[1] - r/2, effective[1] + r/2])
        costs = [profile_cost(data, covariance, np.array([[m], [1.], [1.]])) for m in (2., 1.)]
        by_ancillary.append(np.subtract(*(-np.expm1(-0.5*np.array(costs)))))
    assert by_ancillary[1] == pytest.approx(math.exp(-2) * by_ancillary[0])

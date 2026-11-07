"""Independent residual and stationary-line checks for ours32's numerics.

These verify the existing Student factors; they do not establish consistency
of that factor or resolve the native probability grid's narrow-spike limits.
"""
import math

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from met.full import Reading, isotropic
from met.ours31 import _per_reading
from met.ours32 import _l_terms, bartlett_eta, ours32


def _states(readings, c, nu=10.0):
    return _per_reading(readings, c, [nu] * len(readings))


def _fit_residual(reading, c, phi):
    """Independent full-space least squares minus the outside-span baseline."""
    n = reading.subspace
    zero = np.zeros_like(n)
    span = np.block([[n, zero], [zero, n]])
    chol = np.linalg.cholesky(reading.cov)
    data = np.linalg.solve(chol, np.r_[reading.force, reading.acceleration])
    white_span = np.linalg.solve(chol, span)
    outside = data - white_span @ np.linalg.lstsq(white_span, data, rcond=None)[0]
    candidate = np.vstack([c * math.sin(phi) * n, math.cos(phi) * n])
    candidate = np.linalg.solve(chol, candidate)
    residual = data - candidate @ np.linalg.lstsq(candidate, data, rcond=None)[0]
    return residual @ residual - outside @ outside


def _general_series(seed=3203, scale=1.0, ks=(1, 2, 3, 2)):
    rng = np.random.default_rng(seed)
    readings = []
    for k in ks:
        a = rng.normal(size=(6, 6))
        cov = a @ a.T + 0.4 * np.eye(6)
        n, _ = np.linalg.qr(rng.normal(size=(3, k)))
        readings.append(Reading(scale * rng.normal(size=3), scale * rng.normal(size=3), cov, n))
    return readings


def test_large_antiparallel_reading_has_finite_whole_line_factors():
    """The old 1 + (H+ - H)/(nu + Q+) rounds nonpositive at its negative line."""
    r, nu = 1e8, 0.5
    reading = isotropic([2 * r], [-r], 1.0, 1.0)
    phi = np.array([-math.atan(2.0), -0.9, 0.3, 1.2])
    with np.errstate(all="raise"):
        terms = _l_terms(_states([reading], 1.0, nu), 1.0, phi, 2)[0]
    assert all(np.all(np.isfinite(x)) for x in terms)
    residual = (2 * r * np.cos(phi) + r * np.sin(phi)) ** 2
    expected = -0.5 * nu * np.log1p(residual / nu)
    np.testing.assert_allclose(terms[0], expected, atol=1e-14, rtol=1e-14)
    result = ours32([reading], 1.0, nu)
    assert 0 < result["interval"][0] < result["median"] < result["interval"][1]


def test_general_covariance_restricted_span_matches_independent_fits():
    phi, c, nu = np.array([-1.1, -0.2, 0.3, 1.2]), 1.7, 4.0
    readings = _general_series()
    for reading, terms in zip(readings, _l_terms(_states(readings, c, nu), c, phi, 2)):
        residual = np.array([_fit_residual(reading, c, x) for x in phi])
        expected = -0.5 * (nu + reading.k - 1) * np.log1p(residual / nu)
        np.testing.assert_allclose(terms[0], expected, atol=1e-12, rtol=1e-12)


def test_general_residual_derivatives_match_finite_differences():
    readings, c = _general_series(), 1.7
    states = _states(readings, c, 4.0)
    phi, delta = np.array([-1.1, 0.3, 1.2]), 2e-5
    center = _l_terms(states, c, phi, 2)
    high = _l_terms(states, c, phi + delta, 0)
    low = _l_terms(states, c, phi - delta, 0)
    for middle, up, down in zip(center, high, low):
        np.testing.assert_allclose(middle[1], (up[0] - down[0]) / (2 * delta), atol=2e-8)
        np.testing.assert_allclose(middle[2], (up[0] - 2 * middle[0] + down[0]) / delta**2,
                                   atol=2e-5)


def test_eta_resolves_a_narrow_peak_independently_of_angle_chart():
    rng = np.random.default_rng(90210)
    noise = rng.normal(size=(20, 2, 3))
    theta, strength, nu = math.atan(2.0), 1e4, 10.0
    force, acceleration = noise[:, 0].copy(), noise[:, 1].copy()
    force[:, 0] += strength * math.sin(theta)
    acceleration[:, 0] += strength * math.cos(theta)
    readings = [isotropic(f, a, 1.0, 1.0) for f, a in zip(force, acceleration)]

    def objective(offset):
        angle = theta + offset / strength
        residual = force * math.cos(angle) - acceleration * math.sin(angle)
        return np.log1p(np.sum(residual**2, axis=1) / nu).sum()

    fit = minimize_scalar(objective, bounds=(-10, 10), method="bounded",
                          options={"xatol": 1e-10})
    expected_mass = math.tan(theta + fit.x / strength)
    for c in (0.1, 1.0, 10.0):
        states = _states(readings, c, nu)
        eta, curvature, j, fitted = bartlett_eta(states, c)
        assert c * math.tan(fitted) == pytest.approx(expected_mass, rel=2e-9)
        assert curvature > 0 and eta == pytest.approx(1.0)
        score = sum(x[1][0] for x in _l_terms(states, c, np.array([fitted]), 2))
        assert abs(score) / math.sqrt(j) < 2e-6


def test_eta_matches_exact_symmetric_two_reading_curvature():
    """Check eta's computation, including a small value; not its consistency."""
    theta, strength, delta, nu = math.atan(2.0), 1000.0, 3.0, 10.0
    along = np.array([math.sin(theta), math.cos(theta)])
    normal = np.array([math.cos(theta), -math.sin(theta)])
    data = [strength * along + sign * delta * normal for sign in (-1, 1)]
    readings = [isotropic([y[0]], [y[1]], 1.0, 1.0) for y in data]
    eta, curvature, j, fitted = bartlett_eta(_states(readings, 1.0, nu), 1.0)
    expected = ((nu - delta**2) / (nu * delta**2)
                - (nu + delta**2) / (nu * strength**2))
    assert math.tan(fitted) == pytest.approx(2.0, rel=1e-11)
    assert curvature > 0 and j > 0
    assert eta == pytest.approx(expected, rel=1e-10)


def test_stationary_line_at_the_projective_seam():
    strength = 1e4
    reading = isotropic([strength], [0.0], 1.0, 1.0)
    eta, curvature, _, fitted = bartlett_eta(_states([reading], 1.0), 1.0)
    assert abs(math.cos(fitted)) < 1e-12
    assert curvature == pytest.approx(strength**2, rel=1e-12)
    assert eta == pytest.approx(1.0)


def test_general_covariance_transforms_and_exact_zero_padding():
    readings, c = _general_series(seed=0, scale=4.0, ks=(1, 2, 3, 2, 1, 3)), 1.7
    baseline = bartlett_eta(_states(readings, c), c)[0]
    assert 0 < baseline < 0.9  # Exercise the ratio itself, rather than the eta=1 cap.
    rng = np.random.default_rng(3204)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    frame = np.kron(np.eye(2), rotation)
    changed = [Reading(rotation @ r.force, rotation @ r.acceleration,
                       frame @ r.cov @ frame.T, rotation @ r.subspace) for r in readings]
    assert bartlett_eta(_states(changed, c), c)[0] == pytest.approx(baseline, rel=1e-7)
    scaling = np.diag([3.0] * 3 + [0.4] * 3)
    changed = [Reading(3 * r.force, 0.4 * r.acceleration, scaling @ r.cov @ scaling,
                       r.subspace) for r in readings]
    assert bartlett_eta(_states(changed, c * 7.5), c * 7.5)[0] == pytest.approx(baseline, rel=1e-7)
    swap = np.r_[3:6, 0:3]
    changed = [Reading(r.acceleration, r.force, r.cov[np.ix_(swap, swap)], r.subspace)
               for r in readings]
    assert bartlett_eta(_states(changed, 1 / c), 1 / c)[0] == pytest.approx(baseline, rel=1e-7)
    r = readings[0]
    zero = Reading(np.zeros(3), np.zeros(3), r.cov, r.subspace)
    assert bartlett_eta(_states(readings + [zero], c), c)[0] == pytest.approx(baseline, rel=1e-7)

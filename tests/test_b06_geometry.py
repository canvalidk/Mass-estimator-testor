"""B06: zero-reading geometry and the conditions on zero-axis/1D reduction."""
import math

import numpy as np
import pytest

from met.full import Reading, _h_terms, isotropic, ours30, reading_state


def _h(reading, mass):
    return float(_h_terms(reading.pre(), 1.0, np.array([math.atan(mass)]))[0])


def _block_covariance(retained, zero, retained_indices, zero_indices):
    covariance = np.zeros((len(retained_indices) + len(zero_indices),) * 2)
    covariance[np.ix_(retained_indices, retained_indices)] = retained
    covariance[np.ix_(zero_indices, zero_indices)] = zero
    return covariance


@pytest.mark.parametrize("d,k", [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])
def test_t24_zero_reading_has_zero_projection_and_derivatives(d, k):
    rng = np.random.default_rng(2400 + 10 * d + k)
    matrix = rng.normal(size=(2 * d, 2 * d))
    covariance = matrix @ matrix.T + np.eye(2 * d)
    subspace = np.linalg.qr(rng.normal(size=(d, d)))[0][:, :k]
    zero = Reading(np.zeros(d), np.zeros(d), covariance, subspace)
    phi = np.linspace(-math.pi / 2, math.pi / 2, 19)
    for derivative in _h_terms(zero.pre(), 1.7, phi, 2):
        np.testing.assert_array_equal(derivative, np.zeros_like(phi))
    state = reading_state(zero, 1.7)
    assert zero.pre()["energy"] == state["Q_plus"] == 0
    measured = ours30([zero], 1.7, readout=False)
    assert measured["g"] == [1]
    assert measured["eta"] == 1
    assert measured["median"] == pytest.approx(1.7, rel=1e-12)


def test_t24_zero_reading_leaves_current_ours30_law_and_bartlett_terms_unchanged():
    rng = np.random.default_rng(2401)
    readings = [isotropic(rng.normal(size=3), rng.normal(size=3), 1, 0.7)
                for _ in range(5)]
    matrix = rng.normal(size=(6, 6))
    zero = Reading(np.zeros(3), np.zeros(3), matrix @ matrix.T + np.eye(6),
                   np.linalg.qr(rng.normal(size=(3, 3)))[0][:, :2])
    before = ours30(readings, 2, grid=(12, 2401), readout=False)
    after = ours30(readings + [zero], 2, grid=(12, 2401), readout=False)
    assert after["g"] == before["g"] + [1]
    for key in ("eta", "J", "curv", "median"):
        assert after[key] == pytest.approx(before[key], abs=3e-13)
    np.testing.assert_allclose(after["law"], before["law"], rtol=1e-12, atol=1e-300)


@pytest.mark.parametrize("zero_cov", [np.diag([0.2, 7]), np.array([[4., 1.], [1., 2.]])])
def test_t25_independent_zero_axes_leave_projection_and_misfit_unchanged(zero_cov):
    retained = isotropic([4., 2.], [2., -2.], 1, 1)
    covariance = _block_covariance(np.eye(4), zero_cov, [0, 1, 3, 4], [2, 5])
    padded = Reading([4., 2., 0.], [2., -2., 0.], covariance)
    phi = np.linspace(-1.4, 1.4, 29)
    for left, right in zip(_h_terms(retained.pre(), 1.3, phi, 2),
                           _h_terms(padded.pre(), 1.3, phi, 2)):
        np.testing.assert_allclose(left, right, atol=3e-13)
    first, second = reading_state(retained, 1.3), reading_state(padded, 1.3)
    assert first["Q_plus"] == pytest.approx(second["Q_plus"], abs=3e-13)
    assert retained.pre()["energy"] == padded.pre()["energy"]
    assert (retained.k, padded.k) == (2, 3)
    before, after = (ours30([r], 1.3, grid=(9, 1801), readout=False)
                     for r in (retained, padded))
    assert before["g"][0] < after["g"][0] < 1
    assert np.max(np.abs(before["law"] - after["law"])) > 1e-5


def _correlated_zero_case(subspace=None):
    covariance = np.eye(4)
    covariance[0, 1] = covariance[1, 0] = 0.5
    return Reading([1., 0.], [0., 0.], covariance, subspace)


def test_t25_correlated_zero_axes_need_not_leave_projection_unchanged():
    full = _correlated_zero_case()
    marginal = isotropic([1.], [0.], 1, 1)
    assert _h(full, 1) == pytest.approx(4 / 5)
    assert _h(marginal, 1) == pytest.approx(1 / 2)
    for mass in (0.2, 1, 4):
        residual = (1 + mass * mass) / ((1 + mass * mass) ** 2 - 0.25)
        assert _h(full, mass) == pytest.approx(4 / 3 - residual)


def test_t25_k1_reduction_uses_conditional_covariance():
    restricted = _correlated_zero_case(np.array([[1.], [0.]]))
    conditional = Reading([1.], [0.], np.diag([3 / 4, 1]))
    marginal = isotropic([1.], [0.], 1, 1)
    assert _h(restricted, 1) == pytest.approx(16 / 21)
    for mass in (0.2, 1, 4):
        assert _h(restricted, mass) == pytest.approx(_h(conditional, mass))
        assert _h(restricted, mass) != pytest.approx(_h(marginal, mass))


@pytest.mark.parametrize("d", [2, 3])
def test_t25_general_k1_matches_its_conditional_one_dimensional_reading(d):
    rng = np.random.default_rng(2500 + d)
    rotation = np.linalg.qr(rng.normal(size=(d, d)))[0]
    transform = np.block([[rotation, np.zeros((d, d))],
                          [np.zeros((d, d)), rotation]])
    matrix = rng.normal(size=(2 * d, 2 * d))
    covariance = matrix @ matrix.T + np.eye(2 * d)
    y = rng.normal(size=2 * d)
    ambient = Reading(y[:d], y[d:], covariance, rotation[:, :1])
    rotated_y, rotated_cov = transform.T @ y, transform.T @ covariance @ transform
    retained = [0, d]
    outside = [i for i in range(2 * d) if i not in retained]
    a = rotated_cov[np.ix_(retained, retained)]
    b = rotated_cov[np.ix_(retained, outside)]
    c = rotated_cov[np.ix_(outside, outside)]
    effective_cov = a - b @ np.linalg.solve(c, b.T)
    effective_y = rotated_y[retained] - b @ np.linalg.solve(c, rotated_y[outside])
    reduced = Reading([effective_y[0]], [effective_y[1]], effective_cov)
    phi = np.linspace(-1.4, 1.4, 29)
    for left, right in zip(_h_terms(ambient.pre(), 0.7, phi, 2),
                           _h_terms(reduced.pre(), 0.7, phi, 2)):
        np.testing.assert_allclose(left, right, rtol=3e-12, atol=2e-13)
    assert ambient.pre()["energy"] == pytest.approx(reduced.pre()["energy"], abs=2e-13)
    for key in ("Q_plus", "kappa"):
        assert reading_state(ambient, 0.7)[key] == pytest.approx(reading_state(reduced, 0.7)[key], abs=2e-13)
    before, after = (ours30([r], 0.7, grid=(10, 2001), readout=False)
                     for r in (ambient, reduced))
    assert before["g"] == pytest.approx(after["g"], abs=2e-13)
    assert before["eta"] == after["eta"] == 1
    np.testing.assert_allclose(before["law"], after["law"], rtol=3e-12, atol=1e-300)

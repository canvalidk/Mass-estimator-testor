"""B01's calibration-side declaration: synthetic checks and preserved defaults."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.signal import lfilter
from scipy.spatial.transform import Rotation

from met import robot


def _direct_tau(values):
    values = values - values.mean()
    ac = np.correlate(values, values, "full")[len(values) - 1:] / (values @ values)
    tau = 1.0
    for lag in range(1, len(values) // 4):
        if ac[lag] <= 0:
            break
        tau += 2 * ac[lag]
    return tau


@pytest.mark.parametrize("rho", (0.2, 0.6, 0.85))
def test_tau_on_stationary_ar1(rho):
    rng = np.random.default_rng(20261030)
    values = lfilter([1.0], [1.0, -rho], rng.normal(size=210000))[10000:]
    theoretical = (1 + rho) / (1 - rho)
    assert robot.autocorrelation_time(values) == pytest.approx(theoretical, rel=0.07)


def test_tau_matches_the_existing_direct_cutoff():
    rng = np.random.default_rng(81)
    for rho in (0.0, .5, .9):
        values = lfilter([1.0], [1.0, -rho], rng.normal(size=400))
        assert robot.autocorrelation_time(values) == pytest.approx(_direct_tau(values), abs=2e-13)
    assert robot.autocorrelation_time(np.ones(20)) == 1


def test_similarity_fit_recovers_an_injected_rotation():
    rng = np.random.default_rng(422)
    acceleration = rng.normal(size=(60, 3)) * [1., 2., .5] + [2., 0., -3.]
    rotation = Rotation.from_rotvec([.1, -.07, .03]).as_matrix()
    scale, intercept = .93, np.array([.4, -.7, .2])
    force = scale * acceleration @ rotation.T + intercept
    fitted_scale, fitted_rotation, fitted_intercept = robot.fit_similarity(acceleration, force)
    assert fitted_scale == pytest.approx(scale, abs=1e-14)
    np.testing.assert_allclose(fitted_rotation, rotation, atol=1e-14)
    np.testing.assert_allclose(fitted_intercept, intercept, atol=1e-14)
    assert np.linalg.det(fitted_rotation) == pytest.approx(1., abs=1e-14)


@pytest.fixture(scope="module")
def declared():
    cal = robot.calibrate()
    declaration = robot.calibration_declaration(cal, bootstrap_replicates=80)
    return cal, declaration


@pytest.mark.parametrize("floor", (True, False))
def test_tau_one_reproduces_the_old_white_budget(declared, floor):
    cal, declaration = declared
    declaration = replace(declaration, tau_force=np.ones(3), tau_acceleration=np.ones(3))
    vf, va, ff, fa = robot.budget(cal)
    for nf, na in ((1, 1), (140, 50), (700, 250)):
        expected = np.diag(np.r_[vf / nf + (ff if floor else 0 * ff),
                                 va / na + (fa if floor else 0 * fa)])
        np.testing.assert_array_equal(
            declaration.covariance(nf, na, floor=floor, white="autocorrelated"), expected)


def test_new_modes_transport_covariance_with_the_acceleration(declared):
    _, declaration = declared
    base = declaration.covariance(140, 50, floor=True, white="autocorrelated")
    transformed = declaration.covariance(140, 50, floor=True, white="autocorrelated", rotation="static")
    joint = np.eye(6)
    joint[3:, 3:] = declaration.rotation
    np.testing.assert_allclose(transformed, joint @ base @ joint.T, rtol=1e-14, atol=1e-18)
    assert np.linalg.eigvalsh(transformed).min() > 0
    np.testing.assert_array_equal(transformed[:3, :3], base[:3, :3])


def test_old_window_defaults_and_the_explicit_independent_nominal_mode_agree(declared):
    cal, declaration = declared
    old, old_info = robot.windows("1-baseline", .2, cal)
    new, new_info = robot.windows("1-baseline", .2, cal, declaration=declaration)
    assert old_info == new_info
    for left, right in zip(old, new):
        assert left.window == right.window
        np.testing.assert_array_equal(left.force, right.force)
        np.testing.assert_array_equal(left.acceleration, right.acceleration)
        np.testing.assert_array_equal(left.cov, right.cov)


def test_static_window_rotation_changes_only_the_declared_acceleration_channel(declared):
    cal, declaration = declared
    nominal, info = robot.windows("1-baseline", .2, cal, declaration=declaration,
                                  white="autocorrelated")
    rotated, rotated_info = robot.windows("1-baseline", .2, cal, declaration=declaration,
                                         white="autocorrelated", rotation="static")
    assert info == rotated_info
    for left, right in zip(nominal, rotated):
        np.testing.assert_array_equal(left.force, right.force)
        np.testing.assert_allclose(right.acceleration, declaration.rotation @ left.acceleration,
                                   rtol=1e-14, atol=1e-14)
        assert right.window == left.window


def test_bootstrap_reporting_is_reproducible_and_uses_no_dynamic_files(monkeypatch):
    names = []
    load = robot._load

    def calibration_only(name):
        names.append(name)
        assert name.startswith("0-")
        return load(name)

    monkeypatch.setattr(robot, "_load", calibration_only)
    cal = robot.calibrate()
    a = robot.calibration_declaration(cal, bootstrap_replicates=30, seed=12)
    b = robot.calibration_declaration(cal, bootstrap_replicates=30, seed=12)
    for name in ("tau_force_bootstrap", "tau_acceleration_bootstrap", "rotation_bootstrap",
                 "rotation_vectors_bootstrap", "rotation_covariance", "rotation_matrix_covariance"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    assert a.metadata["force_block_length"] == int(np.ceil(len(cal.static_force) / 5))
    assert a.metadata["acceleration_block_length"] == int(np.ceil(len(cal.static_acceleration) / 5))
    assert a.metadata["rotation_uncertainty_propagated"] is False
    assert a.rotation_covariance.shape == (3, 3)
    assert a.rotation_matrix_covariance.shape == (9, 9)
    np.testing.assert_allclose(a.rotation_covariance,
                               np.cov(a.rotation_vectors_bootstrap, rowvar=False), rtol=1e-14)
    assert names and all(name.startswith("0-") for name in names)

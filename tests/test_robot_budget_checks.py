"""O16, 28 Oct (Claude, pressure): what the robot's shared error is, and is not.

Checks made before briefing the shared-error slot, cited in problems/O16.md (28 Oct). They
use calibration-side data (the static rows, the rest recordings) and the dynamic runs only to
test a candidate cause, never to fit the instrument that is then judged on them, except in
the cross-run offset check, where one run is fitted and the other tested.
"""
import math

import numpy as np
import pytest

from met import robot
from met.full import Reading, ours29

M_CAL = 0.9323


def _rotvec(r):
    angle = math.acos(max(-1.0, min(1.0, (np.trace(r) - 1) / 2)))
    if angle < 1e-12:
        return np.zeros(3)
    return angle / (2 * math.sin(angle)) * np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])


def _procrustes(x, y):
    """y ~ s R x + b by least squares, R a proper rotation."""
    xc, yc = x - x.mean(0), y - y.mean(0)
    u, sv, vt = np.linalg.svd(yc.T @ xc)
    d = np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))])
    r = u @ d @ vt
    s = float(np.trace(np.diag(sv) @ d) / (xc ** 2).sum())
    return s, r, y.mean(0) - s * r @ x.mean(0)


def _tau(v):
    """Integrated autocorrelation time, initial-positive-sequence cutoff."""
    v = v - v.mean()
    n = len(v)
    ac = np.correlate(v, v, "full")[n - 1:] / (v @ v)
    tau = 1.0
    for k in range(1, n // 4):
        if ac[k] <= 0:
            break
        tau += 2 * ac[k]
    return tau


@pytest.fixture(scope="module")
def static():
    cal = robot.calibrate()
    rows = robot._load("0-calibration_fts-accel")
    return cal, rows, cal.acceleration(rows[:, 6:9])


def test_the_floor_is_not_a_shared_offset(static):
    """The floor is the static fit's residual; the fit has a per-axis intercept, so the residual
    is zero-mean over the orientations: it changes with orientation, it is not one error."""
    cal, _, _ = static
    assert np.abs(cal.force_residual.mean(0)).max() < 1e-10
    assert np.abs(cal.acceleration_residual.mean(0)).max() < 1e-10


def test_static_rows_show_a_force_imu_misalignment(static):
    """F_fts ~ s R A + b on the 24 static rows: 1.153 deg, far outside its bootstrap spread."""
    _, rows, acc = static
    s, r, _ = _procrustes(acc, rows[:, :3])
    assert np.degrees(np.linalg.norm(_rotvec(r))) == pytest.approx(1.153, abs=0.005)
    assert s == pytest.approx(0.93126, abs=1e-4)
    rng = np.random.default_rng(1)
    boot = []
    for _ in range(500):
        idx = rng.integers(0, len(rows), len(rows))
        boot.append(_rotvec(_procrustes(acc[idx], rows[idx, :3])[1]))
    assert np.degrees(np.std(boot, axis=0)).max() < 0.2


def test_the_misalignment_does_not_close_the_budget_gap(static):
    """Rotating the IMU channel by the static rotation: baseline 0.2 s medians 0.9191 (floor)
    and 0.9335 (no floor), against 0.9213 and 0.9347 now. The 1.4% gap stays."""
    cal, rows, acc = static
    r = _procrustes(acc, rows[:, :3])[1]
    med = {}
    for floor in (True, False):
        rd, _ = robot.windows("1-baseline", 0.2, cal, floor=floor)
        rot = [Reading(x.force, r @ x.acceleration, x.cov, window=x.window) for x in rd]
        med[floor] = ours29(rot, robot.series_ratio(rd), readout=False)["median"]
    assert med[True] == pytest.approx(0.9191, abs=2e-4)
    assert med[False] == pytest.approx(0.9335, abs=2e-4)
    assert med[False] - med[True] > 0.012


def _omega(run):
    o = robot._load(run + "_orientations")
    t = (o[:, 0] - o[0, 0]) / 1e6
    rot = o[:, 1:].reshape(-1, 3, 3)
    w = np.einsum("nji,njk->nik", rot, np.gradient(rot, t, axis=0))
    return np.stack([w[:, 2, 1], w[:, 0, 2], w[:, 1, 0]], 1)


def _fit_offset(rd):
    f = np.array([x.force for x in rd])
    a = np.array([x.acceleration for x in rd])
    design = np.zeros((3 * len(rd), 4))
    design[:, 0] = a.ravel()
    design[:, 1:] = np.tile(np.eye(3), (len(rd), 1))
    beta = np.linalg.lstsq(design, f.ravel(), rcond=None)[0]
    return beta[0], beta[1:] / beta[0]


def test_a_constant_offset_is_not_a_lever_arm(static):
    """Steady rotation about the FTS y-axis makes an IMU off the load's centre read a constant
    offset -omega^2 r_perp, perpendicular to y. The offset fitted on run 2 is 0.22 m/s^2 (a 42 cm
    arm at 0.72 rad/s) with a y part, and on run 1 it does not bring the two budgets together."""
    cal, _, _ = static
    w = _omega("1-baseline")
    assert np.median(np.linalg.norm(w, axis=1)) == pytest.approx(0.724, abs=0.005)
    assert abs(w.mean(0)[1]) > 0.45 and np.abs(w.mean(0)[[0, 2]]).max() < 0.01
    _, b = _fit_offset(robot.windows("2-vibrations", 0.2, cal, floor=False)[0])
    assert np.linalg.norm(b) / 0.724 ** 2 > 0.3 and abs(b[1]) > 0.04
    out = {}
    for floor in (True, False):
        rd, _ = robot.windows("1-baseline", 0.2, cal, floor=floor)
        moved = [Reading(x.force, x.acceleration + b, x.cov, window=x.window) for x in rd]
        out[floor] = ours29(moved, robot.series_ratio(rd), readout=False)
    assert out[False]["median"] - out[True]["median"] > 0.005
    lo, hi = out[False]["interval"]
    assert not lo <= M_CAL <= hi


def test_rest_noise_is_autocorrelated(static):
    """The white budget var/n assumes independent samples. One second of rest data per sensor:
    tau 3.9, 2.7, 18.8 (force), 1.5, 1.6, 5.1 (acceleration), so window means' variance is
    understated about 2 to 19 times, axis by axis."""
    cal, _, _ = static
    tf = [_tau(cal.static_force[:, j]) for j in range(3)]
    ta = [_tau(cal.static_acceleration[:, j]) for j in range(3)]
    assert tf[2] > 15 and ta[2] > 4
    assert min(tf + ta) > 1.4

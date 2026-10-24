"""The robot experiment: the loader against the lab's numbers, and ours29 on the runs.

The lab's numbers are in legacy/vd-newton-lab/results/robot_cart.json (read, never run).
The ours29 values are pinned to what met/full.py gives today (24 Oct); the one place they
disagree with the experiment is a strict xfail that names both numbers.
"""
import json

import numpy as np
import pytest

from met import robot
from met.full import ours29

LAB = json.loads((robot.ROOT / "legacy" / "vd-newton-lab" / "results" / "robot_cart.json")
                 .read_text())["robot"]


@pytest.fixture(scope="module")
def cal():
    return robot.calibrate()


@pytest.fixture(scope="module")
def baseline(cal):
    rd, info = robot.windows("1-baseline", 0.2, cal)
    return rd, info, ours29(rd, robot.series_ratio(rd))


def test_calibration_matches_lab(cal):
    assert cal.mass == pytest.approx(LAB["reference_fit_mass_kg"], rel=1e-12)
    np.testing.assert_allclose(cal.force_bias, LAB["force_bias_N"], rtol=1e-12)
    np.testing.assert_allclose(cal.imu_bias, LAB["imu_bias_g"], rtol=1e-12)
    np.testing.assert_allclose(cal.static_force.std(0, ddof=1), LAB["stationary_force_sigma_N"],
                               rtol=1e-12)
    np.testing.assert_allclose(cal.static_acceleration.std(0, ddof=1),
                               LAB["stationary_acceleration_sigma_m_s2"], rtol=1e-12)
    vf, va, ff, fa = robot.budget(cal)
    np.testing.assert_allclose(np.sqrt(vf + ff), LAB["total_force_sigma_N"], rtol=1e-12)
    np.testing.assert_allclose(np.sqrt(va + fa), LAB["total_acceleration_sigma_m_s2"], rtol=1e-12)


def test_pairing_at_rest(cal):
    """At rest the pair is (m g, g): both support sign, the same positive m."""
    c = robot._load("0-calibration_fts-accel")
    F, A, g = cal.force(c[:, :3]), cal.acceleration(c[:, 6:9]), c[:, 9:12]
    assert np.all(np.sum(F * g, 1) > 0) and np.all(np.sum(A * g, 1) > 0)
    assert np.sum(F * A) / np.sum(A * A) == pytest.approx(cal.mass, rel=5e-3)


def test_window_is_the_same_interval_for_both_channels(cal, baseline):
    """T33: each reading averages both channels' raw samples over one interval."""
    rd, info, _ = baseline
    tf, F, ta, A = robot.samples("1-baseline", cal)
    for r, (start, nf, na) in zip(rd[:5], info[:5]):
        assert r.window == (start, start + 0.2)
        jf = (tf >= start) & (tf < start + 0.2)
        ia = (ta >= start) & (ta < start + 0.2)
        assert (jf.sum(), ia.sum()) == (nf, na)
        np.testing.assert_allclose(r.force, F[jf].mean(0), rtol=1e-14)
        np.testing.assert_allclose(r.acceleration, A[ia].mean(0), rtol=1e-14)
        # about 700 Hz and 254 Hz (record.json)
        assert 130 <= nf <= 150 and 45 <= na <= 56


def test_recorded_lag_fits_best(cal):
    """The recorded IMU lag (8416 us) gives a smaller Newton II residual than 0 or twice it."""
    def residual(lag):
        tf, F, ta, A = robot.samples("1-baseline", cal, lag)
        ok = (ta > tf[0] + 0.05) & (ta < tf[-1] - 0.05)
        Fi = np.column_stack([np.interp(ta[ok], tf, F[:, j]) for j in range(3)])
        return np.mean((Fi - cal.mass * A[ok]) ** 2)
    assert residual(robot.IMU_LAG_US) < 0.6 * min(residual(0), residual(2 * robot.IMU_LAG_US))


def test_baseline_values_pinned(baseline):
    """Pinned 24 Oct (tools/robot_run.py, width 0.2 s): median 0.9213, [0.9137, 0.9290]."""
    rd, _, r = baseline
    assert len(rd) == 31
    assert r["median"] == pytest.approx(0.92132, abs=2e-4)
    assert r["interval"][0] == pytest.approx(0.91373, abs=2e-4)
    assert r["interval"][1] == pytest.approx(0.92898, abs=2e-4)
    assert r["eta"] == 1.0
    assert r["readout"] == pytest.approx(r["median"], abs=1e-4)


def test_static_heldout_covers_calibration_mass(cal):
    """Static orientations read by ours29 agree with the static fit (0.9302 and 0.9320)."""
    for fit_rows, test_rows in ((slice(0, 12), slice(12, None)), (slice(12, None), slice(0, 12))):
        rd, _ = robot.static_heldout(fit_rows, test_rows)
        r = ours29(rd, robot.series_ratio(rd))
        assert r["interval"][0] < cal.mass < r["interval"][1]


@pytest.mark.xfail(strict=True, reason=(
    "Experiment vs ours29 (24 Oct): 1-baseline, 0.2 s windows, floor budget gives "
    "[0.9137, 0.9290]; the static calibration mass is 0.9323. Dynamic readings carry an "
    "error the budget does not declare (min g 0.15) and the shared floor averages down "
    "across readings (problems/O16.md)."))
def test_baseline_interval_covers_calibration_mass(cal, baseline):
    lo, hi = baseline[2]["interval"]
    assert lo <= cal.mass <= hi

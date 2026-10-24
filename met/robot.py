"""The robot experiment as readings for met.full.

The data are VD-Newton's robot runs, kept in ``legacy/vd-newton-lab/raw/robot``: a
force-torque sensor (FTS) carrying a load on a robot arm, and an IMU on the same
link. Nothing here imports or runs legacy code; the calibration is re-derived and
``tests/test_robot.py`` checks it against the lab's ``results/robot_cart.json``.

The pairing (the lab's, re-derived):

* force: the FTS reading minus its fitted bias, ``F = F_fts - b_F``. At rest this is
  the load's weight, ``m g``, in the sensor frame.
* acceleration: the IMU reading in g, bias removed, rotated into the FTS frame by
  ``R_FA`` and scaled, ``A = -R_FA (a_imu - b_a) G``. At rest this is ``g``.

Both are "support sign": ``F = m (g - a_body)`` and ``A = g - a_body``, so the same
positive m relates them and gravity is part of the push. The frame is the sensor's
own; m is frame free (L16), so the orientation files are not needed.

The window (T33). A reading averages both channels over one time interval
``[t, t + width)``, each channel from its own raw samples, after shifting the IMU
clock by the recorded phase lag (8416 us, ``record.json``). No interpolation.

The noise budget is an input (TASK.md: Sigma_i is declared with the reading). The
one declared here is Claude's, not ruled:

* white part: the steady-state files' per-axis variance, divided by the number of
  samples of that channel in the window;
* floor: the calibration fit's per-axis residual variance, not divided, because it
  is systematic (between-orientation model error), as the lab also treats it.

The floor is the same error in every reading. L13's model makes readings
independent, so the type cannot carry it; treated per reading it averages down
across readings, and the law's width then understates it. That is a finding about
the type, recorded in ``problems/O16.md``, not something this loader fixes.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from met.full import Reading

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "legacy" / "vd-newton-lab" / "raw" / "robot"
G = 9.82085                      # the lab's local g, m/s^2
R_FA = np.array([[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]])
IMU_LAG_US = 8416                # record.json: IMU phase lag
RUNS = ("1-baseline", "2-vibrations", "3-vibrations-contact")
REFERENCE_MASS_KG = 0.932        # the lab's published estimate (static calibration fit)


def _load(name):
    return np.genfromtxt(RAW / (name + ".csv"), delimiter=",", skip_header=1)


@dataclass(frozen=True)
class Calibration:
    mass: float                  # the static fit's mass, kg
    force_bias: np.ndarray       # N, FTS frame
    imu_bias: np.ndarray         # g, IMU frame
    force_residual: np.ndarray   # per orientation, N
    acceleration_residual: np.ndarray  # per orientation, m/s^2
    static_force: np.ndarray     # steady-state samples, N (bias not removed)
    static_acceleration: np.ndarray    # steady-state samples, m/s^2, FTS frame

    def force(self, raw):
        return raw - self.force_bias

    def acceleration(self, raw_g):
        return -(raw_g - self.imu_bias) @ R_FA.T * G


def calibrate():
    """Fit ``F_fts = m g + b_F`` over the 23 static orientations, and the IMU bias.

    ``g`` is the gravity vector in the FTS frame recorded with each orientation.
    """
    c = _load("0-calibration_fts-accel")
    f, a, g = c[:, :3], c[:, 6:9], c[:, 9:12]
    design = np.column_stack([g.ravel(), np.tile(np.eye(3), (len(c), 1))])
    fit = np.linalg.lstsq(design, f.ravel(), rcond=None)[0]
    mass, fb = float(fit[0]), fit[1:]
    ab = np.mean(a + g @ R_FA / G, axis=0)
    ar = -(a - ab) @ R_FA.T * G - g
    return Calibration(mass, fb, ab, f - fb - mass * g, ar,
                       _load("0-steady-state_wrench")[:, :3],
                       -_load("0-steady-state_accel") @ R_FA.T * G)


def budget(cal):
    """Per-axis variances: (white force, white acceleration, floor force, floor acceleration)."""
    return (cal.static_force.var(0, ddof=1), cal.static_acceleration.var(0, ddof=1),
            cal.force_residual.var(0, ddof=1), cal.acceleration_residual.var(0, ddof=1))


def samples(run, cal=None, lag_us=IMU_LAG_US):
    """Both channels of one run, calibrated, on their own clocks (seconds)."""
    cal = cal or calibrate()
    aa, ff = _load(run + "_accel"), _load(run + "_wrench")
    origin = min(aa[0, 0], ff[0, 0])
    ta = (aa[:, 0] - origin - lag_us) / 1e6
    tf = (ff[:, 0] - origin) / 1e6
    return tf, cal.force(ff[:, 1:4]), ta, cal.acceleration(aa[:, 1:4])


def windows(run, width, cal=None, lag_us=IMU_LAG_US, floor=True):
    """One reading per window of ``width`` seconds over the run's common time span.

    Returns the readings and, for each, ``(start, n_force, n_acceleration)``.
    """
    cal = cal or calibrate()
    tf, F, ta, A = samples(run, cal, lag_us)
    vf, va, ff, fa = budget(cal)
    if not floor:
        ff, fa = 0 * ff, 0 * fa
    t0, t1 = max(tf[0], ta[0]), min(tf[-1], ta[-1])
    readings, info = [], []
    start = t0
    while start + width <= t1 + 1e-12:
        jf = (tf >= start) & (tf < start + width)
        ia = (ta >= start) & (ta < start + width)
        nf, na = int(jf.sum()), int(ia.sum())
        if nf and na:
            cov = np.diag(np.r_[vf / nf + ff, va / na + fa])
            readings.append(Reading(F[jf].mean(0), A[ia].mean(0), cov,
                                    window=(start, start + width)))
            info.append((start, nf, na))
        start += width
    return readings, info


def series_ratio(readings):
    """The noise ratio s for the series: rms force sd over rms acceleration sd (L3')."""
    vf = np.mean([np.trace(r.cov[:3, :3]) for r in readings])
    va = np.mean([np.trace(r.cov[3:, 3:]) for r in readings])
    return float(np.sqrt(vf / va))


def static_heldout(fit_rows, test_rows, floor=True):
    """Static readings: fit the biases on ``fit_rows`` of the calibration file, read the
    other orientations as readings. Each orientation is one sample of each channel, so
    the white part is the single-sample variance. Returns (readings, fitted mass)."""
    c = _load("0-calibration_fts-accel")
    ct = c[fit_rows]
    design = np.column_stack([ct[:, 9:12].ravel(), np.tile(np.eye(3), (len(ct), 1))])
    fit = np.linalg.lstsq(design, ct[:, :3].ravel(), rcond=None)[0]
    ab = np.mean(ct[:, 6:9] + ct[:, 9:12] @ R_FA / G, axis=0)
    vf, va, ff, fa = budget(calibrate())
    if not floor:
        ff, fa = 0 * ff, 0 * fa
    cov = np.diag(np.r_[vf + ff, va + fa])
    held = c[test_rows]
    F = held[:, :3] - fit[1:]
    A = -(held[:, 6:9] - ab) @ R_FA.T * G
    return [Reading(F[i], A[i], cov) for i in range(len(held))], float(fit[0])

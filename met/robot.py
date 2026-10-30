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

The floor is the static fit's zero-mean residual over orientations, not one constant
offset shared by every reading. It remains a declared per-reading contribution;
the type has no cross-reading covariance slot (``problems/O16.md``).

``calibration_declaration`` supplies a second, explicitly selected budget. It uses
only the rest recordings and static calibration rows: the white part is multiplied
by the integrated autocorrelation time before division by the window sample count.
Its optional static rotation transports acceleration means and their declared
covariance together. The old defaults and their recorded values are unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

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


def autocorrelation_time(values):
    """Integrated autocorrelation time with the existing initial-positive cutoff.

    Match ``_tau`` in ``tests/test_robot_budget_checks.py``: demean the series,
    normalize its biased linear autocorrelation by the sum of squares, then add
    twice each lag from 1 to ``n//4 - 1``, stopping at the first nonpositive lag.
    The FFT computes the same linear correlations without a quadratic-time cost.
    A constant series uses tau=1; its sample variance is zero either way.
    This is the first-nonpositive-lag rule, not Geyer's paired-lag variant.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("autocorrelation time needs a finite vector with at least two samples")
    centred = values - values.mean()
    norm = float(centred @ centred)
    if norm == 0:
        return 1.0
    fft_length = 1 << (2 * len(values) - 1).bit_length()
    spectrum = np.fft.rfft(centred, n=fft_length)
    ac = np.fft.irfft(spectrum * spectrum.conj(), n=fft_length)[:len(values)] / norm
    tau = 1.0
    for lag in range(1, len(values) // 4):
        if ac[lag] <= 0:
            break
        tau += 2 * ac[lag]
    return float(tau)


def fit_similarity(acceleration, force):
    """Least-squares ``force = scale * R * acceleration + intercept``, R proper.

    Rows are observations. This is the centred orthogonal Procrustes similarity
    fit used by the existing static-row check, with no dynamic readings involved.
    """
    acceleration, force = np.asarray(acceleration, dtype=float), np.asarray(force, dtype=float)
    if (acceleration.ndim != 2 or acceleration.shape[1:] != (3,)
            or force.shape != acceleration.shape or len(force) < 3
            or not np.isfinite(acceleration).all() or not np.isfinite(force).all()):
        raise ValueError("similarity fit needs matching finite arrays with shape (n, 3), n >= 3")
    x, y = acceleration - acceleration.mean(0), force - force.mean(0)
    if np.linalg.matrix_rank(x) < 2:
        raise ValueError("rotation fit needs noncollinear acceleration rows")
    u, singular, vt = np.linalg.svd(y.T @ x)
    correction = np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))])
    rotation = u @ correction @ vt
    scale = float(np.trace(np.diag(singular) @ correction) / (x * x).sum())
    intercept = force.mean(0) - scale * rotation @ acceleration.mean(0)
    return scale, rotation, intercept


def _circular_blocks(values, rng, block_length):
    """One circular moving-block resample, with all sensor axes kept together."""
    count = len(values)
    starts = rng.integers(0, count, size=math.ceil(count / block_length))
    index = (starts[:, None] + np.arange(block_length)) % count
    return values[index.ravel()[:count]]


def _modes(white, rotation):
    if white not in ("independent", "autocorrelated"):
        raise ValueError("white must be independent or autocorrelated")
    if rotation not in ("nominal", "static"):
        raise ValueError("rotation must be nominal or static")


@dataclass(frozen=True)
class NoiseDeclaration:
    """Calibration-side declaration, fixed before any dynamic run is read.

    The old nominal floor vectors are retained. Selecting the static rotation
    transports the whole acceleration covariance, including that fixed floor;
    it does not refit a floor, change a bias or scale, or propagate rotation-fit
    uncertainty. The latter is reported separately in two equivalent coordinates.
    """
    white_force: np.ndarray
    white_acceleration: np.ndarray
    floor_force: np.ndarray
    floor_acceleration: np.ndarray
    tau_force: np.ndarray
    tau_acceleration: np.ndarray
    tau_force_bootstrap: np.ndarray
    tau_acceleration_bootstrap: np.ndarray
    rotation: np.ndarray
    rotation_bootstrap: np.ndarray
    rotation_vectors_bootstrap: np.ndarray
    rotation_covariance: np.ndarray
    rotation_matrix_covariance: np.ndarray
    similarity_scale: float
    similarity_intercept: np.ndarray
    metadata: dict

    @property
    def tau_force_spread(self):
        return self.tau_force_bootstrap.std(0, ddof=1)

    @property
    def tau_acceleration_spread(self):
        return self.tau_acceleration_bootstrap.std(0, ddof=1)

    def covariance(self, n_force, n_acceleration, *, floor=True,
                   white="independent", rotation="nominal"):
        """Rest variance * tau / n, plus the existing optional static-fit floor.

        ``independent`` sets tau=1. ``static`` rotates the acceleration block by
        the fitted R; the uncertainty of R itself is deliberately not propagated.
        """
        _modes(white, rotation)
        if (isinstance(n_force, (bool, np.bool_)) or isinstance(n_acceleration, (bool, np.bool_))
                or n_force <= 0 or n_acceleration <= 0):
            raise ValueError("both channel sample counts must be positive")
        tf = self.tau_force if white == "autocorrelated" else 1.0
        ta = self.tau_acceleration if white == "autocorrelated" else 1.0
        ff = self.floor_force if floor else 0 * self.floor_force
        fa = self.floor_acceleration if floor else 0 * self.floor_acceleration
        covariance = np.diag(np.r_[self.white_force * tf / n_force + ff,
                                   self.white_acceleration * ta / n_acceleration + fa])
        if rotation == "static":
            covariance[3:, 3:] = self.rotation @ covariance[3:, 3:] @ self.rotation.T
        return covariance


def calibration_declaration(cal=None, *, bootstrap_replicates=500, seed=20261030,
                            block_fraction=0.2):
    """Declare the B01 budget using calibration-side data only.

    Tau uses :func:`autocorrelation_time`, exactly the existing first-nonpositive
    cutoff. Its spread is a circular moving-block bootstrap: block length is
    ``ceil(block_fraction * n)`` for each rest recording, default n/5, with all
    three axes sampled together. Independent random streams resample force,
    acceleration and static rows. The method and default seed are fixed before
    reading any dynamic run; one-second recordings give only rough spreads.

    Static rows fit ``F_fts ~= scale * R * A_nominal + intercept``. The bootstrap
    covariance uses left relative rotation vectors ``log(R_boot @ R_hat.T)`` in
    FTS axes, radians squared. A 9x9 covariance of flattened rotation matrices is
    reported as well. Neither covariance enters any Reading.
    """
    if type(bootstrap_replicates) is not int or bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be an integer at least two")
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if (isinstance(block_fraction, bool) or not math.isfinite(block_fraction)
            or not 0 < block_fraction <= 1):
        raise ValueError("block_fraction must be finite and in (0, 1]")
    cal = cal or calibrate()
    rows = _load("0-calibration_fts-accel")
    acceleration = cal.acceleration(rows[:, 6:9])
    scale, rotation, intercept = fit_similarity(acceleration, rows[:, :3])
    streams = [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(3)]
    force_block = max(1, math.ceil(block_fraction * len(cal.static_force)))
    acceleration_block = max(1, math.ceil(block_fraction * len(cal.static_acceleration)))
    tau_force = np.array([autocorrelation_time(v) for v in cal.static_force.T])
    tau_acceleration = np.array([autocorrelation_time(v) for v in cal.static_acceleration.T])
    force_boot, acceleration_boot, rotation_boot = [], [], []
    for _ in range(bootstrap_replicates):
        force_sample = _circular_blocks(cal.static_force, streams[0], force_block)
        acceleration_sample = _circular_blocks(cal.static_acceleration, streams[1], acceleration_block)
        force_boot.append([autocorrelation_time(v) for v in force_sample.T])
        acceleration_boot.append([autocorrelation_time(v) for v in acceleration_sample.T])
        index = streams[2].integers(0, len(rows), len(rows))
        rotation_boot.append(fit_similarity(acceleration[index], rows[index, :3])[1])
    rotation_boot = np.asarray(rotation_boot)
    rotation_vectors = Rotation.from_matrix(rotation_boot @ rotation.T).as_rotvec()
    vf, va, ff, fa = budget(cal)
    metadata = dict(
        seed=seed, bootstrap_replicates=bootstrap_replicates, block_fraction=block_fraction,
        bootstrap="circular moving blocks; axes of each sensor resampled jointly",
        tau="biased normalized ACF; first nonpositive lag; lags 1..n//4-1",
        force_samples=len(cal.static_force), acceleration_samples=len(cal.static_acceleration),
        force_block_length=force_block, acceleration_block_length=acceleration_block,
        static_rows=len(rows), rotation_covariance_coordinates="log(R_boot @ R_hat.T), FTS axes",
        rotation_covariance_units="radian^2", rotation_uncertainty_propagated=False,
        floor="existing nominal static-fit residual variances; transported, not refitted",
        sources=["0-steady-state_wrench.csv", "0-steady-state_accel.csv",
                 "0-calibration_fts-accel.csv"],
    )
    return NoiseDeclaration(vf, va, ff, fa, tau_force, tau_acceleration,
                            np.asarray(force_boot), np.asarray(acceleration_boot),
                            rotation, rotation_boot, rotation_vectors,
                            np.cov(rotation_vectors, rowvar=False),
                            np.cov(rotation_boot.reshape(bootstrap_replicates, 9), rowvar=False),
                            scale, intercept, metadata)


def samples(run, cal=None, lag_us=IMU_LAG_US):
    """Both channels of one run, calibrated, on their own clocks (seconds)."""
    cal = cal or calibrate()
    aa, ff = _load(run + "_accel"), _load(run + "_wrench")
    origin = min(aa[0, 0], ff[0, 0])
    ta = (aa[:, 0] - origin - lag_us) / 1e6
    tf = (ff[:, 0] - origin) / 1e6
    return tf, cal.force(ff[:, 1:4]), ta, cal.acceleration(aa[:, 1:4])


def windows(run, width, cal=None, lag_us=IMU_LAG_US, floor=True, *,
            declaration=None, white="independent", rotation="nominal"):
    """One reading per window of ``width`` seconds over the run's common time span.

    Returns the readings and, for each, ``(start, n_force, n_acceleration)``.
    Defaults retain the original declaration exactly. Explicit new settings need
    a calibration-side ``declaration``; no calibration method is fit on this run.
    """
    cal = cal or calibrate()
    _modes(white, rotation)
    if declaration is None and (white != "independent" or rotation != "nominal"):
        raise ValueError("new white/rotation settings require a calibration declaration")
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
            cov = (np.diag(np.r_[vf / nf + ff, va / na + fa]) if declaration is None
                   else declaration.covariance(nf, na, floor=floor, white=white, rotation=rotation))
            acceleration = A[ia].mean(0)
            if rotation == "static":
                acceleration = declaration.rotation @ acceleration
            readings.append(Reading(F[jf].mean(0), acceleration, cov,
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

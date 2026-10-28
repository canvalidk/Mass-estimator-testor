"""O16's common-frame robot excitation and local shared-error diagnostics.

Uses the maintained robot loader, never legacy programs. Fitted tilt is descriptive:
noisy regressors, shared calibration uncertainty and motion-model error can imitate it.
"""
import hashlib
import math

import numpy as np

from . import robot
from .full import Reading, misfit_weight, reading_state
from .shared_error_geometry import excitation_spectrum, score_tilt, shared_tilt_diagnostic


ATTRIBUTION = {
    "creator": "Aleksander Skrede",
    "title": "Accelerometer and Force/Torque Sensor Measurements for Parameter and State Estimation of an Unknown Robot End Effector",
    "doi": "10.5281/zenodo.11096791",
    "license": "CC-BY-4.0",
}


def describe_readings(readings, reference_mass, alpha=0.05):
    """Observe permitted measured data plus an external calibration mass.

    The reference mass is used for investigator diagnostics, not sent to ours29.
    The same sensor frame is retained across readings. No latent vectors are fitted
    or recovered by the excitation-spectrum calculation.
    """
    readings = list(readings)
    if len(readings) < 2:
        raise ValueError("robot shared-error check needs at least two readings")
    acceleration = np.array([r.acceleration for r in readings])
    covariance = np.array([r.cov[3:, 3:] for r in readings])
    spectrum = excitation_spectrum(acceleration, covariance, alpha=alpha)
    tilt = shared_tilt_diagnostic(readings, reference_mass)
    training = shared_tilt_diagnostic(readings[::2], reference_mass)
    heldout = score_tilt(readings[1::2], reference_mass,
                        training["delta_mass"], training["rotation_vector"])
    c = robot.series_ratio(readings)
    states = [reading_state(r, c) for r in readings]
    q = np.array([s["Q"] for s in states])
    g = np.array([misfit_weight(s["Q"], s["nu"]) for s in states])
    return dict(
        readings=len(readings), reference_mass=float(reference_mass),
        acceleration=acceleration, force=np.array([r.force for r in readings]),
        acceleration_covariance=covariance, spectrum=spectrum,
        tilt=tilt, odd_row_training=training, even_row_score=heldout,
        whole_line_Q=q, misfit_weights=g, summary=dict(
            singular_values=spectrum["singular_values"].tolist(),
            centered_singular_values=spectrum["centered_singular_values"].tolist(),
            noise_radius=spectrum["noise_radius"],
            rank_lower=spectrum["certified_uncentered_rank_lower_bound"],
            centered_rank_lower=spectrum["certified_centered_rank_lower_bound"],
            observed_rotation_singular_values=tilt["rotation_singular_values"].tolist(),
            observed_rotation_rank=tilt["rotation_rank"],
            local_rotation_degrees=float(np.linalg.norm(tilt["rotation_vector"]) * 180 / math.pi),
            local_delta_mass=tilt["delta_mass"],
            raw_cost=tilt["raw_cost"], fitted_cost=tilt["fitted_cost"],
            mass_only_cost=tilt["mass_only_cost"],
            heldout_cost_improvement=heldout["improvement"],
            min_g=float(g.min()), median_Q=float(np.median(q))),
        assumptions=[
            "The physical acceleration includes gravity and stays in the FTS sensor frame.",
            "Uncentered excitation and centered variation are different diagnostics.",
            "Rank lower bounds are conditional on zero-mean correct covariance declarations, not model validation.",
            "Noise-subtracted raw Gram eigenvalues are signed moment diagnostics, not certainty about latent rank.",
            "The floor may be shared; odd/even robot scores are descriptive and are not independent significance tests.",
            "The local tilt fit uses noisy accelerations; biases and other model errors may imitate a shared field.",
            "The type and production estimator are not corrected or extended by these diagnostics.",
        ])


def _manifest(names):
    return {name: hashlib.sha256((robot.RAW / name).read_bytes()).hexdigest() for name in names}


def analyze_robot(widths=(0.05, 0.2, 1.0), floors=(True, False)):
    """Return dynamic cases and a descriptive broad-direction static control.

    Both raw-data hashes and attribution travel with every output case. The static
    control uses all stored calibration rows and is explicitly not held-out validation.
    """
    if (not isinstance(widths, (list, tuple)) or not widths
            or any(isinstance(w, bool) or not isinstance(w, (int, float))
                   or not math.isfinite(w) or w <= 0 for w in widths)):
        raise ValueError("robot widths must be finite positive seconds")
    if not isinstance(floors, (list, tuple)) or not floors or any(type(f) is not bool for f in floors):
        raise ValueError("robot floors must select true and/or false")
    cal = robot.calibrate()
    common = ["record.json", "0-calibration_fts-accel.csv", "0-steady-state_wrench.csv", "0-steady-state_accel.csv"]
    common_hashes = _manifest(common)
    rows = []
    for run in robot.RUNS:
        source = dict(common_hashes, **_manifest([run + "_accel.csv", run + "_wrench.csv"]))
        for width in widths:
            for floor in floors:
                readings, windows = robot.windows(run, width, cal, floor=floor)
                row = describe_readings(readings, cal.mass)
                row.update(name=f"{run}-window-{width:g}-floor-{int(floor)}", run=run,
                           width=float(width), floor=floor, window_samples=windows,
                           attribution=ATTRIBUTION, source_sha256=source,
                           calibration_status=f"Dynamic run separate from the {len(cal.force_residual)} calibration rows; shared calibration uncertainty remains.",
                           pairing_status=("Known contact force: FTS is not the net force; tilt fit is not a sensor-calibration diagnosis."
                                           if run == "3-vibrations-contact" else
                                           "Support-sign pairing from met.robot; unmodeled motion/IMU-offset effects remain possible."))
                rows.append(row)
    raw = robot._load("0-calibration_fts-accel")
    force, acceleration = cal.force(raw[:, :3]), cal.acceleration(raw[:, 6:9])
    vf, va, ff, fa = robot.budget(cal)
    for floor in floors:
        covariance = np.diag(np.r_[vf + (ff if floor else 0), va + (fa if floor else 0)])
        readings = [Reading(f, a, covariance) for f, a in zip(force, acceleration)]
        row = describe_readings(readings, cal.mass)
        row.update(name=f"static-calibration-floor-{int(floor)}", run="static-calibration", width=None,
                   floor=floor, attribution=ATTRIBUTION, source_sha256=common_hashes,
                   recorded_gravity_singular_values=np.linalg.svd(raw[:, 9:12], compute_uv=False),
                   calibration_status=f"Descriptive in-calibration control: these same {len(raw)} rows fit mass and biases; not held out.",
                   pairing_status="Static support-sign pairing; recorded gravity spectrum is a reference with no uncertainty budget supplied.")
        rows.append(row)
    return rows

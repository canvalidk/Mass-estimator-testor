"""O16's serial fixed-tilt experiment and robot diagnostic collection.

The synthetic investigator owns latent vectors and the exact matched hypothesis.
The local fit and the unchanged mass estimator receive observed readings only.
Fitted improvements are descriptive; they are not assigned fitted-data p-values.
"""
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
import yaml

from .analyst import bartlett_weight, misfit_weights
from .full import isotropic
from .law import ReadingSet
from .shared_error_geometry import rotation, score_tilt, shared_tilt_diagnostic
from .world import series_rngs

GEOMETRIES = ("broad", "plane", "line_transverse", "line_axis")
_OURS29 = {"misfit": "kink"}


def _number(value, name, positive=True):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (positive and value <= 0)):
        raise ValueError(f"{name} must be a finite {'positive ' if positive else ''}number")
    return float(value)


def validate(data):
    required = {"kind", "name", "question", "cases", "angles_degrees", "readings",
                "replicates", "seed", "mass", "acceleration_length", "force_sd",
                "acceleration_sd", "rotation_axis", "matched_hypothesis_degrees"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"shared_tilt needs {sorted(required)}")
    if set(data) - required - {"robot", "notes"}:
        raise ValueError("unknown shared_tilt setting")
    if data["kind"] != "shared_tilt":
        raise ValueError("expected kind: shared_tilt")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for ch in name)):
        raise ValueError("shared_tilt name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("shared_tilt needs a nonempty question")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("shared_tilt seed must be a nonnegative integer")
    if type(data["replicates"]) is not int or data["replicates"] < 2:
        raise ValueError("shared_tilt replicates must be an integer >= 2")
    for key in ("cases", "angles_degrees", "readings"):
        if not isinstance(data[key], list) or not data[key]:
            raise ValueError(f"shared_tilt {key} must be a nonempty list")
    if any(case not in GEOMETRIES for case in data["cases"]):
        raise ValueError(f"shared_tilt cases must be drawn from {GEOMETRIES}")
    if len(set(data["cases"])) != len(data["cases"]):
        raise ValueError("shared_tilt cases must be unique")
    angles = [_number(value, "angle", positive=False) for value in data["angles_degrees"]]
    if len(set(angles)) != len(angles):
        raise ValueError("shared_tilt angles must be unique")
    if any(type(value) is not int or value < 4 for value in data["readings"]):
        raise ValueError("shared_tilt readings must be integers >= 4 for the fit/score split")
    if len(set(data["readings"])) != len(data["readings"]):
        raise ValueError("shared_tilt reading counts must be unique")
    settings = {key: _number(data[key], key) for key in
                ("mass", "acceleration_length", "force_sd", "acceleration_sd")}
    hypothesis = _number(data["matched_hypothesis_degrees"], "matched_hypothesis_degrees")
    if hypothesis >= 180:
        raise ValueError("matched_hypothesis_degrees must be less than 180")
    settings["matched_hypothesis_degrees"] = hypothesis
    axis = np.asarray(data["rotation_axis"], dtype=float)
    if axis.shape != (3,) or not np.all(np.isfinite(axis)) or np.max(np.abs(axis)) == 0:
        raise ValueError("rotation_axis must be a nonzero finite length-three vector")
    axis = axis / np.max(np.abs(axis))
    settings["rotation_axis"] = (axis / np.linalg.norm(axis)).tolist()
    robot = data.get("robot")
    if robot is not None:
        if not isinstance(robot, dict) or set(robot) != {"widths", "floors"}:
            raise ValueError("robot needs widths and floors, or null to omit the robot collection")
        if not isinstance(robot["widths"], list) or not robot["widths"]:
            raise ValueError("robot widths must be a nonempty list")
        for width in robot["widths"]:
            _number(width, "robot width")
        if (not isinstance(robot["floors"], list) or not robot["floors"]
                or any(type(value) is not bool for value in robot["floors"])):
            raise ValueError("robot floors must be a nonempty list of booleans")
    latent_ranks = {"broad": (3, 3), "plane": (2, 3), "line_transverse": (1, 2), "line_axis": (1, 2)}
    cells = [dict(id=f"c{i:03d}", geometry=case, angle_degrees=angle, readings=n,
                  latent_excitation_rank=latent_ranks[case][0], latent_rotation_rank=latent_ranks[case][1])
             for i, (case, angle, n) in enumerate(itertools.product(
                 data["cases"], angles, data["readings"]))]
    return dict(kind="shared_tilt", raw=data, settings=settings, cells=cells, robot=robot)


def _directions(rng, count, geometry, axis):
    axis = np.asarray(axis)
    coordinate = np.eye(3)[int(np.argmin(np.abs(axis)))]
    first = coordinate - axis * (axis @ coordinate)
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    if geometry == "broad":
        directions = rng.normal(size=(count, 3))
    elif geometry == "plane":
        directions = rng.normal(size=(count, 2)) @ np.vstack((first, second))
    else:
        line = axis if geometry == "line_axis" else first
        directions = rng.choice([-1.0, 1.0], size=(count, 1)) * line
    return directions / np.linalg.norm(directions, axis=1, keepdims=True)


def investigator_matched_score(force, acceleration, latent_acceleration, reference_mass,
                               injected_rotation, hypothesis_rotation, force_sd, acceleration_sd):
    """Exact Gaussian score conditional on investigator-only latent excitation.

    No fit is performed. The hypothesis must be fixed before these observations.
    Returned NaNs mark a zero-information hypothesis, not a failed detection.
    """
    latent = np.asarray(latent_acceleration)
    template = reference_mass * latent @ (hypothesis_rotation - np.eye(3)).T
    true_pattern = reference_mass * latent @ (injected_rotation - np.eye(3)).T
    residual = np.asarray(force) - reference_mass * np.asarray(acceleration)
    variance = force_sd ** 2 + reference_mass ** 2 * acceleration_sd ** 2
    information = np.sum(template ** 2, axis=(-2, -1)) / variance
    numerator = np.sum(template * residual, axis=(-2, -1)) / variance
    mean_numerator = np.sum(template * true_pattern, axis=(-2, -1)) / variance
    root = np.sqrt(information)
    active = information > 0
    score = np.divide(numerator, root, out=np.full_like(root, np.nan), where=active)
    mean = np.divide(mean_numerator, root, out=np.full_like(root, np.nan), where=active)
    return dict(matched_score=score, matched_mean=mean, centered_score=score - mean,
                matched_information=information)


def simulate_cell(cell, settings, replicates, seed):
    """Generate one cell and inspect each series; no latent truth enters a fit."""
    n, m = cell["readings"], settings["mass"]
    sf, sa = settings["force_sd"], settings["acceleration_sd"]
    axis = settings["rotation_axis"]
    injected = rotation(axis, math.radians(cell["angle_degrees"]))
    hypothesis = rotation(axis, math.radians(settings["matched_hypothesis_degrees"]))
    key = {k: v for k, v in cell.items() if k != "id"}
    key["settings"] = settings
    latent, force, acceleration = (np.empty((replicates, n, 3)) for _ in range(3))
    for i, rng in enumerate(series_rngs(seed, key, replicates)):
        latent[i] = settings["acceleration_length"] * _directions(rng, n, cell["geometry"], axis)
        force[i] = m * latent[i] @ injected.T + sf * rng.normal(size=(n, 3))
        acceleration[i] = latent[i] + sa * rng.normal(size=(n, 3))
    observed = ReadingSet(force, acceleration, np.full((replicates, n), sf),
                          np.full((replicates, n), sa))
    gap = np.hypot(observed.P - observed.Q, 2 * observed.D)
    q = np.maximum(0.0, 0.5 * (observed.P + observed.Q - gap))
    values = dict(mean_Q=q.mean(axis=1), min_Q=q.min(axis=1),
                  min_g=misfit_weights(observed, _OURS29).min(axis=1),
                  eta=bartlett_weight(observed, _OURS29).copy())
    values.update(investigator_matched_score(force, acceleration, latent, m, injected, hypothesis, sf, sa))
    if cell["geometry"] == "line_axis":
        # The exact hypothesis is null on this line; suppress roundoff for non-coordinate axes.
        values["matched_information"][:] = 0
        for key in ("matched_score", "matched_mean", "centered_score"):
            values[key][:] = np.nan
    for key in ("delta_mass", "raw_cost", "fitted_cost", "heldout_improvement", "rotation_rank"):
        values[key] = np.empty(replicates)
    values["rotation_vector"] = np.empty((replicates, 3))
    values["observable_rotation_norm"] = np.empty(replicates)
    for i in range(replicates):
        readings = [isotropic(f, a, sf, sa) for f, a in zip(force[i], acceleration[i])]
        fit = shared_tilt_diagnostic(readings, m)
        training = shared_tilt_diagnostic(readings[::2], m)
        held = score_tilt(readings[1::2], m, training["delta_mass"], training["rotation_vector"])
        for key in ("delta_mass", "raw_cost", "fitted_cost", "rotation_rank", "rotation_vector"):
            values[key][i] = fit[key]
        observable = fit["rotation_vector"]
        if cell["geometry"].startswith("line_"):
            # Investigator-only observable component under the declared latent line.
            # No latent vector was used to obtain the fitted parameters above.
            direction = latent[i, 0] / np.linalg.norm(latent[i, 0])
            observable = observable - direction * (direction @ observable)
        values["observable_rotation_norm"][i] = np.linalg.norm(observable)
        values["heldout_improvement"][i] = held["improvement"]
    return values


def _summary(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    return dict(mean=float(finite.mean()) if len(finite) else None,
                sd=float(finite.std(ddof=1)) if len(finite) > 1 else None,
                minimum=float(finite.min()) if len(finite) else None,
                maximum=float(finite.max()) if len(finite) else None,
                finite_count=int(len(finite)), total_count=int(values.size))


def _markdown(record):
    lines = ["# Shared sensor tilt: controlled field and excitation", "", record["question"], "",
        "A fixed rotation predicts a direction-dependent residual field; raw residual directions need not align.",
        "The centered matched score is conditionally N(0,1), using investigator-only latent excitation and a fixed hypothesis.",
        "Local fitted tilt and held-out cost changes are descriptive; noisy regressors preclude fitted-data p-values.",
        "For synthetic lines, the displayed fit norm projects onto the known latent line's observable rotation components; full vectors remain in the arrays.",
        "No calibration correction or mass-estimator input is changed.", "",
        "| Cell | Geometry | Injected degrees | N | Centered score mean / SD | Observable fitted norm (degrees) | Held-out improvement | Mean Q | eta mean |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in record["cells"]:
        s = row["summary"]
        score = s["centered_score"]
        score_text = (f'{score["mean"]:.5g} / {score["sd"]:.5g}' if score["finite_count"] else "inapplicable")
        lines.append(f'| {row["id"]} | {row["geometry"]} | {row["angle_degrees"]:g} | {row["readings"]} | '
            f'{score_text} | {math.degrees(s["observable_rotation_norm"]["mean"]):.5g} | '
            f'{s["heldout_improvement"]["mean"]:.6g} | {s["mean_Q"]["mean"]:.6g} | {s["eta"]["mean"]:.6g} |')
    if record["robot"]:
        lines += ["", "## Robot excitation and shared-field diagnostics", "",
            "Aleksander Skrede's [robot dataset](https://doi.org/10.5281/zenodo.11096791), CC-BY-4.0.",
            "Spectra use physical acceleration in the common FTS frame, including gravity. Units are m/s^2.",
            "Rank entries are conservative lower bounds under correct zero-mean covariance declarations, at alpha=0.05 per case; they are not simultaneous guarantees over this collection.",
            "A weak third direction is not a proof of exact rank two. Noisy observed-design numerical rank is retained separately.",
            "Static rows are an in-calibration control, not held-out validation; contact rows violate the net-force pairing.", "",
            "| Run | Window (s) | Floor | N | Physical singular values | Centered singular values | Noise radius | Rank lower: raw / centered | Local fit degrees | Even-row cost improvement |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for row in record["robot"]:
            s = row["summary"]
            raw = ", ".join(f"{x:.6g}" for x in s["singular_values"])
            centered = ", ".join(f"{x:.6g}" for x in s["centered_singular_values"])
            width = "static" if row["width"] is None else f'{row["width"]:g}'
            lines.append(f'| {row["run"]} | {width} | {row["floor"]} | {row["readings"]} | {raw} | '
                         f'{centered} | {s["noise_radius"]:.6g} | {s["rank_lower"]} / {s["centered_rank_lower"]} | '
                         f'{s["local_rotation_degrees"]:.6g} | {s["heldout_cost_improvement"]:.6g} |')
    lines += ["", "Robot diagnostics and complete arrays are retained in run.json and values.npz.",
              "Robot dimensions and costs are conditional on the declared calibration/noise budgets; no causal conclusion is assigned."]
    return "\n".join(lines) + "\n"


def run(path, out_root="results", replicates=None, progress=print):
    """Run the research recipe serially; retain exact recipes, arrays and provenance."""
    from .inspection.run import _pack
    from .runner import code_version

    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    reps = spec["raw"]["replicates"] if replicates is None else replicates
    if type(reps) is not int or reps < 2:
        raise ValueError("shared_tilt replicates must be an integer >= 2")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base, suffix = Path(out_root) / spec["raw"]["name"], 1
    out = base / stamp
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    if replicates is None:
        (out / "preset.yaml").write_bytes(source)
    else:
        (out / "source_preset.yaml").write_bytes(source)
        effective = dict(spec["raw"], replicates=reps)
        (out / "preset.yaml").write_text(yaml.safe_dump(effective, sort_keys=False), encoding="utf-8")
    start, arrays, results = time.perf_counter(), {}, []
    for cell in spec["cells"]:
        if progress:
            progress(f'cell {cell["id"]}: {cell["geometry"]}, {cell["angle_degrees"]:g} degrees, N={cell["readings"]}')
        values = simulate_cell(cell, spec["settings"], reps, spec["raw"]["seed"])
        for key, value in values.items():
            arrays[f'{cell["id"]}/{key}'] = value
        summaries = {key: _summary(value) for key, value in values.items() if key != "rotation_vector"}
        summaries["rotation_norm"] = _summary(np.linalg.norm(values["rotation_vector"], axis=1))
        results.append({**cell, "summary": summaries,
            "matched_status": "applicable" if summaries["matched_score"]["finite_count"] else "inapplicable"})
    robot = []
    if spec["robot"] is not None:
        from .robot_shared_error import analyze_robot
        if progress:
            progress("robot: common-frame excitation and shared-field diagnostics")
        robot = analyze_robot(spec["robot"]["widths"], spec["robot"]["floors"])
    packed_robot = _pack(robot, arrays, path="robot")
    code = code_version()
    code.update(preset_sha256=hashlib.sha256(source).hexdigest(),
                effective_preset_sha256=hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest())
    record = dict(kind="shared_tilt", preset=spec["raw"]["name"], question=spec["raw"]["question"],
                  seed=spec["raw"]["seed"], replicates=reps, workers=1, settings=spec["settings"],
                  code=code, started=stamp, seconds=time.perf_counter() - start, cells=results,
                  robot=packed_robot,
                  assumptions=["Synthetic reference mass and matched latent pattern are investigator declarations, never estimator truth inputs.",
                               "The exact Gaussian matched-score reference conditions on latent excitation and the predeclared rotation hypothesis.",
                               "Observed-design rotation_rank can be 3 on a true line; declared latent ranks are separate investigator metadata.",
                               "Local fits and held-out cost differences are descriptive, not p-values or recovered instrument errors."])
    np.savez_compressed(out / "values.npz", **arrays)
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    (out / "summary.md").write_text(_markdown(record), encoding="utf-8")
    with (out / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["cell", "geometry", "angle_degrees", "readings", "matched_status", "centered_score_mean",
                  "centered_score_sd", "full_rotation_norm_degrees", "observable_rotation_norm_degrees",
                  "heldout_improvement", "mean_Q", "eta"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in results:
            s = row["summary"]
            writer.writerow(dict(cell=row["id"], geometry=row["geometry"], angle_degrees=row["angle_degrees"],
                readings=row["readings"], matched_status=row["matched_status"],
                centered_score_mean=s["centered_score"]["mean"], centered_score_sd=s["centered_score"]["sd"],
                full_rotation_norm_degrees=math.degrees(s["rotation_norm"]["mean"]),
                observable_rotation_norm_degrees=math.degrees(s["observable_rotation_norm"]["mean"]),
                heldout_improvement=s["heldout_improvement"]["mean"], mean_Q=s["mean_Q"]["mean"], eta=s["eta"]["mean"]))
    if progress:
        progress(f"wrote {out}")
    return out, record

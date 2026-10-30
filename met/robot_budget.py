"""B01: a serial calibration-side declaration and the full robot rerun matrix.

Calibration methods and the four questions were fixed in O16 before the rerun.
No dynamic data fit an input; no covariance across readings or estimator change
is introduced. Bootstrap spreads are descriptive, not confidence intervals.
"""
import copy
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import yaml

from . import robot
from .full import ours29, ours30


REFERENCE_MASS = 0.9323
RUNS = ("1-baseline", "2-vibrations")
FIELDS = ("id", "run", "width", "white", "floor", "rotation", "n", "c",
          "median", "lo", "hi", "covers", "eta", "min_g", "median_g", "median29",
          "median_g29", "raw_lo", "raw_hi")
ATTRIBUTION = dict(creator="Aleksander Skrede", doi="10.5281/zenodo.11096791",
                   license="CC-BY-4.0")


def validate(data):
    required = {"kind", "name", "question", "runs", "widths", "white", "floors",
                "rotations", "bootstrap_replicates", "seed", "block_fraction"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"robot_budget needs {sorted(required)}")
    if set(data) - required - {"notes"}:
        raise ValueError("unknown robot_budget setting")
    if data["kind"] != "robot_budget":
        raise ValueError("expected kind: robot_budget")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for ch in name)):
        raise ValueError("robot_budget name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("robot_budget needs a nonempty question")
    selectors = {"runs": RUNS, "white": ("independent", "autocorrelated"),
                 "floors": (True, False), "rotations": ("nominal", "static")}
    for key in (*selectors, "widths"):
        values = data[key]
        if not isinstance(values, list) or not values:
            raise ValueError(f"robot_budget {key} must be a nonempty list")
        if key == "widths":
            if any(isinstance(v, bool) or not isinstance(v, (int, float))
                   or not math.isfinite(v) or v <= 0 for v in values):
                raise ValueError("robot_budget widths must be finite positive seconds")
        elif key == "floors":
            if any(type(v) is not bool for v in values):
                raise ValueError("robot_budget floors must contain booleans")
        elif any(not isinstance(v, str) or v not in selectors[key] for v in values):
            raise ValueError(f"robot_budget {key} must select from {selectors[key]}")
        if len(set(values)) != len(values):
            raise ValueError(f"robot_budget {key} must be unique")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("robot_budget seed must be a nonnegative integer")
    if type(data["bootstrap_replicates"]) is not int or data["bootstrap_replicates"] < 2:
        raise ValueError("robot_budget bootstrap_replicates must be an integer >= 2")
    fraction = data["block_fraction"]
    if (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
            or not math.isfinite(fraction) or not 0 < fraction <= 1):
        raise ValueError("robot_budget block_fraction must be in (0, 1]")
    cells = [dict(id=f"c{i:03d}", run=run, width=float(width), white=white,
                  floor=floor, rotation=rotation)
             for i, (run, width, white, floor, rotation) in enumerate(itertools.product(
                 data["runs"], data["widths"], data["white"], data["floors"], data["rotations"]))]
    return dict(kind="robot_budget", raw=copy.deepcopy(data), cells=cells)


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def _questions(rows):
    gaps, widths, weights = [], [], []
    for run, width, white, rotation in sorted({(r["run"], r["width"], r["white"], r["rotation"])
                                               for r in rows if r["run"] == RUNS[0]}):
        selected = {r["floor"]: r for r in rows if (r["run"], r["width"], r["white"], r["rotation"])
                    == (run, width, white, rotation)}
        if set(selected) == {True, False}:
            gap = abs(selected[True]["median"] - selected[False]["median"]) / REFERENCE_MASS
            gaps.append(dict(run=run, width=width, white=white, rotation=rotation,
                             gap_fraction=gap, below_half_percent=gap < .005))
    for run, white, floor, rotation in sorted({(r["run"], r["white"], r["floor"], r["rotation"])
                                               for r in rows}):
        selected = [r for r in rows if (r["run"], r["white"], r["floor"], r["rotation"])
                    == (run, white, floor, rotation)]
        mutual = all(a["lo"] <= b["median"] <= a["hi"] for a in selected for b in selected)
        overlap = max(r["lo"] for r in selected) <= min(r["hi"] for r in selected)
        widths.append(dict(run=run, white=white, floor=floor, rotation=rotation,
                           widths=[r["width"] for r in selected],
                           multiple_widths=len(selected) > 1,
                           mutual_median_containment=mutual, common_interval=overlap,
                           all_widths_cover=all(r["covers"] for r in selected)))
    for run, width, floor, rotation in sorted({(r["run"], r["width"], r["floor"], r["rotation"])
                                               for r in rows}):
        selected = {r["white"]: r for r in rows if (r["run"], r["width"], r["floor"], r["rotation"])
                    == (run, width, floor, rotation)}
        if set(selected) == {"independent", "autocorrelated"}:
            old, new = selected["independent"], selected["autocorrelated"]
            weights.append(dict(run=run, width=width, floor=floor, rotation=rotation,
                                min_g_before=old["min_g"], min_g_after=new["min_g"],
                                median_g_before=old["median_g"], median_g_after=new["median_g"],
                                median_g_ratio=new["median_g"] / old["median_g"]))
    return dict(floor_gaps=gaps, width_compatibility_and_coverage=widths, misfit_weight_changes=weights)


def analyze(spec, progress=None):
    raw = spec["raw"]
    # This must precede dynamic window construction: no results choose the declaration.
    cal = robot.calibrate()
    declaration = robot.calibration_declaration(
        cal, bootstrap_replicates=raw["bootstrap_replicates"], seed=raw["seed"],
        block_fraction=raw["block_fraction"])
    calibration = {key: _plain(getattr(declaration, key)) for key in (
        "tau_force", "tau_acceleration", "tau_force_spread", "tau_acceleration_spread",
        "rotation", "rotation_covariance", "rotation_matrix_covariance", "similarity_scale",
        "similarity_intercept", "metadata")}
    for name in ("force", "acceleration"):
        draws = getattr(declaration, f"tau_{name}_bootstrap")
        calibration[f"tau_{name}_percentiles"] = np.quantile(draws, [.025, .975], axis=0).tolist()
    arrays = {key: np.asarray(getattr(declaration, key)) for key in (
        "tau_force_bootstrap", "tau_acceleration_bootstrap", "rotation_bootstrap",
        "rotation_vectors_bootstrap")}
    rows = []
    for cell in spec["cells"]:
        readings, info = robot.windows(cell["run"], cell["width"], cal,
                                       floor=cell["floor"], declaration=declaration,
                                       white=cell["white"], rotation=cell["rotation"])
        if not readings:
            raise ValueError(f"{cell['id']}: width leaves no robot windows")
        c = robot.series_ratio(readings)
        current = ours30(readings, c, readout=False)
        previous = ours29(readings, c, readout=False)
        lo, hi = current["interval"]
        rlo, rhi = current["raw_interval"]
        row = dict(**cell, n=len(readings), c=c, median=current["median"], lo=lo, hi=hi,
                   covers=bool(lo <= REFERENCE_MASS <= hi), eta=current["eta"],
                   min_g=float(np.min(current["g"])), median_g=float(np.median(current["g"])),
                   median29=previous["median"], median_g29=float(np.median(previous["g"])),
                   raw_lo=rlo, raw_hi=rhi)
        rows.append(_plain(row))
        prefix = cell["id"]
        arrays.update({f"{prefix}_force": np.array([r.force for r in readings]),
                       f"{prefix}_acceleration": np.array([r.acceleration for r in readings]),
                       f"{prefix}_covariance": np.array([r.cov for r in readings]),
                       f"{prefix}_window_info": np.asarray(info),
                       f"{prefix}_g30": np.asarray(current["g"]),
                       f"{prefix}_g29": np.asarray(previous["g"])})
        if progress:
            progress(f"{prefix} {cell['run']} {cell['width']:g}s {cell['white']} "
                     f"floor={cell['floor']} rotation={cell['rotation']}: {current['median']:.6f}")
    sources = {"0-calibration_fts-accel.csv", "0-steady-state_wrench.csv", "0-steady-state_accel.csv"}
    for run in raw["runs"]:
        sources.update((run + "_wrench.csv", run + "_accel.csv"))
    from .runner import code_version
    hashes = {name: hashlib.sha256((robot.RAW / name).read_bytes()).hexdigest() for name in sorted(sources)}
    record = dict(kind="robot_budget", workers=1, calibration_mass=cal.mass,
                  coverage_reference_mass=REFERENCE_MASS, calibration=calibration,
                  rows=rows, questions=_questions(rows),
                  provenance=dict(source_sha256=hashes, attribution=ATTRIBUTION, code=code_version()))
    return _plain(record), arrays


def summary(record):
    c = record["calibration"]
    lines = ["# B01 robot budget rerun", "", "Ours30 is the current probe; ours29 medians are comparisons.",
             "Every declaration comes from calibration-side data. No equation is changed.",
             "", f"Calibration mass: {record['calibration_mass']:.10f} kg; coverage reference: {REFERENCE_MASS} kg.",
             "", "Bootstrap SDs and percentiles are descriptive short-record spreads, not confidence intervals.",
             "", "| channel/axis | tau | bootstrap SD | 2.5% | 97.5% |",
             "| --- | ---: | ---: | ---: | ---: |"]
    for channel in ("force", "acceleration"):
        for j, axis in enumerate("xyz"):
            p = c[f"tau_{channel}_percentiles"]
            lines.append(f"| {channel}/{axis} | {c[f'tau_{channel}'][j]:.6f} | "
                         f"{c[f'tau_{channel}_spread'][j]:.6f} | {p[0][j]:.6f} | {p[1][j]:.6f} |")
    lines += ["", "Additional static rotation R, after nominal R_FA:", "", "```text"]
    lines += [" ".join(f"{v: .10f}" for v in row) for row in c["rotation"]]
    lines += ["```", "", "Bootstrap covariance of log(R_boot R_fit.T), FTS axes, radians squared; not propagated:", "", "```text"]
    lines += [" ".join(f"{v: .8e}" for v in row) for row in c["rotation_covariance"]]
    lines += ["```", "", f"Similarity scale {c['similarity_scale']:.10f} and intercept {c['similarity_intercept']} are diagnostics only.",
              "", "B = 1-baseline; V = 2-vibrations; I = independent white; A = autocorrelated white; N/S = nominal/static rotation.",
              "", "| run | s | white | floor | R | N | median30 | 95% interval30 | covers | eta30 | min g30 | median g30 | median29 |",
              "| --- | ---: | --- | --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |"]
    for r in record["rows"]:
        lines.append(f"| {'B' if r['run']==RUNS[0] else 'V'} | {r['width']:g} | "
                     f"{'I' if r['white']=='independent' else 'A'} | {'on' if r['floor'] else 'off'} | "
                     f"{'N' if r['rotation']=='nominal' else 'S'} | {r['n']} | {r['median']:.6f} | "
                     f"[{r['lo']:.6f}, {r['hi']:.6f}] | {'yes' if r['covers'] else 'no'} | "
                     f"{r['eta']:.5f} | {r['min_g']:.6g} | {r['median_g']:.6g} | {r['median29']:.6f} |")
    q = record["questions"]
    lines += ["", "Q1: baseline floor/no-floor gap / 0.9323; threshold is strictly below 0.5%.",
              "", "| s | white | rotation | gap % | below 0.5% |", "| ---: | --- | --- | ---: | --- |"]
    for r in q["floor_gaps"]:
        lines.append(f"| {r['width']:g} | {r['white']} | {r['rotation']} | {100*r['gap_fraction']:.6f} | "
                     f"{'yes' if r['below_half_percent'] else 'no'} |")
    lines += ["", "Q2/Q3: each median lies in every other width's interval; common intersection; all widths cover 0.9323.",
              "The widths reuse data; these are descriptive comparisons.", "",
              "| run | white | floor | rotation | mutual medians | common interval | all cover |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in q["width_compatibility_and_coverage"]:
        mutual = ("yes" if r["mutual_median_containment"] else "no") if r["multiple_widths"] else "one width"
        lines.append(f"| {r['run']} | {r['white']} | {'on' if r['floor'] else 'off'} | {r['rotation']} | "
                     f"{mutual} | {'yes' if r['common_interval'] else 'no'} | {'yes' if r['all_widths_cover'] else 'no'} |")
    lines += ["", "Q4: matching ours30 rows, independent -> autocorrelated; no declaration is chosen by these changes.", "",
              "| run | s | floor | rotation | min g before/after | median g before/after | median ratio |",
              "| --- | ---: | --- | --- | --- | --- | ---: |"]
    for r in q["misfit_weight_changes"]:
        lines.append(f"| {r['run']} | {r['width']:g} | {'on' if r['floor'] else 'off'} | {r['rotation']} | "
                     f"{r['min_g_before']:.6g}/{r['min_g_after']:.6g} | "
                     f"{r['median_g_before']:.6g}/{r['median_g_after']:.6g} | {r['median_g_ratio']:.6f} |")
    return "\n".join(lines) + "\n"


def run(path, out_root="results", progress=print):
    path = Path(path)
    spec = validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    record, values = analyze(spec, progress)
    out = Path(out_root) / spec["raw"]["name"] / datetime.now().strftime("%Y%m%d-%H%M%S")
    counter = 1
    while out.exists():
        out = out.with_name(out.name.rsplit("-r", 1)[0] + f"-r{counter}")
        counter += 1
    out.mkdir(parents=True)
    (out / "preset.yaml").write_bytes(path.read_bytes())
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (out / "summary.md").write_text(summary(record), encoding="utf-8")
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(record["rows"])
    np.savez_compressed(out / "values.npz", **values)
    if progress:
        progress(f"saved {out}")
    return out, record

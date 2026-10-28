"""Deterministic O2 pressure cases on one unchanged law per measured input.

This recipe compares readouts; it neither generates random data nor fits a new
law for a candidate. Composition checks separately propagate finite candidate
laws and make no Gaussian-adapter reachability claim.
"""
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
import time

import numpy as np
import yaml

from .full import Reading
from .inspection.case import Case, parse_cases
from .inspection.paths import set_prior_center
from .readout_candidates import CANDIDATE_IDS, evaluate
from .readout_law_checks import (cancellation_counterexample, counterexamples,
                                 direction_counterexample)


def _positive(value, label):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value <= 0):
        raise ValueError(f"{label} must be finite and positive")
    return float(value)


def validate(raw):
    """Expand explicit priors and exact-null padding into ordinary replay cases."""
    required = {"kind", "name", "question", "seed", "cases", "prior_centers", "null_padding"}
    if not isinstance(raw, dict) or required - raw.keys():
        raise ValueError(f"readout_candidates needs {sorted(required)}")
    if set(raw) - required - {"notes"}:
        raise ValueError("unknown readout_candidates setting")
    if raw["kind"] != "readout_candidates":
        raise ValueError("expected kind: readout_candidates")
    name = raw["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in name)):
        raise ValueError("readout_candidates name must be a simple identifier")
    if not isinstance(raw["question"], str) or not raw["question"].strip():
        raise ValueError("readout_candidates needs a nonempty question")
    if type(raw["seed"]) is not int or raw["seed"] < 0:
        raise ValueError("readout_candidates seed must be a nonnegative integer (no random draws)")
    if not isinstance(raw["prior_centers"], list) or not raw["prior_centers"]:
        raise ValueError("prior_centers must be a nonempty list")
    centers = [_positive(c, "prior center") for c in raw["prior_centers"]]
    if len(set(centers)) != len(centers):
        raise ValueError("prior_centers must be distinct")
    if not isinstance(raw["null_padding"], list) or not raw["null_padding"]:
        raise ValueError("null_padding must be a nonempty list")
    padding = []
    for item in raw["null_padding"]:
        if not isinstance(item, dict) or set(item) != {"count", "force_sd", "acceleration_sd"}:
            raise ValueError("each null_padding entry needs count, force_sd and acceleration_sd")
        if type(item["count"]) is not int or item["count"] < 0:
            raise ValueError("null-padding count must be a nonnegative integer")
        padding.append(dict(count=item["count"], force_sd=_positive(item["force_sd"], "force_sd"),
                            acceleration_sd=_positive(item["acceleration_sd"], "acceleration_sd")))
    if len({tuple(p.values()) for p in padding}) != len(padding):
        raise ValueError("null_padding entries must be distinct")
    source_cases = parse_cases(raw)
    cases, cells = [], []
    for source in source_cases:
        if len(source.inputs.readings) != 1 or source.inputs.readings[0].force.size != 2:
            raise ValueError("this readout-candidate recipe requires one 2D source reading per case")
        for center_index, center in enumerate(centers):
            base = set_prior_center(source, center)
            for padding_index, pad in enumerate(padding):
                identifier = f"c{len(cells):03d}"
                data = base.to_dict(include_truth=True)
                data["name"] = f"{source.name}--c-{center_index:03d}--padding-{padding_index:03d}"
                ids = {r["id"] for r in data["readings"]}
                covariance = np.diag([pad["force_sd"] ** 2] * 2
                                     + [pad["acceleration_sd"] ** 2] * 2).tolist()
                for index in range(pad["count"]):
                    reading_id = f"null-{index:03d}"
                    while reading_id in ids:
                        reading_id = "_" + reading_id
                    ids.add(reading_id)
                    data["readings"].append(dict(id=reading_id, force=[0.0, 0.0],
                        acceleration=[0.0, 0.0], covariance=deepcopy(covariance),
                        metadata={"role": "exact-null padding"}))
                data["metadata"]["readout_study"] = dict(source=source.name, cell=identifier,
                    prior_center=center, null_padding=deepcopy(pad))
                changed = Case.from_dict(data)
                cases.append(changed)
                cells.append(dict(id=identifier, case=changed.name, source_case=source.name,
                    prior_center=center, null_padding=deepcopy(pad), readings=len(changed.inputs.readings)))
    return dict(kind="readout_candidates", raw=deepcopy(raw), cases=cases, cells=cells)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _save_cell(folder, case, cell, values, seed):
    folder.mkdir(parents=True)
    law = values["law"]
    arrays = dict(log_mass=np.asarray(law["log_mass"]),
                  mass=np.exp(law["log_mass"]), probability=np.asarray(law["probability"]),
                  S2=np.asarray(values["second_moment"]),
                  g=np.asarray(law["g"]), Q=np.asarray(law["Q"]))
    arrays.update({f"weights/{name}": np.asarray(weight) for name, weight in values["weights"].items()})
    np.savez_compressed(folder / "arrays.npz", **arrays)
    replay = dict(case.to_dict(include_truth=True), kind="inspection", seed=seed,
                  numerics={"mass_grid": arrays["mass"].tolist()}, checks=[])
    (folder / "case.yaml").write_text(yaml.safe_dump(replay, sort_keys=False), encoding="utf-8")
    row = dict(cell, law=_json_value({key: value for key, value in law.items()
                                    if key not in ("log_mass", "probability")}),
               candidates=_json_value(values["candidates"]),
               assumptions=list(values["assumptions"]),
               files=dict(case=f'{cell["id"]}/case.yaml', points=f'{cell["id"]}/points.json',
                          arrays=f'{cell["id"]}/arrays.npz'),
               arrays=dict(log_mass="log_mass", mass="mass", probability="probability", S2="S2",
                           size_weights={name: f"weights/{name}" for name in values["weights"]},
                           S2_definition="E[stacked acceleration length squared | m]"))
    (folder / "points.json").write_text(json.dumps(row, indent=2, allow_nan=False), encoding="utf-8")
    return row


def _point(candidate):
    return f'{candidate["point"]:.9g}' if candidate["status"] == "defined" else candidate["status"]


def _write_summary(out, record):
    lines = ["# O2 readout candidates", "", record["question"], "",
        "Every candidate in a row uses the same saved ours29 law. No random draws; one worker.",
        "The perpendicular case is a diagnostic, not a declaration of a Z-null input.", "",
        "| Source case | c | Added nulls (force SD, acceleration SD) | " + " | ".join(CANDIDATE_IDS) + " |",
        "|---|---:|---|" + "---:|" * len(CANDIDATE_IDS)]
    for row in record["cells"]:
        pad = row["null_padding"]
        padding = f'{pad["count"]} ({pad["force_sd"]:g}, {pad["acceleration_sd"]:g})'
        lines.append(f'| {row["source_case"]} | {row["prior_center"]:g} | {padding} | '
                     + " | ".join(_point(row["candidates"][name]) for name in CANDIDATE_IDS) + " |")
    lines.extend(["", "## Finite candidate-law composition checks", "",
        "Each check uses one propagated candidate law and a shared physical channel. It does not refit combined measured inputs.",
        "These checks concern the full finite-law domain; Gaussian-adapter reachability is not asserted.", "",
        "| Candidate | L5 left | L5 right | L14 inverse left | L14 inverse right |", "|---|---:|---:|---:|---:|"])
    for row in record["finite_law_checks"]["composition"]:
        lines.append(f'| {row["candidate"]} | {row["L5_left"]:.12g} | {row["L5_right"]:.12g} | '
                     f'{row["L14_inverse_left"]:.12g} | {row["L14_inverse_right"]:.12g} |')
    lines.extend(["", "The plain-sum-cloud row changes the series join; no composition verdict is inferred for it from these finite magnitude-law checks.",
        "The direction and cancellation records in finite_law_checks.json retain their separate domain statements.", ""])
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    fields = ["cell", "source_case", "prior_center", "null_count", "null_force_sd", "null_acceleration_sd",
              "readings", "candidate", "status", "point", "law_median", "eta"]
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in record["cells"]:
            pad = row["null_padding"]
            for name in CANDIDATE_IDS:
                writer.writerow(dict(cell=row["id"], source_case=row["source_case"],
                    prior_center=row["prior_center"], null_count=pad["count"], null_force_sd=pad["force_sd"],
                    null_acceleration_sd=pad["acceleration_sd"], readings=row["readings"], candidate=name,
                    **row["candidates"][name], law_median=row["law"]["median"], eta=row["law"]["eta"]))


def run(path, out_root="results", progress=print):
    """Evaluate the explicit cells serially, preserving native law arrays and inputs."""
    from .runner import code_version

    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    # Capture sources before evaluation, including the input/replay helpers.
    code = code_version()
    inspection = Path(__file__).parent / "inspection"
    code["inspection_file_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                      for p in sorted(inspection.glob("*.py"))}
    code["dependencies"] = {name: metadata.version(name) for name in ("numpy", "scipy", "PyYAML")}
    code["preset_sha256"] = hashlib.sha256(source).hexdigest()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    base, suffix = Path(out_root) / spec["raw"]["name"], 1
    out = base / stamp
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    (out / "source_preset.yaml").write_bytes(source)
    (out / "preset.yaml").write_text(yaml.safe_dump(spec["raw"], sort_keys=False), encoding="utf-8")
    start, rows = time.perf_counter(), []
    for case, cell in zip(spec["cases"], spec["cells"]):
        if progress:
            progress(f'cell {cell["id"]}: {cell["source_case"]}, c={cell["prior_center"]:g}, '
                     f'added nulls={cell["null_padding"]["count"]}')
        readings = [Reading(r.force, r.acceleration, r.covariance, r.subspace, r.window)
                    for r in case.inputs.readings]
        values = evaluate(readings, cell["prior_center"])
        rows.append(_save_cell(out / cell["id"], case, cell, values, spec["raw"]["seed"]))
    checks = dict(composition=counterexamples(), direction=direction_counterexample(),
                  cancellation=cancellation_counterexample())
    (out / "finite_law_checks.json").write_text(json.dumps(checks, indent=2, allow_nan=False), encoding="utf-8")
    code["effective_preset_sha256"] = hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest()
    record = dict(kind="readout_candidates", preset=spec["raw"]["name"], question=spec["raw"]["question"],
        seed=spec["raw"]["seed"], random_draws=0, deterministic=True, workers=1, code=code,
        started=stamp, seconds=time.perf_counter() - start, candidates=list(CANDIDATE_IDS), cells=rows,
        finite_law_checks=checks,
        assumptions=["The seed is recorded but unused: measured inputs and exact-null additions are deterministic.",
            "Within each cell, all readouts share one unchanged ours29 probability law and its declared prior.",
            "Native probability entries are normalized quadrature weights, not a density per unit mass or log mass.",
            "Every cell saves its explicit prior, covariance, subspace and measured vectors for replay.",
            "Different null-padding cells are separate declared measured inputs; their laws are recomputed.",
            "Finite-law composition checks propagate candidate forces or accelerations without recomputing a measurement law.",
            "No physical criterion or production readout is selected by this comparison."])
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    _write_summary(out, record)
    if progress:
        progress(f"wrote {out}")
    return out, record

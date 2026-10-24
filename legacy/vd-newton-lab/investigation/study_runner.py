"""Run a fixed-size, explicit simulation study with experiment-level checkpoints.

Validate: python study_runner.py --study studies/mass_excitation_pilot.json --validate
Run:      python study_runner.py --study studies/mass_excitation_pilot.json --output PATH
Resume:   repeat the run command with --resume (identical specification and code).

Completed experiments are immutable. A failed/interrupted experiment is restarted
in a new attempt directory; its partial files remain available for inspection.
Checkpoints do not yet resume inside an experiment. Precision calculations are
planning information, never instructions to stop or extend a running study.
"""

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re

import numpy as np

from experiment_config import (MAX_DIAGONAL_SAMPLES, MAX_DIAGONAL_TOTAL_TRIALS,
                               MAX_TOTAL_TRIALS, normalize_config, run_config,
                               scenario_rows, uses_measurement_engine)


HERE = Path(__file__).resolve().parent
SCHEMA_VERSION = 1
ADEMP_FIELDS = ("aims", "data_generating_mechanisms", "estimands", "methods", "performance_measures")
METRIC_RANGES = {("point", "capped_squared_log_factor_2"): math.log(2) ** 2,
                 ("point", "outside_factor_2"): 1., ("uncertainty", "coverage_95"): 1.}


def _object(value, fields, label):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError(f"{label} must contain exactly: {', '.join(fields)}")


def _strings(value, label):
    if not isinstance(value, list) or not value or not all(isinstance(v, str) and v.strip() for v in value):
        raise ValueError(f"{label} must be a nonempty list of nonempty strings")
    return [v.strip() for v in value]


def normalize_study(value):
    """Validate a study and normalize every embedded engine configuration."""
    fields = ("schema_version", "name", "phase", "ademp", "provenance", "interpretation",
              "fixed_replicates_per_cell", "primary_metrics", "experiments")
    _object(value, fields, "study")
    result = deepcopy(value)
    if type(result["schema_version"]) is not int or result["schema_version"] != SCHEMA_VERSION:
        raise ValueError("Unsupported study schema_version")
    if not isinstance(result["name"], str) or not result["name"].strip():
        raise ValueError("Study name must be nonempty")
    result["name"] = result["name"].strip()
    if result["phase"] not in ("exploratory_pilot", "main_assessment"):
        raise ValueError("Study phase must be exploratory_pilot or main_assessment")
    _object(result["ademp"], ADEMP_FIELDS, "ademp")
    result["ademp"] = {key: _strings(result["ademp"][key], f"ademp.{key}") for key in ADEMP_FIELDS}
    result["interpretation"] = _strings(result["interpretation"], "interpretation")
    provenance = result["provenance"]
    _object(provenance, ("managed_repository", "managed_commit", "managed_paths", "methodology"), "provenance")
    if not isinstance(provenance["managed_repository"], str) or not provenance["managed_repository"].strip():
        raise ValueError("provenance.managed_repository must be nonempty")
    if not isinstance(provenance["managed_commit"], str) or not re.fullmatch(r"[0-9a-f]{40}", provenance["managed_commit"]):
        raise ValueError("provenance.managed_commit must be a full Git SHA")
    provenance["managed_paths"] = _strings(provenance["managed_paths"], "provenance.managed_paths")
    if not isinstance(provenance["methodology"], list) or not provenance["methodology"]:
        raise ValueError("provenance.methodology must be a nonempty list")
    for source in provenance["methodology"]:
        _object(source, ("title", "url"), "methodology source")
        if not all(isinstance(source[k], str) and source[k].strip() for k in source) or not source["url"].startswith("https://"):
            raise ValueError("Methodology sources need nonempty titles and HTTPS URLs")
    count = result["fixed_replicates_per_cell"]
    if type(count) is not int or count < 2:
        raise ValueError("fixed_replicates_per_cell must be an integer of at least 2")
    metrics = result["primary_metrics"]
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("primary_metrics must be a nonempty list")
    seen_metrics = set()
    for metric in metrics:
        _object(metric, ("scope", "name", "target_mcse"), "primary metric")
        if not isinstance(metric["scope"], str) or not isinstance(metric["name"], str):
            raise ValueError("Metric scope and name must be strings")
        key = (metric["scope"], metric["name"])
        if key not in METRIC_RANGES or key in seen_metrics:
            raise ValueError("Primary metrics must be supported and distinct")
        seen_metrics.add(key)
        target = metric["target_mcse"]
        if isinstance(target, bool) or not isinstance(target, (int, float)) or not math.isfinite(target) or target <= 0:
            raise ValueError("target_mcse must be finite and positive")
        metric["target_mcse"] = float(target)
    experiments = result["experiments"]
    if not isinstance(experiments, list) or not experiments:
        raise ValueError("experiments must be a nonempty list")
    ids, seeds = set(), set()
    for experiment in experiments:
        _object(experiment, ("id", "config"), "experiment")
        identifier = experiment["id"]
        if not isinstance(identifier, str) or not re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", identifier) or identifier in ids:
            raise ValueError("Experiment ids must be unique safe lowercase names")
        ids.add(identifier)
        config = normalize_config(experiment["config"])
        if config["samples"] != count:
            raise ValueError("Every experiment samples count must equal fixed_replicates_per_cell")
        if not config["save_trials"]:
            raise ValueError("Study experiments must save_trials for audit and rescoring")
        experiment_seeds = [config["training_seed"], config["heldout_seed"]]
        calibration = config.get("measurement", {}).get("calibration", {})
        if calibration.get("mode") in ("estimated", "pooled_isotropic"):
            experiment_seeds.append(calibration["seed"])
        if seeds.intersection(experiment_seeds):
            raise ValueError("Training, held-out and calibration seeds must be distinct across experiments")
        seeds.update(experiment_seeds)
        experiment["config"] = config
    return result


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_fingerprint():
    """Lock execution dependencies; independently versioned analysis can evolve.

    The comparison engine already declares its inference/configuration/scoring
    modules. Include this orchestrator and the production estimator imported
    by the legacy vector-check dependency. Plotting and results browsing do
    not execute in a study and should not prevent resumption.
    """
    from compare_estimators import ENGINE_SOURCES
    paths = {name: HERE / name for name in (*ENGINE_SOURCES, "study_runner.py")}
    paths["../estimator.py"] = HERE.parent / "estimator.py"
    return {"python": platform.python_version(), "numpy": np.__version__,
            "source_sha256": {name: _hash(path) for name, path in sorted(paths.items())}}


def _write_json(path, value):
    # Replace only our own checkpoint/summary; completed experiment files are
    # never rewritten. Atomic replacement prevents a half-written checkpoint.
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _now():
    return datetime.now(timezone.utc).isoformat()


def _precision_entry(metric, mean, mcse, samples, paired=False):
    bounded_range = METRIC_RANGES[(metric["scope"], metric["name"])] * (2 if paired else 1)
    target = metric["target_mcse"]
    available = (isinstance(mcse, (int, float)) and not isinstance(mcse, bool)
                 and math.isfinite(mcse) and mcse >= 0)
    return {"metric": metric["name"], "mean": mean, "mcse": mcse if available else None,
            "target_mcse": target, "target_met": mcse <= target if available else None,
            "status": "available" if available else "mcse_unavailable",
            "planning_status": ("zero_pilot_variance_use_conservative_bound" if available and mcse == 0
                                else "pilot_variance_estimate" if available else "mcse_unavailable"),
            "pilot_variance_suggested_replicates": max(2, math.ceil(samples * (mcse / target) ** 2)) if available and mcse > 0 else None,
            "conservative_bounded_variance_replicates": max(2, math.ceil((bounded_range / (2 * target)) ** 2))}


def precision_summary(report, metrics):
    """Extract existing MCSEs; never invent standard errors for missing metrics.

    The estimated future R uses the pilot sample variance. It can be optimistic,
    especially for zero observed failures. A separate bounded-variance value is
    supplied for conservative planning; neither changes this study's fixed R.
    """
    entries = []
    samples = report["samples_per_scenario"]
    for row in report["rows"]:
        for metric in metrics:
            if metric["scope"] == "point":
                for method, summary in row["methods"].items():
                    entries.append({"scenario": row["id"], "scope": "point", "method": method,
                                    **_precision_entry(metric, summary.get("means", {}).get(metric["name"]),
                                                       summary.get("mcse", {}).get(metric["name"]), samples)})
                for comparator, comparisons in row.get("paired_flat_minus_other", {}).items():
                    comparison = comparisons.get(metric["name"], {})
                    entries.append({"scenario": row["id"], "scope": "paired_point",
                                    "method": "flat_joint_minus_" + comparator,
                                    **_precision_entry(metric, comparison.get("mean_difference"),
                                                       comparison.get("mcse"), samples, paired=True)})
            else:
                for method, law in row.get("uncertainty", {}).items():
                    coverage = law.get("intervals", {}).get("95", {}).get("coverage", {})
                    entries.append({"scenario": row["id"], "scope": "uncertainty", "method": method,
                                    **_precision_entry(metric, coverage.get("mean"), coverage.get("empirical_mcse"), samples)})
                for method, comparisons in row.get("paired_uncertainty_flat_minus_other", {}).items():
                    comparison = comparisons.get(metric["name"], {})
                    entries.append({"scenario": row["id"], "scope": "paired_uncertainty",
                                    "method": "flat_joint_minus_" + method,
                                    **_precision_entry(metric, comparison.get("mean_difference"),
                                                       comparison.get("mcse"), samples, paired=True)})
    config = report.get("experiment_config", {})
    diagonal = config.get("noise_model") == "diagonal_gaussian_3d" or uses_measurement_engine(config)
    cell_count = len(report["rows"])
    capacity = min(MAX_DIAGONAL_SAMPLES, MAX_DIAGONAL_TOTAL_TRIALS // cell_count) if diagonal else min(1_048_576, MAX_TOTAL_TRIALS // cell_count)
    for entry in entries:
        for key in ("pilot_variance_suggested_replicates", "conservative_bounded_variance_replicates"):
            requirement = entry[key]
            entry[key + "_within_current_capacity"] = requirement <= capacity if requirement is not None else None
    return {"samples_per_cell": samples, "current_maximum_replicates_per_cell": capacity,
            "interpretation": "Fixed-R assessment. MCSE targets concern Monte Carlo precision, not acceptable estimator error. Future R is variance-based planning advice, not optional stopping. Zero observed variance does not establish zero population variance. Cells share noise within an experiment and must not be pooled as independent trials. Paired targets describe flat-minus-comparator differences on the same observations. No global winner is selected.",
            "entries": entries}


def _verify_completed(output, record):
    attempt = output / record["attempt_directory"]
    for filename, expected in record["output_sha256"].items():
        path = attempt / filename
        if not path.is_file() or _hash(path) != expected:
            raise ValueError(f"Completed experiment output changed or missing: {path}")
    return json.loads((attempt / "results.json").read_text(encoding="utf-8"))


def _check_report(report, config):
    if report.get("samples_per_scenario") != config["samples"]:
        raise ValueError("Experiment returned an unexpected replicate count")
    expected = [row["id"] for row in scenario_rows(config)]
    if [row["id"] for row in report.get("rows", [])] != expected:
        raise ValueError("Experiment report is incomplete or has unexpected scenario rows")
    if report.get("experiment_config") != config:
        raise ValueError("Experiment report configuration differs from the study")


def run_study(specification, output, resume=False, progress=None):
    """Run or resume a study; failures remain explicit and are propagated."""
    spec = normalize_study(specification)
    output = Path(output)
    fingerprint = source_fingerprint()
    spec_hash = hashlib.sha256(_canonical(spec)).hexdigest()
    checkpoint = output / "study_checkpoint.json"
    if resume:
        if not checkpoint.is_file():
            raise ValueError("--resume requires an existing study_checkpoint.json")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("schema_version") != SCHEMA_VERSION or state.get("spec_sha256") != spec_hash or state.get("specification") != spec:
            raise ValueError("Resume requires the identical normalized study specification")
        if state.get("source_fingerprint") != fingerprint:
            raise ValueError("Resume requires identical source hashes and Python/NumPy versions")
        if [record["id"] for record in state.get("experiments", [])] != [e["id"] for e in spec["experiments"]]:
            raise ValueError("Checkpoint experiment inventory differs from the study")
    else:
        if output.exists():
            raise ValueError("Study output already exists; choose a new directory or use --resume")
        output.mkdir(parents=True)
        state = {"schema_version": SCHEMA_VERSION, "spec_sha256": spec_hash, "specification": spec,
                 "source_fingerprint": fingerprint, "status": "running", "started_utc": _now(),
                 "experiments": [{"id": e["id"], "status": "pending", "attempts": []} for e in spec["experiments"]]}
        _write_json(output / "study.json", spec)
        _write_json(checkpoint, state)
    # Validate every completed result before doing any new work.
    reports = {}
    for experiment, record in zip(spec["experiments"], state["experiments"]):
        if record["status"] == "completed":
            report = _verify_completed(output, record)
            _check_report(report, experiment["config"])
            reports[record["id"]] = report
    if state["status"] == "completed":
        return state
    state["status"] = "running"
    _write_json(checkpoint, state)
    for experiment, record in zip(spec["experiments"], state["experiments"]):
        if record["status"] == "completed":
            continue
        if record["status"] == "running" and record["attempts"]:
            record["attempts"][-1]["status"] = "interrupted"
        number = len(record["attempts"]) + 1
        relative = Path("experiments") / record["id"] / f"attempt-{number:03d}"
        destination = output / relative
        if destination.exists():
            raise ValueError(f"Refusing to overwrite existing experiment attempt: {destination}")
        attempt = {"directory": relative.as_posix(), "status": "running", "started_utc": _now()}
        record["attempts"].append(attempt)
        record["status"] = "running"
        _write_json(checkpoint, state)
        try:
            def on_progress(event):
                if progress:
                    progress({"experiment": record["id"], **event})
            report = run_config(experiment["config"], destination, progress=on_progress)
            _check_report(report, experiment["config"])
            if source_fingerprint() != fingerprint:
                raise ValueError("Source files changed during the experiment; start a new study under stable code")
            result_path = destination / "results.json"
            if json.loads(result_path.read_text(encoding="utf-8")) != report:
                raise ValueError("Saved results differ from the returned experiment report")
            hashes = {p.relative_to(destination).as_posix(): _hash(p)
                      for p in sorted(destination.rglob("*")) if p.is_file()}
            record.update(status="completed", attempt_directory=relative.as_posix(), output_sha256=hashes)
            attempt.update(status="completed", finished_utc=_now())
            reports[record["id"]] = report
        except BaseException as error:
            record["status"] = "failed" if isinstance(error, Exception) else "interrupted"
            attempt.update(status=record["status"], finished_utc=_now(),
                           error={"type": type(error).__name__, "message": str(error)})
            state["status"] = record["status"]
            _write_json(checkpoint, state)
            raise
        _write_json(checkpoint, state)
        _write_json(output / "precision_summary.json", {
            "status": "partial", "completed_experiments": list(reports),
            "experiments": {key: precision_summary(value, spec["primary_metrics"]) for key, value in reports.items()}})
    state.update(status="completed", finished_utc=_now())
    _write_json(output / "precision_summary.json", {
        "status": "completed", "completed_experiments": list(reports),
        "experiments": {key: precision_summary(value, spec["primary_metrics"]) for key, value in reports.items()}})
    _write_json(checkpoint, state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.validate and (args.output is not None or args.resume):
        parser.error("--validate cannot be combined with --output or --resume")
    if not args.validate and args.output is None:
        parser.error("--output is required unless --validate is used")
    try:
        spec = normalize_study(json.loads(args.study.read_text(encoding="utf-8-sig")))
        if args.validate:
            print(json.dumps({"status": "valid", "name": spec["name"],
                              "fixed_replicates_per_cell": spec["fixed_replicates_per_cell"],
                              "experiments": [{"id": e["id"], "cells": len(scenario_rows(e["config"]))} for e in spec["experiments"]]}))
        else:
            state = run_study(spec, args.output, args.resume)
            print(json.dumps({"status": state["status"], "output": str(args.output.resolve())}))
    except (ValueError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

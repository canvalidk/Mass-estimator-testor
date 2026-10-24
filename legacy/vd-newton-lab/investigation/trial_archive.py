"""Versioned simulation trials and point rescoring without inference.

The runner owns ``trial_manifest.json`` and its run-level provenance. Each
entry returned by ``write_trial_archive`` describes one immutable NPZ file.
Saved posterior arrays are summaries, not complete posterior representations;
this module only rescores point estimates. Replaying observations through a
changed estimator is a separate experiment.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import zipfile

import numpy as np

from comparison_scores import paired_summary, point_losses, score_point


ARCHIVE_SCHEMA_VERSION = 1
MANIFEST_NAME = "trial_manifest.json"
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_VECTOR_NAMES = {
    "observation__force", "observation__acceleration",
    "heldout__force", "heldout__acceleration",
    "heldout__true_force", "heldout__true_acceleration",
    "noise__true_force_sd", "noise__true_acceleration_sd",
    "noise__supplied_force_sd", "noise__supplied_acceleration_sd",
    "noise__supplied_mean_force_sd", "noise__supplied_mean_acceleration_sd",
}
_SUMMARY_FIELDS = {
    "methods", "paired_flat_minus_other", "uncertainty",
    "paired_uncertainty_flat_minus_tube", "refinement",
}
_ENTRY_FIELDS = {
    "schema_version", "id", "path", "sha256", "samples", "scenario",
    "arrays", "point_methods", "method_risk",
}
_RISK_FIELDS = {"population_unbounded_risk", "population_risk_note"}


def _method_risk(row, methods):
    """Preserve established analytic limitations separately from sample scores."""
    result = {}
    for name in methods:
        source = row.get("methods", {}).get(name, {})
        if not isinstance(source, dict):
            raise ValueError("Archived point method summaries must be mappings")
        risk = {key: source[key] for key in _RISK_FIELDS if key in source}
        source_paired = row.get("paired_flat_minus_other", {}).get(name, {})
        notes = {metric: summary["mcse_interpretation"] for metric, summary in source_paired.items()
                 if isinstance(summary, dict) and "mcse_interpretation" in summary}
        if notes:
            risk["source_pair_reference"] = "flat_joint"
            risk["source_paired_mcse_notes"] = notes
        if risk:
            result[name] = risk
    return _json_value(result)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value):
    """Round-trip metadata, rejecting non-JSON or nonfinite numeric values."""
    try:
        return json.loads(json.dumps(value, allow_nan=False))
    except (TypeError, ValueError) as error:
        raise ValueError("Archive metadata must be finite JSON data") from error


def _positive_truth(scenario):
    value = scenario.get("true_mass") if isinstance(scenario, dict) else None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("Archive scenario requires a finite positive true_mass")
    try:
        valid = math.isfinite(value) and value > 0
    except OverflowError:
        valid = False
    if not valid:
        raise ValueError("Archive scenario requires a finite positive true_mass")


def _contained_path(run_dir, relative):
    if not isinstance(relative, str) or "\\" in relative:
        raise ValueError("Archive path must be a relative POSIX path")
    part = Path(relative)
    if part.is_absolute() or ".." in part.parts or ":" in relative:
        raise ValueError("Archive path must remain inside the run directory")
    root = Path(run_dir).resolve()
    resolved = (root / part).resolve()
    if not resolved.is_relative_to(root) or resolved == root:
        raise ValueError("Archive path must remain inside the run directory")
    return resolved


def _array_kind(name):
    if not isinstance(name, str):
        raise ValueError("Archive array names must be strings")
    if name == "replicate_id":
        return "replicate"
    if name in _VECTOR_NAMES:
        return "vector"
    parts = name.split("__")
    if len(parts) == 2 and parts[0] == "point" and _TOKEN.fullmatch(parts[1]):
        return "point"
    if (len(parts) == 3 and parts[0] == "posterior"
            and all(_TOKEN.fullmatch(part) for part in parts[1:])):
        return "posterior_summary"
    raise ValueError(f"Unsupported archive array name: {name!r}")


def _validate_arrays(arrays, samples):
    required = {"replicate_id", "observation__force", "observation__acceleration"}
    if not required.issubset(arrays):
        raise ValueError("Archive requires replicate IDs and both observation vectors")
    if not any(isinstance(name, str) and name.startswith("point__") for name in arrays):
        raise ValueError("Archive requires at least one point-estimate array")
    for name, values in arrays.items():
        kind = _array_kind(name)
        if not isinstance(values, np.ndarray) or values.dtype.kind not in "fiu":
            raise ValueError(f"Archive array {name!r} must have a real numeric dtype")
        expected = (samples, 3) if kind == "vector" else (samples,)
        if values.shape != expected:
            raise ValueError(f"Archive array {name!r} shape must be {expected}")
        if kind == "replicate":
            if values.dtype.kind not in "iu" or not np.array_equal(values, np.arange(samples)):
                raise ValueError("replicate_id must be the ordered integer sequence 0..samples-1")
        elif kind == "vector":
            if not np.isfinite(values).all():
                raise ValueError(f"Observation, held-out and noise arrays must be finite: {name}")
            if name.startswith("noise__") and not (values > 0).all():
                raise ValueError(f"Noise SDs must be positive: {name}")
        # Nonfinite point estimates and posterior summaries are retained, not
        # silently removed or converted into missing observations.


def write_trial_archive(output_dir, row, arrays):
    """Save one scenario, add ordered replicate IDs, and return its manifest entry.

    ``row`` contains an ID and positive true mass plus scenario metadata. Point
    estimates are ``point__method`` vectors; observations and supplied noise are
    physical vectors. Existing archives are never overwritten.
    """
    if not isinstance(row, dict) or not isinstance(arrays, dict):
        raise ValueError("row and arrays must be mappings")
    identifier = row.get("id")
    if not isinstance(identifier, str) or not _TOKEN.fullmatch(identifier):
        raise ValueError("Scenario id must be a safe nonempty filename token")
    scenario = _json_value({key: value for key, value in row.items() if key not in _SUMMARY_FIELDS})
    _positive_truth(scenario)
    converted = {name: np.asarray(values) for name, values in arrays.items()}
    force = converted.get("observation__force")
    if force is None or force.ndim != 2 or force.shape[0] < 1:
        raise ValueError("Archive requires nonempty observation__force vectors")
    samples = force.shape[0]
    if "replicate_id" not in converted:
        converted["replicate_id"] = np.arange(samples, dtype=np.int64)
    _validate_arrays(converted, samples)
    methods = sorted(name[7:] for name in converted if name.startswith("point__"))
    method_risk = _method_risk(row, methods)
    relative = f"trials/{identifier}.npz"
    target = _contained_path(output_dir, relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("xb") as handle:
        np.savez_compressed(handle, **converted)
    return {
        "schema_version": ARCHIVE_SCHEMA_VERSION,
        "id": identifier, "path": relative, "sha256": _sha256(target),
        "samples": samples, "scenario": scenario,
        "point_methods": methods, "method_risk": method_risk,
        "arrays": {name: {"shape": list(values.shape), "dtype": values.dtype.str}
                   for name, values in sorted(converted.items())},
    }


def _validate_entry(entry):
    if not isinstance(entry, dict) or set(entry) != _ENTRY_FIELDS:
        raise ValueError("Invalid trial archive manifest entry fields")
    if type(entry["schema_version"]) is not int or entry["schema_version"] != ARCHIVE_SCHEMA_VERSION:
        raise ValueError("Unsupported trial archive schema_version")
    identifier = entry["id"]
    if not isinstance(identifier, str) or not _TOKEN.fullmatch(identifier):
        raise ValueError("Invalid trial archive scenario id")
    if entry["path"] != f"trials/{identifier}.npz":
        raise ValueError("Archive path must be trials/{id}.npz inside the run directory")
    if not isinstance(entry["sha256"], str) or re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) is None:
        raise ValueError("Invalid archive SHA-256")
    if type(entry["samples"]) is not int or entry["samples"] < 1:
        raise ValueError("Archive samples must be a positive integer")
    _positive_truth(entry["scenario"])
    _json_value(entry["scenario"])
    if entry["scenario"].get("id") != identifier:
        raise ValueError("Scenario and archive IDs disagree")
    descriptors = entry["arrays"]
    if not isinstance(descriptors, dict) or not descriptors:
        raise ValueError("Archive arrays require shape and dtype descriptors")
    for name, descriptor in descriptors.items():
        kind = _array_kind(name)
        expected = [entry["samples"], 3] if kind == "vector" else [entry["samples"]]
        if (not isinstance(descriptor, dict) or set(descriptor) != {"shape", "dtype"}
                or descriptor["shape"] != expected
                or any(type(size) is not int for size in descriptor["shape"])):
            raise ValueError(f"Invalid archive descriptor shape: {name}")
        try:
            dtype = np.dtype(descriptor["dtype"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid archive dtype: {name}") from error
        if dtype.kind not in "fiu" or dtype.str != descriptor["dtype"]:
            raise ValueError(f"Unsupported archive dtype: {name}")
    methods = sorted(name[7:] for name in descriptors if name.startswith("point__"))
    if not methods or entry["point_methods"] != methods:
        raise ValueError("Archive point_methods must match the point arrays")
    risks = entry["method_risk"]
    if not isinstance(risks, dict) or set(risks) - set(methods):
        raise ValueError("Archive method_risk must refer to archived point methods")
    for name, risk in risks.items():
        if (not isinstance(risk, dict)
                or set(risk) - (_RISK_FIELDS | {"source_pair_reference", "source_paired_mcse_notes"})
                or any(not isinstance(risk[key], str) for key in _RISK_FIELDS if key in risk)):
            raise ValueError(f"Invalid method risk metadata: {name}")
        if "source_paired_mcse_notes" in risk:
            notes = risk["source_paired_mcse_notes"]
            if (risk.get("source_pair_reference") not in methods or not isinstance(notes, dict)
                    or any(not isinstance(metric, str) or not isinstance(note, str)
                           for metric, note in notes.items())):
                raise ValueError(f"Invalid source paired MCSE notes: {name}")
        elif "source_pair_reference" in risk:
            raise ValueError("source_pair_reference requires source_paired_mcse_notes")


def load_trial_manifest(run_dir):
    """Read a versioned run manifest; old aggregate-only runs cannot be rescored."""
    manifest_path = _contained_path(run_dir, MANIFEST_NAME)
    if not manifest_path.is_file():
        raise ValueError("This run has no versioned trial_manifest.json; aggregate-only or legacy runs cannot be rescored. Generate a new run with save_trials enabled.")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ValueError("Cannot read trial_manifest.json") from error
    if (not isinstance(manifest, dict)
            or set(manifest) - {"schema_version", "entries", "provenance"}
            or type(manifest.get("schema_version")) is not int
            or manifest["schema_version"] != ARCHIVE_SCHEMA_VERSION):
        raise ValueError("Unsupported or invalid trial manifest schema_version")
    entries = manifest.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Trial manifest must contain a nonempty entries list")
    if "provenance" in manifest and not isinstance(manifest["provenance"], dict):
        raise ValueError("Trial manifest provenance must be a mapping")
    _json_value(manifest)
    identifiers = set()
    for entry in entries:
        _validate_entry(entry)
        _contained_path(run_dir, entry["path"])
        # Case-folding also protects Windows filename identity.
        identifier = entry["id"].casefold()
        if identifier in identifiers:
            raise ValueError("Duplicate scenario in trial manifest")
        identifiers.add(identifier)
    return manifest


def load_trial_archive(run_dir, entry):
    """Verify hash and schema before loading arrays with pickle disabled."""
    _validate_entry(entry)
    target = _contained_path(run_dir, entry["path"])
    if not target.is_file():
        raise ValueError(f"Missing trial archive: {entry['path']}")
    if _sha256(target) != entry["sha256"]:
        raise ValueError(f"Trial archive SHA-256 mismatch: {entry['path']}")
    try:
        with np.load(target, allow_pickle=False) as saved:
            if len(saved.files) != len(set(saved.files)) or set(saved.files) != set(entry["arrays"]):
                raise ValueError("Archive array names differ from manifest")
            arrays = {name: saved[name] for name in saved.files}
    except (OSError, TypeError, ValueError, zipfile.BadZipFile) as error:
        raise ValueError(f"Cannot safely load trial archive {entry['path']}: {error}") from error
    _validate_arrays(arrays, entry["samples"])
    for name, values in arrays.items():
        descriptor = entry["arrays"][name]
        if list(values.shape) != descriptor["shape"] or values.dtype.str != descriptor["dtype"]:
            raise ValueError(f"Archive array differs from its descriptor: {name}")
    return arrays


def _factor(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 1:
        raise ValueError(f"{name} must be a finite number greater than one")
    return float(value)


def _custom_losses(points, truth, factor_tolerance, cap_factor):
    valid = np.isfinite(points) & (points > 0)
    with np.errstate(over="ignore", invalid="ignore"):
        success = valid & (points >= truth / factor_tolerance) & (points <= truth * factor_tolerance)
    squared_log = point_losses(points, truth)["squared_log"]
    return {"outside_factor": (~success).astype(float),
            "capped_squared_log": np.minimum(squared_log, math.log(cap_factor) ** 2)}


def _annotate_paired_risk(summary, risks, reference, other, metric):
    relevant = {name: risks[name] for name in (reference, other) if name in risks}
    if relevant:
        summary["source_method_risk"] = relevant
    contextual = {}
    for name, risk in relevant.items():
        note = risk.get("source_paired_mcse_notes", {}).get(metric)
        if note is None:
            continue
        original_reference = risk["source_pair_reference"]
        if {reference, other} == {name, original_reference}:
            # Reversing the same pair changes the sign, not its variance.
            summary["mcse_interpretation"] = note
        else:
            contextual[name] = {"original_reference": original_reference, "note": note}
    if contextual:
        summary["source_pair_mcse_notes"] = contextual
        summary.setdefault("mcse_interpretation", "The source run established population-moment limitations for an input method. This changed pairing has not been shown to have finite population variance; source-pair warnings are retained as context.")
    return summary


def rescore_run(run_dir, output, *, factor_tolerance=2.0, cap_factor=5.0,
                reference_method="flat_joint"):
    """Recompute original point scores and configurable bounded paired losses.

    Results explicitly exclude prediction and distribution rescoring. Source
    files are untouched; ``output`` must be a new JSON file. Monte Carlo SEs
    describe finite samples, without establishing population moment existence.
    """
    factor_tolerance = _factor(factor_tolerance, "factor_tolerance")
    cap_factor = _factor(cap_factor, "cap_factor")
    if not isinstance(reference_method, str) or not _TOKEN.fullmatch(reference_method):
        raise ValueError("reference_method must be a method name")
    run_dir = Path(run_dir).resolve()
    destination = Path(output).resolve()
    if destination.suffix.lower() != ".json":
        raise ValueError("Rescore output must be a new .json file")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {destination}")
    if destination.parent == run_dir and destination.name.lower() in {
            "results.json", "config.json", "experiment_config.json", MANIFEST_NAME}:
        raise ValueError("Rescore output cannot be a run results, config, or manifest file")
    manifest = load_trial_manifest(run_dir)
    rows = []
    for entry in manifest["entries"]:
        if reference_method not in entry["point_methods"]:
            raise ValueError(f"Reference method {reference_method!r} not archived in {entry['id']}")
        arrays = load_trial_archive(run_dir, entry)
        truth = entry["scenario"]["true_mass"]
        methods, comparisons = {}, {}
        reference_points = arrays[f"point__{reference_method}"]
        reference_losses = point_losses(reference_points, truth)
        reference_custom = _custom_losses(reference_points, truth, factor_tolerance, cap_factor)
        for name in entry["point_methods"]:
            points = arrays[f"point__{name}"]
            methods[name] = {"original_point_scores": score_point(points, truth)}
            risk = entry["method_risk"].get(name, {})
            methods[name]["original_point_scores"].update({key: risk[key] for key in _RISK_FIELDS if key in risk})
            if risk:
                methods[name]["source_risk_metadata"] = risk
            original_losses = reference_losses if name == reference_method else point_losses(points, truth)
            configurable = (reference_custom if name == reference_method
                            else _custom_losses(points, truth, factor_tolerance, cap_factor))
            methods[name]["configurable_scores"] = {}
            for metric, values in configurable.items():
                summary = paired_summary(values, np.zeros_like(values))
                methods[name]["configurable_scores"][metric] = {
                    "n": summary["n"], "mean": summary["mean_difference"],
                    "mcse": summary["mcse"], "status": summary["status"],
                }
            if name == reference_method:
                continue
            comparisons[name] = {
                "original_point_losses": {
                    metric: _annotate_paired_risk(
                        paired_summary(reference_losses[metric], values),
                        entry["method_risk"], reference_method, name, metric)
                    for metric, values in original_losses.items()},
                "configurable_losses": {
                    metric: paired_summary(reference_custom[metric], values)
                    for metric, values in configurable.items()},
            }
        rows.append({"id": entry["id"], "scenario": entry["scenario"],
                     "archive_sha256": entry["sha256"], "samples": entry["samples"],
                     "methods": methods, "paired_reference_minus_other": comparisons})
    report = {
        "schema_version": 1, "kind": "saved_trial_point_rescore",
        "source_run": str(run_dir), "source_manifest_sha256": _sha256(run_dir / MANIFEST_NAME),
        "source_provenance": manifest.get("provenance", {}),
        "scoring_source_sha256": {path.name: _sha256(path) for path in
                                  (Path(__file__), Path(__file__).with_name("comparison_scores.py"))},
        "numpy": np.__version__, "reference_method": reference_method,
        "factor_tolerance": factor_tolerance, "cap_factor": cap_factor,
        "scope": "Point scores from saved estimates only; no fitting, posterior reconstruction, prediction rescoring, or pooling across scenarios. Posterior arrays are saved summaries, not complete laws.",
        "failure_policy": "Nonpositive or nonfinite masses fail tolerance and receive maximum capped loss; original unbounded summaries are unavailable when any trial fails. No trials are discarded.",
        "mcse_scope": "Independent trials within each scenario; paired differences use identical replicate IDs. Finite empirical MCSEs do not establish finite population variance. No normal confidence intervals are inferred.",
        "paired_difference_convention": "reference minus comparator; negative favors reference",
        "rows": rows,
    }
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as handle:
        handle.write(encoded)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--factor-tolerance", default=2.0, type=float)
    parser.add_argument("--cap-factor", default=5.0, type=float)
    parser.add_argument("--reference-method", default="flat_joint")
    args = parser.parse_args()
    try:
        report = rescore_run(args.run, args.output, factor_tolerance=args.factor_tolerance,
                             cap_factor=args.cap_factor, reference_method=args.reference_method)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Rescoring failed: {error}\n")
    print(json.dumps({"output": str(args.output.resolve()), "scenarios": len(report["rows"]),
                      "scope": report["scope"]}))


if __name__ == "__main__":
    main()

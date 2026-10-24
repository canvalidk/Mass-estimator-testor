"""Rescore a completed study's zero-excitation trials at different fixed masses.

This is valid only because the declared Gaussian observation noise and supplied
calibration are independent of mass when both latent vectors are zero. No data
or estimates are regenerated. Truth-grid rows reuse the identical observations
and must never be pooled as independent simulated experiments.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np

from trial_archive import _contained_path, load_trial_archive, load_trial_manifest


DEFAULT_MASSES = [1 / 64, 1 / 16, .25, 1., 4., 16., 64.]
LEVELS = (50, 95)
LAWS = ("flat_joint", "tube_joint")
LAW_LABELS = {"flat_joint": "Flat", "tube_joint": "Tube",
              "calibrated_joint": "Calibration integrated", "oracle_joint": "Known-noise oracle"}
LAW_COLORS = {"flat_joint": "#17668a", "tube_joint": "#b84e25",
              "calibrated_joint": "#148568", "oracle_joint": "#8654a3"}
LAW_MARKERS = {"flat_joint": "o", "tube_joint": "s", "calibrated_joint": "D", "oracle_joint": "^"}


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _masses(values):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError("The truth grid must contain at least one positive mass")
    result = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("Truth masses must be finite and positive")
        try:
            number = float(value)
        except OverflowError as error:
            raise ValueError("Truth masses must be finite and positive") from error
        if not math.isfinite(number) or number <= 0:
            raise ValueError("Truth masses must be finite and positive")
        result.append(number)
    if result != sorted(set(result)):
        raise ValueError("Truth masses must be distinct and in increasing order")
    return result


def _mean(values):
    values = np.asarray(values, float)
    if not np.all(np.isfinite(values)):
        return {"mean": None, "mcse": None, "status": "nonfinite_values_present"}
    return {"mean": float(np.mean(values)),
            "mcse": float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else None,
            "status": "finite_sample"}


def _archived_laws(arrays, count):
    """Retain the legacy laws and validate every additional archived law."""
    discovered = set(LAWS)
    for name in arrays:
        if name.startswith("posterior__"):
            parts = name.split("__")
            if len(parts) != 3 or not parts[1] or not parts[2]:
                raise ValueError(f"Invalid archived posterior array name: {name}")
            discovered.add(parts[1])
    laws = tuple(law for law in LAW_LABELS if law in discovered) + tuple(sorted(discovered - LAW_LABELS.keys()))
    for law in laws:
        for level in LEVELS:
            for side in ("lower", "upper"):
                key = f"posterior__{law}__log_{side}_{level}"
                if key not in arrays or np.shape(arrays[key]) != (count,):
                    raise ValueError(f"Null diagnostic requires archived interval array {key}")
    return laws


def score_null_trials(arrays, masses):
    """Calculate all-trial bounded point and interval metrics without refitting."""
    masses = _masses(masses)
    points = arrays.get("point__flat_joint")
    if points is None or np.ndim(points) != 1 or not len(points):
        raise ValueError("Null diagnostic requires nonempty flat_joint point estimates")
    count = len(points)
    archived_laws = _archived_laws(arrays, count)
    valid = np.isfinite(points) & (points > 0)
    rows = []
    for mass in masses:
        truth = math.log(mass)
        logerror = np.full(count, np.inf)
        logerror[valid] = np.log(points[valid]) - truth
        with np.errstate(over="ignore"):
            caploss = np.minimum(logerror ** 2, math.log(2) ** 2)
        # Inclusive physical endpoints match the main runner exactly. If a
        # boundary over/underflows, all finite positive points remain on the
        # appropriate side; subtracting two logs can misclassify exact limits.
        with np.errstate(over="ignore", under="ignore"):
            factor_failure = (~valid) | (points < mass / 2) | (points > mass * 2)
        laws = {}
        for law in archived_laws:
            intervals = {}
            for level in LEVELS:
                lower = arrays[f"posterior__{law}__log_lower_{level}"]
                upper = arrays[f"posterior__{law}__log_upper_{level}"]
                good = np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
                width = np.full(count, np.inf)
                width[good] = upper[good] - lower[good]
                intervals[str(level)] = {
                    "nominal_coverage": level / 100,
                    "coverage": _mean(good & (lower <= truth) & (truth <= upper)),
                    "truth_below_interval": _mean(good & (truth < lower)),
                    "truth_above_interval": _mean(good & (truth > upper)),
                    "invalid_interval": _mean(~good), "log_width": _mean(width)}
            laws[law] = intervals
        rows.append({"true_mass": mass, "samples": count,
                     "flat_point": {"capped_squared_log_factor_2": _mean(caploss),
                                    "outside_factor_2": _mean(factor_failure),
                                    "invalid_estimate": _mean(~valid)},
                     "intervals": laws})
    return rows


def _completed_experiments(study_directory):
    """Verify recorded bytes; historical study code need not match today's code."""
    root = Path(study_directory).resolve()
    checkpoint_path = root / "study_checkpoint.json"
    state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if state.get("schema_version") != 1 or state.get("status") != "completed":
        raise ValueError("Null sweep requires a completed version-1 study checkpoint")
    canonical = json.dumps(state["specification"], sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    if hashlib.sha256(canonical).hexdigest() != state.get("spec_sha256"):
        raise ValueError("Study specification hash mismatch")
    specifications = state["specification"]["experiments"]
    records = state.get("experiments", [])
    if not records or [r["id"] for r in records] != [e["id"] for e in specifications]:
        raise ValueError("Study experiment inventory mismatch")
    verified = []
    for specification, record in zip(specifications, records):
        if record.get("status") != "completed":
            raise ValueError("Every study experiment must be completed")
        folder = _contained_path(root, record["attempt_directory"])
        hashes = record.get("output_sha256", {})
        if "results.json" not in hashes or "trial_manifest.json" not in hashes:
            raise ValueError("Checkpoint must pin results and the trial manifest")
        for relative, expected in hashes.items():
            path = _contained_path(folder, relative)
            if not path.is_file() or _hash(path) != expected:
                raise ValueError(f"Completed output hash mismatch: {relative}")
        report = json.loads((folder / "results.json").read_text(encoding="utf-8"))
        if report.get("experiment_config") != specification["config"]:
            raise ValueError("Completed report configuration differs from the study")
        manifest = load_trial_manifest(folder)
        row_ids = [row["id"] for row in report["rows"]]
        if [e["id"] for e in manifest["entries"]] != row_ids:
            raise ValueError("Trial manifest and completed result inventory differ")
        for entry in manifest["entries"]:
            if hashes.get(entry["path"]) != entry["sha256"]:
                raise ValueError("Checkpoint and trial archive hashes disagree")
        verified.append((record["id"], folder, report, manifest))
    return state, verified


def _null_entries(report, manifest):
    config = report["experiment_config"]
    if config.get("noise_model") not in ("isotropic_unit_3d", "diagonal_gaussian_3d"):
        raise ValueError("Null sweep supports only the declared mass-independent Gaussian noise models")
    mode = config.get("measurement", {}).get("calibration", {}).get("mode", "known")
    if mode not in ("known", "scaled", "estimated", "pooled_isotropic"):
        raise ValueError("Unsupported calibration mechanism for the null diagnostic")
    selected = []
    for row, entry in zip(report["rows"], manifest["entries"]):
        if row.get("true_acceleration_magnitude") == 0 and row.get("true_force_magnitude") == 0:
            for key in ("true_acceleration_magnitude", "true_force_magnitude", "true_mass"):
                if entry["scenario"].get(key) != row[key]:
                    raise ValueError("Archive and report null truth differ")
            if entry["samples"] != report["samples_per_scenario"]:
                raise ValueError("Null archive replicate count differs from its completed experiment")
            selected.append(entry)
    if not selected:
        raise ValueError("Every included experiment must have at least one true-zero-excitation cell")
    return selected


def _invariant_arrays(reference, candidate):
    # Truth-dependent CDFs and scores intentionally differ; observations,
    # supplied calibration, point estimates and interval endpoints must not.
    def keys(arrays):
        return {k for k in arrays if k == "replicate_id" or k.startswith(("observation__", "noise__", "point__"))
                or (k.startswith("posterior__") and ("__log_lower_" in k or "__log_upper_" in k))}
    names = keys(reference)
    if keys(candidate) != names:
        raise ValueError("Existing null cells expose different invariant array inventories")
    for name in names:
        if not np.array_equal(reference[name], candidate[name], equal_nan=True):
            raise ValueError(f"Null observations or fits vary across stipulated truth masses: {name}")
    return sorted(names)


def build_null_sweep(study_directory, masses=DEFAULT_MASSES):
    """Verify source trials and return a JSON-safe diagnostic (no filesystem writes)."""
    masses = _masses(masses)
    state, experiments = _completed_experiments(study_directory)
    results = []
    for identifier, folder, report, manifest in experiments:
        nulls = _null_entries(report, manifest)
        reference = load_trial_archive(folder, nulls[0])
        _archived_laws(reference, nulls[0]["samples"])
        checked = _invariant_arrays(reference, reference)
        for entry in nulls[1:]:
            candidate = load_trial_archive(folder, entry)
            _archived_laws(candidate, entry["samples"])
            _invariant_arrays(reference, candidate)
        results.append({"id": identifier, "measurement": report["experiment_config"].get("measurement", {"repeats": 1}),
                        "source_null_scenarios": [e["id"] for e in nulls],
                        "source_null_masses": [e["scenario"]["true_mass"] for e in nulls],
                        "archive_sha256": {e["path"]: e["sha256"] for e in nulls},
                        "invariance": {"status": "verified" if len(nulls) > 1 else "single_null_cell_no_cross_mass_check",
                                       "checked_arrays": checked},
                        "rows": score_null_trials(reference, masses)})
    return {"schema_version": 1, "status": "completed", "truth_grid": masses,
            "diagnostic_runtime": {"python": platform.python_version(), "numpy": np.__version__},
            "source_study": str(Path(study_directory).resolve()), "source_spec_sha256": state["spec_sha256"],
            "source_checkpoint_sha256": _hash(Path(study_directory) / "study_checkpoint.json"),
            "source_study_fingerprint": state.get("source_fingerprint"),
            "diagnostic_source_sha256": {Path(__file__).name: _hash(__file__), "trial_archive.py": _hash(Path(__file__).with_name("trial_archive.py"))},
            "interpretation": [
                "The original zero-excitation observation and calibration law is independent of true positive mass. These are valid fixed-truth assessments using the same stored null observations and fitted outputs; no estimator was refitted.",
                "Replicates within each row are independent original simulated experiments. Different truth-grid rows reuse the same trials and are statistically dependent; there is no extra independent replication. Do not pool grid rows or interpret grid length as extra simulation replication.",
                "Changing n reduces known measurement noise but supplies no mass information at exactly zero excitation. The specified noise ratio and inference measure still set the fitted mass scale.",
                "Coverage, lower misses, upper misses and invalid intervals exhaust attempted trials. Invalid intervals fail coverage and make the all-trial mean width unavailable. Invalid points receive the maximum capped loss and fail factor-two accuracy.",
                "Reported MCSEs describe marginal simulation precision, not simultaneous bands. Zero empirical MCSE is not evidence of zero population error probability. No global winner or uniform calibration conclusion is implied."],
            "experiments": results}


def plot_null_sweep(result, destination):
    """Save one publication-style coverage figure; import plotting only on demand."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    count = len(result["experiments"])
    figure, axes = plt.subplots(count, len(LEVELS), figsize=(10, max(3.4, 2.9 * count)), squeeze=False, sharex=True, sharey=True)
    laws = list(dict.fromkeys(law for experiment in result["experiments"]
                              for row in experiment["rows"] for law in row["intervals"]))
    colors = {law: LAW_COLORS.get(law, matplotlib.colormaps["tab10"](i % 10))
              for i, law in enumerate(laws)}
    legend = {}
    for i, experiment in enumerate(result["experiments"]):
        x = np.log2(result["truth_grid"])
        experiment_laws = tuple(experiment["rows"][0]["intervals"])
        if any(set(row["intervals"]) != set(experiment_laws) for row in experiment["rows"]):
            raise ValueError("Null sweep rows contain different posterior-law inventories")
        for j, level in enumerate(LEVELS):
            axis = axes[i, j]
            for law in experiment_laws:
                summaries = [row["intervals"][law][str(level)]["coverage"] for row in experiment["rows"]]
                artist = axis.errorbar(x, [s["mean"] for s in summaries], yerr=[s["mcse"] or 0 for s in summaries],
                                      color=colors[law], marker=LAW_MARKERS.get(law, "o"),
                                      markersize=4, linewidth=1.5, capsize=2,
                                      label=LAW_LABELS.get(law, law.replace("_", " ").title()))
                legend.setdefault(law, artist)
            axis.axhline(level / 100, color="#777777", linestyle="--", linewidth=1)
            axis.set_title(f"{experiment['id']} · {level}% interval", fontsize=11)
            axis.set_ylim(0, 1.035)
            axis.grid(alpha=.2)
            axis.set_xticks(x, [format(m, ".4g") for m in result["truth_grid"]])
            if j == 0:
                axis.set_ylabel("Fixed-truth coverage")
            if i == count - 1:
                axis.set_xlabel("Stipulated true mass (log scale)")
    figure.legend([legend[law] for law in laws],
                  [LAW_LABELS.get(law, law.replace("_", " ").title()) for law in laws],
                  loc="lower center", bbox_to_anchor=(.5, .035), ncol=min(4, len(laws)), frameon=False)
    figure.suptitle("Zero excitation: coverage across indistinguishable mass truths", fontsize=13)
    figure.text(.5, .015, "Same trials reused at every mass. Dashed lines: nominal coverage. Bars: ±1 Monte Carlo SE; not simultaneous confidence intervals.", ha="center", fontsize=8)
    figure.tight_layout(rect=(0, .10, 1, .96))
    figure.savefig(destination, dpi=170)
    plt.close(figure)


def write_null_sweep(study_directory, output, masses=DEFAULT_MASSES):
    output = Path(output).resolve()
    root = Path(study_directory).resolve()
    if output == root or output.is_relative_to(root / "experiments"):
        raise ValueError("Diagnostic output must not overwrite the study or its experiment area")
    if output.exists():
        raise ValueError("Diagnostic output already exists; choose a new directory")
    result = build_null_sweep(root, masses)
    # Missing optional plotting dependencies fail before creating an output.
    import matplotlib
    result["diagnostic_runtime"]["matplotlib"] = matplotlib.__version__
    output.mkdir(parents=True)
    plot_null_sweep(result, output / "null_truth_coverage.png")
    (output / "null_truth_sweep.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path, help="Completed study output directory")
    parser.add_argument("--output", required=True, type=Path, help="New diagnostic output directory")
    parser.add_argument("--masses", nargs="+", type=float, default=DEFAULT_MASSES)
    args = parser.parse_args()
    try:
        result = write_null_sweep(args.study, args.output, args.masses)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(json.dumps({"status": result["status"], "output": str(args.output.resolve()), "experiments": len(result["experiments"])}))


if __name__ == "__main__":
    main()

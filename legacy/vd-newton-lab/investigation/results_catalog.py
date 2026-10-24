"""Read-only adapters for the mass-estimator laboratory's saved experiments.

The browser consumes one schema without rewriting historical results. A run ID
selects a discovered result; it is never interpreted as an arbitrary path.
Missing, nonfinite and infinite-risk values remain distinct from zero.
"""

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
import re


RUNS_DIRECTORY = Path(".tools/mass_estimator_runs")
BOUNDARY_DIRECTORY = Path(".tools/mass_boundary_20260918")
STUDIES_DIRECTORY = Path(".tools/mass_estimator_studies")
_RUN_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,95}\Z")
_METHODS = {
    "flat_joint": {
        "label": "Our ratio of means",
        "formula": "E[f | F,a] / E[alpha | F,a]",
        "assumptions": "Flat nonnegative latent magnitudes with a shared unknown direction; independent isotropic unit Gaussian noise in each channel.",
    },
    "flat_geometric": {
        "label": "Flat-law geometric",
        "formula": "exp(E[log m | F,a])",
        "assumptions": "Same flat posterior law; minimizes posterior squared log error.",
    },
    "flat_median": {
        "label": "Flat-law median",
        "formula": "median(m | F,a)",
        "assumptions": "Same flat posterior law; minimizes posterior absolute log error.",
    },
    "calibrated_joint": {
        "label": "Calibration-integrated ratio of means",
        "formula": "E[f | F,a,calibration] / E[alpha | F,a,calibration]",
        "assumptions": "Physical flat magnitude measure; separate unknown isotropic channel variances integrated using external calibration and explicit scale priors.",
    },
    "calibrated_geometric": {"label": "Calibration-integrated geometric"},
    "calibrated_median": {"label": "Calibration-integrated median"},
    "oracle_joint": {
        "label": "Known-noise oracle ratio of means",
        "assumptions": "Flat law given the true noise SDs: a simulation benchmark with information unavailable to calibrated methods.",
    },
    "oracle_geometric": {"label": "Known-noise oracle geometric"},
    "oracle_median": {"label": "Known-noise oracle median"},
    "flat_reciprocal_root": {
        "label": "Flat-law root-loss point",
        "formula": "E[sqrt(m) | F,a] / E[1/sqrt(m) | F,a]",
        "assumptions": "Same flat posterior law; minimizes posterior normalized reciprocal-root loss.",
    },
    "tube_joint": {
        "label": "Tube ratio of means",
        "formula": "E_tube[f | F,a] / E_tube[alpha | F,a]",
        "assumptions": "Latent measure weighted by f² + alpha², with a shared unknown direction and the same Gaussian likelihood.",
    },
    "norm_ratio": {"label": "Vector-length ratio", "formula": "||F|| / ||a||",
                   "assumptions": "Plug-in vector magnitudes; no vector-length noise correction."},
    "norm_floor": {"label": "Floored vector-length ratio",
                   "formula": "max(||F||, sqrt(2/pi)) / max(||a||, sqrt(2/pi))",
                   "assumptions": "Explicit fixed floor in unit-noise coordinates."},
    "positive_profile": {"label": "Positive TLS / profile"},
    "forward_ols": {"label": "Forward OLS"},
    "reverse_ols": {"label": "Reverse OLS"},
    "noise_corrected_ols": {"label": "Noise-corrected moment"},
}
_BOUNDARY_TITLES = {
    "coarse": "Boundary study · broad force-8 path",
    "fine": "Boundary study · refined crossings",
    "force4": "Boundary study · force-4 path",
    "high": "High-signal study · independent trials",
    "high_antithetic": "High-signal study · antithetic pairs",
}
_METRICS = {
    "squared_log": ("Squared log error", "Mean [log(estimate / true mass)]²; multiplicative errors are treated symmetrically."),
    "abs_log": ("Absolute log error", "Mean absolute log(estimate / true mass); multiplicative errors are treated symmetrically."),
    "reciprocal_root": ("Reciprocal-root loss", "Mean sqrt(r) + 1/sqrt(r) − 2, where r = estimate / true mass; normalized mass-only loss."),
    "abs_relative": ("Absolute relative error", "Mean |estimate − true mass| / true mass."),
    "squared_relative": ("Squared relative error", "Mean [(estimate − true mass) / true mass]²."),
    "absolute_mass": ("Absolute mass error", "Mean |estimate − true mass|, in the simulation's standardized mass units."),
    "squared_mass": ("Squared mass error", "Mean (estimate − true mass)², in squared standardized mass units."),
    "joint_root_loss": ("Magnitude-weighted root loss", "The normalized reciprocal-root loss multiplied by sqrt(true force × true acceleration); it changes scale across signal cells."),
    "heldout_known_acceleration_relative_squared_force_error": ("Held-out force error · known acceleration", "Independent force prediction at twice the training excitation, using exact new acceleration; squared vector residual divided by true force magnitude squared."),
    "heldout_noisy_acceleration_relative_squared_force_error": ("Held-out force error · noisy acceleration", "Independent force prediction at twice the training excitation, using a new noisy acceleration; squared vector residual divided by true force magnitude squared."),
    "log_crps": ("Log-mass CRPS", "Continuous ranked probability score for the full log-mass distribution; scores location and spread together."),
    "log_density_score": ("Log-mass density score", "Negative log posterior density at true log mass; lower is better and negative values are possible."),
}
_DESCRIPTIVE_MCSE = {
    "squared_mass", "squared_relative",
    "heldout_known_acceleration_relative_squared_force_error",
    "heldout_noisy_acceleration_relative_squared_force_error",
}
_INFINITE_NOTE = (
    "Under the declared infinite penalty for invalid positive masses, this estimator's "
    "population risk is infinite. A finite observed mean is only a sample diagnostic."
)
_DESCRIPTIVE_NOTE = (
    "The norm ratio's squared-loss population variance is infinite in this 3D model. "
    "This empirical MCSE is descriptive and does not justify an ordinary normal confidence interval."
)


def _signature(path):
    try:
        stat = path.stat()
        return (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
    except OSError:
        return None


@lru_cache(maxsize=16)
def _read_cached(signature):
    return json.loads(Path(signature[0]).read_text(encoding="utf-8-sig"))


def _read(path):
    signature = _signature(path)
    if signature is None:
        raise FileNotFoundError(path)
    return _read_cached(signature)


def _inside(workspace, path):
    """Exclude escaped links as well as traversal without following outside data."""
    try:
        path.resolve().relative_to(workspace.resolve())
    except (ValueError, OSError):
        return False
    return True


def _sources(workspace):
    workspace = Path(workspace)
    candidates = {
        "full-comparison": (workspace / ".tools/mass_comparison_20260917/full/results.json", "Full comparison · 46 signal cells", "comparison"),
        "comparison-smoke": (workspace / ".tools/mass_comparison_20260917/smoke/results.json", "Comparison · validation smoke run", "validation"),
    }
    directory = workspace / BOUNDARY_DIRECTORY
    if _inside(workspace, directory) and directory.is_dir():
        for path in sorted(directory.glob("*.json")):
            if path.stem == "known_force_limit":
                continue
            run_id = "boundary-" + path.stem
            if _RUN_ID.fullmatch(run_id):
                title = _BOUNDARY_TITLES.get(path.stem, "Boundary study · " + path.stem.replace("_", " "))
                candidates[run_id] = (path, title, "boundary")
    directory = workspace / RUNS_DIRECTORY
    if _inside(workspace, directory) and directory.is_dir():
        for folder in sorted(directory.iterdir(), reverse=True):
            if folder.is_dir() and _RUN_ID.fullmatch(folder.name) and folder.name not in candidates:
                candidates[folder.name] = (folder / "results.json", "Experiment · " + folder.name, "comparison")
    directory = workspace / STUDIES_DIRECTORY
    if _inside(workspace, directory) and directory.is_dir():
        for checkpoint in sorted(directory.glob("*/study_checkpoint.json")):
            if not _inside(workspace, checkpoint):
                continue
            try:
                study = _read(checkpoint)
                if not isinstance(study, dict) or study.get("schema_version") != 1:
                    continue
                for record in study.get("experiments", []):
                    if not isinstance(record, dict) or record.get("status") != "completed":
                        continue
                    run_id = f"study-{checkpoint.parent.name}-{record['id']}"
                    path = checkpoint.parent / record["attempt_directory"] / "results.json"
                    if (not _RUN_ID.fullmatch(run_id) or not _inside(checkpoint.parent, path)
                            or not path.is_file()):
                        continue
                    # Only expose the checkpoint's completed, unchanged result.
                    # Full archive verification remains the study resume/loader's job.
                    if hashlib.sha256(path.read_bytes()).hexdigest() != record["output_sha256"].get("results.json"):
                        continue
                    candidates[run_id] = (path, f"Study · {record['id']}", "comparison")
            except (OSError, ValueError, TypeError, KeyError):
                continue
    return {key: value for key, value in candidates.items()
            if _inside(workspace, value[0]) and value[0].is_file()}


def _optional_file(workspace, path):
    if _inside(workspace, path) and path.is_file():
        value = _read(path)
        return value if isinstance(value, dict) else None
    return None


def _config_for(workspace, run_id, path, data):
    executed = data.get("experiment_config")
    if isinstance(executed, dict):
        return executed, "embedded_executed_config"
    config = _optional_file(workspace, path.with_name("config.json"))
    if config is not None:
        return config, "saved_config"
    # These runs predate configuration files, but their report records all
    # simulation settings. Reconstruct only when the full preset is present.
    expected = {"all": 46, "original": 16, "fixed_mass": 30}
    fields = ("samples_per_scenario", "training_seed", "heldout_seed", "quadrature_order", "batch_size")
    selection = data.get("scenario_selection")
    if (run_id in ("full-comparison", "comparison-smoke") and all(key in data for key in fields)
            and selection in expected and len(data.get("rows", [])) == expected[selection]):
        return {
            "schema_version": 1,
            "name": "Reproduction of " + ("full comparison" if run_id == "full-comparison" else "validation comparison"),
            "samples": data["samples_per_scenario"], "order": data["quadrature_order"],
            "batch_size": data["batch_size"], "training_seed": data["training_seed"],
            "heldout_seed": data["heldout_seed"], "save_trials": False,
            "noise_model": "isotropic_unit_3d", "scenario": {"type": "preset", "selection": selection},
        }, "reconstructed_from_recorded_settings"
    return None, None


def _summary(workspace, run_id, entry, data):
    path, title, kind = entry
    config, origin = _config_for(workspace, run_id, path, data)
    if config and origin != "reconstructed_from_recorded_settings" and isinstance(config.get("name", config.get("title")), str):
        title = config.get("name", config.get("title"))
    return {
        "id": run_id, "title": title, "kind": kind,
        "samples": data.get("n") if kind == "boundary" else data.get("samples_per_scenario"),
        "scenario_count": len(data.get("rows", [])),
        "source": path.relative_to(Path(workspace)).as_posix(),
        "config_available": config is not None,
    }


def list_runs(workspace):
    """List completed saved results; pending runs have no results.json yet.

    Unreadable/incomplete JSON is omitted, so an interrupted atomic write or a
    damaged optional experiment cannot prevent browsing the remaining library.
    Calling load_run explicitly still exposes the corresponding read error.
    """
    output = []
    for run_id, entry in _sources(Path(workspace)).items():
        try:
            data = _read(entry[0])
            if not isinstance(data, dict) or not isinstance(data.get("rows"), list):
                continue
            output.append(_summary(Path(workspace), run_id, entry, data))
        except (OSError, ValueError, TypeError):
            continue
    return output


def _finite(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        return value
    return None


def _value(mean, mcse=None, status="finite_sample", note=None):
    finite_mean = _finite(mean)
    if finite_mean is None and status.startswith("finite"):
        status = "unavailable"
    result = {"mean": finite_mean, "mcse": _finite(mcse), "status": status}
    if note:
        result["note"] = note
    return result


def _sample_value(summary):
    return _value(summary.get("mean"), summary.get("empirical_mcse", summary.get("mcse")),
                  summary.get("status", "finite_sample"), summary.get("reason"))


def _risk(value, method, metric, method_summary, paired=False):
    if value["status"] == "not_applicable":
        return value
    bounded = metric.startswith(("outside_factor_", "capped_squared_log_factor_"))
    if method_summary.get("population_unbounded_risk") and not bounded:
        value["status"] = "infinite_population_risk"
        value["note"] = method_summary.get("population_risk_note", _INFINITE_NOTE)
        if paired:
            value["note"] += " The displayed finite-sample difference is not a finite population-risk difference."
    elif method == "norm_ratio" and metric in _DESCRIPTIVE_MCSE:
        value["status"] = "descriptive_mcse" if value["mean"] is not None else value["status"]
        value["note"] = _DESCRIPTIVE_NOTE
    return value


def _method_metadata(names, baseline=None):
    result = []
    for name in names:
        item = dict(_METHODS.get(name, {"label": name.replace("_", " ")}))
        item.update((baseline or {}).get(name, {}))
        if "risk_note" in item:
            item["note"] = item["risk_note"]
        result.append(dict(id=name, **item))
    return result


def _metric_metadata(name):
    item = {"id": name, "better": "lower"}
    if name in _METRICS:
        item["label"], item["description"] = _METRICS[name]
    elif name.startswith("capped_squared_log_factor_"):
        factor = name.removeprefix("capped_squared_log_factor_")
        item.update(label=f"Squared log error · cap at factor {factor}",
                    description=f"Mean min(log(estimate / true mass)², log({factor})²); invalid estimates receive the maximum cap.")
    elif name.startswith("outside_factor_"):
        factor = name.removeprefix("outside_factor_")
        item.update(label=f"Failure outside factor {factor}",
                    description=f"Fraction outside [true mass / {factor}, true mass × {factor}], including every invalid estimate. Lower is better; one minus this value is tolerance success.")
    elif name.startswith("interval_score_"):
        level = name.removeprefix("interval_score_")
        item.update(label=f"{level}% log interval score",
                    description=f"Width of the central {level}% log-mass interval plus the appropriate penalty when truth falls outside it.")
    elif name.startswith("coverage_"):
        level = name.removeprefix("coverage_")
        item.update(label=f"{level}% interval coverage", better="target", target=float(level) / 100,
                    description=f"Fraction of central {level}% intervals containing true mass. Compare with the nominal {level}% target and consider interval width and score together.")
    elif name.startswith("log_width_"):
        level = name.removeprefix("log_width_")
        item.update(label=f"{level}% mean log interval width", better="context",
                    description="Average upper minus lower log-mass interval endpoint. Narrowness alone is not accuracy; interpret together with coverage and interval score.")
    elif name.startswith(("below_interval_", "above_interval_")):
        level = name.rsplit("_", 1)[1]
        below = name.startswith("below")
        item.update(label=f"{level}% interval · truth {'below' if below else 'above'}",
                    better="target", target=(1-float(level)/100)/2,
                    description="Fraction of experiments with true mass outside this end of the central interval. The equal-tail sampling target is a calibration reference, not a guarantee from the conditional law.")
    else:
        item.update(label=name.replace("_", " "), description="Saved experiment metric; consult the run's provenance for its definition.")
    return item


def _comparison(data):
    rows, method_names, metric_names = [], {}, {}
    for original in data["rows"]:
        row = {
            "id": original["id"], "f": original["true_force_snr"],
            "a": original["true_acceleration_snr"], "mass": original["true_mass"],
            "strength": original["total_signal_snr"], "family": original["family"],
            "values": {}, "paired": {}, "diagnostics": {},
            "refinement": original.get("refinement"),
            "heldout_reference_values": original.get("heldout_reference_values"),
        }
        for key in ("true_force_magnitude", "true_acceleration_magnitude", "direction", "repeats",
                    "effective_force_snr", "effective_acceleration_snr", "effective_total_signal_snr", "metric_applicability"):
            if key in original:
                row[key] = original[key]
        tolerance = {"factors": data.get("tolerance_factors", []), "methods": {}}
        for method, summary in original.get("methods", {}).items():
            method_names[method] = None
            values = {}
            for metric, mean in summary.get("means", {}).items():
                value = _value(mean, summary.get("mcse", {}).get(metric),
                               summary.get("mean_status", {}).get(metric, "finite_sample"))
                values[metric] = _risk(value, method, metric, summary)
            for metric, extra in summary.get("additional_tasks", {}).items():
                values[metric] = _risk(_sample_value(extra), method, metric, summary)
            row["values"][method] = values
            metric_names.update(dict.fromkeys(values))
            logs = summary.get("log_diagnostics", {})
            row["diagnostics"][method] = {
                "invalid_count": summary.get("invalid_count"),
                "finite_positive_fraction": summary.get("finite_positive_fraction"),
                "bias": logs.get("signed_log_bias"), "variance": logs.get("log_variance"),
                "status": summary.get("status"),
                "absolute_log_quantiles": summary.get("absolute_log_quantiles", {}),
                "absolute_log_quantile_status": summary.get("absolute_log_quantile_status", {}),
            }
            if "tolerance_curve" in summary:
                tolerance["methods"][method] = summary["tolerance_curve"]
        if tolerance["factors"] and tolerance["methods"]:
            row["tolerance"] = tolerance
        for method, summaries in original.get("paired_flat_minus_other", {}).items():
            pair_values = {}
            for metric, summary in summaries.items():
                value = _value(summary.get("mean_difference"), summary.get("mcse"),
                               summary.get("status", "finite_sample_comparison"), summary.get("reason", summary.get("mcse_interpretation")))
                pair_values[metric] = _risk(value, method, metric, original["methods"].get(method, {}), paired=True)
            row["paired"][method] = pair_values
            metric_names.update(dict.fromkeys(pair_values))
        for method, uncertainty in original.get("uncertainty", {}).items():
            method_names[method] = None
            values = row["values"].setdefault(method, {})
            for metric, summary in uncertainty.get("scores", {}).items():
                values[metric] = _sample_value(summary)
            for level, interval in uncertainty.get("intervals", {}).items():
                for source_key, metric in (("coverage", f"coverage_{level}"),
                                           ("below_interval", f"below_interval_{level}"),
                                           ("above_interval", f"above_interval_{level}"),
                                           ("mean_log_width", f"log_width_{level}"),
                                           ("interval_score", f"interval_score_{level}")):
                    if source_key in interval:
                        values[metric] = _sample_value(interval[source_key])
            metric_names.update(dict.fromkeys(values))
        law_pairs = original.get("paired_uncertainty_flat_minus_other", {
            "tube_joint": original.get("paired_uncertainty_flat_minus_tube", {})})
        for method, comparisons in law_pairs.items():
            for metric, summary in comparisons.items():
                row["paired"].setdefault(method, {})[metric] = _value(
                    summary.get("mean_difference"), summary.get("mcse"), summary.get("status", "finite_sample_comparison"))
                metric_names[metric] = None
        rows.append(row)
    metadata = {key: value for key, value in data.items()
                if key not in ("rows", "baseline_metadata", "tolerance_factors")}
    metrics = [_metric_metadata(key) for key in metric_names]
    if data.get("noise_model") == "diagonal_gaussian_3d":
        for metric in metrics:
            metric["description"] = metric["description"].replace("the simulation's standardized mass units", "physical force/acceleration mass units").replace("squared standardized mass units", "squared physical mass units")
    return rows, _method_metadata(method_names, data.get("baseline_metadata")), metrics, metadata


def _boundary(workspace, data):
    rows, metric_names = [], {}
    names = {"ours": "flat_joint", "norm": "norm_ratio"}
    aliases = {"root": "reciprocal_root", "factor_two_failure": "outside_factor_2"}
    for index, original in enumerate(data["rows"]):
        f, a = original["f"], original["a"]
        row = {"id": f"{index:02d}_boundary", "f": f, "a": a, "mass": original["truth"],
               "strength": math.hypot(f, a), "family": "boundary",
               "values": {}, "paired": {}, "diagnostics": {},
               "refinement": original.get("refinement"),
               "high_signal_gap_prediction": original.get("high_signal_gap_prediction")}
        for old, method in names.items():
            summary = original.get(old, {})
            values = {}
            if "squared_log" in summary:
                values["squared_log"] = _sample_value(summary["squared_log"])
            if "success2" in summary:
                success = _finite(summary["success2"])
                values["outside_factor_2"] = _value(
                    None if success is None else 1 - success, None, "finite_sample",
                    "The legacy run saved the fraction but no individual MCSE for this tolerance.")
            row["values"][method] = values
            row["diagnostics"][method] = {"bias": summary.get("bias"), "variance": summary.get("variance")}
            metric_names.update(dict.fromkeys(values))
        for old, summary in original.get("delta", {}).items():
            method, metric = ("norm_floor", "squared_log") if old == "squared_log_vs_floor" else ("norm_ratio", aliases.get(old, old))
            row["paired"].setdefault(method, {})[metric] = _sample_value(summary)
            metric_names[metric] = None
        rows.append(row)
    metadata = {key: value for key, value in data.items() if key != "rows"}
    antithetic = bool(data.get("antithetic", False))
    metadata.update({
        "antithetic": antithetic,
        "independent_units": data["n"] // 2 if antithetic else data["n"],
        "mcse_basis": ("Independent antithetic pair averages; n/2 independent pairs." if antithetic else "Independent per-trial losses and paired estimator differences."),
        "paired_difference_convention": "flat_joint minus comparator; negative means lower loss for flat_joint",
        "assumptions": "One codirectional latent 3D force/acceleration pair with independent isotropic unit Gaussian component noise. The same noise draws are reused across signal cells; cells are not independent replications.",
        "legacy_scope": "Legacy boundary studies save individual squared log losses and factor-two success, and paired differences for four scores. Missing individual scores and MCSEs are unavailable, not zero. No unrecorded scores are reconstructed.",
    })
    reference_path = Path(workspace) / BOUNDARY_DIRECTORY / "known_force_limit.json"
    reference = _optional_file(workspace, reference_path)
    if reference:
        metadata["known_force_limit_reference"] = {
            "source": reference_path.relative_to(Path(workspace)).as_posix(),
            "description": "Deterministic quadrature/root calculations in the known-force limit; these are separate reference values, not sampled finite-force crossings.",
            "values": reference,
        }
    return rows, _method_metadata(["flat_joint", "norm_ratio", "norm_floor"]), [_metric_metadata(key) for key in metric_names], metadata


@lru_cache(maxsize=12)
def _load_cached(workspace_string, run_id, source_signature, config_signature, status_signature, reference_signature):
    workspace = Path(workspace_string)
    entry = _sources(workspace).get(run_id)
    if entry is None:
        raise KeyError("Unknown run")
    path, _, kind = entry
    data = _read(path)
    summary = _summary(workspace, run_id, entry, data)
    rows, methods, metrics, metadata = _boundary(workspace, data) if kind == "boundary" else _comparison(data)
    metadata["source"] = summary["source"]
    config, origin = _config_for(workspace, run_id, path, data)
    metadata["config_provenance"] = origin
    if origin == "reconstructed_from_recorded_settings":
        metadata["config_note"] = (
            "Configuration reconstructed from the recorded preset, sample count, seeds and integration settings. "
            "Trial archiving defaults to off. The current engine may differ from the recorded source hashes."
        )
    status = _optional_file(workspace, path.with_name("status.json"))
    if status is not None:
        metadata["run_status"] = status
    return dict(summary, metadata=metadata, methods=methods, metrics=metrics, rows=rows, config=config)


def load_run(workspace, run_id):
    """Normalize one known result without changing its files or trial counts.

    Parsed files and normalized views are cached by path, nanosecond mtime and
    size. Returned dictionaries are read-only by convention; serialize them or
    copy them before making changes. Discovery and containment are rechecked.
    """
    workspace = Path(workspace).resolve()
    if not isinstance(run_id, str) or not _RUN_ID.fullmatch(run_id):
        raise KeyError("Unknown run")
    entry = _sources(workspace).get(run_id)
    if entry is None:
        raise KeyError("Unknown run")
    path = entry[0]
    return _load_cached(str(workspace), run_id, _signature(path),
                        _signature(path.with_name("config.json")),
                        _signature(path.with_name("status.json")),
                        _signature(workspace / BOUNDARY_DIRECTORY / "known_force_limit.json"))

"""Verify a completed study and export cell-specific summaries and figures.

This is a read-only results consumer: it never imports the fitting engine,
changes a study, recomputes its estimates, pools truth cells, or ranks winners.
Matplotlib is required only for exporting figures.
"""

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import sys


POINT_METRICS = ("capped_squared_log_factor_2", "outside_factor_2")
# The archived scalar metrics, not the plotting tolerance grid. Keep this
# consumer independent of the scoring/fitting imports used to create a study.
ACCURACY_FACTORS = (1.1, 1.25, 1.5, 2., 3., 5., 10., 20.)
DEFAULT_COMPARATORS = (
    "flat_geometric", "flat_median", "norm_ratio", "norm_floor",
    "positive_profile", "covariance_profile", "forward_ols", "weighted_forward_ols",
    "calibrated_joint", "oracle_joint",
)
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
NOTES = [
    "Every value describes a single experiment and truth cell; no pooling or global winner.",
    "Error bars show plus or minus one Monte Carlo standard error (MCSE), not confidence intervals or simultaneous uncertainty bands.",
    "R is the number of independent simulated experiments per cell; n is the number of stable readings averaged inside each experiment.",
    "The horizontal axis is per-reading acceleration SNR. Force SNR also changes with true mass; noise anisotropy and calibration are retained in each experiment's metadata.",
    "Cells within an experiment may share noise. Lines connect descriptive cell estimates and do not treat those cells as independent replicates.",
    "A zero observed failure rate or zero empirical MCSE does not establish zero population risk or uncertainty.",
    "Unbounded-score MCSEs are finite-sample diagnostics unless population moments are established. Source analytic risk warnings remain attached.",
    "Conditional law coverage and interval tails are assessed at fixed simulation truth. Nominal conditional probability need not equal repeated-experiment coverage.",
    "Zero excitation has no mass information under mass-independent measurement noise; apparent recovery at selected masses is not identification.",
]


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _inside(root, relative):
    if not isinstance(relative, str) or "\\" in relative or ":" in relative:
        raise ValueError("Study paths must be relative POSIX paths")
    part = Path(relative)
    root = Path(root).resolve()
    if part.is_absolute() or ".." in part.parts:
        raise ValueError("Study path escapes its directory")
    result = (root / part).resolve()
    if not result.is_relative_to(root) or result == root:
        raise ValueError("Study path escapes its directory")
    return result


def _read_json(path):
    def reject_nonfinite(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=reject_nonfinite)
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot read finite JSON from {path}: {error}") from error


def load_completed_study(study_dir):
    """Verify checkpoint, specification and every recorded completed artifact.

    Verification is against the saved checkpoint, not today's source tree.
    Consequently historical completed studies remain independently readable.
    Returns a checkpoint plus reports with their source locations and hashes.
    """
    root = Path(study_dir).resolve()
    checkpoint_path = _inside(root, "study_checkpoint.json")
    checkpoint = _read_json(checkpoint_path)
    if (not isinstance(checkpoint, dict) or type(checkpoint.get("schema_version")) is not int
            or checkpoint["schema_version"] != 1 or checkpoint.get("status") != "completed"):
        raise ValueError("A completed schema_version 1 study checkpoint is required")
    spec = checkpoint.get("specification")
    if not isinstance(spec, dict) or hashlib.sha256(_canonical(spec)).hexdigest() != checkpoint.get("spec_sha256"):
        raise ValueError("Checkpoint specification hash mismatch")
    saved_spec = _read_json(_inside(root, "study.json"))
    if _canonical(saved_spec) != _canonical(spec):
        raise ValueError("study.json differs from checkpoint specification")
    specs = spec.get("experiments")
    records = checkpoint.get("experiments")
    if not isinstance(specs, list) or not specs or not isinstance(records, list) or len(specs) != len(records):
        raise ValueError("Study experiment list is missing or incomplete")
    expected_ids = [item.get("id") for item in specs]
    if (any(not isinstance(name, str) or not _TOKEN.fullmatch(name) for name in expected_ids)
            or len(set(name.casefold() for name in expected_ids)) != len(expected_ids)
            or [item.get("id") for item in records] != expected_ids):
        raise ValueError("Study experiment IDs disagree or are duplicated")
    reports = []
    for specification, record in zip(specs, records):
        if record.get("status") != "completed":
            raise ValueError(f"Experiment {record['id']} is not completed")
        attempt = _inside(root, record.get("attempt_directory"))
        hashes = record.get("output_sha256")
        if not isinstance(hashes, dict) or "results.json" not in hashes or "config.json" not in hashes:
            raise ValueError("Completed experiment requires results and configuration hashes")
        for relative, expected in hashes.items():
            if not isinstance(expected, str) or re.fullmatch(r"[0-9a-f]{64}", expected) is None:
                raise ValueError("Invalid recorded SHA-256")
            artifact = _inside(attempt, relative)
            if not artifact.is_file() or _hash(artifact) != expected:
                raise ValueError(f"Completed artifact hash mismatch: {artifact}")
        config = _read_json(attempt / "config.json")
        report = _read_json(attempt / "results.json")
        if config != specification.get("config") or report.get("experiment_config") != config:
            raise ValueError("Completed result configuration disagrees with study")
        if (report.get("samples_per_scenario") != config.get("samples")
                or report.get("samples_per_scenario") != spec.get("fixed_replicates_per_cell")):
            raise ValueError("Completed result replicate count disagrees with study")
        rows = report.get("rows")
        if not isinstance(rows, list) or not rows:
            raise ValueError("Completed result has no scenario rows")
        ids = [row.get("id") for row in rows]
        if any(not isinstance(name, str) for name in ids) or len(set(ids)) != len(ids):
            raise ValueError("Duplicate or invalid result scenario IDs")
        for row in rows:
            if not isinstance(row.get("methods"), dict) or "flat_joint" not in row["methods"]:
                raise ValueError("Every scenario must contain flat_joint point results")
            for summary in row["methods"].values():
                if not isinstance(summary, dict) or summary.get("n") != config["samples"]:
                    raise ValueError("Method attempted-trial denominator differs from planned R")
        reports.append({"id": record["id"], "relative_results":
                        (attempt / "results.json").relative_to(root).as_posix(),
                        "sha256": hashes["results.json"], "verified_artifacts": len(hashes),
                        "report": report, "config": config})
    return {"study_dir": str(root), "checkpoint": checkpoint,
            "checkpoint_sha256": _hash(checkpoint_path), "reports": reports}


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _signed_log_bias(summary):
    diagnostic = summary.get("log_diagnostics", {})
    mean, variance, count = diagnostic.get("signed_log_bias"), diagnostic.get("log_variance"), summary.get("n")
    mcse = None
    if (_finite(mean) and _finite(variance) and variance >= 0
            and diagnostic.get("variance_ddof") == 0 and isinstance(count, int) and count > 1):
        mcse = math.sqrt(variance / (count - 1))
    return {"mean": mean, "mcse": mcse, "status": diagnostic.get("status", "unavailable"),
            "mcse_basis": "sqrt(saved_ddof0_log_variance/(R-1)); empirical finite-R diagnostic, not proof of population moments"}


def summarize_study(loaded, comparators=DEFAULT_COMPARATORS):
    """Extract original numbers, preserving unavailable values and source warnings."""
    if (not isinstance(comparators, (list, tuple)) or not comparators
            or any(not isinstance(name, str) or not _TOKEN.fullmatch(name) for name in comparators)
            or len(set(comparators)) != len(comparators) or "flat_joint" in comparators):
        raise ValueError("Comparators must be unique method names other than flat_joint")
    experiments = []
    for saved in loaded["reports"]:
        report = saved["report"]
        rows = []
        for source in report["rows"]:
            methods = deepcopy(source["methods"])
            for summary in methods.values():
                summary.pop("tolerance_curve", None)
            pairs = {}
            for name in comparators:
                source_pair = source.get("paired_flat_minus_other", {}).get(name)
                pairs[name] = ({"status": "unavailable", "reason": "comparator_not_in_this_experiment"}
                               if source_pair is None else {
                                   "status": "available", "metrics": {
                                       metric: deepcopy(source_pair.get(metric, {"status": "metric_unavailable"}))
                                       for metric in POINT_METRICS}})
            laws = {}
            for name, law in source.get("uncertainty", {}).items():
                laws[name] = {"intervals": {level: deepcopy(law.get("intervals", {}).get(level, {"status": "unavailable"}))
                                            for level in ("50", "95")},
                              "scores": deepcopy(law.get("scores", {}))}
            excluded = {"methods", "paired_flat_minus_other", "uncertainty", "paired_uncertainty_flat_minus_tube",
                        "paired_uncertainty_flat_minus_other"}
            rows.append({"scenario": deepcopy({key: value for key, value in source.items() if key not in excluded}),
                         "methods": methods, "flat_signed_log_bias": _signed_log_bias(methods["flat_joint"]),
                         "laws": laws, "paired_flat_minus_other": pairs,
                         "paired_uncertainty_flat_minus_other": deepcopy(source.get("paired_uncertainty_flat_minus_other", {})),
                         "paired_uncertainty_flat_minus_tube": deepcopy(source.get("paired_uncertainty_flat_minus_tube", {}))})
        experiments.append({"id": saved["id"], "source_results": saved["relative_results"],
                            "source_results_sha256": saved["sha256"], "verified_artifacts": saved["verified_artifacts"],
                            "config": deepcopy(saved["config"]), "samples_per_cell": report["samples_per_scenario"],
                            "metadata": deepcopy({key: value for key, value in report.items()
                                                  if key not in {"rows", "tolerance_factors", "experiment_config"}}),
                            "rows": rows})
    spec = loaded["checkpoint"]["specification"]
    return {"schema_version": 1, "kind": "completed_study_assessment", "study_dir": loaded["study_dir"],
            "source_checkpoint_sha256": loaded["checkpoint_sha256"], "source_spec_sha256": loaded["checkpoint"]["spec_sha256"],
            "study_name": spec.get("name"), "phase": spec.get("phase"),
            "study_provenance": deepcopy(spec.get("provenance", {})),
            "source_interpretation": deepcopy(spec.get("interpretation", [])),
            "comparators": list(comparators), "notes": NOTES,
            "analysis_runtime": {"python": platform.python_version(), "matplotlib": None, "numpy": None},
            "summary_source_sha256": _hash(__file__), "experiments": experiments}


def _accuracy_arguments(tolerance_factors, success_targets, alpha):
    for values, name in ((tolerance_factors, "Tolerance factors"), (success_targets, "Success targets")):
        if (not isinstance(values, (list, tuple)) or not values
                or any(not _finite(value) for value in values)
                or len(set(values)) != len(values)):
            raise ValueError(f"{name} must be a nonempty sequence of unique finite numbers")
    if any(value not in ACCURACY_FACTORS for value in tolerance_factors):
        raise ValueError(f"Tolerance factors must use supported archived metrics: {ACCURACY_FACTORS}")
    if any(not 0 < value < 1 for value in success_targets):
        raise ValueError("Success targets must lie strictly between zero and one")
    if not _finite(alpha) or not 0 < alpha < 1:
        raise ValueError("Alpha must lie strictly between zero and one")


def _accuracy_metric(summary, factor, targets, epsilon):
    metric = "outside_factor_" + format(factor, "g")
    failure = summary.get("means", {}).get(metric)
    source_status = summary.get("mean_status", {}).get(metric)
    if _finite(failure) and not 0 <= failure <= 1:
        raise ValueError(f"Archived failure probability outside [0, 1]: {metric}")
    available = _finite(failure) and source_status in (None, "finite_sample_mean", "finite_sample")
    result = {"source_metric": metric, "source_status": source_status,
              "status": "available" if available else (source_status or "metric_unavailable"),
              "failure_probability": failure, "success_probability": None,
              "empirical_mcse": summary.get("mcse", {}).get(metric),
              "simultaneous_success_interval": None,
              "targets": {str(float(target)): "unavailable" for target in targets}}
    if available:
        success = 1 - failure
        lower, upper = max(0., success - epsilon), min(1., success + epsilon)
        result.update(success_probability=success,
                      simultaneous_success_interval={"lower": lower, "upper": upper},
                      targets={str(float(target)): ("supported" if lower >= target else
                               "below_target" if upper < target else "unresolved") for target in targets})
    return result


def summarize_accuracy(loaded, tolerance_factors=(1.25, 1.5, 2.0), success_targets=(.90, .95), alpha=.05):
    """Assess absolute accuracy using archived losses and simultaneous MC bounds.

    DKW bounds each cell-method error CDF, with invalid estimates at infinity.
    A union bound covers all archived cell-method distributions, so selecting
    fewer displayed factors does not shrink the family or change its bound.
    Trials must be independent within each cell; dependence across methods or
    cells is permitted. No fitting, pooling, interpolation, or observed-data
    acceptance rule is introduced here.
    """
    _accuracy_arguments(tolerance_factors, success_targets, alpha)
    family_size = sum(len(row["methods"]) for saved in loaded["reports"] for row in saved["report"]["rows"])
    if family_size < 1:
        raise ValueError("Accuracy assessment requires at least one archived cell-method distribution")
    experiments = []
    for saved in loaded["reports"]:
        report = saved["report"]
        samples = report["samples_per_scenario"]
        if type(samples) is not int or samples < 1:
            raise ValueError("Accuracy assessment requires a positive integer R per cell")
        epsilon = math.sqrt((math.log(2 * family_size) - math.log(alpha)) / (2 * samples))
        rows = []
        for source in report["rows"]:
            methods = {}
            for name, summary in source["methods"].items():
                if summary.get("n") != samples:
                    raise ValueError("Method attempted-trial denominator differs from R")
                methods[name] = {"n": samples, "invalid_count": summary.get("invalid_count"),
                                 "valid_count": summary.get("valid_count"), "source_status": summary.get("status"),
                                 "accuracy": {format(factor, "g"): _accuracy_metric(summary, factor, success_targets, epsilon)
                                              for factor in tolerance_factors}}
            pairs = {}
            for name in methods:
                if name == "flat_joint":
                    continue
                pairs[name] = {}
                for factor in tolerance_factors:
                    key = format(factor, "g")
                    metric = "outside_factor_" + key
                    source_pair = source.get("paired_flat_minus_other", {}).get(name, {}).get(metric)
                    pair = deepcopy(source_pair) if source_pair is not None else {}
                    pair["source_status"] = pair.get("status")
                    pair["source_metric"] = metric
                    pair["simultaneous_failure_difference_interval"] = None
                    difference = pair.get("mean_difference")
                    if _finite(difference) and not -1 <= difference <= 1:
                        raise ValueError("Archived failure-probability difference outside [-1, 1]")
                    if pair.get("n", samples) != samples:
                        raise ValueError("Paired attempted-trial denominator differs from R")
                    marginals_available = all(methods[method]["accuracy"][key]["status"] == "available"
                                              for method in ("flat_joint", name))
                    if (_finite(difference) and marginals_available
                            and pair["source_status"] in (None, "finite_sample_comparison", "finite_sample")):
                        pair["status"] = "available"
                        pair["simultaneous_failure_difference_interval"] = {
                            "lower": max(-1., difference - 2 * epsilon),
                            "upper": min(1., difference + 2 * epsilon)}
                    else:
                        pair["status"] = pair["source_status"] or "metric_unavailable"
                        if not marginals_available:
                            pair["status"] = "marginal_metric_unavailable"
                    pair.setdefault("mean_difference", None)
                    pair.setdefault("mcse", None)
                    pairs[name][key] = pair
            excluded = {"methods", "paired_flat_minus_other", "uncertainty", "paired_uncertainty_flat_minus_tube",
                        "paired_uncertainty_flat_minus_other"}
            rows.append({"scenario": deepcopy({key: value for key, value in source.items() if key not in excluded}),
                         "dkw_epsilon": epsilon, "methods": methods, "paired_flat_minus_other": pairs})
        experiments.append({"id": saved["id"], "source_results": saved["relative_results"],
                            "source_results_sha256": saved["sha256"], "verified_artifacts": saved["verified_artifacts"],
                            "config": deepcopy(saved["config"]), "samples_per_cell": samples,
                            "source_failure_policy": report.get("failure_policy"),
                            "source_assumptions": report.get("assumptions"), "rows": rows})
    spec = loaded["checkpoint"]["specification"]
    return {"schema_version": 1, "kind": "completed_study_accuracy", "study_dir": loaded["study_dir"],
            "source_checkpoint_sha256": loaded["checkpoint_sha256"], "source_spec_sha256": loaded["checkpoint"]["spec_sha256"],
            "study_name": spec.get("name"), "phase": spec.get("phase"),
            "study_provenance": deepcopy(spec.get("provenance", {})),
            "source_interpretation": deepcopy(spec.get("interpretation", [])),
            "tolerance_factors": list(tolerance_factors), "success_targets": list(success_targets),
            "simultaneous_band": {"alpha": alpha, "family_size": family_size,
                                  "formula": "epsilon_cell = sqrt(log(2 * L / alpha) / (2 * R_cell))",
                                  "scope": "L is the sum of all archived method counts over every report cell, including missing requested metrics; DKW is uniform in tolerance within each distribution.",
                                  "pair_formula": "flat-minus-other failure difference plus or minus 2 * epsilon_cell, clipped to [-1, 1]; empirical paired MCSE is retained separately."},
            "invalid_trial_policy": "Nonpositive or nonfinite estimates fail every finite tolerance; all attempted trials remain in R. Success includes both endpoints [true_mass / factor, true_mass * factor].",
            "notes": [
                "Bounds quantify Monte Carlo sampling error under independent trials within each fixed cell; cross-cell and cross-method dependence does not invalidate the union bound.",
                "The bound covers the complete archived cell-method family and every tolerance, not only the selected displays; it does not cover additional studies or future simulations.",
                "Supported means the simultaneous lower success bound reaches the stated target; below_target means the upper bound is strictly below it; unresolved includes equality at the upper bound.",
                "Missing source metrics retain their statuses and receive no classification or interval; missing never means zero failure.",
                "Zero observed failures still yield a lower bound below one. Empirical MCSEs are not the simultaneous bounds.",
                "This is a truth-indexed operating map on tested grid cells only. It establishes no interval between cells, observed-data acceptance rule, deployment adequacy, or robustness beyond the simulated model.",
                "These bounds do not include numerical fitting error, model misspecification, or uncertainty about the real experiment. Chosen success targets and tolerances are assessment criteria, not universal requirements.",
            ],
            "analysis_runtime": {"python": platform.python_version()},
            "summary_source_sha256": _hash(__file__), "experiments": experiments}


def _plotting():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Matplotlib is required for figures; install it or expose its package directory via PYTHONPATH") from error
    return plt


def _experiment_label(experiment):
    config = experiment["config"]
    measurement = config.get("measurement", {})
    calibration = measurement.get("calibration", {}).get("mode", "known")
    return (f"{experiment['id']} | n={measurement.get('repeats', 1)} stable readings | "
            f"R={experiment['samples_per_cell']} per cell | {config.get('noise_model', 'unspecified noise')} | "
            f"calibration={calibration}")


def _groups(experiment):
    groups = {}
    for row in experiment["rows"]:
        mass = row["scenario"].get("true_mass")
        snr = row["scenario"].get("true_acceleration_snr")
        if not _finite(mass) or not _finite(snr) or snr < 0:
            raise ValueError("Figures require finite true mass and nonnegative per-reading acceleration SNR")
        groups.setdefault(mass, []).append(row)
    return [(mass, sorted(rows, key=lambda row: row["scenario"]["true_acceleration_snr"]))
            for mass, rows in sorted(groups.items())]


def _draw(ax, groups, getter, *, scale=1., label_prefix="", linestyle="-", colors=None):
    """Connect cells only for display; missing values remain gaps."""
    for index, (mass, rows) in enumerate(groups):
        triples = [getter(row) for row in rows]
        x = [row["scenario"]["true_acceleration_snr"] for row in rows]
        y = [value * scale if _finite(value) else float("nan") for value, _ in triples]
        # Draw the mean even when no MCSE is available; never substitute zero.
        color = colors[index % len(colors)] if colors else f"C{index % 10}"
        ax.plot(x, y, marker="o", markersize=4, linestyle=linestyle, color=color,
                label=f"{label_prefix}m={mass:g}", linewidth=1.5)
        finite = [(xx, yy, se * scale) for xx, yy, (_, se) in zip(x, y, triples)
                  if _finite(yy) and _finite(se) and se >= 0]
        if finite:
            xx, yy, se = zip(*finite)
            ax.errorbar(xx, yy, yerr=se, fmt="none", capsize=2, color=color, linewidth=1)
    values = sorted({row["scenario"]["true_acceleration_snr"] for _, rows in groups for row in rows})
    ax.set_xscale("symlog", linthresh=max(min((x for x in values if x > 0), default=1.) / 2, 1e-12))
    ax.set_xticks(values, [f"{value:g}" for value in values])
    ax.set_xlabel("Per-reading acceleration SNR")
    ax.grid(alpha=.22)


def _point(row, metric):
    summary = row["methods"]["flat_joint"]
    return summary.get("means", {}).get(metric), summary.get("mcse", {}).get(metric)


def _interval(row, law, metric, level="95"):
    summary = row.get("laws", {}).get(law, {}).get("intervals", {}).get(level, {}).get(metric, {})
    return summary.get("mean"), summary.get("empirical_mcse")


def _probability_limits(groups, getter, nominal):
    values = [nominal]
    for _, rows in groups:
        for row in rows:
            mean, se = getter(row)
            if _finite(mean):
                uncertainty = se if _finite(se) else 0.
                values.extend((100 * (mean - uncertainty), 100 * (mean + uncertainty)))
    padding = max(1., (max(values) - min(values)) * .15)
    return max(-1., min(values) - padding), min(101., max(values) + padding)


def _save_figure(figure, output, name, footer):
    figure.text(.5, .01, footer, ha="center", va="bottom", fontsize=8, color="#39434c")
    figure.tight_layout(rect=(0, .065, 1, .92))
    files = []
    for extension in ("png", "svg"):
        path = output / f"{name}.{extension}"
        figure.savefig(path, dpi=150, facecolor="white")
        files.append(path.name)
    return files


def render_figures(summary, output):
    """Render one point, law-diagnostic and paired-comparison figure per experiment."""
    plt = _plotting()
    files = []
    footer = ("Bars: ±1 MCSE, not confidence intervals or simultaneous bands. Cells are not pooled.\n"
              "Zero observed failure / MCSE does not imply zero population risk. Missing values remain gaps.")
    for experiment in summary["experiments"]:
        groups = _groups(experiment)
        figure, axes = plt.subplots(2, 2, figsize=(12, 8))
        figure.suptitle("Flat point and conditional law: distinct assessments\n" + _experiment_label(experiment), fontsize=12)
        _draw(axes[0, 0], groups, lambda row: _point(row, POINT_METRICS[0]))
        axes[0, 0].set(title="Point: bounded squared-log loss", ylabel="mean min[log(m̂/m)², log(2)²]")
        axes[0, 0].set_ylim(-.01, math.log(2) ** 2 + .025)
        _draw(axes[0, 1], groups, lambda row: _point(row, POINT_METRICS[1]), scale=100)
        axes[0, 1].set(title="Point: factor-two failure", ylabel="outside [m/2, 2m], including invalid (%)", ylim=(-3, 103))
        _draw(axes[1, 0], groups, lambda row: (row["flat_signed_log_bias"]["mean"], row["flat_signed_log_bias"]["mcse"]))
        axes[1, 0].axhline(0, color="black", linewidth=.8, linestyle=":")
        axes[1, 0].set(title="Point: signed log bias (finite-R diagnostic)", ylabel="mean log(m̂/m)")
        _draw(axes[1, 1], groups, lambda row: _interval(row, "flat_joint", "coverage"), scale=100)
        axes[1, 1].axhline(95, color="black", linestyle=":", label="nominal 95%")
        axes[1, 1].set(title="Flat conditional law: central 95% coverage", ylabel="truth inside interval (%)",
                       ylim=_probability_limits(groups, lambda row: _interval(row, "flat_joint", "coverage"), 95))
        axes[0, 0].legend(title="Fixed true mass", fontsize=8)
        files += _save_figure(figure, output, experiment["id"] + "-flat-assessment", footer)
        plt.close(figure)

        laws = list(dict.fromkeys(law for row in experiment["rows"] for law in row["laws"]))
        figure, axes = plt.subplots(len(laws), 3, figsize=(14, 3.4 * len(laws) + 1), squeeze=False)
        figure.suptitle("Conditional laws: 95% coverage and tails at fixed truth\n" + _experiment_label(experiment), fontsize=12)
        for index, law in enumerate(laws):
            for column, (metric, label, nominal) in enumerate((
                    ("coverage", "Coverage", 95), ("below_interval", "True mass below lower endpoint", 2.5),
                    ("above_interval", "True mass above upper endpoint", 2.5))):
                ax = axes[index, column]
                _draw(ax, groups, lambda row, law=law, metric=metric: _interval(row, law, metric), scale=100)
                ax.axhline(nominal, color="black", linestyle=":", linewidth=.9)
                ax.set(title=f"{law.replace('_joint', '').capitalize()} law: {label}", ylabel="experiments (%)",
                       ylim=_probability_limits(groups, lambda row, law=law, metric=metric: _interval(row, law, metric), nominal))
                if column == 0:
                    ax.legend(title="Fixed true mass", fontsize=8)
        files += _save_figure(figure, output, experiment["id"] + "-coverage-tails", footer)
        plt.close(figure)

        if "calibrated_joint" in laws and "oracle_joint" in laws:
            figure, axes = plt.subplots(len(groups), 3, figsize=(13, 3.2 * len(groups) + 1.2), squeeze=False)
            figure.suptitle("Calibration uncertainty: matched data, three information assumptions\n" + _experiment_label(experiment), fontsize=12)
            labels = {"flat_joint": "Calibration plug-in", "calibrated_joint": "Calibration integrated",
                      "oracle_joint": "True noise known (oracle)"}
            for index, (mass, rows) in enumerate(groups):
                for column, (level, metric, title) in enumerate((
                        ("50", "coverage", "50% coverage (%)"), ("95", "coverage", "95% coverage (%)"),
                        ("95", "mean_log_width", "95% mean log width"))):
                    ax = axes[index, column]
                    coordinates = [row["scenario"]["true_acceleration_snr"] for row in rows]
                    for color_index, (law, label) in enumerate(labels.items()):
                        values = [row["laws"][law]["intervals"][level].get(metric, {}) for row in rows]
                        scale = 100 if metric == "coverage" else 1
                        means = [scale * value["mean"] if _finite(value.get("mean")) else float("nan") for value in values]
                        color = f"C{color_index}"
                        ax.plot(coordinates, means, marker="o", markersize=3, label=label, color=color)
                        bars = [(x, mean, scale * value["empirical_mcse"]) for x, mean, value in zip(coordinates, means, values)
                                if _finite(mean) and _finite(value.get("empirical_mcse"))]
                        if bars:
                            xx, yy, se = zip(*bars)
                            ax.errorbar(xx, yy, yerr=se, fmt="none", capsize=2, color=color)
                    ax.set(title=f"Mass {mass:g}: {title}", xlabel="Per-reading acceleration SNR")
                    ax.set_xscale("symlog", linthresh=.5)
                    ax.set_xticks(coordinates, [f"{x:g}" for x in coordinates])
                    ax.grid(alpha=.2)
                    if metric == "coverage":
                        ax.axhline(float(level), color="black", linestyle=":", linewidth=.8)
                        ax.set_ylim(0, 103)
                    if index == 0 and column == 2:
                        ax.legend(fontsize=8)
            files += _save_figure(figure, output, experiment["id"] + "-calibration-comparison", footer)
            plt.close(figure)

        available = [name for name in summary["comparators"]
                     if any(row["paired_flat_minus_other"][name]["status"] == "available" for row in experiment["rows"])]
        if not available:
            continue
        figure, axes = plt.subplots(len(available), 2, figsize=(12, 2.7 * len(available)), squeeze=False)
        figure.suptitle("Paired point losses: flat ratio of means minus comparator\n" + _experiment_label(experiment), fontsize=12)
        for index, name in enumerate(available):
            for column, metric in enumerate(POINT_METRICS):
                def get_pair(row, name=name, metric=metric):
                    value = row["paired_flat_minus_other"][name].get("metrics", {}).get(metric, {})
                    return value.get("mean_difference"), value.get("mcse")
                ax = axes[index, column]
                _draw(ax, groups, get_pair, scale=100 if column else 1.)
                ax.axhline(0, color="black", linewidth=.9, linestyle=":")
                ax.set(title=f"vs {name}", ylabel=("Δ failure (pp)" if column else "Δ capped log² loss"))
        axes[0, 0].legend(title="Fixed true mass", fontsize=8)
        files += _save_figure(figure, output, experiment["id"] + "-paired-comparisons",
                             "Capped loss = min[log(m̂/m)², log(2)²]. Failure means outside [m/2, 2m], including invalid estimates.\n"
                             "Negative difference favors the flat point on this loss and cell only.\n" + footer)
        plt.close(figure)
    return files


def export_assessment(study_dir, output_dir, *, comparators=DEFAULT_COMPARATORS, figures=True):
    """Export into a new directory; never overwrite an earlier analysis or study."""
    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"Assessment output already exists: {output}")
    loaded = load_completed_study(study_dir)
    if output.is_relative_to(Path(loaded["study_dir"])):
        raise ValueError("Assessment output must be outside the source study directory")
    summary = summarize_study(loaded, comparators)
    # Validate plotting inputs/dependency before creating any output.
    if figures:
        plt = _plotting()
        summary["analysis_runtime"].update(matplotlib=plt.matplotlib.__version__,
                                            numpy=sys.modules["numpy"].__version__)
        for experiment in summary["experiments"]:
            _groups(experiment)
    output.mkdir(parents=True, exist_ok=False)
    summary["figures"] = render_figures(summary, output) if figures else []
    with (output / "summary.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def export_accuracy(study_dir, output_dir, *, tolerance_factors=(1.25, 1.5, 2.0), success_targets=(.90, .95), alpha=.05):
    """Verify and export accuracy.json to a new directory outside the study."""
    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"Accuracy output already exists: {output}")
    loaded = load_completed_study(study_dir)
    if output.is_relative_to(Path(loaded["study_dir"])):
        raise ValueError("Accuracy output must be outside the source study directory")
    summary = summarize_accuracy(loaded, tolerance_factors, success_targets, alpha)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "accuracy.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--comparators", nargs="+", default=list(DEFAULT_COMPARATORS))
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--accuracy-map", action="store_true", help="Export archived absolute-accuracy probabilities and simultaneous MC bounds; no figures")
    parser.add_argument("--tolerance-factors", nargs="+", type=float, default=[1.25, 1.5, 2.])
    parser.add_argument("--success-targets", nargs="+", type=float, default=[.90, .95])
    parser.add_argument("--alpha", type=float, default=.05)
    args = parser.parse_args()
    try:
        if args.accuracy_map:
            summary = export_accuracy(args.study, args.output, tolerance_factors=args.tolerance_factors,
                                      success_targets=args.success_targets, alpha=args.alpha)
        else:
            summary = export_assessment(args.study, args.output, comparators=args.comparators, figures=not args.no_figures)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"Assessment failed: {error}\n")
    print(json.dumps({"output": str(args.output.resolve()), "experiments": len(summary["experiments"]),
                      "figures": len(summary.get("figures", []))}))


if __name__ == "__main__":
    main()

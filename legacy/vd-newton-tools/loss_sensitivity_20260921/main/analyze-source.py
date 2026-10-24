"""Explore saved point-estimate error curves without fitting any estimator.

Uses the production integrity loaders and scoring functions. All analysis
outputs are scratch evidence; README.md remains the interpretation surface.
Run from the workspace root. An existing output directory is never reused.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "03_trace_and_evaluator/mass_estimator/investigation"
sys.path.insert(0, str(ENGINE))
import numpy as np
from assessment_summary import load_completed_study
from trial_archive import load_trial_manifest, load_trial_archive, _custom_losses
from comparison_scores import point_losses, paired_summary

METHODS = ("flat_joint", "flat_geometric", "flat_median", "norm_ratio",
           "norm_floor", "positive_profile")
FACTORS = (1.1, 1.25, 1.5, 2., 3., 5., 10., 20.)
ALPHA = .05


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stats(a, b):
    result = paired_summary(a, b)
    return {"reference_mean": float(np.mean(a)), "comparator_mean": float(np.mean(b)),
            "difference": result["mean_difference"], "mcse": result["mcse"]}


def survival_extrema(a, b):
    # Include both right-continuous jumps and the left limit at zero.
    # Invalid points have infinite absolute-log error and fail every finite cut.
    cuts = np.unique(np.concatenate(([0.], a[np.isfinite(a)], b[np.isfinite(b)])))
    fa = np.searchsorted(np.sort(a), cuts, side="right") / len(a)
    fb = np.searchsorted(np.sort(b), cuts, side="right") / len(b)
    difference = fb - fa  # SF_reference minus SF_comparator, lower is better.
    result = {}
    for key, index in (("minimum", int(np.argmin(difference))),
                       ("maximum", int(np.argmax(difference)))):
        cut = float(cuts[index])
        result[key] = {"log_tolerance": cut, "factor": math.exp(cut),
                       **stats((a > cut).astype(float), (b > cut).astype(float))}
    return result


def check_scoring():
    # Independently check the survival-integral identity, including an invalid
    # point, against the production bounded-loss calculation.
    points = np.array([1., .5, 2., 10., -1.])
    errors = point_losses(points, 1.)["abs_log"]
    for factor in FACTORS:
        c = math.log(factor)
        cuts = np.unique(np.concatenate(([0., c], errors[(errors > 0) & (errors < c)])))
        midpoint = (cuts[:-1] + cuts[1:]) / 2
        survival = np.mean(errors[:, None] > midpoint, axis=0)
        integral = float(np.sum(survival * np.diff(cuts ** 2)) / c ** 2)
        direct = float(np.mean(_custom_losses(points, 1., factor, factor)["capped_squared_log"]) / c ** 2)
        assert abs(integral - direct) < 2e-15
    # Check crossing-sign orientation and a tied/invalid tail.
    ex = survival_extrema(np.array([.1, 2., np.inf]), np.array([.5, .6, np.inf]))
    assert abs(ex["minimum"]["difference"] + 1/3) < 1e-15
    assert abs(ex["maximum"]["difference"] - 1/3) < 1e-15


def analyze(study, output):
    if output.exists():
        raise FileExistsError(output)
    check_scoring()
    loaded = load_completed_study(study)
    cell_count = sum(len(item["report"]["rows"]) for item in loaded["reports"])
    cdf_count = cell_count * len(METHODS)
    rows = []
    verified_comparisons = 0
    for experiment in loaded["reports"]:
        run = Path(loaded["study_dir"]) / Path(experiment["relative_results"]).parent
        manifest = load_trial_manifest(run)
        source_rows = {row["id"]: row for row in experiment["report"]["rows"]}
        assert set(source_rows) == {entry["id"] for entry in manifest["entries"]}
        for entry in manifest["entries"]:
            arrays = load_trial_archive(run, entry)
            source = source_rows[entry["id"]]
            scenario = entry["scenario"]
            truth, R = scenario["true_mass"], entry["samples"]
            epsilon = math.sqrt(math.log(2 * cdf_count / ALPHA) / (2 * R))
            points = {name: arrays["point__" + name] for name in METHODS}
            errors = {name: point_losses(values, truth)["abs_log"] for name, values in points.items()}
            reference = points["flat_joint"]
            row = {"experiment": experiment["id"], "id": entry["id"],
                   "scenario": scenario, "R": R, "cdf_epsilon": epsilon,
                   "uniform_difference_halfwidth": 2 * epsilon,
                   "archive_path": str(run / entry["path"]), "archive_sha256": entry["sha256"],
                   "invalid_counts": {name: int(np.count_nonzero(~np.isfinite(values) | (values <= 0)))
                                      for name, values in points.items()}, "comparisons": {}}
            for other in METHODS[1:]:
                extrema = survival_extrema(errors["flat_joint"], errors[other])
                grid = []
                for factor in FACTORS:
                    left = _custom_losses(reference, truth, factor, factor)
                    right = _custom_losses(points[other], truth, factor, factor)
                    cap = math.log(factor) ** 2
                    grid.append({"factor": factor,
                                 "outside_factor": stats(left["outside_factor"], right["outside_factor"]),
                                 "capped_squared_log": stats(left["capped_squared_log"], right["capped_squared_log"]),
                                 "normalized_capped_squared_log": stats(left["capped_squared_log"] / cap,
                                                                      right["capped_squared_log"] / cap)})
                    # Verify against saved canonical aggregates at all pre-existing metrics.
                    label = format(factor, "g")
                    for new_key, old_key in (("outside_factor", "outside_factor_" + label),
                                             ("capped_squared_log", "capped_squared_log_factor_" + label)):
                        saved = source["paired_flat_minus_other"][other].get(old_key)
                        if saved is not None:
                            assert abs(saved["mean_difference"] - grid[-1][new_key]["difference"]) < 1e-12
                            assert abs(saved["mcse"] - grid[-1][new_key]["mcse"]) < 1e-12
                            verified_comparisons += 1
                negative = extrema["minimum"]["difference"] < -2 * epsilon
                positive = extrema["maximum"]["difference"] > 2 * epsilon
                row["comparisons"][other] = {
                    "survival_extrema": extrema,
                    "crossing_resolved_by_uniform_band": bool(negative and positive),
                    "negative_region_resolved": bool(negative), "positive_region_resolved": bool(positive),
                    "grid": grid,
                }
            rows.append(row)
    report = {
        "kind": "exploratory_archived_point_loss_sensitivity", "study_dir": loaded["study_dir"],
        "checkpoint_sha256": loaded["checkpoint_sha256"], "methods": METHODS, "factors": FACTORS,
        "reference": "flat_joint", "cell_count": cell_count, "verified_saved_comparisons": verified_comparisons,
        "source_hashes": {str(path.relative_to(ROOT)): sha(path) for path in
                          [Path(__file__), ENGINE/"assessment_summary.py", ENGINE/"trial_archive.py", ENGINE/"comparison_scores.py"]},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "scope": "Known isotropic Gaussian SDs; fixed true masses and excitation; n1/n4 stable means; R10000. No inference, curve refitting, or pooling across truth cells.",
        "difference_definition": "reference loss minus comparator loss; negative favors reference",
        "failure_policy": "Invalid point outputs remain infinite absolute-log errors; fail every finite tolerance and receive the maximum capped loss.",
        "pointwise_mcse": "Sample SD of paired differences divided by sqrt(R); never an independent-means approximation.",
        "uniform_band": {
            "alpha": ALPHA, "cdf_count": cdf_count,
            "definition": "Union bound across all cell/method marginal empirical CDFs, using DKW epsilon=sqrt(log(2*cdf_count/alpha)/(2R)). Difference halfwidth=2epsilon. Valid simultaneously over all finite error thresholds. Shared trials/methods/cells do not invalidate the union bound.",
            "cap_extension": "E[min((abs(log(mhat/m))/log(K))^2,1)] is integral_0^logK 2t/logK^2 P(abs(log(mhat/m))>t)dt. Thus the same difference bound holds for every K>1. Raw capped-loss bound multiplies halfwidth by log(K)^2.",
            "invalid_extension": "The extended error variable including +infinity can be mapped monotonically to [0,1] via 2/pi*atan(error); DKW applies to that ordinary distribution, retaining the atom at 1.",
            "limitations": "Conservative bounds for Monte Carlo sampling only, conditional on the fixed numerical implementations and design. Do not cover model error, quadrature bias, unlisted scenarios, or generalization. Absence of a resolved crossing is not evidence of dominance. This is exploratory reuse after earlier results, not fresh confirmation.",
        },
        "rows": rows,
    }
    output.mkdir(parents=True)
    (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"cells": cell_count, "R": rows[0]["R"], "cdf_count": cdf_count,
                      "uniform_halfwidth": rows[0]["uniform_difference_halfwidth"],
                      "verified_saved_comparisons": verified_comparisons,
                      "resolved_crossings": sum(c["crossing_resolved_by_uniform_band"] for r in rows for c in r["comparisons"].values())}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.study, args.output)

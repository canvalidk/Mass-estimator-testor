"""Pretest ('switch') point rules, rescored from saved trials without refitting.

A switch rule reports the vector-length ratio |F|/|a| when the smaller
observed standardized norm, min(|F|/sigma_F, |a|/sigma_a), exceeds a threshold
tau, and a fallback point otherwise.  The threshold is a free parameter; the
decision uses the same noisy readings, so the reported point carries the
uncertainty of the branch choice.  Fallbacks are an archived point rule
(default ``flat_median``) or ``noise_ratio`` = sigma_F/sigma_a.  With unit
channel SDs the noise ratio equals the true mass wherever the two true SNRs
are equal, so its apparent gains there are an artefact of such grids.

An *oracle* choice is also reported: in each cell, the better of the two fixed
branches (|F|/|a| or the fallback), chosen with knowledge of the true cell.  It
is a reference point, not an upper bound: a switch chooses per reading, and
that choice can carry information, so a switch can beat both of its branches.

Each switch is compared with a reference rule (default ``flat_joint``) on
identical replicates, with a two-sided Bonferroni normal bound over the
declared family (cells x thresholds x fallbacks x factors).
"""

import argparse
import json
from pathlib import Path

import numpy as np

from archive_rescoring import (bonferroni_z, cell_record, failures, iter_cells, label, paired,
                               provenance, resolved, write_new_json)
from weighted_posterior import standardize


def switch_points(observed, norm_points, fallback_points, threshold):
    """|F|/|a| where the smaller observed standardized norm exceeds the threshold."""
    return np.where(observed > threshold, norm_points, fallback_points)


def _observed(arrays):
    x, y, s = standardize(arrays["observation__force"], arrays["observation__acceleration"],
                          arrays["noise__supplied_force_sd"], arrays["noise__supplied_acceleration_sd"])
    return np.minimum(np.linalg.norm(x, axis=1), np.linalg.norm(y, axis=1)), s


def analyze(run_dir, thresholds, fallbacks=("flat_median",), factors=(1.25, 1.5, 2.),
            reference="flat_joint", alpha=.05):
    cells = []
    for entry, arrays in iter_cells(run_dir):
        truth = float(entry["scenario"]["true_mass"])
        record = cell_record(entry)
        record["factors"] = {}
        observed, s = _observed(arrays)
        for factor in factors:
            ref_fail = failures(arrays["point__" + reference], truth, factor)
            norm_fail = failures(arrays["point__norm_ratio"], truth, factor)
            block = {"reference_success": float(1 - ref_fail.mean()), "switch": {}, "oracle": {}}
            for fallback in fallbacks:
                fallback_points = s if fallback == "noise_ratio" else arrays["point__" + fallback]
                fallback_fail = failures(fallback_points, truth, factor)
                best = min(norm_fail.mean(), fallback_fail.mean())
                block["oracle"][fallback] = {"success": float(1 - best),
                                             "gain_over_reference": float(ref_fail.mean() - best)}
                rows = {}
                for threshold in thresholds:
                    points = switch_points(observed, arrays["point__norm_ratio"], fallback_points, threshold)
                    fail = failures(points, truth, factor)
                    difference, mcse = paired(fail, ref_fail)
                    rows[label(threshold)] = {"success": float(1 - fail.mean()),
                                              "switch_minus_reference_failure": difference, "mcse": mcse}
                block["switch"][fallback] = rows
            record["factors"][label(factor)] = block
        cells.append(record)
    z = bonferroni_z(len(cells) * len(thresholds) * len(fallbacks) * len(factors), alpha)
    summary = {}
    for factor in factors:
        key = label(factor)
        summary[key] = {}
        for fallback in fallbacks:
            rows = {}
            for threshold in thresholds:
                t = label(threshold)
                gains = []  # positive = switch succeeds more often than the reference
                worse = better = 0
                for record in cells:
                    item = record["factors"][key]["switch"][fallback][t]
                    state = resolved(item["switch_minus_reference_failure"], item["mcse"], z)
                    item["resolved"] = state
                    better += state == -1
                    worse += state == 1
                    gains.append(-item["switch_minus_reference_failure"])
                rows[t] = {"cells_switch_clearly_better": int(better), "cells_switch_clearly_worse": int(worse),
                           "best_cell_gain": float(max(gains)), "worst_cell_gain": float(min(gains)),
                           "mean_cell_gain": float(np.mean(gains))}
            oracle = [record["factors"][key]["oracle"][fallback]["gain_over_reference"] for record in cells]
            summary[key][fallback] = {"thresholds": rows, "oracle_max_gain": float(max(oracle)),
                                      "oracle_cells_with_gain_over_half_point": int(sum(g > .005 for g in oracle))}
    return {"reference": reference, "thresholds": [float(t) for t in thresholds],
            "fallbacks": list(fallbacks), "factors": [float(f) for f in factors], "alpha": alpha,
            "family_size": len(cells) * len(thresholds) * len(fallbacks) * len(factors),
            "bonferroni_z": z, "summary": summary, "cells": cells}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--thresholds", nargs="+", type=float,
                        default=[1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.])
    parser.add_argument("--fallbacks", nargs="+", default=["flat_median", "noise_ratio"])
    parser.add_argument("--factors", nargs="+", type=float, default=[1.25, 1.5, 2.])
    parser.add_argument("--reference", default="flat_joint")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new file")
    try:
        result = analyze(args.run, args.thresholds, tuple(args.fallbacks), tuple(args.factors), args.reference)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    result = {"schema_version": 1, "kind": "switch_rules", **result,
              "provenance": provenance([args.run], ["switch_rules.py", "weighted_posterior.py",
                                                    "comparison_scores.py", "trial_archive.py"])}
    write_new_json(args.output, result)
    print(json.dumps({"status": "completed", "output": str(args.output)}))


if __name__ == "__main__":
    main()

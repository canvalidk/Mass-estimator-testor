"""Cells where one saved point rule beats every other saved rule.

For each cell and tolerance factor K the focal rule's failure indicator is
compared with every other archived point rule on identical replicates.  A
comparison is resolved when the mean paired failure difference lies beyond
z times its MCSE, with z the two-sided Bonferroni normal quantile over the
declared family (cells x factors x compared rules).  A cell is *dominated by
the focal rule* when every comparison is resolved in its favour.

This paired normal bound differs from the operating-range protocol's
simultaneous DKW band: that band bounds each marginal success rate and then
doubles the halfwidth for differences, so it cannot use the pairing and is
more conservative.  Neither says anything about unsampled cells, other noise
models or real data.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from archive_rescoring import (bonferroni_z, cell_record, failures, iter_cells, label, paired,
                               point_names, provenance, resolved, write_new_json)


def analyze(run_dir, focal="flat_joint", factors=(1.25, 1.5, 2.), exclude=(), alpha=.05):
    cells = []
    methods = None
    for entry, arrays in iter_cells(run_dir):
        names = point_names(arrays)
        if methods is None:
            methods = names
            if focal not in methods:
                raise ValueError(f"focal rule {focal} is not archived")
            others = [name for name in methods if name != focal and name not in exclude]
        elif names != methods:
            raise ValueError("cells archive different point rules")
        truth = float(entry["scenario"]["true_mass"])
        record = cell_record(entry)
        base = {factor: failures(arrays["point__" + focal], truth, factor) for factor in factors}
        record["factors"] = {}
        for factor in factors:
            comparisons = {}
            for other in others:
                other_fail = failures(arrays["point__" + other], truth, factor)
                difference, mcse = paired(base[factor], other_fail)
                comparisons[other] = {"success": float(1 - other_fail.mean()),
                                      "focal_minus_other_failure": difference, "mcse": mcse}
            record["factors"][label(factor)] = {"focal_success": float(1 - base[factor].mean()),
                                                "comparisons": comparisons}
        cells.append(record)
    z = bonferroni_z(len(cells) * len(factors) * len(others), alpha)
    summary = {}
    for factor in factors:
        key = label(factor)
        dominated = []
        for record in cells:
            block = record["factors"][key]
            states = {name: resolved(c["focal_minus_other_failure"], c["mcse"], z)
                      for name, c in block["comparisons"].items()}
            for name, state in states.items():
                block["comparisons"][name]["resolved"] = state
            block["focal_beats_all"] = all(state == -1 for state in states.values())
            closest = max(block["comparisons"].items(), key=lambda item: item[1]["focal_minus_other_failure"])
            block["closest_rule"] = closest[0]
            block["closest_margin_success"] = -closest[1]["focal_minus_other_failure"]
            block["closest_z"] = (-closest[1]["focal_minus_other_failure"] / closest[1]["mcse"]
                                  if closest[1]["mcse"] else None)
            if block["focal_beats_all"]:
                dominated.append([record["true_force_snr"], record["true_acceleration_snr"]])
        summary[key] = {"cells_where_focal_beats_all": dominated, "count": len(dominated)}
    return {"focal": focal, "compared_rules": others, "excluded_rules": list(exclude),
            "factors": [float(f) for f in factors], "alpha": alpha,
            "family_size": len(cells) * len(factors) * len(others), "bonferroni_z": z,
            "summary": summary, "cells": cells}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, type=Path, help="experiment attempt directory with trial archives")
    parser.add_argument("--output", required=True, type=Path, help="new JSON file")
    parser.add_argument("--focal", default="flat_joint")
    parser.add_argument("--exclude", nargs="*", default=[])
    parser.add_argument("--factors", nargs="+", type=float, default=[1.25, 1.5, 2.])
    parser.add_argument("--alpha", type=float, default=.05)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new file")
    try:
        result = analyze(args.run, args.focal, tuple(args.factors), tuple(args.exclude), args.alpha)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    result = {"schema_version": 1, "kind": "paired_dominance", **result,
              "provenance": provenance([args.run], ["paired_dominance.py", "comparison_scores.py", "trial_archive.py"])}
    write_new_json(args.output, result)
    print(json.dumps({"status": "completed", "output": str(args.output),
                      "summary": {k: v["count"] for k, v in result["summary"].items()}}))


if __name__ == "__main__":
    main()

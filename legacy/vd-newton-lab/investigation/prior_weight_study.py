"""Score prior-weighted flat laws on saved trials; fit and validate a parameter.

The factors of ``weighted_posterior`` are declared priors on latent mass.  A
value chosen by fitting is a parameter of the equation.  Here it is chosen by
an explicit criterion on named data and then evaluated on other data.

``sweep``      Ratio-of-means points for a grid of factors on every saved
               reading of one or more experiment attempts; tolerance success
               compared with all archived point rules.
``fit``        The grid value minimizing the criterion: mean over tolerance
               factors K of the largest regret over cells.
``summarize``  Per-parameter tables from one or more sweep files.

Regret in a cell is the best success among archived point rules minus the
weighted rule's success.  Archived rules include ``flat_joint``, so the
unweighted law has regret >= 0 and a weighted law can have negative regret.
``sech`` centres are standardized masses (units of sigma_F/sigma_a).
``lognormal_factor`` and ``lognormal_prior`` centres are multiples of each
cell's true mass: that oracle placement measures the cost of a prior that is
wrong by a stated factor.  It is not an estimator anyone can run.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from archive_rescoring import (cell_record, failures, iter_cells, label, point_names, provenance,
                               sha256, write_new_json)
from weighted_posterior import LogNormalFactor, LogNormalPrior, SechTilt, standardize, weighted_points


FAMILIES = ("sech", "lognormal_factor", "lognormal_prior")
_CHUNK = 2000


def parameter_grid(family, lams=(), centres=(1.,), widths=(), centre_factors=()):
    if family == "sech":
        return [{"family": family, "lam": float(l), "centre": float(c)} for c in centres for l in lams]
    if family in ("lognormal_factor", "lognormal_prior"):
        return [{"family": family, "width": float(w), "centre_factor": float(k)}
                for w in widths for k in centre_factors]
    raise ValueError(f"family must be one of {FAMILIES}")


def key(parameter):
    return "|".join([parameter["family"]] + [f"{name}={label(parameter[name])}"
                                             for name in sorted(parameter) if name != "family"])


def _factor(parameter, truth, s):
    if parameter["family"] == "sech":
        return SechTilt(parameter["lam"], parameter["centre"])
    centre = parameter["centre_factor"] * truth / s
    kind = LogNormalFactor if parameter["family"] == "lognormal_factor" else LogNormalPrior
    return kind(centre, parameter["width"])


def _subset(count, half):
    if half == "all":
        return slice(0, count)
    if half == "first":
        return slice(0, count // 2)
    if half == "second":
        return slice(count // 2, count)
    raise ValueError("half must be all, first or second")


def sweep(run_dirs, parameters, factors=(1.25, 1.5, 2.), half="all", cells=None):
    output = []
    for run_dir in run_dirs:
        for entry, arrays in iter_cells(run_dir, cells):
            truth = float(entry["scenario"]["true_mass"])
            part = _subset(len(arrays["observation__force"]), half)
            x, y, s = standardize(arrays["observation__force"][part], arrays["observation__acceleration"][part],
                                  arrays["noise__supplied_force_sd"][part],
                                  arrays["noise__supplied_acceleration_sd"][part])
            points = np.empty((len(parameters), len(x)))
            # Parameters sharing a centre share one panel set; different
            # centres get their own panels instead of one very large union.
            groups = {}
            for index, parameter in enumerate(parameters):
                centre = parameter.get("centre", parameter.get("centre_factor"))
                groups.setdefault((parameter["family"], centre), []).append(index)
            for start in range(0, len(x), _CHUNK):
                block = slice(start, start + _CHUNK)
                for indices in groups.values():
                    factors_ = [_factor(parameters[i], truth, s[block]) for i in indices]
                    points[indices, block] = weighted_points(x[block], y[block], factors_) * s[block]
            record = cell_record(entry)
            record.update(run=str(run_dir), half=half, readings=len(x))
            archived = {name: arrays["point__" + name][part] for name in point_names(arrays)}
            flat = archived["flat_joint"]
            ratio = flat / truth
            record["flat_joint_ratio_quartiles"] = np.quantile(ratio, [.25, .5, .75]).tolist()
            record["success"] = {}
            for factor in factors:
                k = label(factor)
                record["success"][k] = {
                    "archived": {name: float(1 - failures(values, truth, factor).mean())
                                 for name, values in archived.items()},
                    "weighted": {key(p): float(1 - failures(points[i], truth, factor).mean())
                                 for i, p in enumerate(parameters)}}
            unweighted = [i for i, p in enumerate(parameters)
                          if p["family"] == "sech" and p["lam"] == 0.]
            if unweighted:
                record["max_relative_deviation_lam0_vs_archived_flat"] = float(
                    np.max(np.abs(points[unweighted[0]] / flat - 1)))
            output.append(record)
    return output


def summarize(cells, parameters, factors):
    table = {}
    for parameter in parameters:
        name = key(parameter)
        table[name] = {"parameter": parameter}
        for factor in factors:
            k = label(factor)
            regrets, gains, beats = [], [], 0
            for cell in cells:
                success = cell["success"][k]
                best = max(success["archived"].values())
                value = success["weighted"][name]
                regrets.append((best - value, [cell["true_force_snr"], cell["true_acceleration_snr"]]))
                gains.append(value - success["archived"]["flat_joint"])
                beats += value > best
            worst = max(regrets, key=lambda item: item[0])
            table[name][k] = {"max_regret": worst[0], "max_regret_cell": worst[1],
                              "mean_success": float(np.mean([c["success"][k]["weighted"][name] for c in cells])),
                              "mean_gain_over_flat": float(np.mean(gains)),
                              "worst_cell_gain_over_flat": float(min(gains)),
                              "best_cell_gain_over_flat": float(max(gains)),
                              "cells_above_every_archived_rule": int(beats)}
        table[name]["criterion_mean_max_regret"] = float(np.mean([table[name][label(f)]["max_regret"]
                                                                  for f in factors]))
    return table


def _load_sweeps(paths, min_snr=None):
    cells, parameters, factors = [], None, None
    for path in paths:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data.get("kind") != "prior_weight_sweep":
            raise ValueError(f"{path} is not a prior-weight sweep")
        if parameters is None:
            parameters, factors = data["parameters"], data["factors"]
        elif data["parameters"] != parameters or data["factors"] != factors:
            raise ValueError("sweeps must share parameters and factors")
        cells.extend(c for c in data["cells"] if min_snr is None
                     or min(c["true_force_snr"], c["true_acceleration_snr"]) >= min_snr)
    if not cells:
        raise ValueError("no cells remain after the SNR filter")
    return cells, parameters, factors


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("sweep")
    run.add_argument("--run", required=True, action="append", type=Path)
    run.add_argument("--family", required=True, choices=FAMILIES)
    run.add_argument("--lams", nargs="*", type=float, default=[])
    run.add_argument("--centres", nargs="*", type=float, default=[1.])
    run.add_argument("--widths", nargs="*", type=float, default=[])
    run.add_argument("--centre-factors", nargs="*", type=float, default=[])
    run.add_argument("--factors", nargs="+", type=float, default=[1.25, 1.5, 2.])
    run.add_argument("--half", default="all", choices=("all", "first", "second"))
    run.add_argument("--cells", nargs="+")
    run.add_argument("--output", required=True, type=Path)
    for name in ("fit", "summarize"):
        sub = commands.add_parser(name)
        sub.add_argument("--sweep", required=True, nargs="+", type=Path)
        sub.add_argument("--output", required=True, type=Path)
        sub.add_argument("--min-snr", type=float,
                         help="keep only cells whose smaller true effective SNR is at least this value")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new file")
    try:
        if args.command == "sweep":
            parameters = parameter_grid(args.family, args.lams, args.centres, args.widths, args.centre_factors)
            if not parameters:
                parser.error("the parameter grid is empty")
            cells = sweep(args.run, parameters, tuple(args.factors), args.half, args.cells)
            result = {"schema_version": 1, "kind": "prior_weight_sweep", "parameters": parameters,
                      "factors": args.factors, "half": args.half, "cells": cells,
                      "provenance": provenance(args.run, ["prior_weight_study.py", "weighted_posterior.py",
                                                          "comparison_posterior.py", "comparison_scores.py",
                                                          "trial_archive.py"])}
        else:
            cells, parameters, factors = _load_sweeps(args.sweep, args.min_snr)
            table = summarize(cells, parameters, factors)
            result = {"schema_version": 1, "kind": f"prior_weight_{args.command}",
                      "sweeps": [str(p) for p in args.sweep], "factors": factors, "cells": len(cells),
                      "min_snr": args.min_snr,
                      "provenance": {"sweep_sha256": {str(p): sha256(p) for p in args.sweep},
                                     **provenance([], ["prior_weight_study.py"])}}
            if args.command == "fit":
                chosen = min(table.values(), key=lambda row: (row["criterion_mean_max_regret"],
                                                              abs(row["parameter"].get("lam", 0.))))
                result.update(criterion="mean over factors of the largest cell regret; ties prefer smaller |lam|",
                              selected=chosen["parameter"], selected_row=chosen)
            result["table"] = table
    except (ValueError, OSError) as error:
        parser.error(str(error))
    write_new_json(args.output, result)
    print(json.dumps({"status": "completed", "output": str(args.output),
                      "selected": result.get("selected")}))


if __name__ == "__main__":
    main()

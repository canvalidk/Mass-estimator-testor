"""Run a preset: generate each cell's data once, pass it to every estimator, score.

All estimators in a cell see the same readings, so their scores are paired.
Series i of a cell is drawn from SeedSequence(seed, spawn_key=(cell_key, i)),
where cell_key is a hash of the cell's world settings. So any series can be
regenerated from the preset; a quick run with fewer replicates sees the first
series of the full run; and adding a grid value does not change other cells.
"""

import datetime as _dt
import hashlib
import json
import math
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import yaml

from . import preset as preset_mod
from .scores import add_regret, score_interval, score_point, worst_regret
from .world import cell_key, generate, series_rngs

PACKAGE = Path(__file__).resolve().parent


def code_version():
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in sorted(PACKAGE.glob("*.py"))}
    commit, dirty = None, None
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PACKAGE, capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=PACKAGE,
                                    capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    return {"git_commit": commit, "uncommitted_changes": dirty, "file_sha256_16": files,
            "python": platform.python_version(), "numpy": np.__version__}


def run_cell(cell, estimators, replicates, seed, factors, batch_elements, progress=None):
    readings, truth = generate(cell["world"], series_rngs(seed, cell["world"], replicates))
    n = readings.readings
    values, seconds = {}, {}
    for est in estimators:
        start = time.perf_counter()
        per_reading = 1 if (est.direct or est.spec.get("reduce") == "pool_pair") else n
        grid = est.numerics["grid_points"] if est.direct is None else 1
        batch = max(1, int(batch_elements // max(1, per_reading * grid)))
        parts = [est.evaluate(readings.subset(slice(i, i + batch))) for i in range(0, replicates, batch)]
        values[est.name] = {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}
        seconds[est.name] = time.perf_counter() - start
        if progress:
            progress(f"  {cell['id']} {est.name}: {seconds[est.name]:.1f}s")
    scores = {"points": {}, "intervals": {}}
    for est in estimators:
        v = values[est.name]
        for name in est.point_readouts():
            scores["points"][f"{est.name}.{name}"] = score_point(v[name], truth["mass"], factors)
        for name in est.interval_readouts():
            scores["intervals"][f"{est.name}.{name}"] = score_interval(v[name], truth["mass"])
        scores.setdefault("unresolved", {})[est.name] = int(np.sum(v["unresolved"]))
    add_regret(scores, factors)
    return {"id": cell["id"], "axes": cell["axes"], "world": cell["world"],
            "cell_key": cell_key(cell["world"]), "scores": scores,
            "seconds": seconds}, values


def run(path, out_root="results", replicates=None, progress=print):
    spec = preset_mod.load(path)
    raw = spec["raw"]
    reps = replicates or raw["replicates"]
    factors = raw["scores"]["within_factor"]
    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    out = Path(out_root) / raw["name"] / stamp
    suffix = 1
    while out.exists():
        suffix += 1
        out = Path(out_root) / raw["name"] / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    (out / "preset.yaml").write_text(Path(path).read_text(encoding="utf-8"), encoding="utf-8")
    cells, arrays = [], {}
    start = time.perf_counter()
    for cell in spec["cells"]:
        if progress:
            progress(f"cell {cell['id']} {cell['axes']}")
        result, values = run_cell(cell, spec["estimators"], reps, raw["seed"], factors,
                                  spec["numerics"]["batch_elements"], progress)
        cells.append(result)
        for est, readouts in values.items():
            for name, arr in readouts.items():
                arrays[f"{cell['id']}/{est}/{name}"] = arr
    record = {"preset": raw["name"], "question": raw["question"], "preset_file": str(path),
              "replicates": reps, "replicates_override": replicates is not None,
              "seed": raw["seed"], "run_numerics": spec["numerics"],
              "estimator_numerics": {e.name: getattr(e, "numerics", None) for e in spec["estimators"]},
              "code": code_version(), "started": stamp,
              "seconds": time.perf_counter() - start, "cells": cells,
              "worst_regret": worst_regret(cells, factors)}
    (out / "run.json").write_text(json.dumps(record, indent=1, default=_json_default), encoding="utf-8")
    np.savez_compressed(out / "values.npz", **arrays)
    (out / "summary.md").write_text(summary_markdown(record, factors), encoding="utf-8")
    if progress:
        progress(f"wrote {out}")
    return out, record


def _json_default(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    raise TypeError(type(value))


def _pct(x):
    return "—" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{100 * x:.1f}"


def summary_markdown(record, factors):
    lines = [f"# {record['preset']}", "", f"**Question.** {record['question']}", "",
             f"Replicates per cell: {record['replicates']}"
             + (" (overridden on the command line)" if record["replicates_override"] else "")
             + f". Seed {record['seed']}. Code {record['code']['git_commit'] or 'uncommitted'}"
             + (" with uncommitted changes" if record["code"]["uncommitted_changes"] else "") + ".", ""]
    for cell in record["cells"]:
        axes = ", ".join(f"{k} = {v}" for k, v in cell["axes"].items()) or "single cell"
        lines += [f"## {cell['id']}: {axes}", ""]
        head = "| readout | " + " | ".join(f"within ×{k} (%)" for k in factors) + \
               " | median m̂/m | IQR m̂/m | invalid |"
        lines += [head, "|" + "---|" * (len(factors) + 4)]
        for name, s in cell["scores"]["points"].items():
            q = s["ratio_quartiles"]
            within = " | ".join(f"{_pct(s['within_factor'][str(k)]['value'])} ± {_pct(s['within_factor'][str(k)]['mcse'])}"
                                for k in factors)
            lines.append(f"| {name} | {within} | {q[1]:.3g} | {q[0]:.3g}–{q[2]:.3g} | {s['invalid']} |")
        if cell["scores"]["intervals"]:
            lines += ["", "| interval | coverage (%) | truth below (%) | truth above (%) | median log width | invalid |",
                      "|---|---|---|---|---|---|"]
            for name, s in cell["scores"]["intervals"].items():
                lines.append(f"| {name} | {_pct(s['coverage']['value'])} ± {_pct(s['coverage']['mcse'])} | "
                             f"{_pct(s['truth_below']['value'])} | {_pct(s['truth_above']['value'])} | "
                             f"{s['median_log_width']:.3g} | {s['invalid']} |")
        lines.append("")
    lines += ["## Worst-case regret across cells", "",
              "Regret = best within-factor success among all point readouts in the cell, minus this readout's.", "",
              "| readout | " + " | ".join(f"×{k}: pp (cell)" for k in factors) + " |",
              "|" + "---|" * (len(factors) + 1)]
    for name, per in record["worst_regret"].items():
        lines.append(f"| {name} | " + " | ".join(
            f"{100 * per[str(k)]['value']:.1f} ({per[str(k)]['cell']})" for k in factors) + " |")
    return "\n".join(lines) + "\n"

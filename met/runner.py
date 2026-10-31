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
import multiprocessing
import os
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import yaml

from . import preset as preset_mod
from .scores import add_regret, agreement, score_interval, score_point, worst_regret
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


# ratio_of_means (ours24's plain sum of sizes) and ratio_of_means_rss (ours25b's stacked
# length): when a reference has only one of them, the other is compared against it.
RATIO_OF_MEANS_KINDS = ("ratio_of_means", "ratio_of_means_rss")


def run_cell(cell, estimators, replicates, seed, factors, batch_elements, progress=None, compare_to=None):
    readings, truth = generate(cell["world"], series_rngs(seed, cell["world"], replicates))
    n = readings.readings
    values, seconds = {}, {}
    for est in estimators:
        start = time.perf_counter()
        per_reading = 1 if (est.direct or est.spec.get("reduce") == "pool_pair") else n
        if est.spec.get("reduce") == "split_axes":
            per_reading = n * readings.dimension
        grid = est.numerics["grid_points"] if est.numerics else 1
        batch = max(1, int(batch_elements // max(1, per_reading * grid)))
        parts = []
        for i in range(0, replicates, batch):
            part = slice(i, i + batch)
            given = {"mass": truth["mass"], "true_acceleration": truth["true_acceleration"][part]}
            parts.append(est.evaluate(readings.subset(part), given, cell["world"]))
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
    if compare_to:
        ref = values[compare_to]
        scores["agreement"] = {}
        for est in estimators:
            if est.name == compare_to:
                continue
            for name in est.point_readouts() + est.interval_readouts():
                key = name if name in ref else None
                if key is None and name in RATIO_OF_MEANS_KINDS:
                    # the two poolings of the ratio of means are one quantity, compared across equations
                    key = next((k for k in RATIO_OF_MEANS_KINDS if k in ref), None)
                if key is not None:
                    scores["agreement"][f"{est.name}.{name}"] = agreement(values[est.name][name], ref[key])
    return {"id": cell["id"], "axes": cell["axes"], "world": cell["world"],
            "cell_key": cell_key(cell["world"]), "scores": scores,
            "seconds": seconds}, values


# ---------------------------------------------------------------- running cells side by side
#
# Cells are independent (each draws its series from its own seed, see the module docstring),
# so they can run in separate processes and the output is the same bit for bit: each worker
# loads the preset itself, runs whole cells, and the parent puts the results back in preset
# order. Each worker's numpy is held to one thread, so N workers use N cores, not N times
# numpy's own threads. "spawn" is used everywhere because it is what Windows has; a script
# that calls run() with workers > 1 needs the usual `if __name__ == "__main__":` guard.

_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def cpu_count():
    """Cores this process may use."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def default_workers():
    """Half the cores, at least one: leaves the machine usable while a study runs."""
    return max(1, cpu_count() // 2)


def resolve_workers(workers):
    """None -> default_workers(); "max" -> every core; an integer >= 1 as given."""
    if workers is None:
        return default_workers()
    if workers == "max":
        return cpu_count()
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError(f"workers must be a positive integer or 'max', got {workers!r}")
    return workers


_WORKER = {}


def _init_worker(path, args):
    _WORKER["spec"] = preset_mod.load(path)
    _WORKER["args"] = args


def _run_cell_in_worker(index):
    spec, (reps, seed, factors, batch_elements, compare_to) = _WORKER["spec"], _WORKER["args"]
    lines = []
    result, values = run_cell(spec["cells"][index], spec["estimators"], reps, seed, factors,
                              batch_elements, lines.append, compare_to)
    return index, result, values, lines


def _run_cells_parallel(path, spec, args, workers, progress):
    n = len(spec["cells"])
    saved = {k: os.environ.get(k) for k in _THREAD_VARS}
    os.environ.update({k: "1" for k in _THREAD_VARS})          # inherited by the workers only:
    try:                                                        # this process's numpy is loaded
        pool = ProcessPoolExecutor(max_workers=min(workers, n), mp_context=multiprocessing.get_context("spawn"),
                                   initializer=_init_worker, initargs=(str(path), args))
        out = [None] * n
        try:
            futures = [pool.submit(_run_cell_in_worker, i) for i in range(n)]
            for done, future in enumerate(as_completed(futures), 1):
                index, result, values, lines = future.result()
                out[index] = (result, values)
                if progress:
                    cell = spec["cells"][index]
                    progress(f"cell {cell['id']} {cell['axes']} ({done}/{n})")
                    for line in lines:
                        progress(line)
        except BaseException:                   # Ctrl-C or a failed cell: stop the rest now
            _stop_workers(pool)                 # first: shutdown() forgets the processes
            pool.shutdown(wait=False, cancel_futures=True)
            raise
        pool.shutdown(wait=True)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    return out


def _stop_workers(pool):
    """Terminate the pool's processes without waiting for the cells they are running."""
    stop = getattr(pool, "terminate_workers", None)              # Python 3.14+
    if stop is not None:
        stop()
        return
    procs = list((getattr(pool, "_processes", None) or {}).values())
    for proc in procs:
        if proc.is_alive():
            proc.terminate()
    for proc in procs:
        proc.join(timeout=5)


def run(path, out_root="results", replicates=None, progress=print, workers=None):
    """Run a preset and write run.json, values.npz and summary.md under out_root.

    workers: processes running cells side by side; None (default) is half the cores,
    "max" every core, 1 the plain serial loop. The output is the same for any value."""
    spec = preset_mod.load(path)
    if spec.get("kind") == "robot_budget":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("robot_budget is a serial calibration collection; omit --replicates and use --workers 1")
        from .robot_budget import run as run_robot_budget
        return run_robot_budget(path, out_root, progress=progress)
    if spec.get("kind") == "composition_coverage":
        from .composition_coverage import run as run_composition_coverage
        return run_composition_coverage(path, out_root, replicates=replicates,
                                        progress=progress, workers=workers)
    if spec.get("kind") == "readout_candidates":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("readout_candidates is a deterministic serial study; omit --replicates and use --workers 1")
        from .readout_study import run as run_readout_candidates
        return run_readout_candidates(path, out_root, progress=progress)
    if spec.get("kind") == "anisotropic_world":
        if not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("anisotropic_world is a serial research study; use --workers 1")
        from .anisotropic_world import run as run_anisotropic_world
        return run_anisotropic_world(path, out_root, replicates=replicates, progress=progress)
    if spec.get("kind") == "shared_tilt":
        if not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("shared_tilt is a serial research study; use --workers 1")
        from .shared_tilt import run as run_shared_tilt
        return run_shared_tilt(path, out_root, replicates=replicates, progress=progress)
    if spec.get("kind") == "pooled_scale":
        if not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("pooled_scale is a serial research study; use --workers 1")
        from .pooled_scale import run as run_pooled_scale
        return run_pooled_scale(path, out_root, replicates=replicates, progress=progress)
    if spec.get("kind") == "prior_path":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("prior_path is a deterministic serial study; omit --replicates and use --workers 1")
        from .prior_path import run as run_prior_path
        return run_prior_path(path, out_root, progress)
    if spec.get("kind") == "noise_path":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("noise_path is a deterministic serial study; omit --replicates and use --workers 1")
        from .noise_path import run as run_noise_path
        return run_noise_path(path, out_root, progress)
    if spec.get("kind") == "inspection":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("inspection is a deterministic serial study; omit --replicates and use --workers 1")
        from .inspection.run import run as run_inspection
        return run_inspection(path, out_root, progress=progress)
    if spec.get("kind") == "misalignment":
        if replicates is not None or not (workers is None or (type(workers) is int and workers == 1)):
            raise ValueError("misalignment is a deterministic serial study; omit --replicates and use --workers 1")
        from .misalignment import run as run_misalignment
        return run_misalignment(path, out_root, progress)
    workers = resolve_workers(workers)
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
    args = (reps, raw["seed"], factors, spec["numerics"]["batch_elements"], raw["scores"].get("compare_to"))
    workers = min(workers, len(spec["cells"]))
    if workers > 1:
        if progress:
            progress(f"{len(spec['cells'])} cells on {workers} workers ({cpu_count()} cores)")
        finished = _run_cells_parallel(path, spec, args, workers, progress)
    else:
        finished = []
        for cell in spec["cells"]:
            if progress:
                progress(f"cell {cell['id']} {cell['axes']}")
            finished.append(run_cell(cell, spec["estimators"], *args[:4], progress, args[4]))
    for cell, (result, values) in zip(spec["cells"], finished):
        cells.append(result)
        for est, readouts in values.items():
            for name, arr in readouts.items():
                arrays[f"{cell['id']}/{est}/{name}"] = arr
    record = {"preset": raw["name"], "question": raw["question"], "preset_file": str(path),
              "replicates": reps, "replicates_override": replicates is not None,
              "seed": raw["seed"], "run_numerics": spec["numerics"],
              "estimator_numerics": {e.name: getattr(e, "numerics", None) for e in spec["estimators"]},
              "code": code_version(), "started": stamp, "workers": workers,
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
        if cell["scores"].get("agreement"):
            lines += ["", "Distance from the reference estimator, series by series (|log ratio|; intervals: worse endpoint):", "",
                      "| readout | median | 90th pct | 99th pct | max | within 1% | within 5% |",
                      "|---|---|---|---|---|---|---|"]
            for name, s in cell["scores"]["agreement"].items():
                if "median" not in s:
                    lines.append(f"| {name} | — | — | — | — | — | — |")
                    continue
                lines.append(f"| {name} | {s['median']:.2g} | {s['p90']:.2g} | {s['p99']:.2g} | {s['max']:.2g} | "
                             f"{_pct(s['within_1pct'])} | {_pct(s['within_5pct'])} |")
        lines.append("")
    lines += ["## Worst-case regret across cells", "",
              "Regret = best within-factor success among all point readouts in the cell, minus this readout's.", "",
              "| readout | " + " | ".join(f"×{k}: pp (cell)" for k in factors) + " |",
              "|" + "---|" * (len(factors) + 1)]
    for name, per in record["worst_regret"].items():
        lines.append(f"| {name} | " + " | ".join(
            f"{100 * per[str(k)]['value']:.1f} ({per[str(k)]['cell']})" for k in factors) + " |")
    return "\n".join(lines) + "\n"

"""O14's pooled residual-scale check; the production mass equation is unchanged.

This recipe compares repeated-world moments with an independent noncentral-Gram
reference. It estimates a variance multiplier relative to supplied covariance;
Bartlett eta is a different quantity and is retained only for paired spread checks.
"""
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy
import yaml

from .analyst import bartlett_weight
from .world import generate, series_rngs

_UNWEIGHTED = {}
_OURS29 = {"misfit": "kink"}


def residual_statistics(readings):
    """Observe supplied readings only; no world or latent truth enters this function.

    Q is profiled separately for every reading on the whole projective line.
    It is not the residual at a common fitted mass or a positive-domain fit.
    """
    _, n, k = readings.force.shape
    if k not in (2, 3) or n < 1 or not readings.isotropic:
        raise ValueError("pooled residual statistics need nonempty isotropic 2D or 3D readings")
    p, q, dot = readings.P, readings.Q, readings.D
    gap = np.hypot(p - q, 2 * dot)
    per_reading = np.maximum(0.0, 0.5 * (p + q - gap))
    q_sum = per_reading.sum(axis=1)
    pooled = q_sum / (n * (k - 1))
    precision = np.divide(1.0, pooled, out=np.full_like(pooled, np.inf), where=pooled > 0)
    # Module-lived dictionaries avoid reusing an id in ReadingSet's analyst cache.
    return {
        "q_sum": q_sum,
        "pooled_variance": pooled,
        "raw_pooled_precision": precision,
        "capped_pooled_precision": np.minimum(1.0, precision),
        "eta_unweighted": bartlett_weight(readings, _UNWEIGHTED).copy(),
        "eta_ours29": bartlett_weight(readings, _OURS29).copy(),
    }


def validate(data):
    from .preset import expand

    required = {"kind", "name", "question", "world", "replicates", "seed"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"pooled_scale needs {sorted(required)}")
    if set(data) - required - {"notes", "numerics"}:
        raise ValueError("unknown pooled_scale setting")
    if data["kind"] != "pooled_scale":
        raise ValueError("expected kind: pooled_scale")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in name)):
        raise ValueError("pooled_scale name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("pooled_scale needs a nonempty question")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("pooled_scale seed must be a nonnegative integer")
    if type(data["replicates"]) is not int or data["replicates"] < 2:
        raise ValueError("pooled_scale replicates must be an integer >= 2")
    settings = {"reference_nodes": 32}
    supplied = data.get("numerics", {})
    if not isinstance(supplied, dict) or set(supplied) - settings.keys():
        raise ValueError("pooled_scale numerics accepts reference_nodes only")
    settings.update(supplied)
    nodes = settings["reference_nodes"]
    if type(nodes) is not int or not 8 <= nodes <= 64:
        raise ValueError("reference_nodes must be an integer from 8 to 64")
    if not isinstance(data["world"], dict):
        raise ValueError("pooled_scale world must be a mapping")
    cells = expand(data["world"])
    for cell in cells:
        world = cell["world"]
        if world["dimension"] not in (2, 3):
            raise ValueError("pooled_scale requires dimension 2 or 3; the 1D denominator is zero")
        if isinstance(world["acceleration_snr"], (bool, dict)):
            raise ValueError("pooled_scale reference needs a fixed scalar acceleration_snr")
        sf, sa = world["noise"]["force_sd"], world["noise"]["acceleration_sd"]
        if isinstance(sf, dict) or isinstance(sa, dict) or sf != sa:
            raise ValueError("pooled_scale reference needs equal isotropic channel SDs")
        if world["supplied_noise"] != "exact":
            raise ValueError("pooled_scale reference compares against exact supplied covariance")
    return {"kind": "pooled_scale", "raw": data, "cells": cells, "numerics": settings}


def _summary(values):
    finite = bool(np.all(np.isfinite(values)))
    return {"mean": float(np.mean(values)) if finite else None,
            "sd": float(np.std(values, ddof=1)) if finite else None,
            "minimum": float(np.min(values)) if finite else None,
            "maximum": float(np.max(values)) if finite else None,
            "finite_count": int(np.isfinite(values).sum())}


def _correlation(left, right):
    if np.std(left) == 0 or np.std(right) == 0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def _markdown(record):
    lines = ["# Pooled residual scale: bias and repeated-world spread", "", record["question"], "",
        "The pooled ratio estimates a variance multiplier relative to the supplied covariance (target 1 here).",
        "Bartlett eta tempers a working mass law; it does not estimate that multiplier. Their spreads have different targets.",
        "Q is minimized separately per reading on the whole projective line. All runs are serial.",
        "Independent reference error estimates are differences between quadrature resolutions, not rigorous bounds.", "",
        "| Cell | k | Mass | SNR | N | Pooled mean | Reference mean | Pooled SD | Reference SD | eta unweighted mean / SD | eta ours29 mean / SD |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in record["cells"]:
        w, s, ref = row["world"], row["summary"], row["reference"]
        lines.append(f'| {row["id"]} | {w["dimension"]} | {w["mass"]:g} | {w["acceleration_snr"]:g} | '
            f'{w["readings"]} | {s["pooled_variance"]["mean"]:.7g} | {ref["pooled_mean"]:.7g} | '
            f'{s["pooled_variance"]["sd"]:.7g} | {ref["pooled_sd"]:.7g} | '
            f'{s["eta_unweighted"]["mean"]:.7g} / {s["eta_unweighted"]["sd"]:.7g} | '
            f'{s["eta_ours29"]["mean"]:.7g} / {s["eta_ours29"]["sd"]:.7g} |')
    return "\n".join(lines) + "\n"


def run(path, out_root="results", replicates=None, progress=print):
    """Run the saved recipe serially, retaining every replicate for later checks."""
    from .pooled_reference import moments
    from .runner import code_version

    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    reps = spec["raw"]["replicates"] if replicates is None else replicates
    if type(reps) is not int or reps < 2:
        raise ValueError("pooled_scale replicates must be an integer >= 2")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base = Path(out_root) / spec["raw"]["name"]
    out, suffix = base / stamp, 1
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    if replicates is None:
        (out / "preset.yaml").write_bytes(source)
    else:
        # The saved effective recipe must replay the override without CLI history.
        effective = dict(spec["raw"], replicates=reps)
        (out / "source_preset.yaml").write_bytes(source)
        (out / "preset.yaml").write_text(yaml.safe_dump(effective, sort_keys=False), encoding="utf-8")
    start, arrays, results = time.perf_counter(), {}, []
    nodes = spec["numerics"]["reference_nodes"]
    for cell in spec["cells"]:
        if progress:
            progress(f'cell {cell["id"]} {cell["axes"]}')
        w = cell["world"]
        readings, _truth = generate(w, series_rngs(spec["raw"]["seed"], w, reps))
        values = residual_statistics(readings)
        for name, value in values.items():
            arrays[f'{cell["id"]}/{name}'] = value
        k, n, r, mass = w["dimension"], w["readings"], w["acceleration_snr"], w["mass"]
        lam = (1 + mass ** 2) * r ** 2
        ref = moments(k, lam, nodes=nodes)
        reference = {"noncentrality": float(lam), "mean_Q": float(ref["mean"]),
            "variance_Q": float(ref["variance"]), "pooled_mean": float(ref["mean"] / (k - 1)),
            "pooled_sd": math.sqrt(max(0, ref["variance"]) / n) / (k - 1),
            "mean_error_estimate": float(ref["mean_error_estimate"] / (k - 1)),
            "variance_error_estimate": float(ref["variance_error_estimate"] / (n * (k - 1) ** 2)),
            "eta_unweighted_asymptotic": float(lam / (lam + k))}
        results.append({**cell, "summary": {name: _summary(value) for name, value in values.items()},
            "reference": reference,
            "paired_correlations": {
                "pooled_variance_eta_unweighted": _correlation(values["pooled_variance"], values["eta_unweighted"]),
                "pooled_variance_eta_ours29": _correlation(values["pooled_variance"], values["eta_ours29"])}})
    code = code_version()
    code.update(scipy=scipy.__version__, preset_sha256=hashlib.sha256(source).hexdigest(),
                effective_preset_sha256=hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest())
    record = {"kind": "pooled_scale", "preset": spec["raw"]["name"], "question": spec["raw"]["question"],
        "seed": spec["raw"]["seed"], "replicates": reps, "replicates_override": replicates is not None,
        "workers": 1, "numerics": spec["numerics"], "code": code, "started": stamp,
        "seconds": time.perf_counter() - start, "cells": results,
        "scope": "Whole-line per-reading Q; independent isotropic channels; pooled variance and eta have different targets."}
    np.savez_compressed(out / "values.npz", **arrays)
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    (out / "summary.md").write_text(_markdown(record), encoding="utf-8")
    fields = ["cell", "dimension", "mass", "snr", "readings", "noncentrality", "pooled_mean", "pooled_sd",
              "reference_mean", "reference_sd", "eta_unweighted_mean", "eta_unweighted_sd", "eta_ours29_mean", "eta_ours29_sd"]
    with (out / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in results:
            w, s, ref = row["world"], row["summary"], row["reference"]
            writer.writerow(dict(cell=row["id"], dimension=w["dimension"], mass=w["mass"],
                snr=w["acceleration_snr"], readings=w["readings"], noncentrality=ref["noncentrality"],
                pooled_mean=s["pooled_variance"]["mean"], pooled_sd=s["pooled_variance"]["sd"],
                reference_mean=ref["pooled_mean"], reference_sd=ref["pooled_sd"],
                eta_unweighted_mean=s["eta_unweighted"]["mean"], eta_unweighted_sd=s["eta_unweighted"]["sd"],
                eta_ours29_mean=s["eta_ours29"]["mean"], eta_ours29_sd=s["eta_ours29"]["sd"]))
    if progress:
        progress(f"wrote {out}")
    return out, record

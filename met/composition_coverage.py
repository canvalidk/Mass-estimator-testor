"""B03's historical ours29 composition study, as an explicit recorded preset.

The per-replicate calculation is copied from the 28 Oct tool without a numerical
change. Its size-major SeedSequence layout is retained, including the dependence
of later sizes' seeds on the replicate count. This is historical reproduction,
not a new assessment or adoption of the working mass equation.
"""
from copy import deepcopy
import csv
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path

import numpy as np
import yaml


SF, SA, S = 0.2, 0.1, 2.0
M1, M2 = 1.0, 2.0
MU = 1 / (1 / M1 + 1 / M2)
SIZES = (2.0, 0.4, 0.1)
N_READINGS = 4
DRAWS = 20000
ORIGINAL_SOURCE_SHA256 = "106c26e7b3c43c35d39e5625142269a0554e6cdeb4bc1d50e2dac1ac4c412409"
FIXED = dict(dimension=3, n_readings=N_READINGS, draws=DRAWS, sigma_f=SF,
             sigma_a=SA, c=S, m1=M1, m2=M2, equation="ours29")


def _series(rng, mass, size):
    from met.full import isotropic
    a = rng.standard_normal((N_READINGS, 3))
    a *= size / np.linalg.norm(a, axis=1, keepdims=True)
    at = a + SA * rng.standard_normal(a.shape)
    f = mass * a + SF * rng.standard_normal(a.shape)
    return [isotropic(fi, ai, SF, SA) for fi, ai in zip(f, at)]


def _sample(result, rng):
    w = result["law"] / result["law"].sum()
    cdf = np.cumsum(w) - 0.5 * w
    return np.exp(np.interp(rng.uniform(size=DRAWS), cdf, result["lm"]))


def _one(args):
    from met.full import ours29
    seed, size = args
    rng = np.random.default_rng(seed)
    p1 = ours29(_series(rng, M1, size), S)
    p2 = ours29(_series(rng, M2, size), S)
    lump = ours29(_series(rng, M1 + M2, size), S)
    comp = ours29(_series(rng, MU, size), S)
    d1, d2 = _sample(p1, rng), _sample(p2, rng)
    out = {}
    for name, truth, computed, own, point in (
            ("L5", M1 + M2, d1 + d2, lump, lambda x, y: x + y),
            ("L14", MU, 1 / (1 / d1 + 1 / d2), comp, lambda x, y: 1 / (1 / x + 1 / y))):
        lo, hi = np.quantile(computed, [0.025, 0.975])
        olo, ohi = own["interval"]
        out[name] = dict(
            computed_covers=bool(lo <= truth <= hi), own_covers=bool(olo <= truth <= ohi),
            overlap=bool(max(lo, olo) <= min(hi, ohi)),
            defect_cloud=abs(own["readout"] - point(p1["readout"], p2["readout"])) / point(p1["readout"], p2["readout"]),
            defect_median=abs(own["median"] - point(p1["median"], p2["median"])) / point(p1["median"], p2["median"]),
            computed_median=float(np.median(computed)), own_median=own["median"])
    out["parts_cover"] = [bool(p1["interval"][0] <= M1 <= p1["interval"][1]),
                          bool(p2["interval"][0] <= M2 <= p2["interval"][1])]
    return out


def validate(data):
    required = {"kind", "name", "question", "replicates", "seed", "sizes", *FIXED}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"composition_coverage needs {sorted(required)}")
    if set(data) - required - {"notes"}:
        raise ValueError("unknown composition_coverage setting")
    if data["kind"] != "composition_coverage":
        raise ValueError("expected kind: composition_coverage")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for ch in name)):
        raise ValueError("composition_coverage name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("composition_coverage needs a nonempty question")
    for key in ("replicates", "seed"):
        minimum = 2 if key == "replicates" else 0
        if type(data[key]) is not int or data[key] < minimum:
            raise ValueError(f"composition_coverage {key} must be an integer >= {minimum}")
    sizes = data["sizes"]
    if (not isinstance(sizes, list) or not sizes or any(
            isinstance(v, bool) or not isinstance(v, (int, float))
            or not math.isfinite(v) or v <= 0 for v in sizes)):
        raise ValueError("composition_coverage sizes must be a nonempty list of finite positive numbers")
    if len(set(sizes)) != len(sizes):
        raise ValueError("composition_coverage sizes must be unique")
    for key, value in FIXED.items():
        if isinstance(data[key], bool) or data[key] != value:
            raise ValueError(f"composition_coverage preserves historical {key}={value!r}")
        if key in ("dimension", "n_readings", "draws") and type(data[key]) is not int:
            raise ValueError(f"composition_coverage {key} must be an integer")
    raw = deepcopy(data)
    raw["sizes"] = [float(v) for v in sizes]
    return dict(kind="composition_coverage", raw=raw, settings=deepcopy(FIXED),
                cells=[dict(id=f"c{i:03d}", size=size) for i, size in enumerate(raw["sizes"])])


def jobs(seed, replicates, sizes):
    """The original flat uint32 stream, sliced in size-major order."""
    seeds = np.random.SeedSequence(seed).generate_state(replicates * len(sizes))
    return [(int(seeds[k * replicates + r]), size)
            for k, size in enumerate(sizes) for r in range(replicates)]


def summarize(results, sizes, replicates):
    """The original np.mean/np.median aggregation, in its original order."""
    summary = {}
    for k, size in enumerate(sizes):
        block = results[k * replicates:(k + 1) * replicates]
        row = {"parts_coverage": float(np.mean([c for b in block for c in b["parts_cover"]]))}
        for name in ("L5", "L14"):
            rows = [b[name] for b in block]
            row[name] = {key: float(np.mean([r[key] for r in rows]))
                         for key in ("computed_covers", "own_covers", "overlap")}
            row[name].update({f"median_{key}": float(np.median([r[key] for r in rows]))
                              for key in ("defect_cloud", "defect_median")})
        summary[str(size)] = row
    return summary


def legacy_payload(record):
    return {key: record[key] for key in ("seed", "replicates", "sizes", "n_readings",
                                       "sigma_f", "sigma_a", "c", "m1", "m2", "summary")}


def print_legacy(record):
    for size, row in record["summary"].items():
        print(f"size {size}: parts cover {row['parts_coverage']:.3f}")
        for name in ("L5", "L14"):
            r = row[name]
            print(f"   {name:3s}: computed covers {r['computed_covers']:.3f}  own covers {r['own_covers']:.3f}  overlap {r['overlap']:.3f}"
                  f"  | point defect: cloud {r['median_defect_cloud']:.3%}  median {r['median_defect_median']:.3%}")


def _arrays(results, work):
    arrays = {"job_seeds": np.array([seed for seed, size in work], dtype=np.uint32),
              "job_sizes": np.array([size for seed, size in work]),
              "parts_cover": np.array([r["parts_cover"] for r in results], dtype=bool)}
    for name in ("L5", "L14"):
        for key in results[0][name]:
            arrays[f"{name}_{key}"] = np.array([r[name][key] for r in results])
    return arrays


def run(path, out_root="results", replicates=None, progress=print, workers=None):
    from .runner import code_version
    path = Path(path)
    source = path.read_bytes()
    spec = validate(yaml.safe_load(source))
    raw = deepcopy(spec["raw"])
    if replicates is not None:
        if type(replicates) is not int or replicates < 2:
            raise ValueError("composition_coverage replicates must be an integer >= 2")
        raw["replicates"] = replicates
    workers = 1 if workers is None else workers
    if type(workers) is not int or workers < 1:
        raise ValueError("composition_coverage workers must be a positive integer")
    # Retain the original tool's OMP policy. Each _one owns its complete RNG.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    work = jobs(raw["seed"], raw["replicates"], raw["sizes"])
    code = code_version()
    code["source_preset_sha256"] = hashlib.sha256(source).hexdigest()
    if workers == 1:
        results = []
        for i, job in enumerate(work):
            results.append(_one(job))
            if progress and (i + 1) % raw["replicates"] == 0:
                progress(f"size {job[1]}: {raw['replicates']} replicates complete")
    else:
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            results = list(pool.map(_one, work, chunksize=4))
    summary = summarize(results, raw["sizes"], raw["replicates"])
    record = dict(seed=raw["seed"], replicates=raw["replicates"], sizes=raw["sizes"],
                  n_readings=N_READINGS, sigma_f=SF, sigma_a=SA, c=S, m1=M1, m2=M2,
                  summary=summary, kind="composition_coverage", workers=workers, code=code,
                  equation="ours29", draws=DRAWS, mu=MU, cells=spec["cells"],
                  original_source_sha256=ORIGINAL_SOURCE_SHA256,
                  seed_layout="SeedSequence.generate_state(replicates * len(sizes)); size-major; later-size seeds change with replicates")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    base = Path(out_root) / raw["name"]
    out, suffix = base / stamp, 1
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    (out / "source_preset.yaml").write_bytes(source)
    if replicates is None:
        (out / "preset.yaml").write_bytes(source)
    else:
        (out / "preset.yaml").write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    code["effective_preset_sha256"] = hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest()
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    np.savez_compressed(out / "values.npz", **_arrays(results, work))
    fields = ("size", "parts_coverage", "composition", "computed_covers", "own_covers", "overlap",
              "median_defect_cloud", "median_defect_median")
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for size, row in summary.items():
            for name in ("L5", "L14"):
                writer.writerow(dict(size=size, parts_coverage=row["parts_coverage"], composition=name, **row[name]))
    lines = ["# Historical composition coverage (ours29)", "",
             "Legacy fraction formatting is retained, including its rounding ties.", "",
             "| size | parts | rule | computed | own | overlap (fraction) | cloud defect | median defect |",
             "| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    for size, row in summary.items():
        for name in ("L5", "L14"):
            r = row[name]
            lines.append(f"| {size} | {row['parts_coverage']:.3f} | {name} | {r['computed_covers']:.3f} | "
                         f"{r['own_covers']:.3f} | {r['overlap']:.3f} | {r['median_defect_cloud']:.3%} | {r['median_defect_median']:.3%} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if progress:
        progress(f"saved {out}")
    return out, record

"""B21's paired coverage study; research only, no production estimator changes.

Run: python tools/maximal_invariant_study.py --workers 2
Configuration: tools/maximal_invariant_study.yaml. All laws receive identical
readings, declared unit noise and the half-Cauchy at c=s=1. The world uses
fixed-length pushes in independent random directions, as O15's fixed-signal
coverage experiment does. Investigator truth is used only to generate and score.
Central 95% coverage is evaluated by the angular CDF at truth, not a normal
approximation or simulated posterior samples. Two quadrature orders are retained.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

# Keep each worker's numerical libraries to one thread before importing NumPy.
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
from scipy.special import logsumexp
import scipy
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools.maximal_invariant import profile_log_kernel, integrated_log_kernel, profile_delta


def series_statistics(force, acceleration):
    """Full ours30's isotropic algebra, independently anchored by a regression test."""
    x, y = np.asarray(force, float), np.asarray(acceleration, float)
    if x.shape != y.shape or x.ndim != 3 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("finite paired arrays must have shape (series, readings, axes)")
    aa, bb, ab = (x*x).sum(-1), (y*y).sum(-1), (x*y).sum(-1)
    own_gap = np.hypot(aa-bb, 2*ab)
    q_line = np.maximum((aa+bb-own_gap)/2, 0.)
    qp = np.where(ab > 0, q_line, np.minimum(aa, bb))
    g = np.minimum(1., (x.shape[-1]-.5)/np.maximum(qp, 1e-300))

    def matrix(weight):
        a, b, c = (weight*aa).sum(-1), (weight*bb).sum(-1), (weight*ab).sum(-1)
        return np.stack((a, c, c, b), axis=-1).reshape(-1, 2, 2)

    gram, weighted = matrix(1.), matrix(g)
    eigen, vectors = np.linalg.eigh(weighted)
    v = vectors[:, :, -1]
    scores = g*((aa-bb)*(v[:, 0]*v[:, 1])[:, None]
                + ab*(v[:, 1]**2-v[:, 0]**2)[:, None])
    j = (scores*scores).sum(-1)
    curvature = eigen[:, 1]-eigen[:, 0]
    eta = np.ones_like(curvature)
    valid = (j > 0) & (curvature > 0)
    eta[valid] = np.minimum(1., curvature[valid]/j[valid])
    return gram, weighted, eta


def angular_cdfs(gram, weighted, eta, dimensions, truth_angle, nodes):
    """CDFs at one angle with uniform-angle prior (half-Cauchy at c=s)."""
    q, w = np.polynomial.legendre.leggauss(nodes)
    ends = ((0., truth_angle), (truth_angle, math.pi/2))
    angles = np.concatenate([(a+b)/2+(b-a)*q/2 for a, b in ends])
    weights = np.concatenate([(b-a)*w/2 for a, b in ends])
    v = np.stack((np.sin(angles), np.cos(angles)), axis=-1)
    h = np.maximum(np.einsum("qi,bij,qj->bq", v, gram, v), 0.)
    hw = np.maximum(np.einsum("qi,bij,qj->bq", v, weighted, v), 0.)
    kernels = {"profile": profile_log_kernel(h, dimensions),
               "integrated": integrated_log_kernel(h, dimensions),
               "flat": h/2, "ours30": eta[:, None]*hw/2}
    out = {}
    for name, values in kernels.items():
        values = np.asarray(values)+np.log(weights)[None, :]
        out[name] = np.exp(logsumexp(values[:, :nodes], axis=1)-logsumexp(values, axis=1))
    return out


def run_cell(job):
    config, signal_index, signal, count = job
    rng = np.random.default_rng(np.random.SeedSequence([config["seed"], signal_index, count]))
    batch, k, mass = config["series"], config["dimension"], config["mass"]
    direction = rng.normal(size=(batch, count, k))
    direction /= np.linalg.norm(direction, axis=-1, keepdims=True)
    acceleration_mean = signal*direction/math.sqrt(1+mass*mass)
    force = mass*acceleration_mean+rng.normal(size=direction.shape)
    acceleration = acceleration_mean+rng.normal(size=direction.shape)
    gram, weighted, eta = series_statistics(force, acceleration)
    d = count*k
    eigen = np.linalg.eigvalsh(gram)
    fitted_delta = profile_delta(eigen[:, -1], d)
    coarse = angular_cdfs(gram, weighted, eta, d, math.atan(mass), config["quadrature_nodes"])
    fine = angular_cdfs(gram, weighted, eta, d, math.atan(mass), 2*config["quadrature_nodes"])
    rho = fitted_delta/eigen[:, -1]
    gram_eta = 1-eigen[:, 0]/eigen[:, -1]
    target = signal*signal/(signal*signal+k)
    summary = {"signal": signal, "readings": count, "dimensions": d, "series": batch,
               "target": target, "mean_profile_rho": float(rho.mean()),
               "mean_gram_eta": float(gram_eta.mean()), "mean_ours30_eta": float(eta.mean()),
               "profile_rho_rmse": float(np.sqrt(np.mean((rho-target)**2))),
               "gram_eta_rmse": float(np.sqrt(np.mean((gram_eta-target)**2)))}
    records = {"gram": gram, "weighted_gram": weighted, "ours30_eta": eta,
               "profile_delta": fitted_delta, "profile_rho": rho, "gram_eta": gram_eta}
    for name in fine:
        cdf, before = fine[name], coarse[name]
        cover = (cdf >= .025) & (cdf <= .975)
        old_cover = (before >= .025) & (before <= .975)
        error = np.abs(cdf-before)
        probability = float(cover.mean())
        summary[name] = probability
        summary[name+"_mcse"] = math.sqrt(probability*(1-probability)/batch)
        summary[name+"_quadrature_max_difference"] = float(error.max())
        summary[name+"_quadrature_verdict_changes"] = int(np.sum(cover != old_cover))
        records[name+"_cdf"] = cdf
        records[name+"_cdf_coarse"] = before
        records[name+"_covered"] = cover
    for name in ("profile", "integrated"):
        for control in ("ours30", "flat"):
            diff = records[name+"_covered"].astype(float)-records[control+"_covered"]
            summary[name+"_minus_"+control] = float(diff.mean())
            summary[name+"_minus_"+control+"_mcse"] = float(diff.std(ddof=1)/math.sqrt(batch))
    return summary, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", type=Path, default=Path(__file__).with_suffix(".yaml"))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--series", type=int)
    args = parser.parse_args()
    config = yaml.safe_load(args.preset.read_text(encoding="utf-8"))
    if args.series is not None:
        config["series"] = args.series
    validate_config(config)
    if args.workers < 1 or config["series"] < 2 or config["quadrature_nodes"] < 16:
        raise ValueError("positive workers, at least two series and 16 quadrature nodes required")
    jobs = [(config, idx, signal, count) for idx, signal in enumerate(config["signals"])
            for count in config["readings"]]
    out = args.out or ROOT/"results"/"b21-maximal-invariant"/datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out.mkdir(parents=True, exist_ok=False)
    (out/"preset.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    summaries, arrays = [], {}
    pool = ProcessPoolExecutor(args.workers) if args.workers > 1 else None
    try:
        results = pool.map(run_cell, jobs) if pool else map(run_cell, jobs)
        for idx, (summary, records) in enumerate(results):
            summaries.append(summary)
            arrays.update({f"cell_{idx}/{key}": value for key, value in records.items()})
            print(json.dumps(summary), flush=True)
    finally:
        if pool:
            pool.shutdown()
    np.savez_compressed(out/"values.npz", **arrays)
    with (out/"summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    files = [Path(__file__), ROOT/"tools"/"maximal_invariant.py"]
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
    metadata = {"config": config, "workers": args.workers, "python": platform.python_version(),
                "numpy": np.__version__, "scipy": scipy.__version__, "code_sha256": hashes,
                "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
                "summaries": summaries,
                "numerical_scope": "Quadrature differences are error estimates, not certified bounds."}
    (out/"run.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    print(f"Saved {out}", flush=True)


def validate_config(config):
    """Refuse a preset that would misdescribe this deliberately locked experiment."""
    if (config["dimension"] != 3 or config["mass"] != 1.
            or config["prior_center"] != 1. or config["noise_sd"] != 1.):
        raise ValueError("this locked B21 study is 3D with m=c=s=1 and unit noise")
    if config["integrated_measure"] != "flat_dDelta":
        raise ValueError("this comparison integrates against flat_dDelta only")


if __name__ == "__main__":
    main()

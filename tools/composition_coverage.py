"""Composition at the level of the law (L5, L14 as restated 28 Oct, R-2026-10-28-1).

For each replicate: two parts (true masses m1, m2) measured on their own; a lump (m1 + m2)
measured anew; a composite under a common force (reduced mass mu = (1/m1 + 1/m2)^-1) measured
anew. All series: N readings in 3D, isotropic noise, the default starting law c = s.

Reported, per signal size:
- clause (b)/(c): coverage of m1 + m2 (and of mu) by the computed law's 95% interval (sampled
  from the parts' laws) and by the composite's own quoted 95% interval; how often the two
  intervals overlap;
- the point-level defects on data, for the cloud readout and for the law's median, for both L5
  and L14 (the check owed from 28 Oct).

Usage: python tools/composition_coverage.py [--replicates 200] [--workers 4] [--seed 2810]
"""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SF, SA, S = 0.2, 0.1, 2.0
M1, M2 = 1.0, 2.0
MU = 1 / (1 / M1 + 1 / M2)
SIZES = (2.0, 0.4, 0.1)
N_READINGS = 4
DRAWS = 20000


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=2810)
    ap.add_argument("--out", default="results/composition_coverage.json")
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    seeds = np.random.SeedSequence(args.seed).generate_state(args.replicates * len(SIZES))
    jobs = [(int(seeds[k * args.replicates + r]), size) for k, size in enumerate(SIZES) for r in range(args.replicates)]
    with ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(_one, jobs, chunksize=4))
    summary = {}
    for k, size in enumerate(SIZES):
        block = results[k * args.replicates:(k + 1) * args.replicates]
        row = {"parts_coverage": float(np.mean([c for b in block for c in b["parts_cover"]]))}
        for name in ("L5", "L14"):
            rows = [b[name] for b in block]
            row[name] = {key: float(np.mean([r[key] for r in rows])) for key in ("computed_covers", "own_covers", "overlap")}
            row[name].update({f"median_{key}": float(np.median([r[key] for r in rows]))
                              for key in ("defect_cloud", "defect_median")})
        summary[str(size)] = row
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(dict(seed=args.seed, replicates=args.replicates, sizes=SIZES, n_readings=N_READINGS,
                       sigma_f=SF, sigma_a=SA, c=S, m1=M1, m2=M2, summary=summary), fh, indent=1)
    for size, row in summary.items():
        print(f"size {size}: parts cover {row['parts_coverage']:.3f}")
        for name in ("L5", "L14"):
            r = row[name]
            print(f"   {name:3s}: computed covers {r['computed_covers']:.3f}  own covers {r['own_covers']:.3f}  overlap {r['overlap']:.3f}"
                  f"  | point defect: cloud {r['median_defect_cloud']:.3%}  median {r['median_defect_median']:.3%}")


if __name__ == "__main__":
    main()

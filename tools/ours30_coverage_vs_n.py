"""ours30's coverage as readings are added at fixed signal: the Neyman–Scott regime (O15, 30 Oct).

3D, isotropic unit noise, m = c = 1, push of noise-unit size `signal` in a uniform direction per reading.
Compares ours30 with ours30 without eta (the flat factor alone). Seed [3112, 100*i + N].
Run: python tools/ours30_coverage_vs_n.py [--series 100]   (about 7 minutes on 2 cores)
"""
import argparse, math, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from met.full import isotropic, ours30  # noqa: E402


def cell(args):
    idx, r, n, series = args
    rng = np.random.default_rng([3112, 100 * idx + n])
    with_eta, without = [], []
    for _ in range(series):
        rd = []
        for _ in range(n):
            u = rng.normal(size=3); u /= np.linalg.norm(u); x = r * u / math.sqrt(2)
            rd.append(isotropic(x + rng.normal(size=3), x + rng.normal(size=3), 1.0, 1.0))
        lo, hi = ours30(rd, 1.0, readout=False)["interval"]; with_eta.append(lo <= 1 <= hi)
        lo, hi = ours30(rd, 1.0, calibrate=False, readout=False)["interval"]; without.append(lo <= 1 <= hi)
    return r, n, float(np.mean(with_eta)), float(np.mean(without))


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--series", type=int, default=100); p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    jobs = [(i, r, n, a.series) for i, r in enumerate([1.0, 2.0]) for n in (5, 20, 80, 200)]
    with ProcessPoolExecutor(a.workers) as ex:
        for r, n, w, wo in ex.map(cell, jobs):
            print(f"3D signal {r} N={n:4d}: ours30 {w:.2f} | without eta {wo:.2f}", flush=True)

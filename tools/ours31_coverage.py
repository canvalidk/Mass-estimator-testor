"""ours31 against ours30: coverage and width across signal strength (problems/O18.md, hole 3).

    python tools/ours31_coverage.py [--series 150] [--workers 4]

Seeded (base seed 3110, one stream per cell). World: true m = 1, c = s = 1, whitened unit
channel noise scaled per reading by sqrt(v_i); signal length r on a random direction.
Worlds: 'gaussian' (v_i = 1) and 'nu10' (v_i ~ InvGamma(5, 5): sigma stated from 10 dof).
ours31 declares nu_eff = 10 in both. A script, not a preset kind (as tools/composition_coverage.py;
STATE.md Gaps)."""
import argparse, math, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from met.full import isotropic, ours30  # noqa: E402
from met.ours31 import ours31  # noqa: E402

WORLDS = ("gaussian", "nu10")
KS = (1, 3)
SIGNALS = (0.0, 0.5, 1.0, 2.0, 5.0)
NS = (5, 20)
NU_EFF = 10.0


def cell(args):
    idx, world, k, r, n, n_series = args
    rng = np.random.default_rng([3110, idx])
    out = {"ours30": [], "ours31": []}
    for _ in range(n_series):
        rd = []
        for _ in range(n):
            u = rng.normal(size=k); u /= np.linalg.norm(u)
            x = r * u / math.sqrt(2)
            v = 1.0 if world == "gaussian" else (NU_EFF / 2) / rng.gamma(NU_EFF / 2)
            rd.append(isotropic(x + math.sqrt(v) * rng.normal(size=k), x + math.sqrt(v) * rng.normal(size=k), 1.0, 1.0))
        for name, res in (("ours30", ours30(rd, 1.0, readout=False)), ("ours31", ours31(rd, 1.0, NU_EFF))):
            lo, hi = res["interval"]
            out[name].append((lo <= 1.0 <= hi, math.log(hi / lo)))
    return (world, k, r, n, {key: (float(np.mean([c for c, _ in v])), float(np.median([w for _, w in v]))) for key, v in out.items()})


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--series", type=int, default=150)
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args(argv)
    jobs, i = [], 0
    for world in WORLDS:
        for k in KS:
            for n in NS:
                for r in SIGNALS:
                    jobs.append((i, world, k, r, n, a.series)); i += 1
    print(f"series per cell {a.series}; ours31 nu_eff = {NU_EFF}")
    print("world     k   N  signal | ours30 cover  width | ours31 cover  width")
    with ProcessPoolExecutor(a.workers) as ex:
        for world, k, r, n, res in ex.map(cell, jobs):
            (c0, w0), (c1, w1) = res["ours30"], res["ours31"]
            print(f"{world:8s} {k:2d} {n:3d}  {r:5.1f}  |   {c0:.3f}     {w0:6.3f} |   {c1:.3f}     {w1:6.3f}", flush=True)


if __name__ == "__main__":
    main()

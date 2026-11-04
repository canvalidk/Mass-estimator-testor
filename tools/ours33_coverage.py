"""ours30, ours32 and ours33: coverage and width across signal strength (O18; ours33 = T48's exponent + T47's warp).

    python tools/ours33_coverage.py [--series 100] [--workers 2] [--truth 1.0] [--cells all|short]

Seeded (base seed 3110, one stream per cell; with --truth != 1 the stream is [3110, 500 + idx]). World as
tools/ours32_coverage.py: c = s = 1, whitened unit channel noise scaled per reading by sqrt(v_i), signal
length r on a random direction; worlds 'gaussian' (v_i = 1) and 'nu10' (v_i ~ InvGamma(5, 5)); ours32 and
ours33 declare nu_eff = 10 in both. --truth m0 puts the true mass off the starting law's centre (the
28 Oct composition study and every O18 cell so far had m0 = c, which rewards a pull toward c). With
--truth 1.0 the ours30 and ours32 columns reproduce tools/ours32_coverage.py's cells to the digit.
Reported per cell: 95% coverage, median log-width ln(hi/lo), median |ln(m_hat/m0)|. A script, not a
preset kind (STATE.md Gaps)."""
import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from met.full import isotropic, ours30  # noqa: E402
from met.ours32 import ours32  # noqa: E402
from met.ours33 import ours33  # noqa: E402

WORLDS = ("gaussian", "nu10")
KS = (1, 3)
SIGNALS = (0.0, 0.5, 1.0, 2.0, 5.0)
NS = (5, 20)
NU_EFF = 10.0
SHORT = [("gaussian", 3, 20, 0.0), ("gaussian", 3, 20, 1.0), ("gaussian", 3, 20, 2.0), ("gaussian", 3, 20, 5.0),
         ("nu10", 3, 20, 0.0), ("nu10", 3, 20, 1.0), ("nu10", 3, 5, 1.0),
         ("gaussian", 1, 20, 0.0), ("gaussian", 1, 20, 1.0)]


def index(world, k, n, r):
    return WORLDS.index(world) * 20 + KS.index(k) * 10 + NS.index(n) * 5 + SIGNALS.index(r)


def cell(args):
    world, k, r, n, n_series, m0 = args
    idx = index(world, k, n, r)
    rng = np.random.default_rng([3110, idx] if m0 == 1.0 else [3110, 500 + idx])
    out = {"ours30": [], "ours32": [], "ours33": []}
    for _ in range(n_series):
        rd = []
        for _ in range(n):
            u = rng.normal(size=k)
            u /= np.linalg.norm(u)
            a = r * u / math.sqrt(1.0 + m0 * m0)
            v = 1.0 if world == "gaussian" else (NU_EFF / 2) / rng.gamma(NU_EFF / 2)
            rd.append(isotropic(m0 * a + math.sqrt(v) * rng.normal(size=k), a + math.sqrt(v) * rng.normal(size=k), 1.0, 1.0))
        for name, res in (("ours30", ours30(rd, 1.0, readout=False)), ("ours32", ours32(rd, 1.0, NU_EFF)),
                          ("ours33", ours33(rd, 1.0, NU_EFF))):
            lo, hi = res["interval"]
            out[name].append((lo <= m0 <= hi, math.log(hi / lo), abs(math.log(res["median"] / m0))))
    summary = {key: (float(np.mean([c for c, _, _ in v])), float(np.median([w for _, w, _ in v])),
                     float(np.median([e for _, _, e in v]))) for key, v in out.items()}
    return world, k, r, n, summary


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--series", type=int, default=100)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--truth", type=float, default=1.0)
    p.add_argument("--cells", choices=("all", "short"), default="short")
    a = p.parse_args(argv)
    cells = SHORT if a.cells == "short" else [(w, k, n, r) for w in WORLDS for k in KS for n in NS for r in SIGNALS]
    jobs = [(w, k, r, n, a.series, a.truth) for (w, k, n, r) in cells]
    print(f"series per cell {a.series}; true mass {a.truth}, c = s = 1; nu_eff declared {NU_EFF}")
    print("world     k   N  signal | cover  width  |err| : ours30 | ours32 | ours33")
    with ProcessPoolExecutor(a.workers) as ex:
        for world, k, r, n, res in ex.map(cell, jobs):
            cols = " | ".join(f"{c:.3f} {w:6.2f} {e:5.2f}" for c, w, e in (res["ours30"], res["ours32"], res["ours33"]))
            print(f"{world:8s} {k:2d} {n:3d}  {r:5.1f}  | {cols}", flush=True)


if __name__ == "__main__":
    main()

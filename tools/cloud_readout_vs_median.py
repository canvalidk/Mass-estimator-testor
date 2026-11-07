"""The cloud readout (T10) against the median and the geometric mean, on ours30's law and ours33's (O2, G6).

Run on 5 Nov in a scratch folder (Claude) and not committed then; rebuilt here 7 Nov with the same world and
seeds so the 5 Nov table can be regenerated. The cloud readout is T10's ratio of means over the candidate
cloud on the stacked length, unchanged; only the law under it changes. On ours30 it is ours30's own `readout`;
on ours33 it is the same computation on ours33's raw law (`cloud`). The geometric mean, exp E[ln m], is the
comparison readout added 6 Nov.

World: tools/ours33_coverage.py's Gaussian world, c = s = 1 (unit noise on both channels), N = 20, k in (3, 1),
signal in (0, 0.5, 1, 2, 5) in noise units, truth m0 = c or 4c (--truth). Seeds: base 3110, one stream per cell,
[3110, idx] at m0 = 1 and [3110, 500 + idx] otherwise, idx as in ours33_coverage.py, so the series are the same
ones that tool scores. 100 series per cell by default.

Each entry: median |ln(m_hat/m0)|, with the mean ln(m_hat/m0) (the bias) in brackets. A study, not a premise
check: by the standing rule it cannot make a readout better or worse as a point.

    python tools/cloud_readout_vs_median.py --truth 1 --workers 2
    python tools/cloud_readout_vs_median.py --truth 4 --k 3
"""
import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from met.full import _stacked_mean_length, isotropic, ours30  # noqa: E402
from met.ours33 import ours33  # noqa: E402

NU = 10.0
WORLDS = ("gaussian", "nu10")
KS = (1, 3)
SIGNALS = (0.0, 0.5, 1.0, 2.0, 5.0)
NS = (5, 20)
KEYS = ("30 median", "30 geo", "30 cloud", "33 median", "33 geo", "33 cloud")


def index(world, k, n, r):
    """The cell index of tools/ours33_coverage.py, so the seeds match."""
    return WORLDS.index(world) * 20 + KS.index(k) * 10 + NS.index(n) * 5 + SIGNALS.index(r)


def cloud(lm, law, readings):
    """T10's cloud readout on a law given on the grid lm: E[m * a(m)] / E[a(m)], a the stacked mean length."""
    lm = np.asarray(lm)
    law = np.asarray(law)
    keep = law > 1e-14 * law.max()
    a = np.zeros_like(lm)
    a[keep] = _stacked_mean_length(np.exp(lm[keep]), readings)
    m = np.exp(lm)
    return float((m * a * law).sum() / (a * law).sum())


def series(rng, m0, k, r, n):
    out = []
    for _ in range(n):
        u = rng.normal(size=k)
        u /= np.linalg.norm(u)
        a = r * u / math.sqrt(1 + m0 * m0)
        out.append(isotropic(m0 * a + rng.normal(size=k), a + rng.normal(size=k), 1.0, 1.0))
    return out


def cell(args):
    k, r, n, ns, m0 = args
    idx = index("gaussian", k, n, r)
    rng = np.random.default_rng([3110, idx] if m0 == 1.0 else [3110, 500 + idx])
    out = {key: [] for key in KEYS}
    for _ in range(ns):
        rd = series(rng, m0, k, r, n)
        r30 = ours30(rd, 1.0, readout=True)
        r33 = ours33(rd, 1.0, NU)
        vals = {"30 median": r30["median"], "30 geo": r30["geometric_mean"], "30 cloud": r30["readout"],
                "33 median": r33["median"], "33 geo": r33["geometric_mean"],
                "33 cloud": cloud(r33["lm"], r33["law"], rd)}
        for key, v in vals.items():
            out[key].append(math.log(v / m0))
    return k, r, n, m0, {key: (float(np.median(np.abs(v))), float(np.mean(v))) for key, v in out.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--truth", type=float, default=1.0, help="m0 in units of c (1 or 4 in the 5 Nov table)")
    ap.add_argument("--k", type=int, choices=KS, nargs="*", default=[3, 1])
    ap.add_argument("--series", type=int, default=100)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args(argv)
    jobs = [(k, r, 20, args.series, args.truth) for k in args.k for r in SIGNALS]
    print("m0  k  signal  " + "  ".join(f"{key:>14}" for key in KEYS))
    run = map if args.workers <= 1 else ProcessPoolExecutor(args.workers).map
    for k, r, n, m0, res in run(cell, jobs):
        print(f"{m0:<3g} {k}  {r:<6g}  " + "  ".join(f"{e:5.2f} ({b:+.2f})" for e, b in (res[key] for key in KEYS)),
              flush=True)


if __name__ == "__main__":
    main()

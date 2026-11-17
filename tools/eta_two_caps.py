"""Eta's two jobs with a cap on each (Claude, 17 Nov 2026; problems/O14.md, "B22 read").

On B22's worlds and seeds (tools/ours30_eta_gram_coverage.py; unit Gaussian errors, truth m = c = s = 1), the
Bartlett eta is gap/J and the Gram eta is gap/H_max, so

    eta_bartlett = eta_gram * (H_max / J):   the weak-signal share times expected-over-observed squared scores.

The Bartlett form caps the product at 1; capping each factor at 1 gives eta_gram * min(1, H_max/J), which is the
smaller of the two existing etas. This script reports the coverage of that form beside B22's two, per cell, paired
on the same series. It changes nothing in met/.

    python tools/eta_two_caps.py            # 28 cells, 300 series each, about a minute
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ours30_eta_gram_coverage as b22  # noqa: E402

BASES = ("ours30", "ours29")


def cell(k, n, signal, series=300):
    """Per base: coverage of Bartlett, Gram and two-caps, the two-caps median log width, and how often
    the Bartlett eta exceeds the Gram eta (the readings agree better than the declared noise predicts)."""
    rows = {base: [] for base in BASES}
    for replicate in range(series):
        force, acceleration = b22.world(k, n, signal, replicate)
        results = b22.fast_results(force, acceleration)
        for base in BASES:
            qb, eta_b, _ = results[base + "_bartlett"]
            qg, eta_g, _ = results[base + "_gram"]
            qt = qb if eta_b <= eta_g else qg
            rows[base].append((qb[0] <= 1.0 <= qb[2], qg[0] <= 1.0 <= qg[2], qt[0] <= 1.0 <= qt[2],
                               math.log(qt[2] / qt[0]), eta_b > eta_g))
    out = {}
    for base, values in rows.items():
        arr = np.asarray(values, float)
        out[base] = dict(bartlett=arr[:, 0].mean(), gram=arr[:, 1].mean(), two_caps=arr[:, 2].mean(),
                         median_log_width=float(np.median(arr[:, 3])), bartlett_above_gram=arr[:, 4].mean())
    return out


def main():
    print("| k | N | r | 30 B | 30 G | 30 two caps | 29 B | 29 G | 29 two caps | Bartlett above Gram (30 / 29) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for k, n, r in b22.cells():
        res = cell(k, n, r)
        a, b = res["ours30"], res["ours29"]
        print(f"| {k} | {n} | {r:g} | {100*a['bartlett']:.1f} | {100*a['gram']:.1f} | {100*a['two_caps']:.1f} "
              f"| {100*b['bartlett']:.1f} | {100*b['gram']:.1f} | {100*b['two_caps']:.1f} "
              f"| {a['bartlett_above_gram']:.2f} / {b['bartlett_above_gram']:.2f} |")


if __name__ == "__main__":
    main()

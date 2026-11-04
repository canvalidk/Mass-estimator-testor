"""eta against the shared size law's shrinkage, cell by cell (O14, 4 Nov).

    python tools/eta_vs_size_law.py [--series 300] [--workers 2]

Seeded (base seed 2811, one stream per cell). World: true m = 1, c = s = 1, unit noise per component,
pushes Gaussian per component with SD = snr (the ours28 record's world, section 5), k spatial axes, N
readings. For each series: eta of the flat law (met.full.bartlett_eta with g = 1), the shared size
law's profile shrinkage beta_hat = 1 - M/H^2 and arcsine posterior mean beta_bar (met.size_law), and
the oracle Delta/(Delta + M). Prints cell means and the per-series correlation of eta with beta_hat.
A script, not a preset kind (as tools/composition_coverage.py; STATE.md Gaps)."""
import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from met.full import isotropic  # noqa: E402
from met.size_law import beta_bar, beta_hat, explained_energy, flat_eta, latent_count, oracle  # noqa: E402

KS = (1, 3)
NS = (5, 20, 100)
SNRS = (0.5, 1.0, 2.0, 5.0)


def world(rng, k, n, snr):
    """One series: true pairs (a, a) with a ~ N(0, snr^2/2 I_k) so |(F*, a*)|^2 has mean k snr^2."""
    a = rng.normal(scale=snr, size=(n, k)) / math.sqrt(2.0)
    readings = [isotropic(a[i] + rng.normal(size=k), a[i] + rng.normal(size=k), 1.0, 1.0) for i in range(n)]
    return readings, 2.0 * float((a * a).sum())


def cell(args):
    idx, k, n, snr, n_series = args
    rng = np.random.default_rng([2811, idx])
    etas, hats, bars, oracles = [], [], [], []
    for _ in range(n_series):
        readings, delta = world(rng, k, n, snr)
        eta, _, _, phi_hat = flat_eta(readings, 1.0)
        h2, m = explained_energy(readings, 1.0, phi_hat), latent_count(readings)
        etas.append(eta)
        hats.append(beta_hat(readings, 1.0, phi_hat))
        bars.append(beta_bar(h2, m))
        oracles.append(oracle(delta, m))
    etas, hats = np.array(etas), np.array(hats)
    corr = float(np.corrcoef(etas, hats)[0, 1]) if etas.std() > 0 and hats.std() > 0 else float("nan")
    return k, n, snr, float(etas.mean()), float(hats.mean()), float(np.mean(bars)), float(np.mean(oracles)), corr


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--series", type=int, default=300)
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args(argv)
    jobs = [(i, k, n, snr, a.series) for i, (k, n, snr) in enumerate((k, n, s) for k in KS for n in NS for s in SNRS)]
    print(f"series per cell {a.series}; m = 1, c = s = 1, pushes Gaussian with per-component SD = snr")
    print(" k    N   snr | mean eta  mean beta_hat  mean beta_bar  Delta/(Delta+M) | corr(eta, beta_hat)")
    with ProcessPoolExecutor(a.workers) as ex:
        for k, n, snr, e, h, b, o, corr in ex.map(cell, jobs):
            print(f"{k:2d} {n:4d}  {snr:4.1f} |  {e:.3f}      {h:.3f}          {b:.3f}          {o:.3f}       |  {corr:5.2f}", flush=True)


if __name__ == "__main__":
    main()

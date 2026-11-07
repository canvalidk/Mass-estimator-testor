"""eta_G = 1 - Q/H (T70 (7)) against the Bartlett eta and the target Delta/(Delta + M), cell by cell (O14, 7 Nov).

    python tools/eta_from_the_gram.py [--series 200]

Seeded (base seed 2811, one stream per cell). World: tools/eta_vs_size_law.py's (true m = 1, c = s = 1,
unit noise, Gaussian pushes with per-component SD = snr), flat law (g = 1). Prints cell means and each
estimate's RMS distance from the per-series target. A script, as tools/eta_vs_size_law.py."""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from met.full import gram_eta  # noqa: E402
from met.size_law import flat_eta, latent_count, oracle  # noqa: E402
from eta_vs_size_law import KS, NS, SNRS, world  # noqa: E402


def main(n_series):
    print("k    N  snr | mean: Bartlett  eta_G  target | RMS to target: Bartlett  eta_G")
    for idx, (k, n, snr) in enumerate((k, n, s) for k in KS for n in NS for s in SNRS):
        rng = np.random.default_rng([2811, 70, idx])
        b, g, t = [], [], []
        for _ in range(n_series):
            rs, delta = world(rng, k, n, snr)
            b.append(flat_eta(rs, 1.0)[0])
            g.append(gram_eta(rs, 1.0, [1.0] * len(rs))[0])
            t.append(oracle(delta, latent_count(rs)))
        b, g, t = map(np.array, (b, g, t))
        print(f"{k}  {n:3d}  {snr:3.1f} |  {b.mean():.3f}   {g.mean():.3f}  {t.mean():.3f} |  "
              f"{np.sqrt(((b - t) ** 2).mean()):.3f}   {np.sqrt(((g - t) ** 2).mean()):.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", type=int, default=200)
    main(ap.parse_args().series)

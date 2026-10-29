"""T31 under R-2026-09-07-1 (29 Oct, Claude): does eta, read at the series' best line on the whole
projective line, misbehave when that line is at negative mass?

Declared model correct; true m = 1, sigma_F = sigma_a = 1, c = 1; k = 2; random push directions of
whitened size `size`; 250 series per cell; seed 2914. Coverage of the quoted 95% interval of
ours29(..., positivity=True), with eta and with eta = 1, split by the sign of phi_hat.
Recorded in definition/theorems/T31.md (29 Oct). About 3 minutes.

Usage: python tools/eta_sign_check.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import met.full as f  # noqa: E402

CELLS = ((0.0, 20), (0.5, 20), (1.0, 5))


def main(seed=2914, series=250, m=1.0, k=2):
    rng = np.random.default_rng(seed)
    for size in (0.0, 0.5, 1.0, 2.0):
        for n in (5, 20):
            rows = []
            for _ in range(series):
                rd = []
                for _ in range(n):
                    u = rng.standard_normal(k)
                    u /= np.linalg.norm(u)
                    a = size * u
                    rd.append(f.isotropic(m * a + rng.standard_normal(k), a + rng.standard_normal(k), 1.0, 1.0))
                if (size, n) not in CELLS:
                    continue
                with_eta = f.ours29(rd, 1.0, readout=False, positivity=True)
                raw = f.ours29(rd, 1.0, readout=False, positivity=True, calibrate=False)
                neg = with_eta["phi_hat"] < 0
                for name, o in (("eta", with_eta), ("eta=1", raw)):
                    lo, hi = o["interval"]
                    rows.append((name, neg, lo <= m <= hi))
            for name in ("eta", "eta=1"):
                x = [(neg, cov) for nm, neg, cov in rows if nm == name]
                if not x:
                    continue
                neg = np.array([a for a, _ in x])
                cov = np.array([b for _, b in x], dtype=float)
                print(f"size {size} N {n} {name:5s}: covers when phi_hat < 0 {cov[neg].mean():.3f} (n={neg.sum()})"
                      f" | >= 0 {cov[~neg].mean():.3f} | all {cov.mean():.3f}")


if __name__ == "__main__":
    main()

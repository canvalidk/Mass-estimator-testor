"""Median vs geometric mean (exp E[ln m]) under ours30, ours31, ours32: the size of G6's fork
measured, by signal size (O2, 6 Nov). Seeded; python tools/readout_gap.py."""
import math
import sys

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))
from met.full import isotropic, ours30
from met.ours31 import ours31
from met.ours32 import ours32

C, NU, M, D = 2.0, 10.0, 3.0, 3


def series(rng, n, signal, sf=0.2, sa=0.1):
    out = []
    for _ in range(n):
        a = rng.normal(size=D)
        a *= signal * sa / np.linalg.norm(a)
        f = M * a
        out.append(isotropic(f + rng.normal(scale=sf, size=D), a + rng.normal(scale=sa, size=D), sf, sa))
    return out


def run(name, r):
    if name == "ours30":
        return ours30(r, C, readout=False)
    if name == "ours31":
        return ours31(r, C, NU)
    return ours32(r, C, NU)


def main(reps=200, seed=1106):
    print(f"truth {M}, c {C}, 3D, N = 5, {reps} series a cell; gap = |ln(gmean/median)|, "
          "and the RMS log error of each from the truth")
    print("equation signal  median_gap  max_gap   rms_med  rms_gmean")
    for name in ("ours30", "ours31", "ours32"):
        for signal in (0.0, 0.5, 1.0, 2.0, 5.0):
            rng = np.random.default_rng(seed)
            gaps, em, eg = [], [], []
            for _ in range(reps):
                out = run(name, series(rng, 5, signal))
                gaps.append(abs(math.log(out["geometric_mean"] / out["median"])))
                em.append(math.log(out["median"] / M))
                eg.append(math.log(out["geometric_mean"] / M))
            g = np.array(gaps)
            print(f"{name:7s} {signal:6.1f}  {np.median(g):9.4f}  {g.max():8.4f}  "
                  f"{np.sqrt(np.mean(np.square(em))):7.3f}  {np.sqrt(np.mean(np.square(eg))):8.3f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)

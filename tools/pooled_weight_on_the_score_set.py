"""The Gaussian score set (T45) widened by the pooled weight from pairs (18 Nov): scenes, and what it costs on honest
series (problems/O14.md, 23 Nov). Measurements, not a rule. Run: python tools/pooled_weight_on_the_score_set.py"""
import math, sys
import numpy as np
from scipy import stats
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
from met.conditional_score import gaussian_score_test
from met.full import isotropic
from contradiction_in_pairs import pair_weight
CUT = stats.chi2.ppf(0.95, 1)
def weight(pairs, sf, sa):
    return pair_weight([(np.atleast_1d(np.asarray(F, float)) / sf, np.atleast_1d(np.asarray(a, float)) / sa) for F, a in pairs])
def kept(pairs, sf, sa, w, grid):
    rs = [isotropic(np.atleast_1d(F), np.atleast_1d(a), sf, sa) for F, a in pairs]
    st = np.array([gaussian_score_test(rs, m)["statistic"] for m in grid])
    k = grid[st <= CUT / w]
    pieces = 1 + int(np.sum(np.diff(np.where(st <= CUT / w)[0]) > 1)) if len(k) else 0
    return (k.min(), k.max(), pieces) if len(k) else (math.nan, math.nan, 0)
def show(name, pairs, sf, sa, grid):
    w = weight(pairs, sf, sa); a = kept(pairs, sf, sa, 1.0, grid); b = kept(pairs, sf, sa, w, grid)
    print(f"  {name:44s} weight {w:7.4f} | bell alone [{a[0]:.3f}, {a[1]:.3f}]  widened [{b[0]:.3f}, {b[1]:.3f}]{'' if b[2]==1 else '  pieces: '+str(b[2])}")
g = np.linspace(3.0, 8.0, 5001)
print("1. two scalar readings (acceleration 1); with sigma_a = sigma_F each reading own mass is uncertain by about 10 sigma_F")
for s in (0.1, 0.01):
    show(f"near 5.5, 5.6  sigma={s}", [(5.5, 1.0), (5.6, 1.0)], s, s, g); show(f"far  5.0, 6.0  sigma={s}", [(5.0, 1.0), (6.0, 1.0)], s, s, g)
show("far, nearly exact acceleration (sigma_a=1e-4)", [(5.0, 1.0), (6.0, 1.0)], 0.1, 1e-4, g)
show("5, 5.05, 6", [(5.0, 1.0), (5.05, 1.0), (6.0, 1.0)], 0.1, 0.1, g)
far = [(5.0, 1.0), (6.0, 1.0)]
show("far + 40 exact zeros", far + [(0.0, 0.0)] * 40, 0.1, 0.1, g); show("far + 40 readings of size 1e-9", far + [(1e-9, 1e-9)] * 40, 0.1, 0.1, g)
rng = np.random.default_rng(99)
a99 = 1.0 + 0.02 * rng.normal(size=99); many = [(5.0 * 1.0 + 0.1 * rng.normal(), x) for x in a99] 
show("99 honest readings at 5", many, 0.1, 0.02, np.linspace(4.5, 5.5, 4001)); show("  ... and one stray at 10", many + [(10.0, 1.0)], 0.1, 0.02, np.linspace(4.5, 5.8, 5201))
print("2. 3D readings"); e = np.array([1.0, 0.0, 0.0])
show("near 5.5, 5.6 (aligned 3D)", [(5.5 * e, e), (5.6 * e, e)], 0.1, 0.1, g); show("far 5, 6 (aligned 3D)", [(5.0 * e, e), (6.0 * e, e)], 0.1, 0.1, g)
tilt = np.array([math.cos(0.5), math.sin(0.5), 0.0]); show("one reading, force 29 deg off acceleration", [(5.0 * e, tilt)], 0.1, 0.1, np.linspace(0.5, 40, 8001))
print("3. odd scenes"); gg = np.exp(np.linspace(math.log(1e-3), math.log(1e3), 6001))
show("force-only beside acceleration-only", [(2.0, 0.0), (0.0, 2.0)], 1.0, 1.0, gg); show("force against acceleration (10, -10)", [(10.0, -10.0)], 1.0, 1.0, gg)
print("4. honest series (true mass 2; r = signal per reading in noise units, each channel): share of series containing the truth")
for k in (1, 3):
    for n in (2, 5, 20):
        row = []
        for r in (1.0, 3.0, 10.0):
            rg = np.random.default_rng(1000 * k + 10 * n + int(r)); hit0 = hit1 = low = 0; fac = []
            for _ in range(1500):
                a = rg.normal(size=(n, k)); a *= r * 0.1 / np.linalg.norm(a, axis=1, keepdims=True)
                F = 2.0 * a + 0.2 * rg.normal(size=a.shape); A = a + 0.1 * rg.normal(size=a.shape)
                st = gaussian_score_test([isotropic(f, x, 0.2, 0.1) for f, x in zip(F, A)], 2.0)["statistic"]
                w = weight(list(zip(F, A)), 0.2, 0.1); hit0 += st <= CUT; hit1 += st <= CUT / w
                if w < 1: low += 1; fac.append(1 / math.sqrt(w))
            row.append(f"r={r:g}: {hit0/1500:.3f} -> {hit1/1500:.3f} (widened in {low/1500:.0%}, median x{np.median(fac) if fac else 1:.2f})")
        print(f"  k={k} N={n:2d}  " + "   ".join(row))

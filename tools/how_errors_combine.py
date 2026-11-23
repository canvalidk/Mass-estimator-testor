"""How separate readings errors combine, and what that does to two readings that disagree (problems/O14.md, 23 Nov).
Squares: the declared Gaussian. Lengths: an error law exponential in each reading error size. Measurements, not a rule.
Run: python tools/how_errors_combine.py"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scipy import stats
from met.evidence import Curve
from met.full import isotropic
rng = np.random.default_rng(231123)

# ---------- A. direct readings of one mass (acceleration exact): x_i = m + e_i, sd sigma ----------
def law(xs, sig, kind, grid):
    xs = np.asarray(xs, float)[:, None]
    if kind == "gauss":   lg = -0.5 * np.sum((xs - grid) ** 2, 0) / sig ** 2
    if kind == "add":     lg = -np.sum(np.abs(xs - grid), 0) / (sig / math.sqrt(2))     # Laplace with sd sigma
    w = np.exp(lg - lg.max()); return w / w.sum()
def central(grid, w, level=0.95):
    c = np.cumsum(w); return float(grid[np.searchsorted(c, (1 - level) / 2)]), float(grid[np.searchsorted(c, 1 - (1 - level) / 2)])
print("A. two direct readings, each +- sigma; the law of the mass given the readings (exact under each error law)")
for name, xs in (("near 5.5, 5.6", (5.5, 5.6)), ("far  5.0, 6.0", (5.0, 6.0))):
    for sig in (0.1, 0.01):
        g = np.linspace(4.0, 7.0, 600001)
        a, b = central(g, law(xs, sig, "gauss", g)); c, d = central(g, law(xs, sig, "add", g))
        print(f"   {name}, sigma={sig}:  squares [{a:.3f}, {b:.3f}] (+-{(b-a)/2:.3f})   lengths [{c:.3f}, {d:.3f}] (+-{(d-c)/2:.3f})")

# the same additive statistic judged under the DECLARED GAUSSIAN (nothing but sigma and the bell believed):
# drop D(m) = sum|x_i - m| - min; for two readings D = 0 between them, 2*distance outside. Exact 95% cut:
q = stats.norm.isf(math.sqrt(0.025))            # P(both errors on one side and the nearer beyond q) = 2 Phi(-q)^2 = 0.05
e = rng.normal(size=(400000, 2)); lo, hi = e.min(1) - q, e.max(1) + q      # truth 0, sigma 1
cov = (lo <= 0) & (0 <= hi); gap = np.abs(e[:, 0] - e[:, 1])
print(f"B. under the declared Gaussian, the set [lower reading - {q:.3f} sigma, upper reading + {q:.3f} sigma]:")
print(f"   covers the truth in {cov.mean():.4f}; among pairs closer than 0.5 sigma: {cov[gap<0.5].mean():.3f}; 0.5-2: {cov[(gap>=.5)&(gap<2)].mean():.3f}; further than 2 sigma: {cov[gap>2].mean():.3f}")
print(f"   near pair -> [{5.5-q*0.1:.3f}, {5.6+q*0.1:.3f}]   far pair -> [{5-q*0.1:.3f}, {6+q*0.1:.3f}]   (the Gaussian's own: +-0.139 both)")

# ---------- C. E14-6's two series, which have the same Gaussian curve and need a count to tell apart ----------
N, D = 50, 1.0
sA = D / math.sqrt(2 * (N - 1)); A = np.array([5.0, 6.0])
sB = D / 2 * math.sqrt(N / (N - 1)); z = rng.normal(size=N); z = (z - z.mean()) / z.std(ddof=1); B = 5.5 + sB * z
g = np.linspace(4.5, 6.5, 400001)
for nm, xs, s in (("A: two readings far apart, small sigma", A, sA), ("B: fifty honest readings, sigma 0.5 ", B, sB)):
    a, b = central(g, law(xs, s, "gauss", g)); c, d = central(g, law(xs, s, "add", g))
    print(f"C. {nm}: squares +-{(b-a)/2:.3f}   lengths +-{(d-c)/2:.3f}")

# ---------- D. ninety-nine honest readings and one stray ----------
xs = np.concatenate([5.0 + 0.1 * rng.normal(size=99), [10.0]]); g = np.linspace(4.8, 5.3, 500001)
a, b = central(g, law(xs, 0.1, "gauss", g)); c, d = central(g, law(xs, 0.1, "add", g))
print(f"D. 99 readings at 5 +- 0.1 and one at 10: squares [{a:.3f}, {b:.3f}]   lengths [{c:.3f}, {d:.3f}]   (99 alone, squares: +-{1.96*0.1/math.sqrt(99):.3f})")

# ---------- E. what it costs on honest data ----------
for n in (2, 5, 20, 100):
    x = rng.normal(size=(4000, n)); g = np.linspace(-4, 4, 4001); hit = 0; wid = 0.0
    for row in x:
        lo_, hi_ = central(g, law(row, 1.0, "add", g)); hit += lo_ <= 0 <= hi_; wid += hi_ - lo_
    print(f"E. n={n}: honest Gaussian readings; the lengths-law's 95% interval covers {hit/4000:.3f}, mean width {wid/4000/ (2*1.96/math.sqrt(n)):.2f} x the Gaussian's")
def lengths(readings, masses, scale=1.0):          # sum over readings of the leftover's LENGTH (not its square)
    return sum(np.sqrt(np.maximum(Curve([r], scale=scale)(masses), 0.0)) for r in readings)
def squares(readings, masses, scale=1.0):
    return Curve(readings, scale=scale)(masses)
ms = np.array([4.8, 5.0, 5.25, 5.5, 5.75, 6.0, 6.2])
for sig in (0.1, 0.001):
    far = [isotropic([5.0], [1.0], sig, sig), isotropic([6.0], [1.0], sig, sig)]        # force and acceleration both noisy
    L = lengths(far, ms); print(f"F. far pair as (F, a) readings, sigma_F = sigma_a = {sig}: lengths curve above its least at m = {list(ms)}:\n   ", np.round(L - L.min(), 4), "| squares:", np.round(squares(far, ms) - squares(far, ms).min(), 1))
far = [isotropic([5.0], [1.0], 0.1, 0.1), isotropic([6.0], [1.0], 0.1, 0.1)]
base = lengths(far, ms)
zeros = [isotropic([0.0], [0.0], 3.0, 0.02)] * 20; tiny = [isotropic([1e-9], [1e-9], 0.1, 0.1)] * 20
print("G. twenty exact zeros beside the far pair change the lengths curve by", float(np.max(np.abs(lengths(far + zeros, ms) - base))), "; twenty readings of size 1e-9 by at most", f"{float(np.max(np.abs((lengths(far + tiny, ms) - lengths(far+tiny, ms).min()) - (base - base.min())))):.1e}")
right = [isotropic([2.0], [0.0], 1.0, 1.0), isotropic([0.0], [2.0], 1.0, 1.0)]
mm = np.array([0.0, 0.3, 1.0, 3.0, math.inf]); print("H. force-only beside acceleration-only reading: squares", np.round(squares(right, mm), 3), " lengths", np.round(lengths(right, mm), 3))
# one 3D reading whose force and acceleration are 30 degrees apart, sigma shrinking: width of the set within 2 of the least
th = np.linspace(0.6, 1.5, 200001); m = np.tan(th)
for sig in (0.1, 0.01, 0.001):
    r = [isotropic([5.0, 0.0, 0.0], [math.cos(math.pi / 6), math.sin(math.pi / 6), 0.0], sig, sig)]
    out = []
    for f, cut in ((squares, 3.84), (lengths, 3.0)):
        v = f(r, m); keep = m[v - v.min() <= cut]; out.append(keep.max() - keep.min())
    print(f"I. one misaligned 3D reading, sigma={sig}: width kept by squares {out[0]:.4f}, by lengths {out[1]:.4f}")

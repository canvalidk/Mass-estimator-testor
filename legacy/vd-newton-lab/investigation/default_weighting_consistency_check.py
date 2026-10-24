"""Exploratory check (2026-09-23): which repeated-reading weightings converge to the true mass?

N readings of one mass, each with its own unknown excitation. One dimension, known positive
direction, sigma_F = sigma_a = 1 (so s = sigma_F/sigma_a = 1). For each candidate weighting the
large-N limit of the posterior mode is found by maximising the per-reading expected log posterior,
estimated from R simulated readings. Output: limit / true mass (1.00 = converges to the truth).

Candidates:
  symmetric      per reading (1 + m^2/s^2) alpha d(alpha) (theta-uniform prior counted once)
  k=...          instrument-free, channel-symmetric power family m^{N(k+1)-1} prod alpha_i^{2k+1}
                 (k = 0 is contribution 12's eq. 20)
  eq18 (flat m)  prod alpha_i d(alpha_i), flat in m (contribution 12's eq. 18)

Exploratory only: 1D, fixed excitation per reading, product-form weightings only, numerical
(not proved). Usage: python default_weighting_consistency_check.py [R] [seed]
"""
import sys
import numpy as np
from scipy import optimize

T = np.linspace(-12, 12, 1201)
DT = T[1] - T[0]


def log_ell(x, y, m, b):
    """log int_0^inf a^b exp(-(x - m a)^2/2 - (y - a)^2/2) da, per reading."""
    q = 1 + m * m
    c = (m * x + y) / q
    const = -(x * x + y * y) / 2 + (m * x + y) ** 2 / (2 * q)
    a = c[:, None] + T / np.sqrt(q)
    with np.errstate(divide="ignore", invalid="ignore"):
        la = np.where(a > 0, b * np.log(np.where(a > 0, a, 1.0)), -np.inf)
    e = la - T ** 2 / 2
    mx = np.max(e, 1, keepdims=True)
    return const + mx[:, 0] + np.log(np.exp(e - mx).sum(1) * DT / np.sqrt(q))


def limits(m0, a0, R, rng):
    a = np.full(R, float(a0))
    x = m0 * a + rng.normal(size=R)
    y = a + rng.normal(size=R)
    rules = {"symmetric": lambda m: np.mean(log_ell(x, y, m, 1)) + np.log(1 + m * m)}
    for k in (-0.75, -0.5, -0.25, 0.0, 0.5, 1.0):
        rules[f"k={k}"] = (lambda kk: (lambda m: np.mean(log_ell(x, y, m, 2 * kk + 1)) + (kk + 1) * np.log(m)))(k)
    rules["eq18 (flat m)"] = lambda m: np.mean(log_ell(x, y, m, 1))
    out = {}
    for name, f in rules.items():
        r = optimize.minimize_scalar(lambda lm: -f(np.exp(lm)),
                                     bounds=(np.log(m0) - 4, np.log(m0) + 4), method="bounded")
        out[name] = np.exp(r.x) / m0
    return out


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    for m0, a0 in [(1, 1), (1, 3), (4, 1), (4, 3), (0.25, 1), (0.25, 4)]:
        o = limits(m0, a0, R, rng)
        print(f"m/s={m0:4} accSNR={a0}: " + " ".join(f"{k}:{v:.2f}" for k, v in o.items()), flush=True)

"""Judging each candidate mass on its own: three checks on the forms that avoid the fit
(Claude, 18 Nov 2026; problems/O14.md, 18 Nov, evening). Measurements; nothing in met/ changes, no rule adopted.

The standard name. One instrument, noise units: the stacked readings are two n-vectors (n = sum k_i) with means
p sin(phi) and p cos(phi) and unit noise, the pushes p unknown. That is linear instrumental variables with n
instruments and a known reduced-form covariance: the mass is the structural coefficient, the pushes the first stage.
At a candidate line, given the along-line parts t (complete and sufficient for the pushes there, T45), the leftover r
is N(0, I_n) exactly, whatever the pushes, and

    K = sum |r|^2 = S^2/T + J,      S = sum t.r,   T = sum |t|^2,
    S^2/T ~ chi2_1 (T45's score; Kleibergen's K statistic),    J ~ chi2_{n-1}, independent of it (the rest).

K ~ chi2_n is the Anderson-Rubin statistic and exp(-K/2) the Gaussian law. In the Gram's terms J = lam+ lam- / T(phi).
The score statistic is zero at the worst line as well as at the best one: T56's false peak, known for the K statistic.
The textbook repair is Moreira's conditional likelihood ratio (2003): the Gaussian law's own statistic
LR = K - K_min, referred to its exact law given T,

    LR >= x    iff    (S^2/T) / x + J / (x + T) >= 1,

which is exact for any pushes and needs no fitted share. It needs the count n: an exactly zero reading leaves the Gram
alone and raises n by its k.

Five deviances d (a law is pi * exp(-d/2)), each a function of the Gram and the candidate alone, the last also of n:

    gaussian              LR itself
    gaussian, Gram share  LR times the Gram share |Z|/lam+ read at the best line  (7 Nov; ours33's weak-signal part)
    score                 S^2/T                                        (T45, T56)
    signal energy >= 0    the afternoon's phasor quasi-likelihood      (tools/share_and_the_fit.py)
    conditional LR        the chi2_1 value with the same tail probability as LR has given T

The last line's conversion of a tail probability to a deviance is this file's choice: it is the one that returns the
Gaussian law when LR is chi2_1.

    python tools/candidate_by_candidate.py            # about two minutes
"""
import functools
import math
import os
import sys

import numpy as np
from scipy import special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discount_lab as lab                                      # noqa: E402

LAWS = ("gaussian", "gaussian, Gram share", "score", "signal energy >= 0", "conditional LR")
CHI2_95 = 3.841458820694124
SIGMA56 = np.array([[1.0, 0.9], [0.9, 1.0]])                   # T56's example
M0_56 = 1.2


@functools.lru_cache(maxsize=None)
def _nodes(n_rest):
    """Nodes and weights for an expectation over J ~ chi2_{n_rest}: Gauss-Legendre in the probability up to
    1 - 1e-3, and in the log of the upper tail probability from 1e-3 down to 1e-300."""
    x, w = np.polynomial.legendre.leggauss(240)
    u = 0.5 * (x + 1)
    s, ws = u * (1 - 1e-3), 0.5 * w * (1 - 1e-3)
    a = 3 + u * 297
    wa = 0.5 * w * 297 * math.log(10) * 10.0 ** (-a)
    return (np.concatenate([stats.chi2.ppf(s, n_rest), stats.chi2.isf(10.0 ** (-a), n_rest)]),
            np.concatenate([ws, wa]))


def clr_p(x, T, n):
    """P(LR >= x | T) at a candidate: the exact conditional tail probability, any pushes. n = sum k_i."""
    x, T = np.broadcast_arrays(np.asarray(x, float), np.asarray(T, float))
    if n <= 1:
        return special.erfc(np.sqrt(x / 2))
    J, w = _nodes(n - 1)
    out = np.empty(x.shape)
    flat_x, flat_T, flat_o = x.ravel(), T.ravel(), out.ravel()
    for i in range(0, flat_x.size, 20000):
        xs, Ts = flat_x[i:i + 20000, None], flat_T[i:i + 20000, None]
        arg = xs * np.clip(1 - J / np.maximum(xs + Ts, 1e-300), 0, None)
        flat_o[i:i + 20000] = special.erfc(np.sqrt(arg / 2)) @ w
    return out


def clr_deviance(x, T, n):
    return 2 * special.erfcinv(np.clip(clr_p(x, T, n), 1e-300, 1.0)) ** 2


def deviances(lam_plus, lam_minus, s2, n):
    """The five deviances at candidates whose line makes an angle with the best line of squared sine s2."""
    s2 = np.asarray(s2, float)
    H, E = lam_plus - lam_minus, lam_plus + lam_minus
    cos_t = 1 - 2 * s2
    T = lam_plus * (1 - s2) + lam_minus * s2
    D = -E + np.sqrt(np.maximum(E * E + 2 * H * E * cos_t + H * H, 0.0))
    if s2.size > 6000:                                      # tabulate the conditional LR, dense near the best line
        tab = np.linspace(0, 1, 4001) ** 2
        clr = np.interp(s2, tab, clr_deviance(H * tab, lam_plus * (1 - tab) + lam_minus * tab, n))
    else:
        clr = clr_deviance(H * s2, T, n)
    return {"gaussian": H * s2,
            "gaussian, Gram share": H * H * s2 / max(lam_plus, 1e-300),
            "score": H * H * s2 * (1 - s2) / np.maximum(T, 1e-300),
            "signal energy >= 0": np.where(D >= 0, D - H * cos_t, H * H / (2 * max(E, 1e-300))),
            "conditional LR": clr}


def law_on_masses(G, n, W, c=1.0, grid=20001, weight=1.0):
    """Per law: median, central 95% interval, and the grid, for a whitened 2x2 Gram G of n components. W whitens (F, a);
    the prior is the half-Cauchy of centre c (uniform in its own angle). `weight` multiplies the Gram (one weight)."""
    th = (np.arange(grid) + 0.5) / grid * (math.pi / 2)
    m = c * np.tan(th)
    U = W @ np.vstack([m, np.ones_like(m)])
    U /= np.linalg.norm(U, axis=0)
    vals, vecs = np.linalg.eigh(weight * np.asarray(G, float))
    s2 = np.clip(1 - (U[0] * vecs[0, 1] + U[1] * vecs[1, 1]) ** 2, 0, 1)
    out = {}
    for name, dev in deviances(vals[1], max(vals[0], 0.0), s2, n).items():
        w = np.exp(-0.5 * (dev - dev.min()))
        cdf = np.cumsum(w)
        cdf /= cdf[-1]
        lo, med, hi = np.interp([0.025, 0.5, 0.975], cdf, th)
        out[name] = dict(median=c * math.tan(med), lo=c * math.tan(lo), hi=c * math.tan(hi), angle=hi - lo, cdf=cdf, m=m)
    return out


# ---- 1. is the candidate's own verdict exact? ------------------------------------------------------------------------

def size_at_truth(n, amp, reps=100000, seed=18):
    """Scalar readings, identity noise, true line at 45 degrees, every push = amp. Per law: how often the true mass's
    deviance exceeds the chi2_1 95% point. Exact is 0.05."""
    rng = np.random.default_rng([seed, n, int(round(100 * amp))])
    F, A = amp + rng.normal(size=(reps, n)), amp + rng.normal(size=(reps, n))
    t, r = (F + A) / math.sqrt(2), (F - A) / math.sqrt(2)
    T, K, S = (t * t).sum(1), (r * r).sum(1), (t * r).sum(1)
    H = np.sqrt((T - K) ** 2 + 4 * S * S)
    s2 = np.clip(0.5 * (1 - (T - K) / np.maximum(H, 1e-300)), 0, 1)
    lam_plus, lam_minus = 0.5 * (T + K + H), np.maximum(0.5 * (T + K - H), 0.0)
    out = {}
    for name in LAWS:
        if name == "conditional LR":
            dev = clr_deviance(H * s2, T, n)
        elif name == "gaussian":
            dev = H * s2
        elif name == "gaussian, Gram share":
            dev = H * H * s2 / np.maximum(lam_plus, 1e-300)
        elif name == "score":
            dev = S * S / np.maximum(T, 1e-300)
        else:
            E, cos_t = T + K, 1 - 2 * s2
            D = -E + np.sqrt(E * E + 2 * H * E * cos_t + H * H)
            dev = np.where(D >= 0, D - H * cos_t, H * H / (2 * E))
        out[name] = float((dev > CHI2_95).mean())
    return out


def one_reading(energy=9.0, degrees=(10, 30, 60, 85)):
    """One scalar reading: the Gaussian factor is exact. Each law's deviance over the Gaussian's, by the angle between
    the candidate and the reading."""
    s2 = np.sin(np.radians(degrees)) ** 2
    d = deviances(energy, 0.0, s2, 1)
    return {name: d[name] / d["gaussian"] for name in LAWS}


# ---- 2. noise only, and a history that identifies ever more slowly ---------------------------------------------------

def noise_only(n, n_series=200, seed=19):
    """n scalar readings of pure noise. Per law: the median, over series, of the deviance's range over all lines, and
    of the law's central 95% width in the prior's own angle over the prior's."""
    rng = np.random.default_rng([seed, n])
    tab = np.linspace(0, 1, 801) ** 2
    acc = {name: ([], []) for name in LAWS}
    for _ in range(n_series):
        Y = rng.normal(size=(n, 2))
        G = Y.T @ Y
        vals = np.linalg.eigvalsh(G)
        dev = deviances(vals[1], vals[0], tab, n)
        laws = law_on_masses(G, n, np.eye(2), grid=4001)
        for name in LAWS:
            acc[name][0].append(float(np.ptp(dev[name])))
            acc[name][1].append(laws[name]["angle"] / (0.95 * math.pi / 2))
    return {name: (float(np.median(a)), float(np.median(b))) for name, (a, b) in acc.items()}


def history(n, beta=1.0, n_series=200, seed=20, zeros=0, grid=20001):
    """T56's example (covariance SIGMA56, true mass 1.2, false peak at 4/15) with accelerations a_i = i^((beta-1)/2):
    beta = 1 is T56's own, beta < 1 the fading history of T58. `zeros` exactly zero readings are added to the count.
    Per law: mass below m = 0.6, how often the truth is inside the 95% interval, the median log width."""
    rng = np.random.default_rng([seed, n, int(round(100 * beta))])
    W = np.linalg.cholesky(np.linalg.inv(SIGMA56)).T
    C = np.linalg.cholesky(SIGMA56)
    amp = np.arange(1, n + 1) ** (0.5 * (beta - 1))
    acc = {name: ([], [], []) for name in LAWS}
    for _ in range(n_series):
        Y = (amp[:, None] * np.array([M0_56, 1.0])[None, :] + rng.normal(size=(n, 2)) @ C.T) @ W.T
        for name, r in law_on_masses(Y.T @ Y, n + zeros, W, grid=grid).items():
            acc[name][0].append(float(r["cdf"][np.searchsorted(r["m"], 0.6)]))
            acc[name][1].append(r["lo"] <= M0_56 <= r["hi"])
            acc[name][2].append(math.log(r["hi"] / r["lo"]))
    return {name: dict(false_side=float(np.mean(a)), inside=float(np.mean(b)), log_width=float(np.median(c)))
            for name, (a, b, c) in acc.items()}


def weighted_leftover(n, n_series=200, seed=21):
    """T56's example. Q = sum t_i^2 (r_i^2 - 1): given the along-line parts its mean is 0 and its variance
    2 sum t_i^4 exactly, for any pushes, and a reading with zero along-line part is not in it at all. Returns the
    median of Q over its exact standard deviation at the true mass and at the score law's false peak, and the same
    for the unweighted leftover (K - n) / sqrt(2 n), which the count enters."""
    rng = np.random.default_rng([seed, n])
    W = np.linalg.cholesky(np.linalg.inv(SIGMA56)).T
    C = np.linalg.cholesky(SIGMA56)
    out = {"truth": ([], []), "false peak": ([], [])}
    for _ in range(n_series):
        Y = (np.array([M0_56, 1.0])[None, :] + rng.normal(size=(n, 2)) @ C.T) @ W.T
        vals, vecs = np.linalg.eigh(Y.T @ Y)
        u0 = W @ np.array([M0_56, 1.0])
        u0 /= np.linalg.norm(u0)
        for key, u in (("truth", u0), ("false peak", vecs[:, 0])):          # the false peak: the worst line
            t, r = Y @ u, Y @ np.array([-u[1], u[0]])
            out[key][0].append(float((t * t * (r * r - 1)).sum() / math.sqrt(2 * (t ** 4).sum())))
            out[key][1].append(float(((r * r).sum() - n) / math.sqrt(2 * n)))
    return {key: (float(np.median(a)), float(np.median(b))) for key, (a, b) in out.items()}


# ---- 3. the right-angle scene: a no-fit form carries no discount ------------------------------------------------------

def right_angle_scene(sigma, discounted, body="one weight, energy less 2, own best line"):
    """Scene C (force 2 with no acceleration; acceleration 1 with no force; each +- sigma). Per law the median, with
    the Gram as declared or multiplied by the morning's one weight (tools/discount_lab.py)."""
    readings = [([2.0], [0.0]), ([0.0], [1.0])]
    ws = lab.whiten(readings, sigma, sigma)
    w = float(lab.weights(ws, body)[0][0]) if discounted else 1.0
    G = sum(lab.gram(f, al) for f, al in ws)
    laws = law_on_masses(G, 2, np.diag([1 / sigma, 1 / sigma]), weight=w)
    return w, {name: laws[name]["median"] for name in LAWS}


# ---- 4. is the discount the along-line-weighted leftover read at its least? ---------------------------------------

def leftover_at_its_least(readings, sF, sA, grid=200001):
    """Observed over expected, least over positive masses, three ways of reading the leftover:
      squared scores       sum (t_i.r_i)^2          over  sum |t_i|^2          (17 Nov; White's information comparison)
      along-line weighted  sum |t_i|^2 |r_i|^2      over  sum k_i |t_i|^2      (exact given t; no zero reading in it)
      energy weighted      sum (C_i - 2) K_i        over  sum k_i C_i          (the morning's; exact in the mean only)
    For scalar readings the first two are the same number."""
    ws = lab.whiten(readings, sF, sA)
    phi = (np.arange(grid) + 0.5) / grid * (math.pi / 2)
    s, c = np.sin(phi), np.cos(phi)
    B, T, Qo, Qe = (np.zeros(grid) for _ in range(4))
    for f, al in ws:
        t = s[:, None] * f[None, :] + c[:, None] * al[None, :]
        r = c[:, None] * f[None, :] - s[:, None] * al[None, :]
        tt, rr = (t * t).sum(1), (r * r).sum(1)
        B += (t * r).sum(1) ** 2
        T += tt
        Qo += tt * rr
        Qe += len(f) * tt
    least, E0, _ = lab.energy_less_two(ws, corrected=False)
    return {"squared scores": float((B / T).min()), "along-line weighted": float((Qo / Qe).min()),
            "energy weighted": float(least / E0)}


def nearly_zero_readings(n=20, extra=20, size=1e-6, n_series=200, seed=20):
    """T56's example: the conditional LR's median log width with nothing added; with `extra` exactly zero readings
    that the count leaves out (a rule that restores Z at exact zero); and with `extra` readings of `size` noise
    units, which any count must take in."""
    rng = np.random.default_rng([seed, n, 100])
    W = np.linalg.cholesky(np.linalg.inv(SIGMA56)).T
    C = np.linalg.cholesky(SIGMA56)
    out = {"nothing added": [], "exact zeros, left out of the count": [], "nearly zero, counted": []}
    for _ in range(n_series):
        Y = (np.array([M0_56, 1.0])[None, :] + rng.normal(size=(n, 2)) @ C.T) @ W.T
        tiny = np.vstack([Y, size * rng.normal(size=(extra, 2))])
        for key, (G, count) in (("nothing added", (Y.T @ Y, n)), ("exact zeros, left out of the count", (Y.T @ Y, n)),
                                ("nearly zero, counted", (tiny.T @ tiny, n + extra))):
            r = law_on_masses(G, count, W)["conditional LR"]
            out[key].append(math.log(r["hi"] / r["lo"]))
    return {key: float(np.median(v)) for key, v in out.items()}


def main():
    print("1. How often the true mass is rejected at the 5% level (exact = 0.050); 100000 series a cell")
    print("| readings | push | " + " | ".join(LAWS) + " |")
    print("|---|---|" + "---|" * len(LAWS))
    for n in (1, 2, 5, 20, 200):
        for amp in (0.0, 0.5, 2.0):
            r = size_at_truth(n, amp)
            print(f"| {n} | {amp} | " + " | ".join(f"{r[name]:.3f}" for name in LAWS) + " |")
    r = one_reading()
    print("\n   one scalar reading, deviance over the Gaussian's at 10, 30, 60, 85 degrees from the reading:")
    for name in LAWS:
        print(f"   {name:20s} " + "  ".join(f"{v:.3f}" for v in r[name]))

    print("\n2a. Noise only (the tilt stays bounded; O15, 9 Nov): median range of the deviance over all lines | median 95% width over the prior's; 200 series")
    print("| readings | " + " | ".join(LAWS) + " |")
    print("|---|" + "---|" * len(LAWS))
    for n in (5, 50, 500, 5000):
        r = noise_only(n)
        print(f"| {n} | " + " | ".join(f"{r[name][0]:.2f}, {r[name][1]:.2f}" for name in LAWS) + " |")

    for beta, label in ((1.0, "2b. T56's example"), (0.75, "2c. the fading history, signal energy ~ N^0.75")):
        print(f"\n{label}: mass below m = 0.6 | truth inside the 95% interval | median log width; 200 series")
        print("| readings | " + " | ".join(LAWS) + " |")
        print("|---|" + "---|" * len(LAWS))
        for n in (20, 200, 2000, 20000) if beta == 1.0 else (100, 1000, 10000, 100000):
            r = history(n, beta)
            print(f"| {n} | " + " | ".join(f"{r[name]['false_side']:.3f}, {r[name]['inside']:.2f}, {r[name]['log_width']:.4f}"
                                           for name in LAWS) + " |")

    print("\n2d. T56's example, 20 readings: median log width with exactly zero readings added (200 series)")
    for z in (0, 20, 200):
        r = history(20, 1.0, zeros=z)
        print(f"   {z:3d} zeros: " + "  ".join(f"{name} {r[name]['log_width']:.3f}" for name in LAWS))

    print("\n2e. T56's example: an exact statistic no zero reading enters, sum t^2 (r^2 - 1), in its own standard deviations"
          " | the unweighted leftover, which the count enters")
    for n in (20, 200, 2000):
        r = weighted_leftover(n)
        print(f"   {n:5d} readings: at the truth {r['truth'][0]:+.2f} | {r['truth'][1]:+.2f};  "
              f"at the false peak {r['false peak'][0]:+.2f} | {r['false peak'][1]:+.2f}")

    print("\n3. Scene C: median as sigma shrinks, Gram as declared | times the morning's one weight")
    print("| sigma | weight | " + " | ".join(LAWS) + " |")
    print("|---|---|" + "---|" * len(LAWS))
    for sigma in (1.0, 0.2, 0.05, 0.01):
        _, plain = right_angle_scene(sigma, False)
        w, disc = right_angle_scene(sigma, True)
        print(f"| {sigma} | {w:.3g} | " + " | ".join(f"{plain[name]:.3g}, {disc[name]:.3g}" for name in LAWS) + " |")


def main_least():
    print("\n4. Observed over expected leftover at its least over positive masses")
    print("| scene | squared scores | along-line weighted | energy weighted |")
    print("|---|---|---|---|")
    far = [([5.0], [1.0]), ([6.0], [1.0])]
    rows = [("B: 5 and 6, each +- 0.1", far, 0.1, 0.002)]
    for sigma in (0.2, 0.05, 0.01):
        rows.append((f"C as two scalar readings, sigma = {sigma}", [([2.0], [0.0]), ([0.0], [1.0])], sigma, sigma))
    for sigma in (0.2, 0.05, 0.01):
        rows.append((f"the same numbers as one 2D reading, sigma = {sigma}", [([2.0, 0.0], [0.0, 1.0])], sigma, sigma))
    for name, readings, sF, sA in rows:
        r = leftover_at_its_least(readings, sF, sA)
        print(f"| {name} | {r['squared scores']:.3g} | {r['along-line weighted']:.3g} | {r['energy weighted']:.3g} |")
    r = nearly_zero_readings()
    print("\n   T56's example, 20 readings, the conditional LR's median log width: "
          + "; ".join(f"{k} {v:.3f}" for k, v in r.items()))


if __name__ == "__main__":
    main()
    main_least()

"""Discount lab (Claude, 18 Nov 2026; problems/O14.md, 18 Nov). Measurements on ways to read the leftover, on one
instrument. Nothing in met/ changes and no rule is adopted here.

One instrument: every reading has force sd sF and acceleration sd sA (isotropic, uncorrelated). In noise units
f = F/sF, al = a/sA, a mass m is the line at angle phi with m = (sF/sA) tan(phi) (T36, T70). For reading i (any k_i):

    t_i(phi) = sin(phi) f_i + cos(phi) al_i     along the line
    r_i(phi) = cos(phi) f_i - sin(phi) al_i     left over;   K_i(phi) = |r_i|^2  (T36's residual, the leftover)
    C_i      = |f_i|^2 + |al_i|^2               the reading's energy (equation.md's E_i)

The law is pi(m) exp(-1/2 sum_i w_i K_i(phi)) on m > 0, pi half-Cauchy with centre c. The bodies differ in the weights,
all read once, at the best positive line of the summed Gaussian leftover (the anchor):

    gauss                  w = 1
    eta as on the pages    curvature over squared scores at the whole-line best line, one cap (T31), as a power
    one weight, by count   (sum k_i - 1) / total leftover: T30's rule on the series stacked as one reading (Birge 1932)
    one weight, by leverage     sum I_i (k_i - h_i) / sum I_i K_i, I_i = |t_i|^2 at the anchor: eta's second factor
                                with the fit correction (O14, 17 Nov, B22 read; MacKinnon-White HC2)
    one weight, by energy       sum C_i (k_i - h_i) / sum C_i K_i
    one weight, by signal energy   sum (C_i - 2k_i + 2)(k_i - h_i) / sum (C_i - 2k_i) K_i   (tools/signal_energy_identity.py)
    one weight per reading      (k_i - h_i) / K_i for each reading (O14, 17 Nov evening, with the fit-corrected count)

h_i is reading i's share of the along-line energy at the anchor (0 when positivity binds): the expected leftover of a
reading at the fitted line is k_i - h_i, not k_i. Each weight is capped at 1.

    python tools/discount_lab.py            # the scenes of O14, 18 Nov, under every body
"""
import math

import numpy as np

BODIES = ("gauss", "eta as on the pages", "one weight, by count", "one weight, by leverage", "one weight, by energy",
          "one weight, by signal energy", "one weight per reading")
TINY = 1e-9


def whiten(readings, sF, sA):
    return [(np.atleast_1d(np.asarray(F, float)) / sF, np.atleast_1d(np.asarray(a, float)) / sA) for F, a in readings]


def gram(f, al):
    return np.array([[f @ f, f @ al], [f @ al, al @ al]])


def K_of(M, phi):
    s, c = np.sin(phi), np.cos(phi)
    return M[0, 0] * c * c - 2 * M[0, 1] * s * c + M[1, 1] * s * s


def H_of(M, phi):
    s, c = np.sin(phi), np.cos(phi)
    return M[0, 0] * s * s + 2 * M[0, 1] * s * c + M[1, 1] * c * c


def anchor_positive(Ms):
    """Best line of the summed Gaussian leftover on the closed positive quadrant, and whether it is interior."""
    S = sum(Ms)
    vals, vecs = np.linalg.eigh(S)
    if vals[1] - vals[0] <= 1e-12 * max(vals[1], 1e-300):
        return math.pi / 4, False, S                       # every line alike
    v = vecs[:, 1] * (1 if vecs[1, 1] >= 0 else -1)
    phi = math.atan2(v[0], v[1])
    if -1e-12 <= phi <= math.pi / 2 + 1e-12:
        return min(max(phi, 0.0), math.pi / 2), True, S
    return (0.0 if K_of(S, 0.0) <= K_of(S, math.pi / 2) else math.pi / 2), False, S


def _ratio(num, den, scale):
    """Expected over observed, capped at 1; nan where both vanish (0/0)."""
    if abs(num) <= TINY * scale and abs(den) <= TINY * scale:
        return math.nan
    return 1.0 if den <= max(num, 0.0) else max(num, 0.0) / den


def weights(ws, body):
    """Per-reading weights (equal for the one-weight bodies) and what was read at the anchor."""
    Ms = [gram(f, al) for f, al in ws]
    ks = np.array([len(f) for f, _ in ws], float)
    n = len(ws)
    phi, interior, S = anchor_positive(Ms)
    K = np.array([K_of(M, phi) for M in Ms])
    I = np.array([H_of(M, phi) for M in Ms])
    C = np.array([np.trace(M) for M in Ms])
    h = I / I.sum() if (interior and I.sum() > 0) else np.zeros(n)
    kappa = ks - h
    scale = max(float(C.sum()), 1.0)
    info = dict(anchor=phi, interior=interior, leftover=float(K.sum()), K=K, C=C, kappa=kappa)
    if body == "gauss":
        return np.ones(n), info
    if body == "eta as on the pages":
        vals, vecs = np.linalg.eigh(S)
        ph = math.atan2(vecs[0, 1], vecs[1, 1])
        s, c = math.sin(ph), math.cos(ph)
        J = sum(((c * f - s * al) @ (s * f + c * al)) ** 2 for f, al in ws)
        curv = vals[1] - vals[0]
        return np.full(n, 1.0 if (J <= 1e-300 or curv <= 0) else min(1.0, curv / J)), info
    if body == "one weight, by count":
        return np.full(n, _ratio(kappa.sum(), K.sum(), 1.0)), info
    if body == "one weight, by leverage":
        return np.full(n, _ratio((I * kappa).sum(), (I * K).sum(), scale ** 2)), info
    if body == "one weight, by energy":
        return np.full(n, _ratio((C * kappa).sum(), (C * K).sum(), scale ** 2)), info
    if body == "one weight, by signal energy":
        return np.full(n, _ratio(((C - 2 * ks + 2) * kappa).sum(), ((C - 2 * ks) * K).sum(), scale ** 2)), info
    if body == "one weight per reading":
        return np.array([_ratio(kp, kk, scale) for kp, kk in zip(kappa, K)]), info
    raise ValueError(body)


def swing(ws, w):
    """How far -1/2 sum w_i K_i swings over all lines (twice the log-law's range): the eigengap of the weighted Gram."""
    return float(np.ptp(np.linalg.eigvalsh(sum(wi * gram(f, al) for (f, al), wi in zip(ws, w)))))


def law(readings, sF, sA, c, body, grid=200001):
    """Median and central 95% interval of the law of m; the weights; whether a weight was 0/0 (then taken as 1)."""
    ws = whiten(readings, sF, sA)
    w, info = weights(ws, body)
    undefined = bool(np.isnan(w).any())
    w = np.where(np.isnan(w), 1.0, w)
    th = (np.arange(grid) + 0.5) / grid * (math.pi / 2)      # uniform in the half-Cauchy's own angle
    m = c * np.tan(th)
    phi = np.arctan(m / (sF / sA))
    lp = np.zeros(grid)
    for (f, al), wi in zip(ws, w):
        lp -= 0.5 * wi * K_of(gram(f, al), phi)
    cdf = np.cumsum(np.exp(lp - lp.max()))
    lo, med, hi = np.interp([0.025, 0.5, 0.975], cdf / cdf[-1], th)
    return dict(point=c * math.tan(med), lo=c * math.tan(lo), hi=c * math.tan(hi), w=w, undefined=undefined, **info)


def ordinary_share(body, k, n, signal, reps=100, seed=5, truth=1.0):
    """Ordinary data (the declared sigma true, n readings with k components at one push size, noise units): the share
    of series the body discounts and the median weight when it does. One-weight bodies only."""
    rng = np.random.default_rng(seed)
    phi0 = math.atan(truth)
    ws_all = []
    for _ in range(reps):
        d = rng.normal(size=(n, k))
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        f = signal * math.sin(phi0) * d + rng.normal(size=(n, k))
        al = signal * math.cos(phi0) * d + rng.normal(size=(n, k))
        w, _ = weights(list(zip(f, al)), body)
        ws_all.append(1.0 if np.isnan(w[0]) else float(w[0]))
    ws_all = np.array(ws_all)
    low = ws_all < 1 - 1e-12
    return float(low.mean()), float(np.median(ws_all[low])) if low.any() else 1.0


def scenes(seed=3):
    """name -> (readings as (F, a) pairs, sigmas to run, prior centre, sA/sF)."""
    rng = np.random.default_rng(seed)
    s60, c60 = math.sin(math.radians(60)), math.cos(math.radians(60))
    shrink = (1.0, 0.2, 0.05, 0.01)
    far = [([5.0], [1.0]), ([6.0], [1.0])]
    nulls = [([0.1 * rng.normal()], [0.002 * rng.normal()]) for _ in range(40)]
    ordinary = [([5 + 0.1 * rng.normal()], [1.0]) for _ in range(99)]
    return {
        "C: force 2, no acceleration; acceleration 1, no force": ([([2.0], [0.0]), ([0.0], [1.0])], shrink, 1.0, 1.0),
        "C': two readings (2, 0) and one (0, 1)": ([([2.0], [0.0])] * 2 + [([0.0], [1.0])], shrink, 1.0, 1.0),
        "C'': two readings on m = 1.73 and one at right angles": ([([2 * s60], [2 * c60])] * 2 + [([-c60], [s60])], shrink, 1.0, 1.0),
        "one 2D reading F = (2, 0), a = (0, 1)": ([([2.0, 0.0], [0.0, 1.0])], shrink, 1.0, 1.0),
        "the same with an exactly zero third axis": ([([2.0, 0.0, 0.0], [0.0, 1.0, 0.0])], shrink, 1.0, 1.0),
        "anti-parallel scalar reading F = 1, a = -1": ([([1.0], [-1.0])], shrink, 1.0, 1.0),
        "aligned pair (2, 1) and (4, 2)": ([([2.0], [1.0]), ([4.0], [2.0])], shrink, 1.0, 1.0),
        "A: 5.5, 5.6 (each +- 0.1)": ([([5.5], [1.0]), ([5.6], [1.0])], (0.1,), 5.0, 0.02),
        "B: 5, 6": (far, (0.1,), 5.0, 0.02),
        "5, 5.05, 6": ([([5.0], [1.0]), ([5.05], [1.0]), ([6.0], [1.0])], (0.1,), 5.0, 0.02),
        "B and 40 noise-only readings": (far + nulls, (0.1,), 5.0, 0.02),
        "B and 40 exactly zero readings": (far + [([0.0], [0.0])] * 40, (0.1,), 5.0, 0.02),
        "99 ordinary near 5": (ordinary, (0.1,), 5.0, 0.02),
        "99 ordinary and one at 9": (ordinary + [([9.0], [1.0])], (0.1,), 5.0, 0.02),
    }


def _fmt(x):
    return f"{x:.4g}" if x < 1e6 else f"{x:.2e}"


def main():
    for name, (readings, sigmas, c, ratio) in scenes().items():
        print(f"\n### {name}   (sigma = {', '.join(str(s) for s in sigmas)})")
        for body in BODIES:
            cells = []
            for s in sigmas:
                r = law(readings, s, s * ratio, c, body)
                w = r["w"]
                wtxt = "0/0" if r["undefined"] else (f"{w[0]:.3g}" if np.allclose(w, w[0]) else "(" + ", ".join(f"{x:.2g}" for x in w[:3]) + (", ..." if len(w) > 3 else "") + ")")
                cells.append(f"{_fmt(r['point'])} [{_fmt(r['lo'])}, {_fmt(r['hi'])}] w = {wtxt}")
            print(f"  {body:30s} | " + " | ".join(cells))


if __name__ == "__main__":
    main()

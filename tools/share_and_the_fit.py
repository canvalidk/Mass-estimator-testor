"""The weak-signal share and the fit; the score law's false peak as negative signal energy
(Claude, 18 Nov 2026; problems/O14.md, 18 Nov, afternoon). Measurements; nothing in met/ changes, no rule adopted.

One instrument, noise units, c = s. A series has energy E = sum C_i and phasor Z (T70): H(phi), the explained energy at
the line phi, is E/2 + (|Z|/2) cos 2(phi - phi_hat). At a candidate line

    A(phi) = |Z| cos 2(phi - phi_hat) = sum (|t_i|^2 - |r_i|^2)     mean Delta (the signal energy) at the true line
    T(phi) = (E + A(phi)) / 2      = sum |t_i|^2                    mean Delta + M there

for any pushes, with no count, and an exactly zero reading adds nothing to either. So the share Delta/(Delta + M) has
no count in it at the true line; the Gram eta of 7 Nov is A/T at the fitted line, where A = |Z| is at its largest.

1. `share_table`: on the 7 Nov world and seeds (tools/eta_from_the_gram.py), the Gram share, the same with the
   first-order fit correction |Z| - 2 T/|Z| (the leverage correction again), and the share with A and T averaged over
   the law pi(phi) exp(eta |Z| cos 2(phi - phi_hat) / 4), solved for its own power eta.
2. `false_peak`: T56's example (Codex, 4 Nov; exploration branch). The law pi exp(-T^2/2) built from the exact
   conditional score has a second peak at the profile minimum. There A = -|Z|: the candidate explains the readings
   only with negative signal energy. The Gaussian quasi-likelihood of the phasor sum (mean D e^{i psi}, variance
   4 (E + D), both exact) profiled over D >= 0 has no such peak:
       2 g(theta) = D* - |Z| cos theta,  D* = -E + sqrt(E^2 + 2 |Z| E cos theta + |Z|^2)  where that is >= 0,
       2 g(theta) = |Z|^2 / (2 E)  elsewhere;   theta = 2 (phi - phi_hat).
   Its form beyond "D >= 0" is a choice (it is not T45's pivot); what is measured is that the peak goes and the local
   width stays the Godambe one.

    python tools/share_and_the_fit.py            # about a minute
"""
import math

import numpy as np

KS, NS, SNRS = (1, 3), (5, 20, 100), (0.5, 1.0, 2.0, 5.0)
_TH = (np.arange(4001) + 0.5) / 4001 * (math.pi / 2)


def series(rng, k, n, snr):
    """tools/eta_vs_size_law.py's world, draw for draw: true m = 1, pushes N(0, snr^2/2 I_k) per component."""
    a = rng.normal(scale=snr, size=(n, k)) / math.sqrt(2.0)
    F, A = np.empty((n, k)), np.empty((n, k))
    for i in range(n):
        F[i] = a[i] + rng.normal(size=k)
        A[i] = a[i] + rng.normal(size=k)
    return F, A, 2.0 * float((a * a).sum())


def phasor(F, A):
    """|Z|, the best line's angle, and the series' energy."""
    Sff, Sfa, Saa = float((F * F).sum()), float((F * A).sum()), float((A * A).sum())
    return math.hypot(Saa - Sff, 2 * Sfa), 0.5 * math.atan2(2 * Sfa, Saa - Sff), Sff + Saa


def gram_share(absZ, energy):
    return absZ / (0.5 * (energy + absZ))


def true_line_share(F, A, phi0=math.pi / 4):
    t, r = math.sin(phi0) * F + math.cos(phi0) * A, math.cos(phi0) * F - math.sin(phi0) * A
    return float(((t * t).sum() - (r * r).sum()) / (t * t).sum())


def first_order_share(absZ, energy):
    lam = 0.5 * (energy + absZ)
    return max(0.0, absZ - 2 * lam / absZ) / lam


def _rho(eta, absZ, phi_hat):
    lw = eta * 0.25 * absZ * np.cos(2 * (_TH - phi_hat))
    w = np.exp(lw - lw.max())
    return float((w * np.cos(2 * (_TH - phi_hat))).sum() / w.sum())


def share_under_the_law(absZ, phi_hat, energy):
    """eta = <A>/<T> under the law with power eta, iterated down from 1; and the number of solutions on a grid."""
    f = lambda e: 2 * absZ * _rho(e, absZ, phi_hat) / (energy + absZ * _rho(e, absZ, phi_hat))
    e = 1.0
    for _ in range(200):
        new = min(1.0, max(0.0, f(e)))
        if abs(new - e) < 1e-10:
            break
        e = 0.5 * e + 0.5 * new
    grid = np.linspace(0, 1, 201)
    g = np.array([f(x) - x for x in grid])
    return e, int((np.sign(g[:-1]) != np.sign(g[1:])).sum())


def share_cell(idx, k, n, snr, n_series=200):
    rng = np.random.default_rng([2811, 70, idx])
    rows, several = [], 0
    for _ in range(n_series):
        F, A, delta = series(rng, k, n, snr)
        absZ, ph, energy = phasor(F, A)
        law, crossings = share_under_the_law(absZ, ph, energy)
        several += crossings > 1
        rows.append((delta / (delta + n * k), gram_share(absZ, energy), first_order_share(absZ, energy), law))
    r = np.array(rows)
    rms = lambda c: float(np.sqrt(((r[:, c] - r[:, 0]) ** 2).mean()))
    return dict(target=r[:, 0].mean(), gram=r[:, 1].mean(), first_order=r[:, 2].mean(), under_the_law=r[:, 3].mean(),
                rms_gram=rms(1), rms_first_order=rms(2), rms_under_the_law=rms(3), several_solutions=several)


def share_table():
    print("| k | N | snr | target | Gram | first-order | under the law | RMS to target: Gram | first-order | under the law | series with several solutions |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for idx, (k, n, snr) in enumerate((k, n, s) for k in KS for n in NS for s in SNRS):
        c = share_cell(idx, k, n, snr)
        print(f"| {k} | {n} | {snr} | {c['target']:.3f} | {c['gram']:.3f} | {c['first_order']:.3f} | {c['under_the_law']:.3f} | "
              f"{c['rms_gram']:.3f} | {c['rms_first_order']:.3f} | {c['rms_under_the_law']:.3f} | {c['several_solutions']} |")


# ---- T56's example ------------------------------------------------------------------------------------------------

SIGMA = np.array([[1.0, 0.9], [0.9, 1.0]])
M0, M_FALSE = 1.2, 4 / 15


def false_peak_laws(F, lm):
    """For a whitened Gram F and a grid of log masses: per law, the mass below m = 0.6 and the median and 95% interval.
    Laws: the Gaussian flat-push law, the score law pi exp(-T^2/2) (T56), and the phasor quasi-likelihood with D >= 0."""
    W = np.linalg.cholesky(np.linalg.inv(SIGMA)).T
    m = np.exp(lm)
    U = W @ np.vstack([m, np.ones_like(m)])
    U /= np.linalg.norm(U, axis=0)
    log_prior = -np.log1p(m * m) + lm
    vals, vecs = np.linalg.eigh(F)
    H, E = vals[1] - vals[0], vals[1] + vals[0]
    cos_d = U[0] * vecs[0, 1] + U[1] * vecs[1, 1]
    cos_t, sin2 = 2 * cos_d ** 2 - 1, 1 - cos_d ** 2
    D = -E + np.sqrt(E * E + 2 * H * E * cos_t + H * H)
    stats = {"gaussian": H * sin2,
             "score law": H * H * sin2 * cos_d ** 2 / (vals[1] * cos_d ** 2 + vals[0] * sin2),
             "signal energy >= 0": np.where(D >= 0, D - H * cos_t, H * H / (2 * E))}
    out = {}
    for name, stat in stats.items():
        lp = log_prior - 0.5 * stat
        cdf = np.cumsum(np.exp(lp - lp.max()))
        cdf /= cdf[-1]
        lo, med, hi = np.exp(np.interp([0.025, 0.5, 0.975], cdf, lm))
        out[name] = dict(false_side=float(cdf[np.searchsorted(m, 0.6)]), median=float(med), lo=float(lo), hi=float(hi))
    return out


def false_peak(n, n_series=20, seed=56, grid=200001):
    """T56's world: one body, a = 1, m0 = 6/5, covariance SIGMA, half-Cauchy c = 1. Returns per-law lists over series."""
    rng = np.random.default_rng([seed, n])
    W = np.linalg.cholesky(np.linalg.inv(SIGMA)).T
    C = np.linalg.cholesky(SIGMA)
    lm = np.linspace(math.log(1e-4), math.log(1e4), grid)
    acc = {}
    for _ in range(n_series):
        Y = (np.array([M0, 1.0])[None, :] + rng.normal(size=(n, 2)) @ C.T) @ W.T
        for name, r in false_peak_laws(Y.T @ Y, lm).items():
            acc.setdefault(name, []).append(r)
    return acc


def main():
    share_table()
    print("\nT56's example: mass of the law below m = 0.6 (false peak at 4/15, truth at 1.2); median [95% interval]; 20 series")
    for n in (20, 200, 2000, 20000):
        print(f"N = {n}")
        for name, rows in false_peak(n).items():
            fs = np.mean([r["false_side"] for r in rows])
            med, lo, hi = (np.median([r[key] for r in rows]) for key in ("median", "lo", "hi"))
            cover = np.mean([r["lo"] <= M0 <= r["hi"] for r in rows])
            print(f"   {name:20s} false-side mass {fs:.4f} | {med:.3f} [{lo:.3f}, {hi:.3f}] | truth inside 95%: {cover:.2f}")


if __name__ == "__main__":
    main()

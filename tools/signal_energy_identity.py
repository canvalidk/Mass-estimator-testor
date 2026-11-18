"""What the declared noise says about a reading's signal energy, and an identity it gives at the true mass
(Claude, 18 Nov 2026; problems/O14.md, 18 Nov). One instrument, noise units; nothing in met/ changes, no rule adopted.

A reading with k components has energy C = |f|^2 + |al|^2 (equation.md's E_i). Under the declared Gaussian noise,
whatever the push rho,   E[C] = rho^2 + 2k,   so C - 2k is an unbiased estimate of the reading's signal energy at every
signal size; among functions of C it is the only one (the noncentral chi-square family is complete in rho^2).

At the true mass the leftover K = |r|^2 is chi-square_k and independent of the along-line part t (T36 (4)), and
C = |t|^2 + |r|^2. Hence, for every configuration of pushes,

    E[ (C - 2k) K | t ] = k (|t|^2 - k + 2) = E[ k (C - 2k + 2) | t ],

so   O = sum_i (C_i - 2k_i) K_i   and   E0 = sum_i k_i (C_i - 2k_i + 2)   have the same mean at the true mass.
O is the squared-score sum of T31 with each reading's information read from its own energy, not from where the
reading sits along a fitted line.

This script checks the identity by simulation at the true line, and reports the same sums at the FITTED line, where
the fit removes part of the leftover. The correction used there, E1 = sum (C - 2k + 2)(k - h), h the reading's share of
along-line energy at the fitted line, is a guess and is not derived.

**Later the same day: the covariance form, and the fit** (`fit_cell`, `fit_table`). The reading's energy contains its
own leftover, C = |t|^2 + K with Var(K) = 2k, so

    E[ (C - 2) K | t ] = k (|t|^2 + k) = k E[ C | t ]        at the true mass, for every push and every k.

O(m) = sum_i (C_i - 2) K_i(m) and E0 = sum_i k_i C_i. For k = 1 this is the identity above. For k > 1 it differs from it:
an exactly zero reading adds nothing to either side here, while above it adds k (2 - 2k) to the expected side.
Read at the line where the weighted leftover is least, O* = min_m O(m) <= O(true mass) for every data set, so
E[O*] <= E[E0] for every configuration of pushes: an inequality, exact. To first order in the fit,
E[O*] ~ E0 - beta with beta = sum a_i^2 C_i / sum a_j C_j, a = C - 2. Read at the Gaussian best line instead, the mean
can exceed E0 (noise-only readings pull that line).

    python tools/signal_energy_identity.py          # both tables, about two minutes
"""
import math

import numpy as np


def series(rng, k, signals, reps, truth=1.0):
    """reps series; signals: push size per reading. Returns f, al of shape (reps, n, k) and the true angle."""
    signals = np.asarray(signals, float)
    phi0 = math.atan(truth)
    d = rng.normal(size=(reps, len(signals), k))
    d /= np.linalg.norm(d, axis=2, keepdims=True)
    f = signals[None, :, None] * math.sin(phi0) * d + rng.normal(size=d.shape)
    al = signals[None, :, None] * math.cos(phi0) * d + rng.normal(size=d.shape)
    return f, al, phi0


def parts(f, al, phi):
    """Leftover K, along-line energy T and energy C per reading, at angle phi (a number or one per series)."""
    s, c = np.sin(phi), np.cos(phi)
    if np.ndim(phi):
        s, c = s[:, None, None], c[:, None, None]
    t, r = s * f + c * al, c * f - s * al
    return (r * r).sum(2), (t * t).sum(2), (f * f).sum(2) + (al * al).sum(2)


def best_line(f, al):
    """The whole-line best line of the summed Gaussian leftover, per series."""
    return 0.5 * np.arctan2(2 * (f * al).sum((1, 2)), (al * al).sum((1, 2)) - (f * f).sum((1, 2)))


def cell(rng, k, signals, reps):
    """Mean O over mean E at the true line (E0) and at the fitted line (E0 and E1); how often O exceeds E1 there."""
    f, al, phi0 = series(rng, k, signals, reps)
    K0, _, C = parts(f, al, phi0)
    a = C - 2 * k
    O0, E0 = (a * K0).sum(1), (k * (C - 2 * k + 2)).sum(1)
    K1, T1, _ = parts(f, al, best_line(f, al))
    h = T1 / T1.sum(1, keepdims=True)
    O1, E1 = (a * K1).sum(1), ((C - 2 * k + 2) * (k - h)).sum(1)
    over = O1 > np.maximum(E1, 0.0)
    w = np.where(over, np.maximum(E1, 0.0) / np.where(O1 > 0, O1, 1.0), 1.0)
    return dict(true=O0.mean() / E0.mean(), fitted_uncorrected=O1.mean() / E0.mean(), fitted=O1.mean() / E1.mean(),
                discounted=float(over.mean()), median_weight=float(np.median(w[over])) if over.any() else 1.0)


def fit_cell(rng, k, signals, reps, positive=True):
    """Covariance form: mean O over mean E0 at the true line, at the Gaussian best line and at the weighted leftover's
    own least line (on the closed positive quadrant if positive); the last also over E0 - beta; and how often
    O* > E0 - beta on these ordinary series, with the weight then."""
    f, al, phi0 = series(rng, k, signals, reps)
    ff, fa, aa = (f * f).sum(2), (f * al).sum(2), (al * al).sum(2)
    C = ff + aa
    a = C - 2

    def K_at(phi):
        s, c = np.sin(phi), np.cos(phi)
        if np.ndim(phi):
            s, c = s[:, None], c[:, None]
        return ff * c * c - 2 * fa * s * c + aa * s * s

    E0 = (k * C).sum(1)
    O_true = (a * K_at(phi0)).sum(1)
    O_gauss = (a * K_at(0.5 * np.arctan2(2 * fa.sum(1), aa.sum(1) - ff.sum(1)))).sum(1)
    Sff, Sfa, Saa = (a * ff).sum(1), (a * fa).sum(1), (a * aa).sum(1)
    least = 0.5 * (Sff + Saa - np.hypot(Sff - Saa, 2 * Sfa))
    if positive:
        phi = 0.5 * np.arctan2(2 * Sfa, Saa - Sff)
        least = np.where((phi >= 0) & (phi <= math.pi / 2), least, np.minimum(Sff, Saa))
    den = (a * C).sum(1)
    E1 = E0 - np.where(den > 0, (a * a * C).sum(1) / np.where(den > 0, den, 1.0), 0.0)
    over = least > np.maximum(E1, 0.0)
    w = np.where(over, np.maximum(E1, 0.0) / np.where(least > 0, least, 1.0), 1.0)
    return dict(true=O_true.mean() / E0.mean(), gauss=O_gauss.mean() / E0.mean(), own=least.mean() / E0.mean(),
                own_corrected=least.mean() / E1.mean(), never_above_true=bool((least <= O_true + 1e-9 * np.abs(O_true) + 1e-9).all()),
                discounted=float(over.mean()), median_weight=float(np.median(w[over])) if over.any() else 1.0)


def fit_table(seed=20261118):
    rng = np.random.default_rng(seed)
    print("| k | series | O/E0 true line | Gaussian best line | own least line | own, over E0 - beta | O* > E0 - beta | median weight then |")
    print("|---|---|---|---|---|---|---|---|")
    for k in (1, 3):
        for name, signals in cases().items():
            c = fit_cell(rng, k, signals, 20000 if len(signals) <= 60 else 4000)
            print(f"| {k} | {name} | {c['true']:.3f} | {c['gauss']:.3f} | {c['own']:.3f} | {c['own_corrected']:.3f} | {100*c['discounted']:.0f}% | {c['median_weight']:.2f} |")


def cases():
    out = {}
    for n in (2, 5, 20, 200):
        for sig in (0.0, 1.0, 5.0):
            out[f"{n} readings at signal {sig:g}"] = [sig] * n
    out["5 strong (10) and 50 weak (0.5)"] = [10] * 5 + [0.5] * 50
    out["one strong (30) among 60 nulls"] = [30] + [0.0] * 60
    out["two (5) among 200 nulls"] = [5] * 2 + [0.0] * 200
    out["40 readings, signals spread 0.1 to 8"] = list(np.linspace(0.1, 8, 40))
    return out


def main(seed=20261118):
    rng = np.random.default_rng(seed)
    print("| k | series | true line O/E0 | fitted line O/E0 | fitted line O/E1 | O > E1 | median E1/O then |")
    print("|---|---|---|---|---|---|---|")
    for k in (1, 3):
        for name, signals in cases().items():
            c = cell(rng, k, signals, 20000 if len(signals) <= 60 else 4000)
            print(f"| {k} | {name} | {c['true']:.3f} | {c['fitted_uncorrected']:.3f} | {c['fitted']:.3f} | {100*c['discounted']:.0f}% | {c['median_weight']:.2f} |")


if __name__ == "__main__":
    main()
    print()
    fit_table()

"""What no single mass explains, read with no mass: readings taken in pairs
(Claude, 18 Nov 2026; problems/O14.md, 18 Nov, night). Measurements and three identities; nothing in met/ changes,
no rule adopted.

One instrument, noise units. Reading i has k_i components, f_i and alpha_i, energy C_i = |f_i|^2 + |alpha_i|^2 and its
own 2x2 Gram G_i. Two readings lie on one line through the origin (one mass explains both) exactly when

    X_ij = |f_i|^2 |alpha_j|^2 + |alpha_i|^2 |f_j|^2 - 2 (f_i.alpha_i)(f_j.alpha_j)

is zero: for scalar readings X_ij = (f_i alpha_j - alpha_i f_j)^2, the squared area the two readings span. A vector
reading is on one line by itself exactly when det G_i = |f_i|^2 |alpha_i|^2 - (f_i.alpha_i)^2 is zero (its force and
acceleration are parallel). The series' Gram has det(sum G_i) = sum_i det G_i + sum_{i<j} X_ij (Cauchy-Binet): what no
single mass explains is this sum over pairs, and no candidate mass appears in it.

Under the declared Gaussian noise, for readings whose signals share one line and are otherwise arbitrary:

    between two readings    E[(C_i - 2)(C_j - 2) X_ij] = E[C_i C_j (k_j (C_i - k_i - 2) + k_i (C_j - k_j - 2))]
    inside one reading      E[(C_i - 4) det G_i]       = (k_i - 1) E[C_i (C_i - k_i - 2)]
    the series' signal      E[|Z|^2 - 4 E]             = Delta (Delta + 4)

No true mass, no fitted mass, no count of readings, no size of signal; an exactly zero reading adds nothing to either
side of any of them. Among weights C - c the constants 2 (between) and 4 (inside) are the ones for which the expected
side vanishes at a zero reading (L20). The first two are compared as "expected over observed, capped at 1" below: a
measurement of what that reading of them gives, not a rule.

    python tools/contradiction_in_pairs.py            # about a minute
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discount_lab as lab                                      # noqa: E402

FITTED = "one weight, energy less 2, own best line"            # the morning's body, read at a fitted line


def parts(ws):
    """Per reading: |f|^2, |alpha|^2, f.alpha, k."""
    a = np.array([f @ f for f, _ in ws])
    b = np.array([al @ al for _, al in ws])
    c = np.array([f @ al for f, al in ws])
    k = np.array([len(f) for f, _ in ws], float)
    return a, b, c, k


def pair_sums(ws):
    """Observed and expected sides, between readings and inside readings."""
    a, b, c, k = parts(ws)
    C = a + b
    n = len(ws)
    obs = exp = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            X = a[i] * b[j] + b[i] * a[j] - 2 * c[i] * c[j]
            obs += (C[i] - 2) * (C[j] - 2) * X
            exp += C[i] * C[j] * (k[j] * (C[i] - k[i] - 2) + k[i] * (C[j] - k[j] - 2))
    inside_obs = float(((C - 4) * (a * b - c * c)).sum())
    inside_exp = float(((k - 1) * C * (C - k - 2)).sum())
    return dict(between_observed=float(obs), between_expected=float(exp),
                inside_observed=inside_obs, inside_expected=inside_exp)


def expected_between(Ci, Cj, c=2.0):
    """Scalar readings, weights C - c: the function of the two energies whose mean equals E[(C_i - c)(C_j - c) X_ij]
    for every pair of signals. At C_i = 0 it is (2 - c) (C_j^2 - (2 + c) C_j + 2 c - 4): zero for every C_j only at c = 2."""
    A = lambda C: C * C - (2 + c) * C + 2 * c - 4
    B = lambda C: C + 2 - c
    return A(Ci) * B(Cj) + B(Ci) * A(Cj) + 2 * B(Ci) * B(Cj)


def expected_inside(C, k, c=4.0):
    """A k-component reading, weight C - c: the function of its energy whose mean equals E[(C - c) det G]. At C = 0
    it is -(k - 1) k (4 - c): zero only at c = 4."""
    return (k - 1) * (C * C + (2 - c - k) * (C - k) - k * (k + 2))


def pair_weight(ws):
    """Expected over observed, capped at 1, both sides summed over all pairs and all readings' insides. With a
    non-positive expected side (a quiet series) the weight is 1."""
    s = pair_sums(ws)
    obs = s["between_observed"] + s["inside_observed"]
    exp = s["between_expected"] + s["inside_expected"]
    return 1.0 if (exp <= 0 or obs <= exp) else exp / obs


def draw(rng, k, pushes, reps, phi0=0.7):
    """reps series; reading i has k components and a signal of length pushes[i] in a random direction, on the line phi0."""
    n = len(pushes)
    d = rng.normal(size=(reps, n, k))
    d /= np.linalg.norm(d, axis=2, keepdims=True)
    mu = d * np.asarray(pushes, float)[None, :, None]
    return (math.sin(phi0) * mu + rng.normal(size=(reps, n, k)),
            math.cos(phi0) * mu + rng.normal(size=(reps, n, k)))


def identity_cell(k, pushes, reps, seed=22):
    """Mean observed over mean expected for the three identities (the inside one only for k > 1), with the standard
    error of the first."""
    rng = np.random.default_rng([seed, k, len(pushes), int(round(100 * sum(pushes)))])
    f, al = draw(rng, k, pushes, reps)
    a, b, c = (f * f).sum(2), (al * al).sum(2), (f * al).sum(2)
    C = a + b
    n = C.shape[1]
    D, V = np.zeros(reps), np.zeros(reps)
    for i in range(n):
        for j in range(i + 1, n):
            X = a[:, i] * b[:, j] + b[:, i] * a[:, j] - 2 * c[:, i] * c[:, j]
            D += (C[:, i] - 2) * (C[:, j] - 2) * X
            V += C[:, i] * C[:, j] * k * (C[:, i] + C[:, j] - 2 * k - 4)
    Dw, Vw = ((C - 4) * (a * b - c * c)).sum(1), ((k - 1) * C * (C - k - 2)).sum(1)
    energy = C.sum(1)
    Z2 = (b - a).sum(1) ** 2 + 4 * c.sum(1) ** 2
    delta = float(np.sum(np.square(pushes)))
    return dict(between=float(D.mean() / V.mean()), between_se=float(D.std() / math.sqrt(reps) / V.mean()),
                inside=float(Dw.mean() / Vw.mean()) if k > 1 else float("nan"),
                signal=float(np.mean(Z2 - 4 * energy)), signal_exact=delta * (delta + 4),
                signal_se=float(np.std(Z2 - 4 * energy) / math.sqrt(reps)))


CASES = (("2 readings, no signal", [0.0, 0.0]), ("2 readings at signal 1", [1.0, 1.0]),
         ("2 readings at signal 5", [5.0, 5.0]), ("5 readings at signal 0.5", [0.5] * 5),
         ("signals 0, 2 and 7", [0.0, 2.0, 7.0]), ("20 readings, no signal", [0.0] * 20),
         ("20 readings at signal 1", [1.0] * 20), ("two at 5 among 40 at none", [5.0, 5.0] + [0.0] * 40))


def ordinary(k, pushes, reps=400, seed=23):
    """Ordinary series (the declaration true): how often each weight is below 1, and its median when it is."""
    rng = np.random.default_rng([seed, k, len(pushes), int(round(100 * sum(pushes)))])
    f, al = draw(rng, k, pushes, reps)
    out = {}
    for name, fn in (("pairs", pair_weight), ("fitted", lambda ws: float(lab.weights(ws, FITTED)[0][0]))):
        w = np.array([fn([(f[r, i], al[r, i]) for i in range(len(pushes))]) for r in range(reps)])
        low = w < 1
        out[name] = (float(low.mean()), float(np.median(w[low])) if low.any() else 1.0)
    return out


def scene_table():
    rows = []
    for name, (readings, sigmas, c, ratio) in lab.scenes().items():
        for sigma in sigmas:
            ws = lab.whiten(readings, sigma, sigma * ratio)
            rows.append((name, sigma, pair_weight(ws), float(lab.weights(ws, FITTED)[0][0])))
    return rows


def main():
    print("1. The identities: mean observed over mean expected (1 = exact), no mass anywhere")
    print("| series | k | between readings (se) | inside a reading | E[|Z|^2 - 4E] | Delta (Delta + 4) |")
    print("|---|---|---|---|---|---|")
    for k in (1, 3):
        for name, pushes in CASES:
            r = identity_cell(k, pushes, 400000 if len(pushes) <= 5 else 40000)
            inside = "" if k == 1 else f"{r['inside']:.3f}"
            print(f"| {name} | {k} | {r['between']:.3f} ({r['between_se']:.3f}) | {inside} | "
                  f"{r['signal']:.4g} ({r['signal_se']:.2g}) | {r['signal_exact']:.4g} |")

    print("\n2. Scenes: the weight from pairs (no fit) | the morning's weight at its own least line")
    print("| scene | sigma | pairs | fitted |")
    print("|---|---|---|---|")
    for name, sigma, wp, wf in scene_table():
        print(f"| {name} | {sigma} | {wp:.4g} | {wf:.4g} |")

    print("\n3. Ordinary series: share with weight below 1, and the median weight then (400 series)")
    print("| series | k | pairs | fitted |")
    print("|---|---|---|---|")
    for k in (1, 3):
        for name, pushes in (("5 readings at signal 5", [5.0] * 5), ("20 readings at signal 5", [5.0] * 20),
                             ("20 readings at signal 1", [1.0] * 20), ("two at 5 among 40 at none", [5.0, 5.0] + [0.0] * 40)):
            r = ordinary(k, pushes)
            print(f"| {name} | {k} | {r['pairs'][0]:.2f}, {r['pairs'][1]:.2f} | {r['fitted'][0]:.2f}, {r['fitted'][1]:.2f} |")


if __name__ == "__main__":
    main()

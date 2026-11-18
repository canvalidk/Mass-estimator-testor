"""The discount on the scalar slice: readings x_i +- s_i of one mass (Claude, 18 Nov 2026; problems/O14.md, 18 Nov).
Measurements; nothing in met/ changes and no rule is adopted.

K_i(m) = (x_i - m)^2 / s_i^2,  I_i = 1 / s_i^2,  h_i = I_i / sum I (leverage). The expected leftover of reading i at the
fitted mass is 1 - h_i, not 1 (simulated below; the textbook hat-matrix identity). With that count:

  - two readings: one weight and one weight per reading are the same number, (s_1^2 + s_2^2) / gap^2, at any sigma ratio;
  - lumping: two bodies through common exact accelerations, forces measured separately, the glued body's force the sum.
    With one weight the glued estimate is the sum of the parts' estimates exactly; with one weight per reading it is not.

    python tools/discount_scalar_slice.py
"""
import math

import numpy as np

Z95 = 1.959964


def gauss(x, s):
    I = 1 / s ** 2
    return float((I * x).sum() / I.sum()), 1 / math.sqrt(I.sum())


def one_weight(x, s, count="leverage"):
    """Estimate, sd and weight: expected over observed leftover at the Gaussian estimate, capped at 1."""
    I = 1 / s ** 2
    m, sd = gauss(x, s)
    K = (x - m) ** 2 * I
    h = I / I.sum()
    if count == "leverage":
        num, den = (I * (1 - h)).sum(), (I * K).sum()
    elif count == "N":                      # the 17 Nov scratch (tools/leftover_at_series_best_mass.py)
        num, den = float(len(x)), K.sum()
    else:
        raise ValueError(count)
    w = 1.0 if den <= num else num / den
    return m, sd / math.sqrt(w), w


def per_reading(x, s, steps=1, start=None):
    """One weight per reading, (1 - h_i) / K_i capped at 1, read at the Gaussian estimate (steps = 1) or iterated."""
    I = 1 / s ** 2
    kappa = 1 - I / I.sum()
    m = gauss(x, s)[0] if start is None else start
    w = np.ones(len(x))
    for _ in range(steps):
        K = (x - m) ** 2 * I
        w = np.where(K > kappa, kappa / np.maximum(K, 1e-300), 1.0)
        m = float((w * I * x).sum() / (w * I).sum())
    return m, 1 / math.sqrt((w * I).sum()), w


def glue(F1, F2, s1, s2, rule):
    """Sum of the parts' estimates and the glued body's estimate under a rule (x, s) -> (estimate, sd, ...)."""
    F1, F2 = np.asarray(F1, float), np.asarray(F2, float)
    p1, p2 = rule(F1, np.full(len(F1), s1)), rule(F2, np.full(len(F2), s2))
    pc = rule(F1 + F2, np.full(len(F1), math.hypot(s1, s2)))
    return dict(parts=p1[0] + p2[0], parts_sd=math.hypot(p1[1], p2[1]), glued=pc[0], glued_sd=pc[1])


def expected_leftover_at_the_fit(seed=3, reps=400000):
    rng = np.random.default_rng(seed)
    s = np.array([0.05, 0.1, 0.1, 0.3, 1.0]); I = 1 / s ** 2
    X = 5 + s * rng.normal(size=(reps, 5))
    mhat = (X * I).sum(1) / I.sum()
    return (((X - mhat[:, None]) / s) ** 2).mean(0), 1 - I / I.sum()


def lumping_on_ordinary_data(rule, seed=1811, reps=20000, n=10):
    """Glued estimate minus the sum of the parts' estimates, in units of the glued sd."""
    rng = np.random.default_rng(seed)
    off = []
    for _ in range(reps):
        F1, F2 = 5 + 0.1 * rng.normal(size=n), 3 + 0.1 * rng.normal(size=n)
        g = glue(F1, F2, 0.1, 0.1, rule)
        off.append((g["glued"] - g["parts"]) / (0.1 * math.sqrt(2) / math.sqrt(n)))
    return np.array(off)


def main():
    sim, formula = expected_leftover_at_the_fit()
    print("expected leftover at the fitted mass, simulated", np.round(sim, 4), " 1 - h_i", np.round(formula, 4))
    for s1, s2 in ((0.1, 0.1), (0.05, 0.1), (0.01, 0.3)):
        x, s = np.array([5.0, 6.0]), np.array([s1, s2])
        print(f"two readings at 5 and 6, sd {s1}, {s2}: one weight {one_weight(x, s)[2]:.6f}, per reading "
              f"{np.round(per_reading(x, s)[2], 6)}, (s1^2 + s2^2)/gap^2 {s1**2 + s2**2:.6f}")
    x, s = np.array([5.0, 5.05, 6.0]), np.full(3, 0.1)
    for start in (None, 6.0):
        m, sd, w = per_reading(x, s, steps=500, start=start)
        print(f"5, 5.05, 6 per reading, iterated from {'the Gaussian estimate' if start is None else start}: {m:.4f} +- {Z95*sd:.4f}, weights {np.round(w, 4)}")
    rules = {"gauss": gauss, "one weight": one_weight, "per reading": per_reading}
    for label, F1, F2 in (("stray in body 1", [5.0, 5.05, 6.0], [3.1, 2.9, 3.0]),
                          ("a stray in each, different trials", [5.0, 5.05, 6.0], [3.0, 3.9, 3.05]),
                          ("a stray in each, same trial", [5.0, 5.05, 6.0], [3.0, 3.05, 3.9])):
        for name, rule in rules.items():
            g = glue(F1, F2, 0.1, 0.1, rule)
            print(f"  {label:34s} {name:12s} parts {g['parts']:.3f} +- {Z95*g['parts_sd']:.3f} | glued {g['glued']:.3f} +- {Z95*g['glued_sd']:.3f} | off by {g['glued']-g['parts']:+.3f}")
    for name in ("one weight", "per reading"):
        off = lumping_on_ordinary_data(rules[name])
        print(f"  ordinary data, 10 trials, {name}: glued minus parts in glued sd: rms {np.sqrt((off**2).mean()):.3f}, largest {np.abs(off).max():.2f}")


if __name__ == "__main__":
    main()

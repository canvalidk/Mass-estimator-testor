"""Where a reading's leftover is read (Claude, 17 Nov 2026; problems/O14.md, "the sign problem and where the
leftover is read"). Scratch on one-number readings; nothing in met/ changes.

Three ways to apply "weight = expected leftover over observed leftover, capped at 1":
  every candidate   the weight re-evaluated at each candidate mass (Gaussian core, log tail)
  pooled            one weight for the series, from the total leftover at the series' best mass
  per reading       one weight per reading, from its own leftover at the series' best mass
and, for the sign problem, ours30's way: the reading's own misfit at its own best positive mass.

    python tools/leftover_at_series_best_mass.py
"""
import math

import numpy as np

KINDS = ("gauss", "every candidate", "pooled", "per reading")


def law(forces, sigma_f=0.1, sigma_a=0.002, c=5.0, kind="gauss", kappa=1.0):
    """Median and central 95% interval of pi times the weighted evidence, for readings (F_i, a = 1)."""
    lm = np.linspace(math.log(c) - 12, math.log(c) + 12, 400001)
    m = np.exp(lm)
    K = np.array([(f - m) ** 2 / (sigma_f ** 2 + (m * sigma_a) ** 2) for f in forces])
    log_pi = -np.log1p((m / c) ** 2) + lm
    at_best = K[:, np.argmin(K.sum(0))]                 # each reading's leftover at the series' best mass
    if kind == "every candidate":
        rho = np.where(K <= kappa, K, kappa * (1 + np.log(np.maximum(K, 1e-300) / kappa)))
        log_evidence = -0.5 * rho.sum(0)
    else:
        if kind == "gauss":
            w = np.ones(len(forces))
        elif kind == "pooled":
            w = np.full(len(forces), min(1.0, kappa * len(forces) / at_best.sum()))
        elif kind == "per reading":
            w = np.minimum(1.0, kappa / np.maximum(at_best, 1e-300))
        else:
            raise ValueError(kind)
        log_evidence = -0.5 * (w[:, None] * K).sum(0)
    lp = log_pi + log_evidence
    cdf = np.cumsum(np.exp(lp - lp.max()))
    lo, med, hi = np.exp(np.interp([0.025, 0.5, 0.975], cdf / cdf[-1], lm))
    return med, lo, hi


def sign_check(signal, n=2_000_000, seed=11, truth=2.0):
    """Many scalar readings at one fixed weak signal (c = 1, unit noise): the weighted best line's mass with no
    weights, with weights read at the series' best mass, and with ours30's (own best positive mass, kappa 1/2)."""
    rng = np.random.default_rng(seed)
    th0 = math.atan(truth)
    Y = signal * np.array([math.sin(th0), math.cos(th0)]) + rng.normal(size=(n, 2))

    def best(w):
        vals, vecs = np.linalg.eigh((Y * w[:, None]).T @ Y)
        v = vecs[:, 1] if vecs[1, 1] > 0 else -vecs[:, 1]
        return math.atan2(v[0], v[1])

    th = best(np.ones(n))
    at_series = (Y @ np.array([math.cos(th), -math.sin(th)])) ** 2
    own = np.where(Y[:, 0] * Y[:, 1] < 0, np.minimum(Y[:, 0] ** 2, Y[:, 1] ** 2), 0.0)
    w_series = np.minimum(1.0, 1.0 / np.maximum(at_series, 1e-300))
    w_own = np.where(own > 0, np.minimum(1.0, 0.5 / np.maximum(own, 1e-300)), 1.0)
    return math.tan(th), math.tan(best(w_series)), math.tan(best(w_own))


def cases(seed=3):
    ordinary = list(5 + 0.1 * np.random.default_rng(seed).normal(size=99))
    return {"5.5, 5.6": [5.5, 5.6], "5, 6": [5.0, 6.0], "5, 5.05, 6": [5.0, 5.05, 6.0],
            "99 ordinary near 5": ordinary, "the same 99 and one at 9": ordinary + [9.0]}


def main():
    print("| readings (each ± 0.1) | " + " | ".join(KINDS) + " |")
    print("|---|" + "---|" * len(KINDS))
    for name, forces in cases().items():
        cells = []
        for kind in KINDS:
            med, lo, hi = law(forces, kind=kind)
            cells.append(f"{med:.3f} ± {0.5 * (hi - lo):.3f}")
        print(f"| {name} | " + " | ".join(cells) + " |")
    print()
    print("| signal per reading | no weights | read at the series' best mass | read at the reading's own best positive mass |")
    print("|---|---|---|---|")
    for signal in (1.0, 2.0):
        a, b, c = sign_check(signal)
        print(f"| {signal:g} | {a:.3f} | {b:.3f} | {c:.3f} |")


if __name__ == "__main__":
    main()

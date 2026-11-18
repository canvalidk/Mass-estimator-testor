"""Readings at right angles in the force-acceleration plane (Claude, 18 Nov 2026; problems/O14.md, 18 Nov).

Scene C: one reading with force 2 and no acceleration, one with acceleration 1 and no force, each +- sigma in both
channels (sigma_F = sigma_a = sigma, c = 1). No mass explains both; the leftover at the best line is 1/sigma^2. The two
readings' own lines are at right angles in noise units, so each has zero score at the series' best line: eta's squared
scores vanish (J = 0) and eta = 1. What the equations in met/ return as sigma shrinks, and the same four numbers
declared as one 2D reading, where ours29 and ours30's g acts.

    python tools/right_angle_scenes.py
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from met.full import Reading, ours29, ours30  # noqa: E402
from met.ours33 import ours33  # noqa: E402
import discount_lab as lab  # noqa: E402

SIGMAS = (1.0, 0.5, 0.2, 0.1, 0.05)
EQUATIONS = {"ours29": lambda rs: ours29(rs, 1.0), "ours30": lambda rs: ours30(rs, 1.0),
             "ours33, nu = inf": lambda rs: ours33(rs, 1.0, math.inf), "ours33, nu = 10": lambda rs: ours33(rs, 1.0, 10)}


def scene_c(sigma, split=True):
    if split:
        cov = np.diag([sigma ** 2] * 2)
        return [Reading(np.array([2.0]), np.array([0.0]), cov), Reading(np.array([0.0]), np.array([1.0]), cov)]
    return [Reading(np.array([2.0, 0.0]), np.array([0.0, 1.0]), np.diag([sigma ** 2] * 4))]


def run(name, sigma, split=True):
    r = EQUATIONS[name](scene_c(sigma, split))
    interval = r.get("interval_95", r.get("interval"))
    return dict(median=float(r.get("median", r.get("point"))), lo=float(interval[0]), hi=float(interval[1]),
                eta=float(r["eta"]), g=None if r.get("g") is None else [float(x) for x in np.atleast_1d(r["g"])])


def near_right_angle(eps_deg, sigma, body, grid=200001):
    """Two readings on m = 1.73 and one eps degrees off the right angle to them; the law's 95% interval."""
    a, b = math.radians(60), math.radians(-30 + eps_deg)
    readings = [([2 * math.sin(a)], [2 * math.cos(a)])] * 2 + [([math.sin(b)], [math.cos(b)])]
    r = lab.law(readings, sigma, sigma, 1.0, body, grid=grid)
    return r["lo"], r["hi"]


def main():
    for split in (True, False):
        print("\ntwo scalar readings (2, 0) and (0, 1)" if split else "\nthe same numbers as one 2D reading F = (2, 0), a = (0, 1)")
        print("| sigma | " + " | ".join(EQUATIONS) + " |")
        print("|---|" + "---|" * len(EQUATIONS))
        for sigma in SIGMAS:
            cells = []
            for name in EQUATIONS:
                r = run(name, sigma, split)
                g = "" if r["g"] is None else ", g " + "/".join(f"{x:.3g}" for x in r["g"])
                cells.append(f"{r['median']:.3g} [{r['lo']:.3g}, {r['hi']:.3g}], eta {r['eta']:.3g}{g}")
            print(f"| {sigma} | " + " | ".join(cells) + " |")
    print("\ntwo readings on m = 1.73, one eps degrees off the right angle: 95% interval at sigma = 0.002")
    for eps in (0.0, 1.0, 5.0, 20.0):
        row = [f"eps {eps:g}"]
        for body in ("one weight, by count", "one weight per reading"):
            lo, hi = near_right_angle(eps, 0.002, body)
            row.append(f"{body}: [{lo:.3f}, {hi:.3f}]")
        print("  " + " | ".join(row))


if __name__ == "__main__":
    main()

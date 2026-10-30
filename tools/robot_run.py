"""Run the definition (ours30, met/full.py) on the robot experiment and print one table.

    python tools/robot_run.py                 # all runs, the default widths
    python tools/robot_run.py --widths 0.2 --runs 1-baseline
    python tools/robot_run.py --equation ours29   # the table as it stood before 29 Oct
    python tools/robot_run.py --white autocorrelated --rotation static
    python tools/robot_run.py --matrix            # B01's complete seeded preset

The original default mode uses no random numbers. The new declaration modes use
the fixed calibration-side bootstrap seed, 20261030. Each row is one run cut into windows of one width;
the readings are those of ``met.robot.windows`` and c is the series' noise ratio
(``met.robot.series_ratio``, Claude's declaration for several noise pairs, O11).
Columns: the law's median and 95% interval, the raw law's interval (eta = 1), the
readout, eta, the smallest misfit weight, and whether the interval covers the lab's
static calibration mass. Two last rows read the static calibration orientations
themselves, biases fitted on the other half. ``--no-floor`` drops the systematic floor from the budget.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from met import robot  # noqa: E402
from met.full import ours29, ours30  # noqa: E402

EQUATIONS = {"ours30": ours30, "ours29": ours29}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--runs", nargs="*", default=None)
    p.add_argument("--widths", nargs="*", type=float, default=None)
    p.add_argument("--no-floor", action="store_true")
    p.add_argument("--equation", choices=sorted(EQUATIONS), default="ours30")
    p.add_argument("--white", choices=("independent", "autocorrelated"), default=None)
    p.add_argument("--rotation", choices=("nominal", "static"), default=None)
    p.add_argument("--matrix", action="store_true", help="run all B01 preset rows and save their provenance")
    p.add_argument("--out", type=Path, help="matrix output root (default results/)")
    a = p.parse_args(argv)
    if a.matrix:
        if (a.runs is not None or a.widths is not None or a.no_floor or a.equation != "ours30"
                or a.white is not None or a.rotation is not None):
            p.error("--matrix uses all selectors in presets/robot_budget.yaml; edit a copied preset for a subset")
        from met.robot_budget import run, summary
        root = Path(__file__).resolve().parent.parent
        out, record = run(root / "presets" / "robot_budget.yaml", a.out or root / "results", progress=None)
        print(summary(record), end="")
        print(f"saved {out}")
        return 0
    if a.out is not None:
        p.error("--out is for --matrix; individual tables are printed to stdout")
    a.runs = list(robot.RUNS) if a.runs is None else a.runs
    a.widths = [0.05, 0.2, 1.0] if a.widths is None else a.widths
    a.white = "independent" if a.white is None else a.white
    a.rotation = "nominal" if a.rotation is None else a.rotation
    equation = EQUATIONS[a.equation]
    cal = robot.calibrate()
    declaration = (robot.calibration_declaration(cal)
                   if a.white != "independent" or a.rotation != "nominal" else None)
    print(f"{a.equation}; calibration mass {cal.mass:.4f} kg (published {robot.REFERENCE_MASS_KG})")
    if declaration is not None:
        print(f"white={a.white}, rotation={a.rotation}; calibration-side seed={declaration.metadata['seed']}")
    print("run                   width   n     c  median  interval          raw interval      "
          "readout   eta   min g  covers")
    for run in a.runs:
        for w in a.widths:
            rd, _ = robot.windows(run, w, cal, floor=not a.no_floor,
                                  declaration=declaration, white=a.white, rotation=a.rotation)
            c = robot.series_ratio(rd)
            r = equation(rd, c)
            lo, hi = r["interval"]
            rlo, rhi = r["raw_interval"]
            cov = "yes" if lo <= cal.mass <= hi else "no"
            print(f"{run:21s} {w:5.2f} {len(rd):4d} {c:5.3f} {r['median']:.4f}  "
                  f"[{lo:.4f}, {hi:.4f}]  [{rlo:.4f}, {rhi:.4f}]  {r['readout']:.4f}  "
                  f"{r['eta']:.3f}  {np.min(r['g']):.1e}  {cov}")
    if declaration is not None:
        print("static controls retain their original nominal single-sample budget")
    for name, fit_rows, test_rows in (("static, fit 1-12", slice(0, 12), slice(12, None)),
                                      ("static, fit 13-23", slice(12, None), slice(0, 12))):
        rd, _ = robot.static_heldout(fit_rows, test_rows, floor=not a.no_floor)
        c = robot.series_ratio(rd)
        r = equation(rd, c)
        lo, hi = r["interval"]
        rlo, rhi = r["raw_interval"]
        cov = "yes" if lo <= cal.mass <= hi else "no"
        print(f"{name:21s}     - {len(rd):4d} {c:5.3f} {r['median']:.4f}  "
              f"[{lo:.4f}, {hi:.4f}]  [{rlo:.4f}, {rhi:.4f}]  {r['readout']:.4f}  "
              f"{r['eta']:.3f}  {np.min(r['g']):.1e}  {cov}")
    return 0


if __name__ == "__main__":
    main()

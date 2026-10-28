"""Command line.

  python -m met list                       presets and the question each answers
  python -m met validate PRESET            check a preset without running it
  python -m met run PRESET [--replicates N] [--out DIR] [--workers N|max]
  python -m met inspect CASE [--case NAME] [--check ID] [--out DIR]

--workers: cells run side by side in this many processes; the default is half the cores
(the machine stays usable), "max" uses every core, 1 runs one cell at a time. The results
are the same for any value.
"""

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _workers(text):
    if text == "max":
        return "max"
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a positive integer or 'max', got {text!r}")
    if value < 1:
        raise argparse.ArgumentTypeError(f"expected a positive integer or 'max', got {text!r}")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m met", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list presets")
    v = sub.add_parser("validate", help="validate a preset")
    v.add_argument("preset", type=Path)
    r = sub.add_parser("run", help="run a preset")
    r.add_argument("preset", type=Path)
    r.add_argument("--replicates", type=int, help="override the preset's replicates (recorded in the output)")
    r.add_argument("--out", type=Path, default=ROOT / "results")
    r.add_argument("--workers", type=_workers, default=None,
                   help="processes running cells side by side: a number, or 'max' (default: half the cores)")
    i = sub.add_parser("inspect", help="inspect a saved case or case collection, serially")
    i.add_argument("preset", type=Path)
    i.add_argument("--case", dest="cases", action="append", help="case name (repeat to select several)")
    i.add_argument("--check", dest="checks", action="append", help="check id (repeat to select several)")
    i.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args(argv)
    serial_inspection = args.command == "inspect"
    if args.command == "run":
        # Inspect the lightweight YAML envelope before any numerical imports so
        # the legacy run alias honours the same single-thread resource policy.
        import yaml
        with args.preset.open(encoding="utf-8") as handle:
            envelope = yaml.safe_load(handle)
        serial_inspection = isinstance(envelope, dict) and envelope.get("kind") in ("inspection", "noise_path", "prior_path", "pooled_scale", "shared_tilt", "anisotropic_world", "readout_candidates")
    if serial_inspection:
        # Do this before importing NumPy on the standalone CLI path. Library users
        # control their own numerical runtime; inspection never spawns processes.
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[key] = "1"
    from . import preset as preset_mod
    if args.command == "list":
        for name, question in preset_mod.list_presets(ROOT / "presets"):
            print(f"{name}\n    {question}")
        return 0
    try:
        if args.command == "inspect":
            from .inspection.run import run as inspect_run
            inspect_run(args.preset, args.out, case_names=args.cases, check_ids=args.checks)
            return 0
        if args.command == "validate":
            spec = preset_mod.load(args.preset)
            if spec.get("kind") == "readout_candidates":
                print(f"ok: {len(spec['cells'])} readout-candidate diagnostic cell(s); serial, no random draws")
            elif spec.get("kind") == "anisotropic_world":
                print(f"ok: {len(spec['cells'])} anisotropic-world research cell(s); serial")
            elif spec.get("kind") == "shared_tilt":
                print(f"ok: {len(spec['cells'])} shared-tilt research cell(s); serial")
            elif spec.get("kind") == "pooled_scale":
                print(f"ok: {len(spec['cells'])} pooled-scale research cell(s); serial")
            elif spec.get("kind") in ("inspection", "noise_path", "prior_path"):
                print(f"ok: {len(spec['cases'])} saved inspection case(s); no random draws")
            elif spec.get("kind") == "misalignment":
                print(f"ok: {len(spec['cells'])} fixed-reading diagnostic cell(s); no random draws")
            else:
                print(f"ok: {len(spec['cells'])} cell(s), {len(spec['estimators'])} estimator(s)")
            return 0
        from .runner import run
        run(args.preset, args.out, args.replicates, workers=args.workers)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("interrupted: no results written", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())

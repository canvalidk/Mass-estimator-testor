"""Command line.

  python -m met list                       presets and the question each answers
  python -m met validate PRESET            check a preset without running it
  python -m met run PRESET [--replicates N] [--out DIR]
"""

import argparse
import sys
from pathlib import Path

from . import preset as preset_mod
from .runner import run

ROOT = Path(__file__).resolve().parent.parent


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
    args = parser.parse_args(argv)
    if args.command == "list":
        for name, question in preset_mod.list_presets(ROOT / "presets"):
            print(f"{name}\n    {question}")
        return 0
    try:
        if args.command == "validate":
            spec = preset_mod.load(args.preset)
            print(f"ok: {len(spec['cells'])} cell(s), {len(spec['estimators'])} estimator(s)")
            return 0
        run(args.preset, args.out, args.replicates)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

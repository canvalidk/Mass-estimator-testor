"""Compatibility CLI for the historical ours29 composition-coverage preset.

    python -m met run presets/composition_coverage.yaml --workers 4
    python tools/composition_coverage.py --replicates 400 --workers 4

The legacy CLI keeps its default 200 replicates, four workers, seed 2810 and JSON --out path.
The shipped preset uses the recorded 400 replicates. The numerical results are
independent of worker count. Native preset artifacts
are saved under results/ as well as the legacy JSON output.
"""
import argparse
import json
from pathlib import Path
import sys
import tempfile

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from met.composition_coverage import (  # noqa: E402,F401 -- compatibility imports
    SF, SA, S, M1, M2, MU, SIZES, N_READINGS, DRAWS,
    _one, _series, _sample, legacy_payload, print_legacy, run,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=2810)
    parser.add_argument("--out", type=Path, default=Path("results/composition_coverage.json"))
    args = parser.parse_args(argv)
    raw = yaml.safe_load((ROOT / "presets" / "composition_coverage.yaml").read_text(encoding="utf-8"))
    raw["replicates"], raw["seed"] = args.replicates, args.seed
    with tempfile.TemporaryDirectory(prefix="met-composition-preset-") as temporary:
        path = Path(temporary) / "preset.yaml"
        path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
        _, record = run(path, ROOT / "results", workers=args.workers, progress=None)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(legacy_payload(record), indent=1), encoding="utf-8")
    print_legacy(record)
    return 0


if __name__ == "__main__":
    main()

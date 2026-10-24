"""Shared helpers for analyses that rescore saved trial archives.

Nothing here refits an archived method.  Success means a positive finite
estimate within [m/K, K*m] with inclusive endpoints; every invalid estimate
is a failure (``comparison_scores.point_losses``).  Paired comparisons use
identical replicate IDs and ``comparison_scores.paired_summary``.
"""

import hashlib
import json
import math
from pathlib import Path
import platform
from statistics import NormalDist

import numpy as np

from comparison_scores import paired_summary, point_losses
from trial_archive import load_trial_archive, load_trial_manifest


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def label(value):
    return format(float(value), "g")


def failures(points, truth, factor):
    """Per-trial failure indicator (1 = outside factor or invalid)."""
    key = "outside_factor_" + label(factor)
    losses = point_losses(points, truth)
    if key not in losses:
        raise ValueError(f"factor {factor} is not one of the scored tolerance factors")
    return losses[key]


def paired(failures_a, failures_b):
    """Mean (A - B) failure difference and MCSE; negative favours A."""
    summary = paired_summary(failures_a, failures_b)
    return summary["mean_difference"], summary["mcse"]


def bonferroni_z(family_size, alpha=.05):
    if family_size < 1:
        raise ValueError("family must contain at least one comparison")
    return NormalDist().inv_cdf(1 - alpha / (2 * family_size))


def resolved(difference, mcse, z):
    """-1: A significantly better, 1: significantly worse, 0: unresolved."""
    if difference is None or mcse is None or mcse <= 0:
        return 0
    if difference < -z * mcse:
        return -1
    if difference > z * mcse:
        return 1
    return 0


def iter_cells(run_dir, cells=None):
    """Yield (entry, arrays) for each saved scenario of one experiment attempt."""
    run_dir = Path(run_dir)
    manifest = load_trial_manifest(run_dir)
    found = False
    for entry in manifest["entries"]:
        if cells is not None and entry["id"] not in cells:
            continue
        found = True
        yield entry, load_trial_archive(run_dir, entry)
    if not found:
        raise ValueError("no requested cells were found")


def point_names(arrays):
    return sorted(name[len("point__"):] for name in arrays if name.startswith("point__"))


def cell_record(entry):
    scenario = entry["scenario"]
    return {"cell": entry["id"], "archive_sha256": entry["sha256"],
            "true_mass": float(scenario["true_mass"]),
            "true_force_snr": float(scenario["true_force_snr"]),
            "true_acceleration_snr": float(scenario["true_acceleration_snr"])}


def provenance(run_dirs, sources):
    here = Path(__file__).resolve().parent
    return {"runs": [{"run": str(run), "manifest_sha256": sha256(Path(run) / "trial_manifest.json")}
                     for run in run_dirs],
            "runtime": {"python": platform.python_version(), "numpy": np.__version__},
            "source_sha256": {name: sha256(here / name) for name in sorted(set(sources) | {"archive_rescoring.py"})}}


def write_new_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x", encoding="utf-8") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def finite(value):
    return None if value is None or not math.isfinite(value) else float(value)

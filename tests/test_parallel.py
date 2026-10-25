"""Cells run side by side (met.runner, workers > 1) give the same output as one at a time.

Each cell draws its series from its own seed, so the only way the worker count could
show is through state a process carries from one cell to the next. Found on 25 Oct:
stacked_mean_radius widened its spline table on demand, so ratio_of_means_rss in the
200-reading cells moved by ~1e-13 with which cells had run earlier in the process; the
table is now a fixed ladder (tests/test_stacked.py).
"""

import json

import numpy as np
import pytest

from met import runner
from met.__main__ import main

PRESETS = ["presets/ours29_noise_vs_signal.yaml", "presets/repetition_new_excitation.yaml"]


def _load(out):
    run_dir = next(out.glob("*/*/"))
    values = np.load(run_dir / "values.npz")
    record = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    return {k: values[k] for k in values.files}, record, (run_dir / "summary.md").read_text(encoding="utf-8")


@pytest.mark.parametrize("preset", PRESETS)
def test_two_workers_give_the_serial_output_bit_for_bit(preset, tmp_path):
    reps = 3 if "ours29" in preset else 4
    runner.run(preset, tmp_path / "one", reps, progress=None, workers=1)
    runner.run(preset, tmp_path / "two", reps, progress=None, workers=2)
    v1, r1, s1 = _load(tmp_path / "one")
    v2, r2, s2 = _load(tmp_path / "two")
    assert sorted(v1) == sorted(v2)
    for key in v1:
        assert np.array_equal(v1[key], v2[key], equal_nan=True), key
    assert s1 == s2
    assert (r1["workers"], r2["workers"]) == (1, 2)
    strip = lambda cells: [{k: v for k, v in c.items() if k != "seconds"} for c in cells]   # noqa: E731
    assert strip(r1["cells"]) == strip(r2["cells"])


def test_worker_count_resolution(monkeypatch):
    monkeypatch.setattr(runner, "cpu_count", lambda: 8)
    assert runner.resolve_workers(None) == 4                     # half, by default
    assert runner.resolve_workers("max") == 8
    assert runner.resolve_workers(3) == 3
    monkeypatch.setattr(runner, "cpu_count", lambda: 1)
    assert runner.resolve_workers(None) == 1
    for bad in (0, -2, 1.5, "all", True):
        with pytest.raises(ValueError):
            runner.resolve_workers(bad)


def test_command_line_workers_flag(capsys):
    for bad in ("0", "two", "-1"):
        with pytest.raises(SystemExit):
            main(["run", PRESETS[0], "--workers", bad])
    assert "positive integer or 'max'" in capsys.readouterr().err


def test_a_failed_cell_stops_the_run_and_its_workers(tmp_path, monkeypatch):
    """An error in one worker reaches the caller, and no worker is left running."""
    calls = []
    real = runner._stop_workers
    monkeypatch.setattr(runner, "_stop_workers", lambda pool: (calls.append(1), real(pool)))
    monkeypatch.setattr(runner, "_run_cell_in_worker", _fails_on_cell_3)
    with pytest.raises(RuntimeError, match="cell 3"):
        runner.run(PRESETS[1], tmp_path, 2, progress=None, workers=2)
    assert calls == [1]


def _fails_on_cell_3(index):
    if index == 3:
        raise RuntimeError("cell 3 failed on purpose")
    return _real_cell(index)


def _real_cell(index):
    """The worker's own cell (this function runs in the worker, which loaded the preset)."""
    from met.runner import _WORKER, run_cell
    spec, (reps, seed, factors, batch_elements, compare_to) = _WORKER["spec"], _WORKER["args"]
    result, values = run_cell(spec["cells"][index], spec["estimators"], reps, seed, factors, batch_elements,
                              None, compare_to)
    return index, result, values, []

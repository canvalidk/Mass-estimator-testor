"""B01's fixed factorial declaration, validation and maintained runner route."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met import preset, runner


ROOT = Path(__file__).resolve().parents[1]


def recipe():
    return dict(kind="robot_budget", name="robot-budget-test", question="Check the robot budget.",
                runs=["1-baseline"], widths=[1.0], white=["autocorrelated"],
                floors=[False], rotations=["static"], bootstrap_replicates=12,
                seed=20261030, block_fraction=0.2)


def test_shipped_recipe_has_all_forty_eight_declared_cases():
    spec = preset.load(ROOT / "presets" / "robot_budget.yaml")
    assert spec["kind"] == "robot_budget"
    assert len(spec["cells"]) == 48
    assert spec["raw"]["seed"] == 20261030
    assert spec["raw"]["bootstrap_replicates"] == 500
    assert spec["raw"]["block_fraction"] == 0.2
    keys = ("run", "width", "white", "floor", "rotation")
    selections = {tuple(cell[key] for key in keys) for cell in spec["cells"]}
    assert len(selections) == 48
    assert {row[0] for row in selections} == {"1-baseline", "2-vibrations"}
    assert {row[1] for row in selections} == {.05, .2, 1.}
    assert {row[2] for row in selections} == {"independent", "autocorrelated"}
    assert {row[3] for row in selections} == {False, True}
    assert {row[4] for row in selections} == {"nominal", "static"}


@pytest.mark.parametrize("key", tuple(recipe()))
def test_recipe_requires_each_declared_setting(key):
    from met.robot_budget import validate

    data = recipe()
    del data[key]
    with pytest.raises(ValueError):
        validate(data)


@pytest.mark.parametrize("change", [
    lambda d: d.update(extra="not declared"),
    lambda d: d.update(name="../escape"),
    lambda d: d.update(question=""),
    lambda d: d.update(runs=[]),
    lambda d: d.update(runs=["3-vibrations-contact"]),
    lambda d: d.update(runs=["1-baseline", "1-baseline"]),
    lambda d: d.update(widths=[0.]),
    lambda d: d.update(widths=[float("nan")]),
    lambda d: d.update(widths=[True]),
    lambda d: d.update(widths=[.2, .2]),
    lambda d: d.update(white=["unknown"]),
    lambda d: d.update(white=["independent", "independent"]),
    lambda d: d.update(floors=[1]),
    lambda d: d.update(floors=[True, True]),
    lambda d: d.update(rotations=["fitted-on-dynamic"]),
    lambda d: d.update(seed=-1),
    lambda d: d.update(seed=True),
    lambda d: d.update(bootstrap_replicates=1),
    lambda d: d.update(bootstrap_replicates=True),
    lambda d: d.update(block_fraction=0.),
    lambda d: d.update(block_fraction=1.1),
])
def test_recipe_scope_and_selectors_are_validated(change):
    from met.robot_budget import validate

    data = deepcopy(recipe())
    change(data)
    with pytest.raises(ValueError):
        validate(data)


@pytest.mark.parametrize("workers", [2, "max", True, 1.0])
def test_runner_refuses_parallel_collection(tmp_path, workers):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    with pytest.raises(ValueError, match="serial"):
        runner.run(path, tmp_path, workers=workers, progress=None)


def test_runner_refuses_a_dynamic_replicate_override(tmp_path):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    with pytest.raises(ValueError, match="replicates"):
        runner.run(path, tmp_path, replicates=2, workers=1, progress=None)


def test_validate_cli_reports_robot_cases(capsys):
    from met.__main__ import main

    assert main(["validate", str(ROOT / "presets" / "robot_budget.yaml")]) == 0
    assert "48 robot-budget calibration case(s); serial" in capsys.readouterr().out


def test_list_cli_includes_robot_budget_question(capsys):
    from met.__main__ import main

    assert main(["list"]) == 0
    output = capsys.readouterr().out
    assert "robot_budget.yaml" in output
    assert "declared white-noise autocorrelation" in output


def test_recipe_replay_keeps_inputs_results_and_arrays(tmp_path):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    first, record = runner.run(path, tmp_path, workers=1, progress=None)
    second, replay = runner.run(path, tmp_path, workers=1, progress=None)
    assert first != second
    assert (first / "preset.yaml").read_bytes() == path.read_bytes()
    assert (second / "preset.yaml").read_bytes() == path.read_bytes()
    assert json.loads((first / "run.json").read_text(encoding="utf-8")) == record
    assert record["workers"] == replay["workers"] == 1
    assert record["coverage_reference_mass"] == 0.9323
    assert record["provenance"]
    assert record["rows"] == replay["rows"]
    with np.load(first / "values.npz") as a, np.load(second / "values.npz") as b:
        assert a.files
        assert sorted(a.files) == sorted(b.files)
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])
    assert (first / "summary.csv").is_file()
    assert (first / "summary.md").is_file()

"""B03: preserve the original composition study while adding preset replay.

The small fixture was captured before refactoring from the original script.
Frozen unprinted floats allow only rounding-scale platform variation; wrapper
and preset comparisons in the same process are exact. The slow check compares
the historical record's printed digits and does not select another equation.
"""
from copy import deepcopy
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest
import yaml

from met import preset, runner


ROOT = Path(__file__).resolve().parent.parent
PRESET = ROOT / "presets" / "composition_coverage.yaml"
ORIGINAL = json.loads((ROOT / "tests" / "data" / "composition_coverage_legacy_small.json")
                      .read_text(encoding="utf-8"))


def recipe():
    data = yaml.safe_load(PRESET.read_text(encoding="utf-8"))
    data["replicates"] = ORIGINAL["replicates"]
    return data


def _same_original(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            _same_original(actual[key], expected[key])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for got, want in zip(actual, expected):
            _same_original(got, want)
    elif isinstance(expected, bool):
        assert type(actual) is bool and actual is expected
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, rel=2e-12, abs=2e-12)
    else:
        assert actual == expected


def _wrapper():
    path = ROOT / "tools" / "composition_coverage.py"
    spec = importlib.util.spec_from_file_location("composition_coverage_wrapper", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixed_historical_preset_is_explicit():
    spec = preset.load(PRESET)
    assert spec["kind"] == "composition_coverage"
    assert len(spec["cells"]) == 3
    assert spec["raw"]["replicates"] == 400
    assert spec["raw"]["seed"] == 2810
    for key in ("m1", "m2", "sigma_f", "sigma_a", "c", "n_readings", "draws", "sizes"):
        assert spec["raw"][key] == ORIGINAL[key]
    assert spec["raw"]["dimension"] == 3 and spec["raw"]["equation"] == "ours29"


@pytest.mark.parametrize("job", ORIGINAL["jobs"])
def test_relocated_core_matches_the_pre_refactor_oracle(job):
    from met.composition_coverage import _one

    _same_original(_one((job["seed"], job["size"])), job["result"])


@pytest.fixture(scope="module")
def small_run(tmp_path_factory):
    folder = tmp_path_factory.mktemp("composition-coverage")
    path = folder / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe(), sort_keys=False), encoding="utf-8")
    out, record = runner.run(path, folder, workers=1, progress=None)
    return folder, path, out, record


def test_preset_preserves_the_original_seed_layout_and_summary(small_run):
    _, _, out, record = small_run
    _same_original(record["summary"], ORIGINAL["summary"])
    with np.load(out / "values.npz") as arrays:
        np.testing.assert_array_equal(arrays["job_seeds"].reshape(-1),
                                      [job["seed"] for job in ORIGINAL["jobs"]])


def test_wrapper_and_preset_agree_exactly_in_the_same_environment(small_run):
    folder, _, _, record = small_run
    target = folder / "wrapper.json"
    assert _wrapper().main(["--replicates", "2", "--workers", "1", "--seed", "2810",
                            "--out", str(target)]) in (None, 0)
    wrapper = json.loads(target.read_text(encoding="utf-8"))
    assert wrapper["summary"] == record["summary"]
    for key in ("seed", "replicates", "sizes", "n_readings", "sigma_f", "sigma_a", "c", "m1", "m2"):
        assert wrapper[key] == record[key]


def test_replay_preserves_effective_inputs_and_numerical_arrays(small_run):
    folder, source, out, record = small_run
    replay_path, replay = runner.run(source, folder, workers=1, progress=None)
    assert replay_path != out
    assert record["summary"] == replay["summary"]
    assert record["workers"] == replay["workers"] == 1
    assert (out / "source_preset.yaml").read_bytes() == source.read_bytes()
    assert yaml.safe_load((out / "preset.yaml").read_text(encoding="utf-8")) == recipe()
    assert json.loads((out / "run.json").read_text(encoding="utf-8")) == record
    assert (out / "summary.csv").is_file() and (out / "summary.md").is_file()
    with np.load(out / "values.npz") as one, np.load(replay_path / "values.npz") as two:
        assert one.files and sorted(one.files) == sorted(two.files)
        for key in one.files:
            np.testing.assert_array_equal(one[key], two[key])


def test_replicate_override_is_recorded_without_changing_the_source(small_run):
    folder, source, _, _ = small_run
    out, record = runner.run(source, folder, replicates=3, workers=1, progress=None)
    assert record["replicates"] == 3
    assert yaml.safe_load((out / "preset.yaml").read_text(encoding="utf-8"))["replicates"] == 3
    assert yaml.safe_load((out / "source_preset.yaml").read_text(encoding="utf-8"))["replicates"] == 2


@pytest.mark.parametrize("key", ["kind", "name", "question", "seed", "replicates", "sizes",
                                  "m1", "m2", "sigma_f", "sigma_a", "c", "n_readings",
                                  "draws", "dimension", "equation"])
def test_preset_requires_all_scientific_inputs(key):
    from met.composition_coverage import validate

    data = recipe()
    del data[key]
    with pytest.raises(ValueError):
        validate(data)


@pytest.mark.parametrize("change", [
    {"unexpected": 1}, {"name": "../escape"}, {"question": ""},
    {"seed": True}, {"seed": -1}, {"replicates": True}, {"replicates": 1},
    {"sizes": []}, {"sizes": [1., 1.]}, {"sizes": [-1.]}, {"sizes": [float("nan")]},
    {"sizes": [True]}, {"m1": 2.}, {"m2": 3.}, {"sigma_f": .3}, {"sigma_a": .2},
    {"c": 1.}, {"n_readings": 3}, {"n_readings": True}, {"draws": 100},
    {"dimension": 2}, {"equation": "ours30"},
])
def test_preset_refuses_changes_to_the_historical_scope(change):
    from met.composition_coverage import validate

    data = deepcopy(recipe())
    data.update(change)
    with pytest.raises(ValueError):
        validate(data)


@pytest.mark.parametrize("workers", [0, -1, "max", True, 1.0])
def test_runner_refuses_invalid_worker_count(tmp_path, workers):
    source = tmp_path / "recipe.yaml"
    source.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    with pytest.raises(ValueError, match="workers"):
        runner.run(source, tmp_path, workers=workers, progress=None)


def test_two_workers_preserve_all_numerical_results_exactly(small_run):
    folder, source, out, record = small_run
    parallel_out, parallel = runner.run(source, folder, workers=2, progress=None)
    assert parallel["workers"] == 2
    assert parallel["summary"] == record["summary"]
    with np.load(out / "values.npz") as serial_values, np.load(parallel_out / "values.npz") as parallel_values:
        assert sorted(serial_values.files) == sorted(parallel_values.files)
        for key in serial_values.files:
            np.testing.assert_array_equal(serial_values[key], parallel_values[key])


@pytest.mark.parametrize("replicates", [True, 1, 2.5])
def test_runner_refuses_invalid_replicate_override(tmp_path, replicates):
    source = tmp_path / "recipe.yaml"
    source.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    with pytest.raises(ValueError, match="replicates"):
        runner.run(source, tmp_path, replicates=replicates, workers=1, progress=None)


def test_validate_and_list_cli_expose_the_new_kind(capsys):
    from met.__main__ import main

    assert main(["validate", str(PRESET)]) == 0
    assert "3 composition-coverage signal-size cell(s); ours29, workers supported" in capsys.readouterr().out
    assert main(["list"]) == 0
    assert "composition_coverage.yaml" in capsys.readouterr().out


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="400 replicates per signal size")
def test_full_preset_reproduces_o2_printed_entries(tmp_path):
    _, record = runner.run(PRESET, tmp_path, workers=4, progress=None)
    # Printed fractions are the original script's .3f convention. In particular,
    # .9925 -> .993 corresponds to O2's 99.3%; .1% would expose a different tie.
    expected = {
        "2.0": ("0.954", ("0.948", "0.953", "0.993"), ("0.958", "0.970", "0.988")),
        "0.4": ("0.963", ("0.950", "0.955", "0.995"), ("0.948", "0.960", "1.000")),
        "0.1": ("0.926", ("0.927", "0.902", "0.988"), ("0.958", "0.920", "0.990")),
    }
    defects = {"2.0": (2, "5.84%", "5.85%"), "0.4": (1, "31.3%", "32.2%"),
               "0.1": (1, "91.4%", "90.5%")}
    for size, (parts, lump, reduced) in expected.items():
        row = record["summary"][size]
        assert f"{row['parts_coverage']:.3f}" == parts
        for name, wanted in (("L5", lump), ("L14", reduced)):
            got = tuple(f"{row[name][key]:.3f}" for key in ("computed_covers", "own_covers", "overlap"))
            assert got == wanted
        places, cloud, median = defects[size]
        assert format(row["L14"]["median_defect_cloud"], f".{places}%") == cloud
        assert format(row["L14"]["median_defect_median"], f".{places}%") == median

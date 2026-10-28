"""Deterministic recipe expansion, input replay and explicit unavailable points."""
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met.inspection.case import load_cases
from met.readout_candidates import CANDIDATE_IDS
from met.readout_study import run, validate


def recipe():
    return dict(kind="readout_candidates", name="readout-test", question="Compare one fixed law.",
        seed=0, prior_centers=[0.1, 1.0], null_padding=[
            dict(count=0, force_sd=1, acceleration_sd=1),
            dict(count=1, force_sd=3, acceleration_sd=1)],
        cases=[dict(name="perpendicular", question="A flat relative-evidence curve.",
            prior=dict(family="half_cauchy", center=1), truth={"private": [7]},
            metadata={"instrument": "fixed"}, readings=[dict(id="null-000", force=[1, 0],
                acceleration=[0, 1], covariance=np.eye(4).tolist(), window=[0, 1])])])


def test_expansion_preserves_inputs_truth_and_declarations_without_mutating_source():
    raw = recipe()
    original = deepcopy(raw)
    spec = validate(raw)
    assert raw == original
    assert [c["id"] for c in spec["cells"]] == ["c000", "c001", "c002", "c003"]
    assert [c["prior_center"] for c in spec["cells"]] == [0.1, 0.1, 1.0, 1.0]
    assert [c["readings"] for c in spec["cells"]] == [1, 2, 1, 2]
    for case in spec["cases"]:
        first = case.inputs.readings[0]
        np.testing.assert_array_equal(first.force, [1, 0])
        np.testing.assert_array_equal(first.acceleration, [0, 1])
        np.testing.assert_array_equal(first.covariance, np.eye(4))
        assert first.window == (0, 1)
        assert case.investigator_truth == {"private": [7]}
        assert case.metadata["instrument"] == "fixed"
    added = spec["cases"][1].inputs.readings[1]
    assert added.id != spec["cases"][1].inputs.readings[0].id
    np.testing.assert_array_equal(added.force, [0, 0])
    np.testing.assert_array_equal(added.acceleration, [0, 0])
    np.testing.assert_array_equal(added.covariance, np.diag([9, 9, 1, 1]))
    np.testing.assert_array_equal(added.subspace, np.eye(2))
    spec["cases"][1].investigator_truth["private"].append(8)
    assert spec["cases"][0].investigator_truth == {"private": [7]}
    assert raw == original


def test_shipped_recipe_has_84_cells_with_perpendicular_first_and_an_explicit_line():
    path = Path(__file__).resolve().parents[1] / "presets" / "readout_candidates.yaml"
    spec = validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    assert spec["raw"]["seed"] == 0
    assert len(spec["cases"]) == len(spec["cells"]) == 84
    assert [cell["source_case"] for cell in spec["cells"][:12]] == ["perpendicular"] * 12
    assert [spec["cells"][i]["source_case"] for i in range(0, 84, 12)] == [
        "perpendicular", "zero", "weak_aligned", "force_only", "acceleration_only",
        "strong_aligned", "outside_declared_span"]
    assert [cell["null_padding"]["count"] for cell in spec["cells"][:4]] == [0, 1, 1, 8]
    np.testing.assert_array_equal(spec["cases"][72].inputs.readings[0].subspace, [[1], [0]])
    np.testing.assert_allclose(spec["cases"][60].inputs.readings[0].covariance, 0.04 * np.eye(4))


@pytest.mark.parametrize("change", [
    lambda d: d.update(seed=-1), lambda d: d.update(seed=True),
    lambda d: d.update(name="../escape"), lambda d: d.update(question=""),
    lambda d: d.update(prior_centers=[]), lambda d: d.update(prior_centers=[1, 1]),
    lambda d: d.update(prior_centers=[0]), lambda d: d.update(prior_centers=[True]),
    lambda d: d.update(prior_centers=[float("nan")]), lambda d: d.update(null_padding=[]),
    lambda d: d["null_padding"][0].update(count=True),
    lambda d: d["null_padding"][0].update(count=-1),
    lambda d: d["null_padding"][0].update(force_sd=0),
    lambda d: d["null_padding"][0].update(acceleration_sd=float("inf")),
    lambda d: d["null_padding"].append(deepcopy(d["null_padding"][0])),
    lambda d: d["cases"][0]["prior"].update(family="unknown"),
    lambda d: d.update(candidate_module="untrusted.module"),
])
def test_invalid_recipe_is_rejected_before_evaluation(change):
    data = recipe()
    change(data)
    with pytest.raises(ValueError):
        validate(data)


def test_study_saves_native_arrays_undefined_status_and_exact_replay(tmp_path, monkeypatch):
    import met.readout_study as study
    import met.runner as runner

    data = recipe()
    data["prior_centers"] = [1]
    path = tmp_path / "source.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    calls = []

    def fake_evaluate(readings, center):
        calls.append((readings, center))
        return dict(law=dict(log_mass=np.log([0.5, 1, 2]), probability=np.array([0.25, 0.5, 0.25]),
            median=1.0, interval=[0.5, 2], raw_interval=[0.5, 2], eta=1.0,
            g=np.ones(len(readings)), Q=np.zeros(len(readings))),
            candidates={name: dict(status="undefined", point=None) if name == "mean_vector"
                        else dict(status="defined", point=1.0) for name in CANDIDATE_IDS},
            weights={name: np.ones(3) for name in
                     ("cloud_mean", "mean_vector", "conditional_rms", "plain_sum_cloud")},
            second_moment=np.array([4, 2, 1]), assumptions=["A replaceable ordinary evaluation function."])

    def no_rng(*args, **kwargs):
        raise AssertionError("deterministic recipe must not draw random numbers")

    monkeypatch.setattr(study, "evaluate", fake_evaluate)
    monkeypatch.setattr(runner, "code_version", lambda: {"test_version": True})
    monkeypatch.setattr(np.random, "default_rng", no_rng)
    first, record = run(path, tmp_path / "results", progress=None)
    second, repeated = run(path, tmp_path / "results", progress=None)
    assert first != second
    assert len(calls) == 4
    assert [len(readings) for readings, _ in calls] == [1, 2, 1, 2]
    assert all(center == 1 for _, center in calls)
    assert record["workers"] == 1 and record["random_draws"] == 0 and record["deterministic"] is True
    assert (first / "source_preset.yaml").read_bytes() == path.read_bytes()
    assert validate(yaml.safe_load((first / "preset.yaml").read_text(encoding="utf-8")))["cells"] == validate(data)["cells"]
    assert record["code"]["preset_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert {"case.py", "paths.py"} <= record["code"]["inspection_file_sha256"].keys()
    assert all(len(h) == 64 for h in record["code"]["inspection_file_sha256"].values())
    assert set(record["code"]["dependencies"]) == {"numpy", "scipy", "PyYAML"}
    assert record["code"]["effective_preset_sha256"] == hashlib.sha256((first / "preset.yaml").read_bytes()).hexdigest()
    saved = json.loads((first / "run.json").read_text(encoding="utf-8"))
    assert saved["cells"] == record["cells"] == repeated["cells"]
    expected = validate(data)["cases"]
    for row, case in zip(record["cells"], expected):
        folder = first / row["id"]
        document, cases = load_cases(folder / "case.yaml")
        assert document["kind"] == "inspection" and document["seed"] == 0
        assert cases[0].to_dict(include_truth=True) == case.to_dict(include_truth=True)
        assert row["candidates"]["mean_vector"] == {"status": "undefined", "point": None}
        assert json.loads((folder / "points.json").read_text(encoding="utf-8")) == row
        with np.load(folder / "arrays.npz", allow_pickle=False) as values, np.load(second / row["id"] / "arrays.npz", allow_pickle=False) as replayed:
            assert sorted(values.files) == sorted(replayed.files)
            for key in values.files:
                np.testing.assert_array_equal(values[key], replayed[key])
            np.testing.assert_array_equal(values["S2"], [4, 2, 1])
            np.testing.assert_array_equal(values["probability"], [0.25, 0.5, 0.25])
            np.testing.assert_array_equal(values["mass"], document["numerics"]["mass_grid"])
            np.testing.assert_array_equal(values["weights/cloud_mean"], np.ones(3))
    with (first / "summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12
    undefined = [row for row in rows if row["candidate"] == "mean_vector"]
    assert len(undefined) == 2 and all(row["point"] == "" and row["status"] == "undefined" for row in undefined)
    checks = json.loads((first / "finite_law_checks.json").read_text(encoding="utf-8"))
    assert checks == record["finite_law_checks"]
    assert checks["cancellation"]["status"] == "undefined"
    assert {row["candidate"] for row in checks["composition"]} == {
        "cloud_mean", "mean_vector", "conditional_rms", "global_rms", "law_median"}
    assert "undefined" in (first / "summary.md").read_text(encoding="utf-8")


def test_runner_and_cli_preserve_deterministic_resource_policy(tmp_path, monkeypatch, capsys):
    from met import runner
    from met.__main__ import main
    import met.readout_study as study
    import os

    path = tmp_path / "source.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    calls = []

    def fake_run(path, out_root, progress):
        calls.append((path, out_root))
        return out_root, {}

    monkeypatch.setattr(study, "run", fake_run)
    for workers in (2, "max", True):
        with pytest.raises(ValueError, match="deterministic serial"):
            runner.run(path, tmp_path, workers=workers, progress=None)
    with pytest.raises(ValueError, match="omit --replicates"):
        runner.run(path, tmp_path, replicates=2, workers=1, progress=None)
    assert calls == []
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        monkeypatch.setenv(key, "4")
    assert main(["run", str(path), "--out", str(tmp_path), "--workers", "1"]) == 0
    assert calls == [(path, tmp_path)]
    assert all(os.environ[key] == "1" for key in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"))
    assert main(["validate", str(path)]) == 0
    assert "serial, no random draws" in capsys.readouterr().out


def test_real_saved_case_replays_all_candidate_points_and_native_arrays(tmp_path):
    from met.full import Reading
    from met.readout_candidates import evaluate

    data = recipe()
    data["prior_centers"] = [.1]
    data["null_padding"] = [dict(count=1, force_sd=3, acceleration_sd=1)]
    path = tmp_path / "source.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    out, record = run(path, tmp_path / "results", progress=None)
    _, cases = load_cases(out / "c000" / "case.yaml")
    inputs = cases[0].inputs
    readings = [Reading(r.force, r.acceleration, r.covariance, r.subspace, r.window)
                for r in inputs.readings]
    replayed = evaluate(readings, inputs.prior["center"])
    for name in CANDIDATE_IDS:
        saved = record["cells"][0]["candidates"][name]
        assert saved["status"] == replayed["candidates"][name]["status"]
        assert saved["point"] == pytest.approx(replayed["candidates"][name]["point"], rel=1e-12)
    with np.load(out / "c000" / "arrays.npz", allow_pickle=False) as saved:
        np.testing.assert_array_equal(saved["log_mass"], replayed["law"]["log_mass"])
        np.testing.assert_array_equal(saved["probability"], replayed["law"]["probability"])
        np.testing.assert_array_equal(saved["S2"], replayed["second_moment"])
        for name, weights in replayed["weights"].items():
            np.testing.assert_allclose(saved[f"weights/{name}"], weights, rtol=1e-12)

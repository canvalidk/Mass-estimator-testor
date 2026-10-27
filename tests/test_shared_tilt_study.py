"""Independent matched-pattern reference, replay and scope of the O16 recipe."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met.shared_error_geometry import rotation
from met.shared_tilt import investigator_matched_score, run, simulate_cell, validate


def recipe():
    return dict(kind="shared_tilt", name="tilt-test", question="Check a shared tilt field.",
        cases=["broad"], angles_degrees=[0.0, 5.0], readings=[6], replicates=4, seed=20261027,
        mass=1.5, acceleration_length=2.0, force_sd=0.02, acceleration_sd=0.02,
        rotation_axis=[0.0, 0.0, 1.0], matched_hypothesis_degrees=5.0, robot=None)


@pytest.mark.parametrize("degrees", [0.0, 5.0])
def test_investigator_matched_score_obeys_exact_conditional_gaussian_reference(degrees):
    """A predeclared template and latent vectors make this an unfitted N(0,1) check."""
    rng = np.random.default_rng(20261027 + int(degrees))
    b, n, mass, sf, sa = 3000, 8, 1.5, 0.02, 0.02
    latent = np.broadcast_to(np.array([2.0, 0.0, 0.0]), (b, n, 3)).copy()
    injected = rotation([0, 0, 1], np.deg2rad(degrees))
    hypothesis = rotation([0, 0, 1], np.deg2rad(5.0))
    f = mass * latent @ injected.T + sf * rng.normal(size=latent.shape)
    a = latent + sa * rng.normal(size=latent.shape)
    observed = investigator_matched_score(f, a, latent, mass, injected, hypothesis, sf, sa)
    score = observed["centered_score"]
    assert abs(score.mean()) < 6 / np.sqrt(b)
    assert abs(score.var(ddof=1) - 1) < 6 * np.sqrt(2 / (b - 1))
    if degrees == 0:
        np.testing.assert_array_equal(observed["matched_mean"], np.zeros(b))
    else:
        np.testing.assert_allclose(observed["matched_mean"], np.sqrt(observed["matched_information"]), rtol=1e-13)


def test_axis_hypothesis_is_explicitly_inapplicable():
    data = recipe()
    data["cases"] = ["line_axis"]
    spec = validate(data)
    values = simulate_cell(spec["cells"][0], spec["settings"], 4, data["seed"])
    np.testing.assert_array_equal(values["matched_information"], np.zeros(4))
    assert np.all(np.isnan(values["matched_score"]))
    assert np.all(np.isnan(values["centered_score"]))
    # The unavailable axial rotation is excluded from the investigator metric.
    np.testing.assert_allclose(values["observable_rotation_norm"],
                               np.linalg.norm(values["rotation_vector"][:, :2], axis=1))


def test_recipe_replay_retains_prefix_and_exact_effective_settings(tmp_path):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    first, record = run(path, tmp_path, progress=None)
    second, _ = run(path, tmp_path, progress=None)
    longer, _ = run(path, tmp_path, replicates=6, progress=None)
    with np.load(first / "values.npz") as a, np.load(second / "values.npz") as b, np.load(longer / "values.npz") as more:
        assert sorted(a.files) == sorted(b.files)
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])
            np.testing.assert_array_equal(a[key], more[key][:4])
    assert (first / "preset.yaml").read_bytes() == path.read_bytes()
    assert yaml.safe_load((longer / "preset.yaml").read_text(encoding="utf-8"))["replicates"] == 6
    assert (longer / "source_preset.yaml").read_bytes() == path.read_bytes()
    restored = json.loads((first / "run.json").read_text(encoding="utf-8"))
    assert restored["workers"] == record["workers"] == 1
    assert restored["robot"] == []
    assert restored["code"]["effective_preset_sha256"]
    assert (first / "summary.csv").exists() and (first / "summary.md").exists()


@pytest.mark.parametrize("change", [
    lambda d: d.update(seed=-1), lambda d: d.update(replicates=True),
    lambda d: d.update(cases=["unknown"]), lambda d: d.update(readings=[2]),
    lambda d: d.update(rotation_axis=[0, 0, 0]), lambda d: d.update(force_sd=0),
    lambda d: d.update(mass=float("nan")), lambda d: d.update(name="../escape"),
    lambda d: d.update(matched_hypothesis_degrees=0),
    lambda d: d.update(robot={"widths": [0.2], "floors": [1]}),
])
def test_recipe_scope_is_validated(change):
    data = deepcopy(recipe())
    change(data)
    with pytest.raises(ValueError):
        validate(data)


def test_shipped_recipe_has_twenty_four_seeded_cells():
    from met.preset import load

    spec = load(Path(__file__).resolve().parents[1] / "presets" / "shared_sensor_tilt.yaml")
    assert spec["kind"] == "shared_tilt"
    assert len(spec["cells"]) == 24
    assert spec["raw"]["seed"] == 20261027


def test_existing_runner_refuses_parallel_shared_tilt(tmp_path):
    from met.runner import run as runner

    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    for workers in (2, "max", True, 1.0):
        with pytest.raises(ValueError, match="serial"):
            runner(path, tmp_path, workers=workers, progress=None)

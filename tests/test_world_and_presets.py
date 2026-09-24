"""The world generator, preset validation and a small end-to-end run."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met import preset as preset_mod
from met.runner import run
from met.world import generate

ROOT = Path(__file__).resolve().parent.parent

WORLD = {"dimension": 3, "design": "same_pair", "mass": 2.0, "acceleration_snr": 3.0,
         "direction": "random", "readings": 4, "noise": {"force_sd": 1.0, "acceleration_sd": 0.5},
         "supplied_noise": "exact"}


def test_same_pair_shares_the_latent_pair_and_new_excitation_does_not():
    rng = np.random.default_rng(0)
    _, truth = generate(WORLD, 50, rng)
    assert np.allclose(truth["true_acceleration"], truth["true_acceleration"][:, :1])
    other = dict(WORLD, design="new_excitation")
    _, truth = generate(other, 50, rng)
    assert not np.allclose(truth["true_acceleration"], truth["true_acceleration"][:, :1])
    assert np.allclose(truth["true_force"], 2.0 * truth["true_acceleration"])
    norms = np.linalg.norm(truth["true_acceleration"], axis=2)
    assert np.allclose(norms, 3.0 * 0.5)


def test_noise_and_supplied_noise():
    rng = np.random.default_rng(1)
    readings, truth = generate(dict(WORLD, readings=1), 40000, rng)
    err_f = readings.force - truth["true_force"]
    err_a = readings.acceleration - truth["true_acceleration"]
    assert np.std(err_f) == pytest.approx(1.0, rel=0.02)
    assert np.std(err_a) == pytest.approx(0.5, rel=0.02)
    assert np.allclose(readings.force_sd, 1.0) and np.allclose(readings.acceleration_sd, 0.5)
    scaled = dict(WORLD, supplied_noise={"force_scale": 2.0, "acceleration_scale": 0.5})
    readings, _ = generate(scaled, 3, rng)
    assert np.allclose(readings.force_sd, 2.0) and np.allclose(readings.acceleration_sd, 0.25)


def test_zero_excitation_and_one_dimension():
    rng = np.random.default_rng(2)
    _, truth = generate(dict(WORLD, acceleration_snr=0.0), 5, rng)
    assert np.all(truth["true_acceleration"] == 0)
    readings, truth = generate(dict(WORLD, dimension=1), 200, rng)
    assert readings.dimension == 1
    assert set(np.unique(np.sign(truth["true_acceleration"]))) == {-1.0, 1.0}


def _minimal():
    return {"name": "t", "question": "q?", "world": copy.deepcopy(WORLD),
            "estimators": [{"name": "flat", "reduce": "none",
                            "law": {"reference": "flat", "nuisance": "radius"},
                            "combine": {"rule": "likelihood_product",
                                        "prior": [{"uniform_angle": {"center": "first_reading"}}]},
                            "readouts": ["ratio_of_means", "interval_95"]},
                           {"name": "len", "direct": "norm_ratio"}],
            "scores": {"within_factor": [1.5]}, "replicates": 20, "seed": 1}


@pytest.mark.parametrize("path,value,message", [
    (("world", "design"), None, "design"),
    (("world", "supplied_noise"), None, "supplied_noise"),
    (("estimators", 0, "combine", "prior"), None, "prior"),
    (("estimators", 0, "law", "nuisance"), None, "reference and nuisance"),
    (("seed",), None, "seed"),
])
def test_required_settings_have_no_default(path, value, message):
    data = _minimal()
    target = data
    for key in path[:-1]:
        target = target[key]
    del target[path[-1]]
    with pytest.raises(ValueError, match=message):
        preset_mod.validate(data)


def test_prior_is_never_implicit():
    data = _minimal()
    data["estimators"][0]["combine"]["prior"] = []
    with pytest.raises(ValueError, match="never implicit"):
        preset_mod.validate(data)


def test_grid_expansion():
    data = _minimal()
    data["world"]["mass"] = [0.5, 1, 2]
    data["world"]["noise"]["force_sd"] = [1.0, 2.0]
    cells = preset_mod.validate(data)["cells"]
    assert len(cells) == 6
    assert cells[-1]["axes"] == {"mass": 2, "noise.force_sd": 2.0}


def test_all_shipped_presets_validate():
    paths = sorted((ROOT / "presets").glob("*.yaml"))
    assert paths
    for path in paths:
        preset_mod.load(path)


def test_end_to_end_run_is_reproducible(tmp_path):
    path = tmp_path / "p.yaml"
    data = _minimal()
    data["world"]["mass"] = [1.0, 3.0]
    path.write_text(yaml.safe_dump(data))
    out1, rec1 = run(path, tmp_path / "a", progress=None)
    out2, rec2 = run(path, tmp_path / "b", progress=None)
    assert (out1 / "summary.md").exists() and (out1 / "values.npz").exists()
    v1, v2 = np.load(out1 / "values.npz"), np.load(out2 / "values.npz")
    for key in v1.files:
        assert np.array_equal(v1[key], v2[key], equal_nan=True)
    cell = json.loads((out1 / "run.json").read_text())["cells"][0]
    assert set(cell["scores"]["points"]) == {"flat.ratio_of_means", "len.norm_ratio"}
    assert "flat.interval_95" in cell["scores"]["intervals"]

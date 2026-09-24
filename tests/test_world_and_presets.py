"""The world generator, preset validation and a small end-to-end run."""

import copy
import json

import numpy as np
import pytest
import yaml

from pathlib import Path

from met import preset as preset_mod
from met.runner import run
from met.world import generate, series_rngs

ROOT = Path(__file__).resolve().parent.parent

WORLD = {"dimension": 3, "design": "same_pair", "mass": 2.0, "acceleration_snr": 3.0,
         "direction": "random", "readings": 4, "noise": {"force_sd": 1.0, "acceleration_sd": 0.5},
         "supplied_noise": "exact"}


def _gen(world, count, seed=0):
    return generate(world, series_rngs(seed, world, count))


def test_same_pair_shares_the_latent_pair_and_new_excitation_does_not():
    _, truth = _gen(WORLD, 50)
    assert np.allclose(truth["true_acceleration"], truth["true_acceleration"][:, :1])
    _, truth = _gen(dict(WORLD, design="new_excitation"), 50)
    assert not np.allclose(truth["true_acceleration"], truth["true_acceleration"][:, :1])
    assert np.allclose(truth["true_force"], 2.0 * truth["true_acceleration"])
    assert np.allclose(np.linalg.norm(truth["true_acceleration"], axis=2), 3.0 * 0.5)


def test_noise_and_supplied_noise():
    readings, truth = _gen(dict(WORLD, readings=1), 20000)
    assert np.std(readings.force - truth["true_force"]) == pytest.approx(1.0, rel=0.02)
    assert np.std(readings.acceleration - truth["true_acceleration"]) == pytest.approx(0.5, rel=0.02)
    assert np.allclose(readings.force_sd, 1.0) and np.allclose(readings.acceleration_sd, 0.5)
    readings, _ = _gen(dict(WORLD, supplied_noise={"force_scale": 2.0, "acceleration_scale": 0.5}), 3)
    assert np.allclose(readings.force_sd, 2.0) and np.allclose(readings.acceleration_sd, 0.25)


def test_zero_excitation_and_one_dimension():
    _, truth = _gen(dict(WORLD, acceleration_snr=0.0), 5)
    assert np.all(truth["true_acceleration"] == 0)
    readings, truth = _gen(dict(WORLD, dimension=1), 200)
    assert readings.dimension == 1
    assert set(np.unique(np.sign(truth["true_acceleration"]))) == {-1.0, 1.0}


def test_a_quick_run_sees_the_first_series_of_the_full_run():
    small, _ = _gen(WORLD, 5, seed=7)
    large, _ = _gen(WORLD, 50, seed=7)
    assert np.array_equal(small.force, large.force[:5])
    assert np.array_equal(small.acceleration, large.acceleration[:5])


def test_cell_data_do_not_depend_on_grid_position():
    data = _minimal()
    data["world"]["mass"] = [1.0, 3.0]
    cells_a = preset_mod.validate(data)["cells"]
    data["world"]["mass"] = [0.5, 1.0, 3.0]
    cells_b = preset_mod.validate(data)["cells"]
    a, _ = _gen(cells_a[1]["world"], 3, seed=1)
    b, _ = _gen(cells_b[2]["world"], 3, seed=1)
    assert np.array_equal(a.force, b.force)


def _minimal():
    return {"name": "t", "question": "q?", "world": copy.deepcopy(WORLD),
            "estimators": [{"name": "flat", "reduce": "none",
                            "law": {"reference": "flat", "nuisance": "radius"},
                            "combine": {"rule": "likelihood_product",
                                        "prior": [{"uniform_angle": {"center": "first_reading"}}]},
                            "readouts": ["ratio_of_means", "interval_95"]},
                           {"name": "len", "direct": "norm_ratio", "reduce": "pool_pair"}],
            "scores": {"within_factor": [1.5]}, "replicates": 20, "seed": 1}


def _refuses(data, message):
    with pytest.raises(ValueError, match=message):
        preset_mod.validate(data)


REQUIRED = [("world", "dimension"), ("world", "design"), ("world", "mass"), ("world", "acceleration_snr"),
            ("world", "direction"), ("world", "readings"), ("world", "noise"), ("world", "supplied_noise"),
            ("estimators", 0, "reduce"), ("estimators", 0, "law"), ("estimators", 0, "combine"),
            ("estimators", 0, "readouts"), ("estimators", 0, "combine", "prior"),
            ("estimators", 0, "combine", "rule"), ("estimators", 0, "law", "nuisance"),
            ("estimators", 0, "law", "reference"), ("estimators", 1, "reduce"),
            ("scores",), ("replicates",), ("seed",), ("question",)]


@pytest.mark.parametrize("path", REQUIRED, ids=lambda p: ".".join(map(str, p)))
def test_required_settings_have_no_default(path):
    data = _minimal()
    target = data
    for key in path[:-1]:
        target = target[key]
    del target[path[-1]]
    with pytest.raises(ValueError):
        preset_mod.validate(data)


def test_prior_is_never_implicit():
    data = _minimal()
    data["estimators"][0]["combine"]["prior"] = []
    _refuses(data, "never implicit")


@pytest.mark.parametrize("combine,nuisance", [
    ({"rule": "likelihood_product", "prior": ["flat_log_mass"]}, "radius"),
    ({"rule": "likelihood_product", "prior": ["flat_mass"]}, "radius"),
    ({"rule": "likelihood_product", "prior": ["flat_log_mass"]}, "acceleration"),
    ({"rule": "posterior_product", "coordinate": "log_mass", "prior": ["flat_mass"]}, "radius"),
])
def test_improper_laws_are_refused_before_running(combine, nuisance):
    data = _minimal()
    data["world"]["readings"] = 1
    data["estimators"][0]["combine"] = combine
    data["estimators"][0]["law"]["nuisance"] = nuisance
    _refuses(data, "cannot be normalised")


@pytest.mark.parametrize("prior", [["flat_mass"], ["flat_inverse_mass"]])
def test_improper_single_reading_priors_are_refused(prior):
    data = _minimal()
    data["world"]["readings"] = 1
    data["estimators"][0]["combine"] = {"rule": "single", "prior": prior}
    _refuses(data, "cannot be normalised")


def test_single_rule_needs_one_reading_in_every_cell():
    data = _minimal()
    data["world"]["readings"] = [1, 4]
    data["estimators"][0]["combine"] = {"rule": "single", "prior": []}
    _refuses(data, "exactly one reading")


def test_pooling_a_new_excitation_world_must_be_declared():
    data = _minimal()
    data["world"]["design"] = "new_excitation"
    _refuses(data, "declared_mismatch")
    data["estimators"][1]["declared_mismatch"] = "on purpose"
    preset_mod.validate(data)


def test_new_excitation_that_cannot_vary_is_refused():
    data = _minimal()
    data["world"].update(design="new_excitation", direction="fixed")
    _refuses(data, "same latent pair")


@pytest.mark.parametrize("change,message", [
    (lambda d: d["world"].update(dimension=3.0), "dimension"),
    (lambda d: d["estimators"][0].update(numerics={"grid_point": 16001}), "unknown numerics"),
    (lambda d: d["estimators"][0].update(numerics={"cutoff": "big"}), "cutoff"),
    (lambda d: d["estimators"][0].update(numerics={"grid_points": 2001.0}), "grid_points"),
    (lambda d: d["estimators"][0]["combine"].update(
        prior=[{"lognormal": {"center": 1.0, "width": 0}}]), "width"),
    (lambda d: d["estimators"][0]["combine"].update(
        prior=[{"sech_tilt": {"lambda": "x", "center": "first_reading"}}]), "lambda"),
    (lambda d: d["estimators"][0]["combine"].update(
        prior=[{"uniform_angle": {"center": True}}]), "centre"),
    (lambda d: d["estimators"][1].pop("reduce"), "reduce: pool_pair"),
])
def test_validation_holes_are_closed(change, message):
    data = _minimal()
    change(data)
    _refuses(data, message)


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
    out1, _ = run(path, tmp_path / "a", progress=None)
    out2, _ = run(path, tmp_path / "a", progress=None)
    assert out1 != out2
    assert (out1 / "summary.md").exists() and (out1 / "values.npz").exists()
    v1, v2 = np.load(out1 / "values.npz"), np.load(out2 / "values.npz")
    for key in v1.files:
        assert np.array_equal(v1[key], v2[key], equal_nan=True)
    cell = json.loads((out1 / "run.json").read_text())["cells"][0]
    assert set(cell["scores"]["points"]) == {"flat.ratio_of_means", "len.norm_ratio"}
    assert "flat.interval_95" in cell["scores"]["intervals"]

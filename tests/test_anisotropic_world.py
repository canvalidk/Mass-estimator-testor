"""O17 worlds: investigator references, paired inputs and ordinary replay."""
from copy import deepcopy
import inspect
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met.anisotropic_world import (evaluate, generate_cell, investigator_scalar_compatibility,
                                   known_tensor_cost, measured_scalar_diagnostics, run, validate)
from met.law import ReadingSet


def recipe():
    return dict(kind="anisotropic_world", name="tensor-test", question="Check scalar compatibility.",
        tensors=[dict(name="isotropic", matrix=(2 * np.eye(3)).tolist())],
        geometries=["principal_line"], readings=[3], noise_sd=[.1], replicates=2, seed=20261027,
        acceleration_length=1.0, prior_center=1.0,
        numerics={"representative_log_half_span": 6.0, "representative_nodes": 101})


def test_estimator_boundary_has_no_tensor_or_latent_argument():
    assert list(inspect.signature(evaluate).parameters) == ["readings", "prior_center"]


def test_principal_eigenline_hides_tensor_anisotropy_in_paired_measured_inputs():
    data = recipe()
    data["tensors"].append(dict(name="anisotropic", matrix=np.diag([1, 2, 4]).tolist()))
    spec = validate(data)
    first, truth1 = generate_cell(spec["cells"][0], spec["settings"], 3, data["seed"])
    second, truth2 = generate_cell(spec["cells"][1], spec["settings"], 3, data["seed"])
    np.testing.assert_array_equal(first.force, second.force)
    np.testing.assert_array_equal(first.acceleration, second.acceleration)
    one, two = evaluate(first, 1), evaluate(second, 1)
    for key in one:
        np.testing.assert_array_equal(one[key], two[key])
    for truth in (truth1, truth2):
        compatible = investigator_scalar_compatibility(truth["true_force"], truth["true_acceleration"])
        np.testing.assert_array_equal(compatible["compatible_scalar"], np.full(3, 2.0))


@pytest.mark.parametrize("geometry,exact", [("principal_line", True), ("principal_mix", False),
                                            ("oblique_line", False), ("broad", False)])
def test_latent_scalar_status_is_kept_outside_the_estimator(geometry, exact):
    data = recipe()
    data["tensors"][0]["matrix"] = np.diag([1, 2, 4]).tolist()
    data["geometries"] = [geometry]
    spec = validate(data)
    _, truth = generate_cell(spec["cells"][0], spec["settings"], 3, data["seed"])
    compatible = investigator_scalar_compatibility(truth["true_force"], truth["true_acceleration"])
    np.testing.assert_array_equal(compatible["scalar_exact"], np.full(3, exact))
    if not exact:
        assert np.all(np.isnan(compatible["compatible_scalar"]))


def test_known_tensor_unfitted_cost_matches_exact_chi_square_moments():
    rng = np.random.default_rng(170027)
    repetitions, n, sd = 3000, 3, .1
    tensor = np.diag([1., 2., 4.])
    latent = np.broadcast_to(np.eye(3), (repetitions, n, 3))
    force = latent @ tensor.T + sd * rng.normal(size=latent.shape)
    acceleration = latent + sd * rng.normal(size=latent.shape)
    cost = known_tensor_cost(force, acceleration, tensor, sd)
    df = 3 * n
    assert abs(cost.mean() - df) < 6 * np.sqrt(2 * df / repetitions)
    assert abs(cost.var(ddof=1) - 2 * df) < .15 * 2 * df


@pytest.mark.parametrize("scale", [1e-9, 1.0, 1e9])
def test_scalar_compatibility_tolerance_is_relative_and_scale_covariant(scale):
    acceleration = scale * np.array([[[1., 1., 0.]]])
    force = acceleration @ np.diag([1., 2., 4.])
    result = investigator_scalar_compatibility(force, acceleration)
    assert result["scalar_status"][0] == "none"
    assert not result["scalar_exact"][0]


def test_exact_zero_latent_experiment_has_every_scalar_but_no_unique_target():
    result = investigator_scalar_compatibility(np.zeros((1, 2, 3)), np.zeros((1, 2, 3)))
    assert result["scalar_status"][0] == "every"
    assert np.isnan(result["compatible_scalar"][0])


def test_scalar_conservative_tail_is_valid_against_a_fixed_true_mass_cost():
    rng = np.random.default_rng(170028)
    f0, a0 = 2 * np.eye(3), np.eye(3)
    force = f0[None] + .1 * rng.normal(size=(10, 3, 3))
    acceleration = a0[None] + .1 * rng.normal(size=(10, 3, 3))
    readings = ReadingSet(force, acceleration, np.full((10, 3), .1), np.full((10, 3), .1))
    diagnostic = measured_scalar_diagnostics(readings)
    fixed_cost = np.sum((force - 2 * acceleration) ** 2, axis=(1, 2)) / (.1 ** 2 * 5)
    assert np.all(diagnostic["common_scalar_cost"] <= fixed_cost + 1e-10)
    from scipy.stats import chi2
    assert np.all(diagnostic["scalar_h0_conservative_p"] >= chi2.sf(fixed_cost, 9) - 1e-14)


def test_recipe_replays_prefix_and_saves_effective_override_and_truth_separately(tmp_path):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    first, record = run(path, tmp_path, progress=None)
    second, _ = run(path, tmp_path, progress=None)
    longer, _ = run(path, tmp_path, replicates=3, progress=None)
    with np.load(first / "values.npz") as a, np.load(second / "values.npz") as b, np.load(longer / "values.npz") as more:
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])
            if "/representative/" not in key:
                np.testing.assert_array_equal(a[key], more[key][:2])
    assert (first / "preset.yaml").read_bytes() == path.read_bytes()
    assert yaml.safe_load((longer / "preset.yaml").read_text(encoding="utf-8"))["replicates"] == 3
    replay = yaml.safe_load((first / "c000" / "case.yaml").read_text(encoding="utf-8"))
    assert replay["truth"]["tensor"] == (2 * np.eye(3)).tolist()
    assert "mass" not in replay["truth"]
    restored = json.loads((first / "run.json").read_text(encoding="utf-8"))
    assert restored["workers"] == record["workers"] == 1
    assert restored["cells"][0]["representative"]["grid_role"].startswith("Display")


@pytest.mark.parametrize("change", [lambda d: d.update(replicates=True), lambda d: d.update(seed=-1),
    lambda d: d.update(name="../escape"), lambda d: d.update(geometries=["unknown"]),
    lambda d: d.update(readings=[0]), lambda d: d.update(noise_sd=[0]),
    lambda d: d["tensors"][0].update(matrix=np.diag([-1, 2, 4]).tolist()),
    lambda d: d["tensors"][0].update(matrix=[[2, .1, 0], [.1, 2, 0], [0, 0, 2]]),
    lambda d: d.update(prior_center=float("nan"))])
def test_scope_validation(change):
    data = deepcopy(recipe())
    change(data)
    with pytest.raises(ValueError):
        validate(data)


def test_shipped_recipe_has_forty_eight_seeded_worlds():
    from met.preset import load
    spec = load(Path(__file__).resolve().parents[1] / "presets" / "anisotropic_world.yaml")
    assert len(spec["cells"]) == 48
    assert spec["raw"]["seed"] == 20261027


def test_existing_runner_rejects_parallel_tensor_recipe(tmp_path):
    from met.runner import run as runner
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    for workers in (2, "max", True, 1.0):
        with pytest.raises(ValueError, match="serial"):
            runner(path, tmp_path, workers=workers, progress=None)

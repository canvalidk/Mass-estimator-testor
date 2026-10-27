"""The O14 recipe observes readings only and retains reproducible paired draws."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from met.law import ReadingSet
from met.pooled_scale import residual_statistics, run, validate


def recipe():
    return dict(kind="pooled_scale", name="pooled-test", question="Check residual scale.", seed=20261027,
        replicates=8, numerics={"reference_nodes": 8}, world=dict(dimension=3, design="new_excitation",
            mass=1.0, acceleration_snr=0.0, direction="random", readings=5,
            noise={"force_sd": 1.0, "acceleration_sd": 1.0}, supplied_noise="exact"))


def readings(force, acceleration):
    force, acceleration = np.asarray(force, dtype=float), np.asarray(acceleration, dtype=float)
    shape = force.shape[:2]
    return ReadingSet(force, acceleration, np.ones(shape), np.ones(shape))


@pytest.mark.parametrize("k", [2, 3])
def test_q_sum_matches_independent_singular_values(k):
    f = np.array([[[3, 0, 1], [0, 2, 0], [1, 2, 3]]], dtype=float)[..., :k]
    a = np.array([[[1, 2, 0], [1, 0, 2], [-1, -2, -3]]], dtype=float)[..., :k]
    q = [np.linalg.svd(np.column_stack((x, y)), compute_uv=False)[-1] ** 2
         for x, y in zip(f[0], a[0])]
    result = residual_statistics(readings(f, a))
    assert result["q_sum"][0] == pytest.approx(sum(q), abs=1e-12)
    assert result["pooled_variance"][0] == pytest.approx(sum(q) / (3 * (k - 1)))


def test_whole_line_q_is_invariant_to_channel_column_rotation():
    f = np.array([[[3., 1., 0.], [1., 4., -2.]]])
    a = np.array([[[1., -1., 2.], [2., 0., 1.]]])
    angle = 0.73
    rotated_f = np.cos(angle) * f + np.sin(angle) * a
    rotated_a = -np.sin(angle) * f + np.cos(angle) * a
    before = residual_statistics(readings(f, a))
    after = residual_statistics(readings(rotated_f, rotated_a))
    np.testing.assert_allclose(before["q_sum"], after["q_sum"], atol=1e-12)


@pytest.mark.parametrize("k", [2, 3])
def test_pooled_variance_is_homogeneous_in_declared_noise_units(k):
    f = np.array([[[3., 1., 0.], [1., 4., -2.]]])[..., :k]
    a = np.array([[[1., -1., 2.], [2., 0., 1.]]])[..., :k]
    before = residual_statistics(readings(f, a))
    after = residual_statistics(readings(2 * f, 2 * a))
    np.testing.assert_allclose(after["pooled_variance"], 4 * before["pooled_variance"], rtol=1e-13)


@pytest.mark.parametrize("k", [2, 3])
def test_null_mean_and_sampling_spread_match_exact_moments(k):
    """A fixed seeded repeated-world check, with six-standard-error tolerances.

    This verifies the moment formulas; reported studies use the saved recipe.
    """
    from met.pooled_reference import moments

    rng = np.random.default_rng(27102026 + k)
    replicates, n = 2000, 20
    result = residual_statistics(readings(rng.normal(size=(replicates, n, k)),
                                          rng.normal(size=(replicates, n, k))))["pooled_variance"]
    ref = moments(k, 0, nodes=8)
    mean = ref["mean"] / (k - 1)
    sd = np.sqrt(ref["variance"] / n) / (k - 1)
    assert abs(result.mean() - mean) < 6 * sd / np.sqrt(replicates)
    assert abs(result.std(ddof=1) - sd) < 0.08 * sd


def test_recipe_replays_identically_and_keeps_replicate_prefix(tmp_path):
    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    first, record = run(path, tmp_path, progress=None)
    second, _ = run(path, tmp_path, progress=None)
    longer, _ = run(path, tmp_path, replicates=12, progress=None)
    assert first != second
    with np.load(first / "values.npz") as one, np.load(second / "values.npz") as two, np.load(longer / "values.npz") as more:
        assert sorted(one.files) == sorted(two.files)
        for key in one.files:
            np.testing.assert_array_equal(one[key], two[key])
            np.testing.assert_array_equal(one[key], more[key][:8])
    assert (first / "preset.yaml").read_bytes() == path.read_bytes()
    assert yaml.safe_load((longer / "preset.yaml").read_text(encoding="utf-8"))["replicates"] == 12
    assert (longer / "source_preset.yaml").read_bytes() == path.read_bytes()
    assert (first / "summary.md").exists() and (first / "summary.csv").exists()
    restored = json.loads((first / "run.json").read_text(encoding="utf-8"))
    assert restored["workers"] == record["workers"] == 1
    assert restored["code"]["preset_sha256"]
    reference = record["cells"][0]["reference"]
    assert all(np.isfinite(v) for v in reference.values())


@pytest.mark.parametrize("change", [
    lambda d: d.update(replicates=True),
    lambda d: d.update(seed=-1),
    lambda d: d.update(name="../escape"),
    lambda d: d["world"].update(dimension=1),
    lambda d: d["world"].update(acceleration_snr={"uniform": [0, 1]}),
    lambda d: d["world"]["noise"].update(force_sd=2),
    lambda d: d["world"].update(supplied_noise={"force_scale": 2, "acceleration_scale": 2}),
    lambda d: d.update(numerics={"reference_nodes": 128}),
])
def test_reference_scope_is_validated(change):
    data = deepcopy(recipe())
    change(data)
    with pytest.raises(ValueError):
        validate(data)


def test_shipped_recipe_has_explicit_seed_and_ninety_cells():
    from met.preset import load

    path = Path(__file__).resolve().parents[1] / "presets" / "pooled_noise_scale.yaml"
    spec = load(path)
    assert spec["raw"]["seed"] == 20261027
    assert len(spec["cells"]) == 90


def test_existing_runner_rejects_parallel_pooled_recipe(tmp_path):
    from met.runner import run as runner

    path = tmp_path / "recipe.yaml"
    path.write_text(yaml.safe_dump(recipe()), encoding="utf-8")
    for workers in (2, "max", True, 1.0):
        with pytest.raises(ValueError, match="serial"):
            runner(path, tmp_path, workers=workers, progress=None)

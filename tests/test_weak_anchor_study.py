"""Exact scatter-generator moments and reproducible fixed-history study."""
from pathlib import Path
import math

import numpy as np
import pytest
from scipy import special
import yaml

from met.weak_anchor_study import scatter_from_components


def test_bartlett_accelerator_has_the_fixed_gaussian_history_scatter_moments():
    count, energy, direction = 5, 3., np.array([.6, .8])
    nodes, normal_weights = np.polynomial.hermite.hermgauss(3)
    nodes, normal_weights = math.sqrt(2)*nodes, normal_weights/math.sqrt(math.pi)
    first, first_weights = special.roots_genlaguerre(3, (count-1)/2-1)
    second, second_weights = special.roots_genlaguerre(3, (count-2)/2-1)
    first, second = 2*first, 2*second
    first_weights /= special.gamma((count-1)/2)
    second_weights /= special.gamma((count-2)/2)
    expected = count*np.eye(2)+energy*np.outer(direction, direction)
    mean, variance = np.zeros((2, 2)), 0.
    for i, j, k, a, b in np.ndindex(3, 3, 3, 3, 3):
        probability = normal_weights[i]*normal_weights[j]*normal_weights[k]*first_weights[a]*second_weights[b]
        scatter = scatter_from_components(energy, direction, [nodes[i], nodes[j]],
                                          [first[a], second[b]], nodes[k])
        mean += probability*scatter
        variance += probability*np.sum((scatter-expected)**2)
    np.testing.assert_allclose(mean, expected, rtol=2e-14, atol=2e-14)
    assert variance == pytest.approx(6*(count+energy), rel=2e-14)


def test_standard_runner_stores_replayable_scatter_worlds_and_scope(tmp_path):
    from met import preset, runner
    config = yaml.safe_load(Path("presets/weak_anchor_sanity.yaml").read_text(encoding="utf-8"))
    config.update(replicates=2, counts=[100])
    config["cases"] = config["cases"][:1]
    path = tmp_path/"recipe.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    assert preset.load(path)["kind"] == "weak_anchor_sanity"
    first, record = runner.run(path, tmp_path/"first", workers=1, progress=None)
    second, again = runner.run(path, tmp_path/"second", workers=1, progress=None)
    assert record == again
    assert len(record["summary"]) == 2
    assert (first/"summary.md").is_file()
    with np.load(first/"values.npz") as a, np.load(second/"values.npz") as b:
        assert a.files == b.files
        assert a[a.files[0]].shape == (2, 2, 2)
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])
    config["counts"] = [100, 100]
    with pytest.raises(ValueError, match="distinct"):
        preset.validate(config)
    with pytest.raises(ValueError, match="serial"):
        runner.run(path, tmp_path/"invalid", workers=2, progress=None)

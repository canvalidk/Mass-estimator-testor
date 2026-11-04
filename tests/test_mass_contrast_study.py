"""Exact central-scatter moments and standard-runner replay."""
from pathlib import Path

import numpy as np
import pytest
from scipy import special
import yaml

from met.mass_contrast_study import central_physical_scatter
from met.runner import run


def test_central_bartlett_sampler_matches_gaussian_scatter_mean_and_variance():
    count = 7
    covariance = np.array([[2., .4], [.4, .8]])
    chol = np.linalg.cholesky(covariance)
    gaussian_nodes, gaussian_weights = np.polynomial.hermite.hermgauss(3)
    components = [special.roots_genlaguerre(3, df/2-1) for df in (count, count-1)]
    mean, variance = np.zeros((2, 2)), 0.

    class Components:
        def __init__(self, chi, normal):
            self.chi, self.value = iter(chi), normal

        def chisquare(self, _df):
            return next(self.chi)

        def normal(self):
            return self.value

    for i, j, k in np.ndindex(3, 3, 3):
        probability = (components[0][1][i]/special.gamma(count/2)
            *components[1][1][j]/special.gamma((count-1)/2)
            *gaussian_weights[k]/np.sqrt(np.pi))
        rng = Components([2*components[0][0][i], 2*components[1][0][j]],
                         np.sqrt(2)*gaussian_nodes[k])
        scatter = central_physical_scatter(rng, count, covariance)
        white = np.linalg.solve(chol, np.linalg.solve(chol, scatter).T).T
        mean += probability*scatter
        variance += probability*np.sum((white-count*np.eye(2))**2)
    np.testing.assert_allclose(mean, count*covariance, rtol=2e-14, atol=2e-14)
    assert variance == pytest.approx(6*count, rel=3e-14)


def test_mass_contrast_study_replays_actual_both_group_scatters_through_standard_runner(tmp_path):
    source = Path(__file__).resolve().parents[1]/"presets/mass_contrast_sanity.yaml"
    config = yaml.safe_load(source.read_text(encoding="utf-8"))
    config.update(replicates=2, informative_counts=[8], prior_nodes=64, epsrel=1e-8)
    preset = tmp_path/"study.yaml"
    preset.write_text(yaml.safe_dump(config), encoding="utf-8")
    first, metadata = run(preset, tmp_path/"results", workers=1, progress=None)
    second, _ = run(preset, tmp_path/"results", workers=1, progress=None)
    assert (first/"worlds.csv").read_bytes() == (second/"worlds.csv").read_bytes()
    with np.load(first/"values.npz") as original, np.load(second/"values.npz") as replay:
        assert original["scatters_8"].shape == (2, 2, 2, 2)
        np.testing.assert_array_equal(original["scatters_8"], replay["scatters_8"])
    assert len(metadata["summary"]) == 3
    assert "Both covariance groups are actually sampled" in " ".join(metadata["limitations"])
    for name in ("summary.md", "population_scatter.csv", "run.json", "preset.yaml"):
        assert (first/name).is_file()
    with pytest.raises(ValueError, match="serial"):
        run(preset, tmp_path/"invalid", workers=2, progress=None)
    config["informative_counts"] = [8, 8]
    preset.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValueError, match="distinct"):
        run(preset, tmp_path/"invalid", workers=1, progress=None)

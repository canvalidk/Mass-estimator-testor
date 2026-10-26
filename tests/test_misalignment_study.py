"""Independent law checks and the reproducible diagnostic-run contract.

These verify calculations and recording, without declaring how a physically correct
definition must respond to a fixed incompatible pair.
"""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.integrate import quad

from met import preset, runner
from met.__main__ import main
from met.full import _stacked_mean_length, isotropic
from met.misalignment import candidate_mean_length, investigate_cell

PRESET = Path("presets/misalignment_shrinking_noise.yaml")


def _small_preset():
    data = yaml.safe_load(PRESET.read_text(encoding="utf-8"))
    data["cases"] = [data["cases"][0], data["cases"][-1]]
    data["noise"]["scales"] = [0.2]
    data["plots"] = False
    return data


@pytest.mark.parametrize("cosine", [1.0, math.sqrt(3) / 2, 0.0, -1.0])
def test_reported_law_quantiles_against_independent_angle_quadrature(cosine):
    """The half-Cauchy(c=1) is uniform in angle: integrate the explicit residual cost.

    This oracle neither uses met.full's explained-energy path nor its log-grid CDF.
    It checks the study's interpretation of raw/quoted output, not physical adequacy.
    """
    sigma = 0.2
    cell = dict(id="c000", case="test", force=[2.0, 0.0],
                acceleration=[cosine, math.sqrt(max(0, 1 - cosine ** 2))],
                noise_scale=sigma, force_sd=sigma, acceleration_sd=sigma)
    result, _ = investigate_cell(cell, 1.0, dict(profile_log_half_span=12, profile_nodes=101))
    gap = math.sqrt(9 + 16 * cosine ** 2)
    projective_cost = (5 - gap) / (2 * sigma ** 2)
    positive_cost = projective_cost if cosine > 0 else 1 / sigma ** 2
    for label, weight in (("raw", 1.0), ("quoted", min(1.0, 1 / projective_cost) if projective_cost > 0 else 1.0)):
        def integrand(phi):
            sine, cosine_phi = math.sin(phi), math.cos(phi)
            cost = (4 * cosine_phi ** 2 + sine ** 2 - 4 * cosine * sine * cosine_phi) / sigma ** 2
            return math.exp(-0.5 * weight * (cost - positive_cost))

        area = quad(integrand, 0, math.pi / 2, epsabs=1e-12, epsrel=1e-12)[0]
        values = [result[label]["interval"][0], result[label]["median"], result[label]["interval"][1]]
        for mass, probability in zip(values, [0.025, 0.5, 0.975]):
            measured = quad(integrand, 0, math.atan(mass), epsabs=1e-12, epsrel=1e-12)[0] / area
            assert measured == pytest.approx(probability, abs=8e-8)


def test_exported_cloud_size_has_no_computational_cutoff_zeros():
    reading = isotropic([2, 0], [1, 0], 0.05, 0.05)
    masses = np.exp(np.linspace(-18, 18, 35))
    got = candidate_mean_length(reading, masses)
    # The Laplace-transform integral is independent of the Rice/Bessel formula.
    expected = _stacked_mean_length(masses, [reading])
    assert np.all(got > 0)
    np.testing.assert_allclose(got, expected, rtol=2e-12)


def test_diagnostic_runner_keeps_fixed_readings_and_reproduces_all_curves(tmp_path):
    data = _small_preset()
    path = tmp_path / "study.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    out1, first = runner.run(path, tmp_path / "runs", progress=None, workers=1)
    out2, second = runner.run(path, tmp_path / "runs", progress=None)
    assert first["cells"] == second["cells"]
    assert first["deterministic"] and first["seed"] == 0
    assert first["prior"] == data["prior"]
    assert "scores" not in first and "truth" not in first
    assert (out1 / "preset.yaml").read_bytes() == path.read_bytes()
    assert json.loads((out1 / "run.json").read_text(encoding="utf-8"))["cells"] == first["cells"]
    assert (out1 / "summary.csv").is_file() and (out1 / "summary.md").is_file()
    with np.load(out1 / "values.npz") as one, np.load(out2 / "values.npz") as two:
        assert sorted(one.files) == sorted(two.files)
        for key in one.files:
            np.testing.assert_array_equal(one[key], two[key])
        for row in first["cells"]:
            for label in ("raw", "quoted"):
                assert np.trapezoid(one[f"{row['id']}/{label}_density"],
                                    one[f"{row['id']}/{label}_lm"]) == pytest.approx(1, abs=2e-15)
                assert np.all(one[f"{row['id']}/{label}_A"] > 0)


@pytest.mark.parametrize("replicates,workers", [(3, None), (None, 2), (None, "max"), (None, True), (None, 1.0)])
def test_stochastic_overrides_are_rejected_for_fixed_observations(tmp_path, replicates, workers):
    with pytest.raises(ValueError, match="deterministic serial"):
        runner.run(PRESET, tmp_path, replicates=replicates, workers=workers, progress=None)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("change", [
    {"prior": {"family": "flat", "center": 1}},
    {"prior": {"family": "half_cauchy", "center": float("nan")}},
    {"noise": {"force_sd": 1, "acceleration_sd": 1, "scales": [1, 1]}},
    {"noise": {"force_sd": 1, "acceleration_sd": 1, "scales": [1e-200]}},
    {"noise": {"force_sd": 1, "acceleration_sd": 1, "scales": [1e200]}},
    {"numerics": {"profile_nodes": 2}},
    {"cases": []},
    {"cases": [{"name": "bad", "force": [1, float("inf")], "acceleration": [1, 0]}]},
    {"name": "../elsewhere"},
    {"unrecognized": 1},
])
def test_fixed_reading_preset_rejects_invalid_declarations(change):
    data = copy.deepcopy(_small_preset())
    data.update(change)
    with pytest.raises(ValueError):
        preset.validate(data)


def test_cli_describes_diagnostic_cells_without_claiming_simulated_worlds(capsys):
    assert main(["validate", str(PRESET)]) == 0
    assert "32 fixed-reading diagnostic cell(s); no random draws" in capsys.readouterr().out


def test_optional_report_writes_local_figures_and_explains_boundaries(tmp_path):
    pytest.importorskip("matplotlib")
    from met.investigation_report import write_report
    spec = preset.validate(_small_preset())
    row, arrays = investigate_cell(spec["cells"][1], 1.0, spec["numerics"])
    data = dict(name="smoke", question="Inspect a fixed antiparallel pair.",
                seed=0, deterministic=True, prior=spec["raw"]["prior"], cells=[row],
                code_version=runner.code_version())
    assert write_report(tmp_path, data, {f"{row['id']}/{k}": v for k, v in arrays.items()})
    page = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<svg" in page and "m → ∞" in page
    assert "not recovered actual instrument errors" in page
    assert (tmp_path / "overview.png").is_file()

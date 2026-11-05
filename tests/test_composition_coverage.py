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


# B18 adds a separate protocol; every B03 test above remains unchanged.
def _b18_recipe(equation="ours33"):
    data = deepcopy(recipe())
    data.update(name="composition-b18-check", protocol="b18", equation=equation,
                sizes=[0.4], truth_cases=["historical", "c", "4c"])
    if equation in ("ours32", "ours33"):
        data["nu_eff"] = 10.0
    return data


@pytest.mark.parametrize("equation", ["ours30", "ours32", "ours33"])
def test_b18_protocol_is_explicit_and_preserves_declared_geometry(equation):
    from met.composition_coverage import validate

    spec = validate(_b18_recipe(equation))
    assert spec["raw"]["protocol"] == "b18"
    assert spec["raw"]["equation"] == equation
    assert spec["raw"]["truth_cases"] == ["historical", "c", "4c"]
    historical = recipe()
    for key in ("sigma_f", "sigma_a", "c", "dimension", "n_readings", "draws", "m1", "m2"):
        assert spec["raw"][key] == historical[key]
    assert len(spec["cells"]) == 6


@pytest.mark.parametrize("change", [
    {"protocol": "other"}, {"equation": "ours29"}, {"equation": "unknown"},
    {"nu_eff": 0.0}, {"nu_eff": float("inf")}, {"nu_eff": float("nan")},
    {"nu_eff": True}, {"nu_eff": "10"}, {"truth_cases": []},
    {"truth_cases": ["historical", "c", "c"]},
    {"truth_cases": ["historical", "c", "unknown"]}, {"sigma_f": 0.3},
])
def test_b18_refuses_undeclared_or_invalid_scientific_scope(change):
    from met.composition_coverage import validate

    data = _b18_recipe()
    data.update(change)
    with pytest.raises(ValueError):
        validate(data)


def test_b18_student_equations_require_noise_degrees_of_freedom():
    from met.composition_coverage import validate

    for equation in ("ours32", "ours33"):
        data = _b18_recipe(equation)
        del data["nu_eff"]
        with pytest.raises(ValueError):
            validate(data)


def _half_cauchy_sum_cdf(value, first_scale, second_scale):
    """Independent density integral, rather than probability-node composition."""
    import math
    from scipy.integrate import quad

    cdf, error = quad(
        lambda x: (2 * first_scale / (math.pi * (first_scale ** 2 + x ** 2))
                   * 2 / math.pi * math.atan((value - x) / second_scale)),
        0, value, epsabs=2e-13, epsrel=2e-13, limit=200)
    assert error < 1e-10
    return cdf


def _b18_null_law(scale):
    from met.composition_laws import analyse

    return analyse([], scale, equation="ours33", nu_eff=10.0)


@pytest.mark.parametrize("value", [0.2, 2.0, 20.0])
def test_b18_composed_cdf_matches_independent_half_cauchy_density_integral(value):
    import math
    from met.composition_laws import composed_logcdf

    first_scale, second_scale = 0.7, 2.3
    first, second = _b18_null_law(first_scale), _b18_null_law(second_scale)
    logcdf, logsf, logbound = composed_logcdf(
        first, second, math.log(value), "L5", nodes=4096)
    expected = _half_cauchy_sum_cdf(value, first_scale, second_scale)
    assert abs(math.exp(logcdf) - expected) <= math.exp(logbound) + 2e-11
    assert abs(math.exp(logsf) - (1 - expected)) <= math.exp(logbound) + 2e-11
    assert np.logaddexp(logcdf, logsf) == pytest.approx(0, abs=2e-12)


@pytest.mark.parametrize("value", [0.2, 0.7, 3.0])
def test_b18_parallel_composition_obeys_reciprocal_duality(value):
    import math
    from met.composition_laws import composed_logcdf

    scales = (0.7, 2.3)
    first, second = (_b18_null_law(c) for c in scales)
    parallel = composed_logcdf(first, second, math.log(value), "L14", nodes=4096)
    inverse_first, inverse_second = (_b18_null_law(1 / c) for c in scales)
    reciprocal_sum = composed_logcdf(
        inverse_first, inverse_second, -math.log(value), "L5", nodes=4096)
    joint_bound = math.exp(parallel[2]) + math.exp(reciprocal_sum[2]) + 2e-11
    assert abs(math.exp(parallel[0]) - math.exp(reciprocal_sum[1])) <= joint_bound
    assert abs(math.exp(parallel[1]) - math.exp(reciprocal_sum[0])) <= joint_bound
    expected = 1 - _half_cauchy_sum_cdf(1 / value, 1 / scales[0], 1 / scales[1])
    assert abs(math.exp(parallel[0]) - expected) <= math.exp(parallel[2]) + 2e-11


class _LognormalCompositionControl:
    """An exact distribution interface, independent of the fitted-law machinery."""
    eta = 0.25
    sigma = 0.2

    def logcdf(self, logmass, quoted=True):
        import math
        from scipy.special import log_ndtr

        z = np.asarray(logmass) / self.sigma
        if quoted:
            z = math.sqrt(self.eta) * z
        return log_ndtr(z), log_ndtr(-z)

    def logquantile(self, logp, logsf, quoted=True):
        import math
        from scipy.special import ndtri_exp

        logp, logsf = np.asarray(logp), np.asarray(logsf)
        z = np.where(logp < logsf, ndtri_exp(logp), -ndtri_exp(logsf))
        return self.sigma * z / (math.sqrt(self.eta) if quoted else 1.0)


def test_b18_equal_eta_still_does_not_commute_with_composition():
    """T51's continuous lognormal example uses the same eta for both parts."""
    import math
    from scipy.special import log_ndtr, ndtri_exp
    from met.composition_laws import composed_logcdf

    law = _LognormalCompositionControl()
    raw = composed_logcdf(law, law, math.log(2), "L5", quoted=False, nodes=4096)
    before = composed_logcdf(law, law, math.log(2), "L5", quoted=True, nodes=4096)
    after = math.exp(log_ndtr(math.sqrt(law.eta) * ndtri_exp(raw[0])))
    assert abs(math.exp(raw[0]) - 0.4723820114869846) <= math.exp(raw[2]) + 2e-10
    assert abs(math.exp(before[0]) - 0.4475999380925859) <= math.exp(before[2]) + 2e-10
    assert after == pytest.approx(0.4861827208446029, abs=2e-4)
    assert abs(math.exp(before[0]) - after) > 0.03


@pytest.mark.parametrize("quoted", [False, True])
def test_b18_student_law_resolves_log_tail_ranks_beyond_float_probability(quoted):
    import math
    from met.composition_laws import analyse
    from met.full import isotropic

    law = analyse([isotropic([2.0], [1.0], 0.2, 0.1)], 2.0,
                  equation="ours33", nu_eff=10.0)
    low = law.logquantile(-900.0, 0.0, quoted=quoted)
    high = law.logquantile(0.0, -900.0, quoted=quoted)
    farther_low = law.logquantile(-1000.0, 0.0, quoted=quoted)
    farther_high = law.logquantile(0.0, -1000.0, quoted=quoted)
    assert all(math.isfinite(x) for x in (low, high, farther_low, farther_high))
    assert farther_low < low < -700 and farther_high > high > 700
    assert law.logcdf(low, quoted=quoted)[0] == pytest.approx(-900.0, abs=2e-6)
    assert law.logcdf(high, quoted=quoted)[1] == pytest.approx(-900.0, abs=2e-6)
    # At either endpoint the finite-nu factor tends to a positive constant,
    # so the half-Cauchy prior supplies a unit exponential log-mass tail.
    assert low - farther_low == pytest.approx(100.0, abs=2e-6)
    assert farther_high - high == pytest.approx(100.0, abs=2e-6)


def test_b18_midpoint_error_bound_controls_convergence():
    import math
    from met.composition_laws import composed_logcdf

    first, second = _b18_null_law(0.7), _b18_null_law(2.3)
    exact = _half_cauchy_sum_cdf(2.0, 0.7, 2.3)
    coarse = composed_logcdf(first, second, math.log(2), "L5", nodes=128)
    fine = composed_logcdf(first, second, math.log(2), "L5", nodes=2048)
    coarse_error, fine_error = abs(math.exp(coarse[0]) - exact), abs(math.exp(fine[0]) - exact)
    assert coarse_error <= math.exp(coarse[2]) + 2e-11
    assert fine_error <= math.exp(fine[2]) + 2e-11
    assert fine[2] < coarse[2] - math.log(8)
    assert fine_error < coarse_error


def test_b18_new_protocol_replays_exactly_across_worker_counts(tmp_path):
    data = _b18_recipe()
    data["truth_cases"] = ["c"]
    source = tmp_path / "b18.yaml"
    source.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    serial_out, serial = runner.run(source, tmp_path, workers=1, progress=None)
    parallel_out, parallel = runner.run(source, tmp_path, workers=2, progress=None)
    assert serial["summary"] == parallel["summary"]
    assert serial["workers"] == 1 and parallel["workers"] == 2
    assert serial["equation"] == "ours33" and serial["nu_eff"] == 10.0
    assert serial["truth_cases"] == ["c"] and len(serial["cells"]) == 2
    assert (serial_out / "source_preset.yaml").read_bytes() == source.read_bytes()
    assert yaml.safe_load((serial_out / "preset.yaml").read_text(encoding="utf-8")) == data
    assert json.loads((serial_out / "run.json").read_text(encoding="utf-8")) == serial
    procedures = {row["procedure"] for row in serial["summary"]}
    assert len(procedures) == 3
    for row in serial["summary"]:
        assert row["truth_case"] == "c" and row["truth"] == pytest.approx(data["c"])
        expected = {"L5": (data["c"] / 3, 2 * data["c"] / 3),
                    "L14": (1.5 * data["c"], 3 * data["c"])}[row["composition"]]
        assert (row["m1"], row["m2"]) == pytest.approx(expected)
    with np.load(serial_out / "values.npz") as one, np.load(parallel_out / "values.npz") as two:
        assert {"seed", "size", "cell_id", "parts_cover", "eta_parts", "eta_own"} <= set(one.files)
        assert sorted(one.files) == sorted(two.files)
        for key in one.files:
            np.testing.assert_array_equal(one[key], two[key])


@pytest.mark.parametrize("equation", ["ours30", "ours32", "ours33"])
@pytest.mark.parametrize("scenario", ["weak", "strong", "negative-cross"])
def test_b18_isotropic_factory_matches_native_full_reference(monkeypatch, equation, scenario):
    import math
    from met import composition_laws
    from met.full import isotropic

    index = {"weak": 0, "strong": 1, "negative-cross": 2}[scenario]
    rng = np.random.default_rng([1818, index])
    signal = rng.normal(size=(4, 3))
    signal *= (0.02 if scenario == "weak" else 1.2) / np.linalg.norm(signal, axis=1)[:, None]
    mass = -3.0 if scenario == "negative-cross" else 3.0
    force = mass * signal + 0.2 * rng.normal(size=signal.shape)
    acceleration = signal + 0.1 * rng.normal(size=signal.shape)
    readings = [isotropic(f, a, 0.2, 0.1) for f, a in zip(force, acceleration)]
    fast = composition_laws.analyse(readings, 2.0, equation, nu_eff=10.0)
    # The fallback calls the native estimator and its full covariance residual
    # geometry. It shares no isotropic Gram/residual shortcut with the fast path.
    with monkeypatch.context() as patch:
        patch.setattr(composition_laws, "_fast_isotropic_law", lambda *args, **kwargs: None)
        reference = composition_laws.analyse(readings, 2.0, equation, nu_eff=10.0)
    assert fast.eta == pytest.approx(reference.eta, rel=2e-9, abs=2e-11)
    # ours33's raw ranks avoid invoking the native tiny-eta quoted interval.
    quoted = equation != "ours33"
    probabilities = np.array([0.05, 0.5, 0.95])
    points = reference.logquantile(np.log(probabilities), np.log1p(-probabilities), quoted=quoted)
    actual = fast.logcdf(points, quoted=quoted)
    expected = reference.logcdf(points, quoted=quoted)
    for got, want in zip(actual, expected):
        np.testing.assert_allclose(np.exp(got), np.exp(want), rtol=0, atol=2e-6)
    assert math.isfinite(fast.median_log(quoted=quoted))


def test_b18_nontrivial_warp_round_trips_extreme_log_ranks():
    import math
    from met.composition_laws import analyse
    from met.full import isotropic

    # T39/T48's opposed pair has a non-capped eta near 1/99.
    along = np.array([2.0, 1.0]) / math.sqrt(5)
    perpendicular = np.array([1.0, -2.0]) / math.sqrt(5)
    readings = []
    for sign in (-1, 1):
        observed = 1000 * along + sign * 3 * perpendicular
        readings.append(isotropic([observed[0]], [observed[1]], 1.0, 1.0))
    law = analyse(readings, 1.0, "ours33", nu_eff=10.0)
    assert law.eta == pytest.approx(1 / 99 - 19 / (11 * 1000 ** 2), rel=2e-7)
    logp = np.array([-950.0, -1050.0, 0.0, 0.0])
    logsf = np.array([0.0, 0.0, -950.0, -1050.0])
    quantiles = law.logquantile(logp, logsf, quoted=True)
    assert np.all(np.isfinite(quantiles))
    assert quantiles[1] < quantiles[0] < -10000 < 10000 < quantiles[2] < quantiles[3]
    cdf, survival = law.logcdf(quantiles, quoted=True)
    np.testing.assert_allclose(cdf[:2], logp[:2], rtol=0, atol=2e-6)
    np.testing.assert_allclose(survival[2:], logsf[2:], rtol=0, atol=2e-6)
    assert law.median_log(quoted=True) == law.median_log(quoted=False)


def test_b18_post_min_coverage_and_overlap_contain_post_geometric():
    from met.composition_probe import _one, post_eta, validate

    data = _b18_recipe()
    data["truth_cases"] = ["c"]
    spec = validate(data)
    # Fixed seed 1 exercises genuinely different post-warp parameters.
    result = _one((1, spec["cells"][0], spec["raw"]))
    first, second = result["eta_parts"]
    assert 0 < post_eta(first, second, "warp_after_min") < post_eta(
        first, second, "warp_after_geometric") < 1
    minimum = result["procedures"]["warp_after_min"]
    geometric = result["procedures"]["warp_after_geometric"]
    for key in ("computed_covers", "computed_lower", "computed_upper",
                "overlap", "overlap_lower", "overlap_upper"):
        assert not geometric[key] or minimum[key]

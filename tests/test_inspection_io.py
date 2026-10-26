"""Replay, extension seams and report provenance for case inspection."""
import hashlib
import json

import numpy as np
import pytest
import yaml

from met.__main__ import main
from met.inspection.case import Case, load_cases, parse_cases
from met.inspection.checks import DEFAULT_CHECKS
from met.inspection.geometry import profile
from met.inspection.run import inspect_case, run
from met.preset import load
from met.runner import run as run_preset


def example():
    return dict(name="small", question="Check a saved declaration", prior=dict(family="half_cauchy", center=1),
                readings=[dict(force=[2, 0], acceleration=[1, 0], covariance=np.eye(4).tolist())])


def lightweight_adapter(inputs, masses, **kwargs):
    # Deliberately a different, minimal adapter. Geometry-only custom checks
    # need no ours29 internals, native grid or mandatory output class.
    return dict(adapter="test-observer", number_of_readings=len(inputs.readings),
                novel_curve=np.arange(len(masses), dtype=float))


def lightweight_reporter(folder, result):
    (folder / "report.html").write_text("<p>Custom presentation</p>", encoding="utf-8")
    return False


def test_unknown_annotations_and_truth_round_trip_without_fixing_prior_family():
    data = example()
    data["prior"] = dict(family="future-declared-family", center=2, shape={"future": 3})
    data["future_annotation"] = {"note": "retained, not a hidden input"}
    data["readings"][0]["future_reading_note"] = [1, 2]
    data["truth"] = {"mass": 777}
    case = Case.from_dict(data)
    saved = case.to_dict(include_truth=True)
    assert saved["future_annotation"] == data["future_annotation"]
    assert saved["readings"][0]["future_reading_note"] == [1, 2]
    assert saved["prior"] == data["prior"]
    assert saved["truth"] == data["truth"]
    assert "truth" not in case.to_dict() and "truth" not in case.inputs.to_dict()
    data["readings"][0]["force"][0] = 999
    assert case.inputs.readings[0].force[0] == 2


@pytest.mark.parametrize("change", [
    {"force": [float("nan"), 0]},
    {"acceleration": [1]},
    {"covariance": [[1, 0], [0, 1]]},
    {"covariance": np.diag([1, 1, 0, 1]).tolist()},
    {"covariance": [[1, 0.1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]},
    {"subspace": [[1, 1], [0, 0]]},
    {"window": [1, 0]},
])
def test_validation_protects_mathematical_inputs(change):
    data = example()
    data["readings"][0].update(change)
    with pytest.raises(ValueError):
        Case.from_dict(data)


def test_case_names_and_reading_ids_are_unambiguous():
    with pytest.raises(ValueError, match="case names"):
        parse_cases(dict(cases=[example(), example()]))
    data = example()
    data["readings"] *= 2
    data["readings"][0]["id"] = "duplicate"
    with pytest.raises(ValueError, match="reading ids"):
        Case.from_dict(data)


def test_a_diagnostic_and_check_can_be_added_without_changing_case_or_renderer():
    case = Case.from_dict(example())

    def with_extra(observation, masses):
        return dict(profile(observation, masses), a_new_named_curve=np.ones(len(masses)))

    def extra_check(context):
        return dict(id="new-question", title="New question", status="unresolved", source="research",
                    expected="Not established", observed={"new_curve_present":
                        "a_new_named_curve" in context.diagnostics["readings"][0]})

    result = inspect_case(case, masses=np.exp(np.linspace(-2, 2, 9)), adapter=lightweight_adapter,
                          diagnostic=with_extra, checks=[extra_check])
    assert result["checks"][0]["status"] == "unresolved"
    assert result["checks"][0]["observed"]["new_curve_present"]
    assert result["evaluation"]["adapter"] == "test-observer"


def test_selected_case_replays_with_generic_arrays_and_provenance(tmp_path):
    document = dict(kind="inspection", name="saved-examples", seed=0, numerics={"nodes": 9},
                    cases=[example(), dict(example(), name="second")])
    path = tmp_path / "examples.yaml"
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    out, bundle = run(path, tmp_path / "runs", case_names=["second"], check_ids=["projection"],
                      adapter=lightweight_adapter, reporter=lightweight_reporter, progress=None)
    assert bundle["selected_cases"] == ["second"] and bundle["selected_checks"] == ["projection"]
    assert bundle["cases"][0]["checks"] == {"pass": 1}
    assert bundle["provenance"]["preset_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert "geometry.py" in bundle["provenance"]["inspection_file_sha256"]
    assert (out / "preset.yaml").read_bytes() == path.read_bytes()
    saved = json.loads((out / "case_000" / "run.json").read_text(encoding="utf-8"))
    reference = saved["evaluation"]["novel_curve"]["npz"]
    with np.load(out / "case_000" / "arrays.npz", allow_pickle=False) as arrays:
        np.testing.assert_array_equal(arrays[reference], np.arange(9))
        assert any("residual_covariance" in key for key in arrays.files)
    _, replay = load_cases(out / "case_000" / "case.yaml")
    assert replay[0].name == "second"
    replay_out, replay_bundle = run(out / "case_000" / "case.yaml", tmp_path / "replay",
        adapter=lightweight_adapter, reporter=lightweight_reporter, progress=None)
    assert replay_bundle["selected_checks"] == ["projection"]
    with np.load(replay_out / "case_000" / "arrays.npz") as arrays:
        np.testing.assert_array_equal(arrays[reference], np.arange(9))
    assert bundle["provenance"]["deterministic"] is None
    assert bundle["provenance"]["callables"]["adapter"]["name"] == "lightweight_adapter"
    assert (out / "report.html").is_file() and (out / "summary.md").is_file()


@pytest.mark.parametrize("selection", [{"case_names": ["missing"]}, {"check_ids": ["missing"]}])
def test_unknown_selections_fail_before_writing_results(tmp_path, selection):
    path = tmp_path / "case.yaml"
    path.write_text(yaml.safe_dump(example()), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown"):
        run(path, tmp_path / "runs", adapter=lightweight_adapter, reporter=lightweight_reporter, **selection)
    assert not (tmp_path / "runs").exists()


def test_inspection_preset_is_visible_to_the_existing_cli_and_refuses_parallel_overrides(capsys, tmp_path):
    assert main(["validate", "presets/inspection_cases.yaml"]) == 0
    assert "8 saved inspection case(s)" in capsys.readouterr().out
    assert len(load("presets/inspection_cases.yaml")["cases"]) == 8
    for overrides in (dict(workers=2), dict(workers="max"), dict(replicates=10)):
        with pytest.raises(ValueError, match="deterministic serial"):
            run_preset("presets/inspection_cases.yaml", tmp_path, **overrides)
    assert not list(tmp_path.iterdir())


def test_report_preserves_unknown_sections_and_escapes_case_text_without_figures(monkeypatch, tmp_path):
    from met.inspection import report
    monkeypatch.setattr(report, "_figures", lambda folder, result: ([], ["Plotting intentionally unavailable"]))
    result = dict(case=dict(name="<script>unsafe</script>", prior={}),
                  checks=[dict(id="novel", title="A check", kind="algebra", source="independent source",
                      status="fail", expected="one", observed={"value": 2}, assumptions=["stated"], tolerance={})],
                  investigator_truth={"mass": 7}, future_section={"new_quantity": [1, 2, 3]},
                  evaluation={"laws": {"quoted": {"support_edge_warning": True,
                      "support_edge_density_fraction": 0.01, "log_grid_bounds": [-4, 4]}}})
    assert report.write_report(tmp_path, result) is False
    html = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<script>unsafe</script>" not in html and "&lt;script&gt;unsafe" in html
    assert "future_section" in html and "new_quantity" in html
    assert "independent source" in html and "Investigator-only truth" in html
    assert "Plotting intentionally unavailable" in html
    assert "Native-grid support requires review" in html and "not omitted probability" in html


@pytest.mark.parametrize("masses", [[1, 1, 2], [3, 2, 1], [0, 1, 2], [1, float("inf"), 3]])
def test_invalid_inspection_grids_are_rejected(masses):
    with pytest.raises(ValueError, match="increasing finite positive"):
        inspect_case(Case.from_dict(example()), masses=masses, adapter=lightweight_adapter, checks=[])


def test_extension_array_names_cannot_collide_in_saved_npz(tmp_path):
    path = tmp_path / "case.yaml"
    path.write_text(yaml.safe_dump(example()), encoding="utf-8")

    def extra_adapter(inputs, masses, **kwargs):
        return {"a/b": np.ones(3), "a": {"b": np.zeros(3)}, "a~b": np.full(3, 7)}

    out, _ = run(path, tmp_path / "runs", adapter=extra_adapter, reporter=lightweight_reporter,
                 checks=(), progress=None)
    record = json.loads((out / "case_000" / "run.json").read_text(encoding="utf-8"))
    a, b, c = (record["evaluation"]["a/b"]["npz"], record["evaluation"]["a"]["b"]["npz"],
               record["evaluation"]["a~b"]["npz"])
    assert len({a, b, c}) == 3
    with np.load(out / "case_000" / "arrays.npz") as arrays:
        np.testing.assert_array_equal(arrays[a], np.ones(3))
        np.testing.assert_array_equal(arrays[b], np.zeros(3))
        np.testing.assert_array_equal(arrays[c], np.full(3, 7))


def test_custom_prior_is_not_assumed_to_have_half_cauchy_symmetry():
    data = example()
    data["prior"] = dict(family="future", center=1, another_dimensional_parameter=3)
    selected = tuple(c for c in DEFAULT_CHECKS if c.id in ("units", "reciprocal_median"))
    result = inspect_case(Case.from_dict(data), masses=np.exp(np.linspace(-2, 2, 9)),
                          adapter=lightweight_adapter, checks=selected)
    assert all(check["status"] == "inapplicable" for check in result["checks"])


def test_experimental_check_errors_are_visible_without_aborting_other_checks():
    def broken_check(context):
        raise RuntimeError("temporary experimental issue")

    def later_check(context):
        return dict(id="later", status="unresolved")

    result = inspect_case(Case.from_dict(example()), masses=np.exp(np.linspace(-2, 2, 9)),
                          adapter=lightweight_adapter, checks=(broken_check, later_check))
    assert [c["status"] for c in result["checks"]] == ["error", "unresolved"]
    assert "temporary experimental issue" in result["checks"][0]["detail"]


def test_generated_series_adapter_uses_supplied_noise_and_keeps_actual_errors_separate():
    from met.inspection.sources import from_reading_set
    from met.world import generate, series_rngs
    world = dict(dimension=2, design="same_pair", mass=3.0, acceleration_snr=2.0,
        direction="random", readings=3, noise=dict(force_sd={"per_axis": [0.4, 0.7]}, acceleration_sd=0.2),
        supplied_noise=dict(force_scale=2, acceleration_scale=0.5))
    observed, truth = generate(world, series_rngs(12, world, 2))
    case = from_reading_set(observed, truth=truth, prior=dict(family="half_cauchy", center=1.5),
                           name="generated", index=1, metadata={"seed": 12, "world": world})
    assert case.inputs.prior["center"] == 1.5 and case.investigator_truth["mass"] == 3.0
    for reading, actual in zip(case.inputs.readings, case.investigator_truth["readings"]):
        np.testing.assert_allclose(np.diag(reading.covariance), [0.8 ** 2, 1.4 ** 2, 0.1 ** 2, 0.1 ** 2])
        np.testing.assert_allclose(reading.force - 3 * reading.acceleration,
            np.array(actual["force_error"]) - 3 * np.array(actual["acceleration_error"]), atol=2e-16)
    assert "truth" not in case.inputs.to_dict() and case.metadata["seed"] == 12


def test_law_check_detects_a_corrupted_native_density_even_when_sampled_kernel_is_unchanged():
    def inconsistent_adapter(inputs, masses, **kwargs):
        density = np.array([1, 1, 1 if len(inputs.readings) == 1 else 3], dtype=float)
        return dict(laws={"quoted": dict(log_kernel=np.zeros(len(masses)), point=1,
                    log_mass=np.array([-1, 0, 1]), density=density)})

    checks = tuple(c for c in DEFAULT_CHECKS if c.id == "null_law")
    result = inspect_case(Case.from_dict(example()), masses=np.exp(np.linspace(-2, 2, 9)),
                          adapter=inconsistent_adapter, checks=checks)
    finding = result["checks"][0]
    assert finding["status"] == "fail"
    assert finding["observed"]["sampled_log_kernel_error_modulo_constant"] == 0
    assert finding["observed"]["native_cdf_maximum_error"] > 0.1

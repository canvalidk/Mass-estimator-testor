"""Study boundaries, interrupted-run recovery and pilot precision calculations."""

from copy import deepcopy
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch
from uuid import uuid4

from experiment_config import scenario_rows
import study_runner as study


def specification(experiments=2):
    return {
        "schema_version": 1, "name": "Study test", "phase": "exploratory_pilot",
        "ademp": {key: ["Explicit test design"] for key in study.ADEMP_FIELDS},
        "provenance": {"managed_repository": "https://github.com/canvalidk/VD-docs",
                       "managed_commit": "a" * 40, "managed_paths": ["source.md"],
                       "methodology": [{"title": "Methods", "url": "https://doi.org/10.1002/sim.8086"}]},
        "interpretation": ["Fixed-size pilot, not a stopping rule"],
        "fixed_replicates_per_cell": 100,
        "primary_metrics": [
            {"scope": "point", "name": "outside_factor_2", "target_mcse": .01},
            {"scope": "uncertainty", "name": "coverage_95", "target_mcse": .005}],
        "experiments": [{"id": f"experiment_{i}", "config": {
            "name": f"Experiment {i}", "samples": 100, "save_trials": True,
            "training_seed": i * 2, "heldout_seed": i * 2 + 1,
            "scenario": {"type": "pairs", "pairs": [[1, 1]]}}} for i in range(experiments)]}


def report_for(config):
    rows = []
    for row in scenario_rows(config):
        rows.append({**row, "methods": {
            "flat_joint": {"means": {"outside_factor_2": .1}, "mcse": {"outside_factor_2": .03}},
            "norm_ratio": {"means": {"outside_factor_2": .2}, "mcse": {"outside_factor_2": .04}}},
            "paired_flat_minus_other": {"norm_ratio": {"outside_factor_2": {"mean_difference": -.1, "mcse": .02}}},
            "uncertainty": {name: {"intervals": {"95": {"coverage": {"mean": .95, "empirical_mcse": .02}}}}
                            for name in ("flat_joint", "tube_joint")}})
    return {"samples_per_scenario": config["samples"], "experiment_config": config, "rows": rows}


def execute_fixture(config, output, progress=None):
    output.mkdir(parents=True)
    report = report_for(config)
    (output / "results.json").write_text(json.dumps(report), encoding="utf-8")
    (output / "trials.npz").write_bytes(b"a saved trial archive")
    return report


class StudySpecificationTests(unittest.TestCase):
    def test_fingerprint_locks_execution_but_not_analysis_dependencies(self):
        baseline = study.source_fingerprint()
        sources = baseline["source_sha256"]
        for name in ("study_runner.py", "experiment_config.py", "compare_estimators.py",
                     "comparison_scores.py", "trial_archive.py", "../estimator.py"):
            self.assertIn(name, sources)
        for name in ("assessment_summary.py", "null_truth_sweep.py", "results_catalog.py"):
            self.assertNotIn(name, sources)
        original = study._hash
        with patch.object(study, "_hash", side_effect=lambda p:
                          "changed" if Path(p).name == "comparison_scores.py" else original(p)):
            self.assertNotEqual(study.source_fingerprint(), baseline)

    def test_main_assessment_phase_preserves_same_validation_rules(self):
        spec = specification(1)
        spec["phase"] = "main_assessment"
        self.assertEqual(study.normalize_study(spec)["phase"], "main_assessment")
        spec["experiments"][0]["config"]["samples"] -= 1
        with self.assertRaisesRegex(ValueError, "fixed_replicates_per_cell"):
            study.normalize_study(spec)
        for phase in ("confirmatory_winner", "", None):
            spec = specification(1)
            spec["phase"] = phase
            with self.assertRaisesRegex(ValueError, "phase"):
                study.normalize_study(spec)

    def test_normalization_is_immutable_and_idempotent(self):
        original = specification()
        before = deepcopy(original)
        normalized = study.normalize_study(original)
        self.assertEqual(original, before)
        self.assertEqual(study.normalize_study(normalized), normalized)
        self.assertIn("order", normalized["experiments"][0]["config"])

    def test_invalid_scientific_and_output_boundaries(self):
        cases = []
        for key, value in (("samples", 99), ("save_trials", False)):
            item = specification()
            item["experiments"][0]["config"][key] = value
            cases.append(item)
        for identifier in ("../escape", "CON", "experiment_1"):
            item = specification()
            item["experiments"][0]["id"] = identifier
            cases.append(item)
        item = specification()
        item["experiments"][1]["config"]["training_seed"] = 0
        cases.append(item)
        for target in (0, -1, float("nan"), True):
            item = specification()
            item["primary_metrics"][0]["target_mcse"] = target
            cases.append(item)
        item = specification()
        item["primary_metrics"][0]["name"] = "unbounded_unknown_loss"
        cases.append(item)
        item = specification()
        item["ademp"]["aims"] = []
        cases.append(item)
        item = specification()
        item["arbitrary_extra_field"] = 1
        cases.append(item)
        for item in cases:
            with self.subTest(item=item), self.assertRaises(ValueError):
                study.normalize_study(item)

    def test_production_pilot_has_fixed_truth_and_independent_experiments(self):
        path = Path(study.__file__).parent / "studies" / "mass_excitation_pilot.json"
        spec = study.normalize_study(json.loads(path.read_text(encoding="utf-8")))
        self.assertEqual(spec["fixed_replicates_per_cell"], 512)
        self.assertEqual([e["config"]["measurement"]["repeats"] for e in spec["experiments"]], [1, 4])
        for experiment in spec["experiments"]:
            rows = scenario_rows(experiment["config"])
            self.assertEqual(len(rows), 12)
            null_rows = [r for r in rows if r["true_acceleration_snr"] == 0]
            self.assertEqual([r["true_mass"] for r in null_rows], [.25, 1., 4.])

    def test_main_and_calibration_protocols_have_fresh_seeds_and_fixed_precision(self):
        directory = Path(study.__file__).parent / "studies"
        protocols = [study.normalize_study(json.loads((directory / name).read_text()))
                     for name in ("mass_excitation_pilot.json", "mass_excitation_main.json",
                                  "calibration_sensitivity_pilot.json", "calibration_integration_pilot.json")]
        seen = set()
        for spec in protocols:
            for experiment in spec["experiments"]:
                config = experiment["config"]
                seeds = [config["training_seed"], config["heldout_seed"]]
                calibration = config["measurement"]["calibration"]
                if calibration["mode"] in ("estimated", "pooled_isotropic"):
                    seeds.append(calibration["seed"])
                self.assertFalse(seen.intersection(seeds))
                seen.update(seeds)
                self.assertEqual(config["samples"], spec["fixed_replicates_per_cell"])
                self.assertEqual(len(scenario_rows(config)), 12)
        main = protocols[1]
        self.assertEqual(main["phase"], "main_assessment")
        self.assertEqual(main["fixed_replicates_per_cell"], 10000)
        for metric in main["primary_metrics"]:
            width = study.METRIC_RANGES[(metric["scope"], metric["name"])]
            if metric["scope"] == "point":
                width *= 2  # Also cover a paired difference, whose range doubles.
            self.assertLessEqual(width/(2*math.sqrt(10000)), metric["target_mcse"])
        robust = protocols[2]
        self.assertEqual(robust["phase"], "exploratory_pilot")
        self.assertEqual(robust["fixed_replicates_per_cell"], 1024)
        self.assertEqual(len(robust["experiments"]), 5)


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        # Keep test artifacts on the authorized scratch surface, as the other
        # integration tests do. They can also help diagnose a failed check.
        scratch = Path(study.__file__).resolve().parents[3] / ".tools" / "study_runner_tests"
        self.output = scratch / uuid4().hex / "study"
        self.fingerprint = patch.object(study, "source_fingerprint", return_value={"source_sha256": {"engine.py": "unchanged"}})
        self.fingerprint.start()
        self.addCleanup(self.fingerprint.stop)

    def checkpoint(self):
        return json.loads((self.output / "study_checkpoint.json").read_text(encoding="utf-8"))

    def test_completed_resume_is_read_only_and_new_run_never_overwrites(self):
        with patch.object(study, "run_config", side_effect=execute_fixture) as execute:
            result = study.run_study(specification(), self.output)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(execute.call_count, 2)
            before = {p: p.read_bytes() for p in self.output.rglob("*") if p.is_file()}
            execute.reset_mock()
            self.assertEqual(study.run_study(specification(), self.output, resume=True), result)
            execute.assert_not_called()
            self.assertEqual({p: p.read_bytes() for p in before}, before)
            with self.assertRaisesRegex(ValueError, "already exists"):
                study.run_study(specification(), self.output)

    def test_partial_failure_retains_attempt_and_resumes_only_incomplete_experiment(self):
        def fail_second(config, output, progress=None):
            if config["training_seed"] == 2:
                output.mkdir(parents=True)
                (output / "partial.txt").write_text("keep me", encoding="utf-8")
                raise RuntimeError("Numerical refinement failed")
            return execute_fixture(config, output, progress)
        with patch.object(study, "run_config", side_effect=fail_second), self.assertRaises(RuntimeError):
            study.run_study(specification(), self.output)
        state = self.checkpoint()
        self.assertEqual(state["status"], "failed")
        self.assertEqual([r["status"] for r in state["experiments"]], ["completed", "failed"])
        with patch.object(study, "run_config", side_effect=execute_fixture) as execute:
            resumed = study.run_study(specification(), self.output, resume=True)
            self.assertEqual(execute.call_count, 1)
        record = resumed["experiments"][1]
        self.assertEqual([a["status"] for a in record["attempts"]], ["failed", "completed"])
        self.assertTrue(record["attempt_directory"].endswith("attempt-002"))
        self.assertEqual((self.output / record["attempts"][0]["directory"] / "partial.txt").read_text(), "keep me")

    def test_keyboard_interrupt_is_never_reported_as_complete(self):
        with patch.object(study, "run_config", side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
            study.run_study(specification(1), self.output)
        self.assertEqual(self.checkpoint()["status"], "interrupted")
        with patch.object(study, "run_config", side_effect=execute_fixture):
            result = study.run_study(specification(1), self.output, resume=True)
        self.assertEqual([a["status"] for a in result["experiments"][0]["attempts"]], ["interrupted", "completed"])

    def test_resume_rejects_changed_spec_source_and_any_saved_trial_bytes(self):
        with patch.object(study, "run_config", side_effect=execute_fixture):
            state = study.run_study(specification(1), self.output)
        different = specification(1)
        different["ademp"]["aims"] = ["Different research question"]
        with self.assertRaisesRegex(ValueError, "identical normalized"):
            study.run_study(different, self.output, resume=True)
        with patch.object(study, "source_fingerprint", return_value={"source_sha256": {"engine.py": "changed"}}):
            with self.assertRaisesRegex(ValueError, "identical source"):
                study.run_study(specification(1), self.output, resume=True)
        archive = self.output / state["experiments"][0]["attempt_directory"] / "trials.npz"
        archive.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "output changed or missing"):
            study.run_study(specification(1), self.output, resume=True)

    def test_incomplete_report_does_not_create_completed_checkpoint(self):
        def incomplete(config, output, progress=None):
            report = execute_fixture(config, output, progress)
            report["rows"] = []
            return report
        with patch.object(study, "run_config", side_effect=incomplete), self.assertRaisesRegex(ValueError, "incomplete"):
            study.run_study(specification(1), self.output)
        self.assertEqual(self.checkpoint()["status"], "failed")
        self.assertNotIn("output_sha256", self.checkpoint()["experiments"][0])


class PrecisionPlanningTests(unittest.TestCase):
    def test_additional_law_and_paired_coverage_keep_their_own_mcse(self):
        spec = study.normalize_study(specification(1))
        report = report_for(spec["experiments"][0]["config"])
        row = report["rows"][0]
        row["uncertainty"]["calibrated_joint"] = {"intervals": {"95": {"coverage": {
            "mean": .94, "empirical_mcse": .011}}}}
        row["paired_uncertainty_flat_minus_other"] = {"calibrated_joint": {"coverage_95": {
            "mean_difference": -.03, "mcse": .007}}}
        entries = study.precision_summary(report, spec["primary_metrics"])["entries"]
        law = next(e for e in entries if e["scope"] == "uncertainty" and e["method"] == "calibrated_joint")
        pair = next(e for e in entries if e["scope"] == "paired_uncertainty")
        self.assertEqual((law["mean"], law["mcse"]), (.94, .011))
        self.assertEqual((pair["mean"], pair["mcse"]), (-.03, .007))

    def test_reads_real_paired_mcse_and_flags_capacity_without_changing_R(self):
        spec = study.normalize_study(specification(1))
        config = spec["experiments"][0]["config"]
        config["measurement"] = {"repeats": 4, "calibration": {"mode": "known"}}
        report = report_for(config)
        summary = study.precision_summary(report, spec["primary_metrics"])
        entries = summary["entries"]
        paired = next(e for e in entries if e["scope"] == "paired_point")
        self.assertEqual(paired["mcse"], .02)
        self.assertEqual(paired["pilot_variance_suggested_replicates"], 400)
        coverage = next(e for e in entries if e["scope"] == "uncertainty")
        self.assertEqual(coverage["pilot_variance_suggested_replicates"], 1600)
        self.assertEqual(coverage["conservative_bounded_variance_replicates"], 10000)
        self.assertEqual(coverage["conservative_bounded_variance_replicates_within_current_capacity"],
                         10000 <= summary["current_maximum_replicates_per_cell"])
        self.assertEqual(report["samples_per_scenario"], 100)

    def test_missing_and_zero_mcse_do_not_become_zero_replicate_recommendations(self):
        spec = study.normalize_study(specification(1))
        report = report_for(spec["experiments"][0]["config"])
        row = report["rows"][0]
        row["methods"]["flat_joint"]["mcse"].clear()
        row["uncertainty"]["flat_joint"]["intervals"]["95"]["coverage"]["empirical_mcse"] = 0.
        entries = study.precision_summary(report, spec["primary_metrics"])["entries"]
        missing = next(e for e in entries if e["scope"] == "point" and e["method"] == "flat_joint")
        self.assertEqual(missing["status"], "mcse_unavailable")
        self.assertIsNone(missing["target_met"])
        zero = next(e for e in entries if e["scope"] == "uncertainty" and e["method"] == "flat_joint")
        self.assertIsNone(zero["pilot_variance_suggested_replicates"])
        self.assertEqual(zero["planning_status"], "zero_pilot_variance_use_conservative_bound")
        self.assertEqual(zero["conservative_bounded_variance_replicates"], 10000)


if __name__ == "__main__":
    unittest.main()

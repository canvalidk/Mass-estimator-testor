"""Integrity, interpretation and non-destructive completed-study analysis."""

from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest
from uuid import uuid4

from assessment_summary import (
    export_accuracy, export_assessment, load_completed_study, summarize_accuracy, summarize_study,
)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AssessmentSummaryTests(unittest.TestCase):
    def setUp(self):
        self.base = (Path(__file__).resolve().parents[3] / ".tools" /
                     "assessment_summary_tests" / uuid4().hex)
        self.study = self.base / "study"
        self.attempt = self.study / "experiments" / "known" / "attempt-001"
        self.attempt.mkdir(parents=True)
        config = {"name": "Fixture", "samples": 4, "noise_model": "diagonal_gaussian_3d",
                  "noise": {"force_sd": [1, 2, 3], "acceleration_sd": [3, 2, 1]},
                  "measurement": {"repeats": 2, "calibration": {"mode": "estimated", "samples": 20}}}
        score = {"n": 4, "invalid_count": 0, "valid_count": 4,
                 "means": {"capped_squared_log_factor_2": .1, "outside_factor_2": .25},
                 "mcse": {"capped_squared_log_factor_2": .01, "outside_factor_2": .25},
                 "mean_status": {"capped_squared_log_factor_2": "finite_sample_mean", "outside_factor_2": "finite_sample_mean"},
                 "tolerance_curve": [0., .5, 1.],
                 "log_diagnostics": {"signed_log_bias": -.1, "log_variance": .12, "variance_ddof": 0, "status": "finite_sample_diagnostics"},
                 "additional_tasks": {"prediction": {"mean": None, "empirical_mcse": None, "status": "not_applicable"}}}
        comparator = deepcopy(score)
        comparator.update(invalid_count=1, valid_count=3, population_unbounded_risk="infinite_under_declared_failure_policy",
                          population_risk_note="Preserve this established risk warning.")
        comparator["log_diagnostics"].update(signed_log_bias=None, log_variance=None, status="invalid_estimates")
        pair = {"mean_difference": -.05, "mcse": .02, "status": "finite_sample_comparison", "n": 4,
                "difference_definition": "A_minus_B", "mcse_interpretation": "Paired source warning"}
        interval = {"coverage": {"mean": .75, "empirical_mcse": .25, "status": "finite_sample"},
                    "below_interval": {"mean": .25, "empirical_mcse": .25, "status": "finite_sample"},
                    "above_interval": {"mean": 0., "empirical_mcse": 0., "status": "finite_sample"},
                    "mean_log_width": {"mean": 2.3, "empirical_mcse": .1, "status": "finite_sample"}}
        row = {"id": "00_null", "true_mass": 1., "true_force_snr": 0., "true_acceleration_snr": 0.,
               "metric_applicability": {"prediction": {"status": "not_applicable", "reason": "zero_true_force"}},
               "methods": {"flat_joint": score, "flat_median": comparator},
               "paired_flat_minus_other": {"flat_median": {name: pair for name in score["means"]}},
               "uncertainty": {"flat_joint": {"intervals": {"50": interval, "95": interval}}}}
        self.report = {"samples_per_scenario": 4, "experiment_config": config,
                       "assumptions": "Fixture source assumption", "failure_policy": "Fixture source failure policy",
                       "rows": [row]}
        self.spec = {"schema_version": 1, "name": "Fixture study", "phase": "exploratory_pilot",
                     "fixed_replicates_per_cell": 4, "interpretation": ["Do not pool cells"],
                     "experiments": [{"id": "known", "config": config}]}
        write_json(self.attempt / "config.json", config)
        write_json(self.attempt / "results.json", self.report)
        (self.attempt / "extra.bin").write_bytes(b"recorded extra artifact")
        self.checkpoint = {"schema_version": 1, "status": "completed", "specification": self.spec,
                           "spec_sha256": hashlib.sha256(json.dumps(self.spec, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                           "experiments": [{"id": "known", "status": "completed", "attempt_directory": "experiments/known/attempt-001",
                                            "output_sha256": {name: digest(self.attempt / name) for name in ("results.json", "config.json", "extra.bin")}}]}
        write_json(self.study / "study.json", self.spec)
        write_json(self.study / "study_checkpoint.json", self.checkpoint)

    def save_checkpoint(self):
        write_json(self.study / "study_checkpoint.json", self.checkpoint)

    def test_verified_summary_preserves_numeric_values_statuses_and_risk(self):
        loaded = load_completed_study(self.study)
        summary = summarize_study(loaded, ("flat_median", "norm_ratio"))
        experiment = summary["experiments"][0]
        row = experiment["rows"][0]
        self.assertEqual(experiment["verified_artifacts"], 3)
        self.assertEqual(experiment["config"], self.spec["experiments"][0]["config"])
        self.assertEqual(row["methods"]["flat_median"]["invalid_count"], 1)
        self.assertEqual(row["methods"]["flat_median"]["population_risk_note"], "Preserve this established risk warning.")
        self.assertEqual(row["methods"]["flat_joint"]["additional_tasks"]["prediction"]["status"], "not_applicable")
        self.assertEqual(row["scenario"]["metric_applicability"], self.report["rows"][0]["metric_applicability"])
        self.assertEqual(row["laws"]["flat_joint"]["intervals"], self.report["rows"][0]["uncertainty"]["flat_joint"]["intervals"])
        self.assertEqual(row["paired_flat_minus_other"]["norm_ratio"]["status"], "unavailable")
        pair = row["paired_flat_minus_other"]["flat_median"]["metrics"]["outside_factor_2"]
        self.assertEqual(pair["mcse"], .02)
        self.assertEqual(pair["mcse_interpretation"], "Paired source warning")
        self.assertAlmostEqual(row["flat_signed_log_bias"]["mcse"], .2)
        self.assertEqual(row["flat_signed_log_bias"]["mean"], -.1)
        self.assertNotIn("tolerance_curve", row["methods"]["flat_joint"])
        self.assertIn("tolerance_curve", self.report["rows"][0]["methods"]["flat_joint"])
        self.assertTrue(any("not confidence intervals" in note for note in summary["notes"]))

    def test_export_json_changes_no_source_and_refuses_overwrite(self):
        originals = {path: digest(path) for path in self.study.rglob("*") if path.is_file()}
        output = self.base / "analysis"
        result = export_assessment(self.study, output, figures=False)
        self.assertEqual(result, json.loads((output / "summary.json").read_text()))
        self.assertEqual(result["figures"], [])
        self.assertEqual(originals, {path: digest(path) for path in originals})
        with self.assertRaises(FileExistsError):
            export_assessment(self.study, output, figures=False)
        with self.assertRaisesRegex(ValueError, "outside"):
            export_assessment(self.study, self.study / "analysis", figures=False)

    def test_hash_verification_checks_every_recorded_artifact(self):
        (self.attempt / "extra.bin").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            export_assessment(self.study, self.base / "analysis", figures=False)
        self.assertFalse((self.base / "analysis").exists())

    def test_refuses_partial_study_and_unsafe_checkpoint_paths(self):
        self.checkpoint["status"] = "running"
        self.save_checkpoint()
        with self.assertRaisesRegex(ValueError, "completed"):
            load_completed_study(self.study)
        self.checkpoint["status"] = "completed"
        self.checkpoint["experiments"][0]["attempt_directory"] = "../elsewhere"
        self.save_checkpoint()
        with self.assertRaisesRegex(ValueError, "escapes"):
            load_completed_study(self.study)

    def test_rejects_specification_tampering(self):
        self.spec["name"] = "Changed study"
        write_json(self.study / "study.json", self.spec)
        with self.assertRaisesRegex(ValueError, "differs"):
            load_completed_study(self.study)
        self.checkpoint["specification"] = self.spec
        self.save_checkpoint()
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            load_completed_study(self.study)

    def test_rejects_wrong_trial_denominator_after_valid_file_hash(self):
        self.report["rows"][0]["methods"]["flat_joint"]["n"] = 3
        write_json(self.attempt / "results.json", self.report)
        self.checkpoint["experiments"][0]["output_sha256"]["results.json"] = digest(self.attempt / "results.json")
        self.save_checkpoint()
        with self.assertRaisesRegex(ValueError, "denominator"):
            load_completed_study(self.study)

    def test_missing_bias_remains_unavailable(self):
        loaded = load_completed_study(self.study)
        loaded["reports"][0]["report"]["rows"][0]["methods"]["flat_joint"]["log_diagnostics"].update(
            signed_log_bias=None, log_variance=None, status="invalid_estimates")
        result = summarize_study(loaded)
        bias = result["experiments"][0]["rows"][0]["flat_signed_log_bias"]
        self.assertIsNone(bias["mean"])
        self.assertIsNone(bias["mcse"])
        self.assertEqual(bias["status"], "invalid_estimates")

    def test_additional_law_and_its_pair_survive_summary(self):
        loaded = load_completed_study(self.study)
        source = loaded["reports"][0]["report"]["rows"][0]
        source["uncertainty"]["calibrated_joint"] = deepcopy(source["uncertainty"]["flat_joint"])
        source["paired_uncertainty_flat_minus_other"] = {"calibrated_joint": {"coverage_95": {
            "mean_difference": -.05, "mcse": .01, "status": "finite_sample_comparison"}}}
        row = summarize_study(loaded)["experiments"][0]["rows"][0]
        self.assertEqual(row["laws"]["calibrated_joint"]["intervals"], source["uncertainty"]["calibrated_joint"]["intervals"])
        self.assertEqual(row["paired_uncertainty_flat_minus_other"], source["paired_uncertainty_flat_minus_other"])
        self.assertNotIn("paired_uncertainty_flat_minus_other", row["scenario"])

    def test_accuracy_family_counts_every_archived_method_once_per_cell(self):
        loaded = load_completed_study(self.study)
        row = loaded["reports"][0]["report"]["rows"][0]
        other = deepcopy(row)
        other["id"] = "another_cell"
        other["methods"]["missing_metrics"] = {"n": 4, "means": {}}
        loaded["reports"][0]["report"]["rows"].append(other)
        summary = summarize_accuracy(loaded)
        self.assertEqual(summary["simultaneous_band"]["family_size"], 5)
        self.assertAlmostEqual(summary["experiments"][0]["rows"][0]["dkw_epsilon"], math.sqrt(math.log(200) / 8))
        single = summarize_accuracy(loaded, tolerance_factors=(2.,))
        self.assertEqual(single["simultaneous_band"], summary["simultaneous_band"])

    def test_accuracy_zero_failures_and_target_classification(self):
        loaded = load_completed_study(self.study)
        report = loaded["reports"][0]["report"]
        report["samples_per_scenario"] = 10000
        row = report["rows"][0]
        for summary in row["methods"].values():
            summary.update(n=10000, valid_count=10000, invalid_count=0)
        row["methods"]["flat_joint"]["means"].update(
            {"outside_factor_1.25": .5, "outside_factor_1.5": .1, "outside_factor_2": 0.})
        row["methods"]["flat_joint"]["mcse"]["outside_factor_2"] = 0.
        row["paired_flat_minus_other"]["flat_median"]["outside_factor_2"].update(n=10000, mean_difference=-.25)
        summary = summarize_accuracy(loaded)
        result = summary["experiments"][0]["rows"][0]
        accuracy = result["methods"]["flat_joint"]["accuracy"]
        self.assertEqual(accuracy["1.25"]["targets"], {"0.9": "below_target", "0.95": "below_target"})
        self.assertEqual(accuracy["1.5"]["targets"]["0.9"], "unresolved")
        self.assertEqual(accuracy["2"]["targets"], {"0.9": "supported", "0.95": "supported"})
        self.assertEqual(accuracy["2"]["empirical_mcse"], 0.)
        bounds = accuracy["2"]["simultaneous_success_interval"]
        self.assertEqual(bounds["upper"], 1.)
        self.assertLess(bounds["lower"], 1.)
        self.assertAlmostEqual(bounds["lower"], 1 - result["dkw_epsilon"])
        pair = result["paired_flat_minus_other"]["flat_median"]["2"]
        self.assertEqual(pair["mcse"], .02)
        self.assertEqual(pair["mcse_interpretation"], "Paired source warning")
        self.assertAlmostEqual(pair["simultaneous_failure_difference_interval"]["upper"], -.25 + 2 * result["dkw_epsilon"])
        self.assertEqual(pair["source_status"], "finite_sample_comparison")
        # Equality at either limit uses the declared asymmetric classification:
        # reaching a target supports it; merely reaching it from below does not.
        lower, upper = accuracy["1.5"]["simultaneous_success_interval"].values()
        boundary = summarize_accuracy(loaded, success_targets=(lower, upper))
        classifications = boundary["experiments"][0]["rows"][0]["methods"]["flat_joint"]["accuracy"]["1.5"]["targets"]
        self.assertEqual(classifications[str(lower)], "supported")
        self.assertEqual(classifications[str(upper)], "unresolved")

    def test_accuracy_missing_metrics_and_statuses_are_not_successes(self):
        loaded = load_completed_study(self.study)
        source = loaded["reports"][0]["report"]["rows"][0]["methods"]["flat_joint"]
        source["means"]["outside_factor_2"] = None
        source["mean_status"]["outside_factor_2"] = "not_applicable"
        result = summarize_accuracy(loaded)["experiments"][0]["rows"][0]
        missing = result["methods"]["flat_joint"]["accuracy"]["1.25"]
        self.assertEqual(missing["status"], "metric_unavailable")
        explicit = result["methods"]["flat_joint"]["accuracy"]["2"]
        self.assertEqual(explicit["source_status"], "not_applicable")
        self.assertIsNone(explicit["success_probability"])
        self.assertIsNone(explicit["simultaneous_success_interval"])
        self.assertEqual(explicit["targets"]["0.9"], "unavailable")
        self.assertEqual(result["paired_flat_minus_other"]["flat_median"]["2"]["status"], "marginal_metric_unavailable")
        self.assertEqual(result["methods"]["flat_median"]["invalid_count"], 1)
        self.assertEqual(result["methods"]["flat_median"]["accuracy"]["2"]["success_probability"], .75)

    def test_accuracy_validates_arguments_and_probabilities(self):
        loaded = load_completed_study(self.study)
        for arguments in ({"tolerance_factors": ()}, {"tolerance_factors": (2., 2.)},
                          {"tolerance_factors": (2.5,)}, {"tolerance_factors": (True,)},
                          {"success_targets": (.9, .9)}, {"success_targets": (1.,)},
                          {"success_targets": (float("nan"),)}, {"alpha": 0.}, {"alpha": True}):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                summarize_accuracy(loaded, **arguments)
        loaded["reports"][0]["report"]["rows"][0]["methods"]["flat_joint"]["means"]["outside_factor_2"] = 1.01
        with self.assertRaisesRegex(ValueError, "outside"):
            summarize_accuracy(loaded)

    def test_accuracy_export_preserves_sources_and_refuses_overwrite(self):
        originals = {path: digest(path) for path in self.study.rglob("*") if path.is_file()}
        output = self.base / "accuracy"
        result = export_accuracy(self.study, output)
        self.assertEqual(result, json.loads((output / "accuracy.json").read_text()))
        self.assertEqual(result["source_checkpoint_sha256"], digest(self.study / "study_checkpoint.json"))
        self.assertEqual(result["experiments"][0]["source_results_sha256"], digest(self.attempt / "results.json"))
        self.assertEqual(originals, {path: digest(path) for path in originals})
        with self.assertRaises(FileExistsError):
            export_accuracy(self.study, output)
        with self.assertRaisesRegex(ValueError, "outside"):
            export_accuracy(self.study, self.study / "accuracy")

    def test_accuracy_cli_exports_without_plotting(self):
        script = Path(__file__).with_name("assessment_summary.py")
        output = self.base / "cli-accuracy"
        result = subprocess.run([sys.executable, "-B", str(script), "--study", str(self.study),
                                 "--output", str(output), "--accuracy-map", "--tolerance-factors", "2",
                                 "--success-targets", ".9", "--alpha", ".1"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        saved = json.loads((output / "accuracy.json").read_text())
        self.assertEqual(saved["tolerance_factors"], [2.])
        self.assertEqual(saved["simultaneous_band"]["alpha"], .1)
        self.assertFalse((output / "summary.json").exists())

    def test_cli_json_export_has_no_engine_import(self):
        script = Path(__file__).with_name("assessment_summary.py")
        result = subprocess.run([sys.executable, "-B", str(script), "--study", str(self.study),
                                 "--output", str(self.base / "cli-analysis"), "--no-figures"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        code = "import sys, assessment_summary; assert 'study_runner' not in sys.modules; assert 'comparison_posterior' not in sys.modules"
        imports = subprocess.run([sys.executable, "-B", "-c", code], cwd=script.parent, capture_output=True, text=True)
        self.assertEqual(imports.returncode, 0, imports.stderr)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is an optional plotting dependency")
    def test_png_svg_export_includes_missing_law_and_method(self):
        output = self.base / "figures"
        result = export_assessment(self.study, output, comparators=("flat_median", "norm_ratio"))
        self.assertEqual(len(result["figures"]), 6)
        for name in result["figures"]:
            self.assertGreater((output / name).stat().st_size, 1000)
        svg = (output / "known-flat-assessment.svg").read_text()
        self.assertIn("Per-reading acceleration SNR", svg)
        self.assertIn("not confidence intervals", svg)


if __name__ == "__main__":
    unittest.main()

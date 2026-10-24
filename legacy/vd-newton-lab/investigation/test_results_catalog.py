"""Check saved-result compatibility and statistical labels used by the UI."""

import hashlib
import json
from pathlib import Path
import unittest
import uuid

from results_catalog import list_runs, load_run


WORKSPACE = Path(__file__).resolve().parents[3]


def comparison_data():
    point = {
        "means": {"squared_log": 0.4, "squared_mass": 2., "outside_factor_2": 0.2},
        "mcse": {"squared_log": 0.01, "squared_mass": 0.4, "outside_factor_2": 0.03},
        "mean_status": {"squared_log": "finite_sample_mean"},
        "additional_tasks": {
            "heldout_known_acceleration_relative_squared_force_error": {
                "mean": 3., "empirical_mcse": 0.5, "status": "finite_sample"}},
        "log_diagnostics": {"signed_log_bias": 0.1, "log_variance": 0.39},
        "invalid_count": 0, "finite_positive_fraction": 1.,
        "tolerance_curve": [0., 0.8],
    }
    failed = dict(point, population_unbounded_risk="infinite_under_declared_failure_policy")
    missing = dict(point, means={"squared_log": None, "outside_factor_2": 0.3},
                   mcse={"squared_log": None}, mean_status={"squared_log": "invalid_estimates"},
                   invalid_count=2, finite_positive_fraction=0.98)
    pair = {"mean_difference": -0.1, "mcse": 0.01, "status": "finite_sample_comparison"}
    sample = {"mean": 0.9, "empirical_mcse": 0.02, "status": "finite_sample"}
    return {
        "samples_per_scenario": 100, "training_seed": 10, "heldout_seed": 11,
        "quadrature_order": 96, "batch_size": 32, "scenario_selection": "original",
        "tolerance_factors": [1., 2.], "source_sha256": {"engine.py": "abc123"},
        "baseline_metadata": {"norm_ratio": {"label": "Recorded norm ratio", "formula": "NF / NA"}},
        "rows": [{
            "id": "00_original", "family": "original", "true_force_snr": 8.,
            "true_acceleration_snr": 2., "true_mass": 4., "total_signal_snr": 68.**0.5,
            "methods": {"flat_joint": point, "norm_ratio": point, "positive_profile": failed, "forward_ols": missing},
            "paired_flat_minus_other": {
                "norm_ratio": {"squared_log": pair, "squared_mass": pair},
                "positive_profile": {"squared_log": pair, "outside_factor_2": pair}},
            "uncertainty": {"flat_joint": {
                "scores": {"log_crps": sample, "log_density_score": dict(sample, mean=-0.3)},
                "intervals": {"95": {"coverage": sample, "mean_log_width": sample, "interval_score": sample}}}},
            "paired_uncertainty_flat_minus_tube": {"log_crps": pair},
            "refinement": {"passed": True},
        }],
    }


class CatalogTests(unittest.TestCase):
    def setUp(self):
        # Normal mkdir keeps inherited workspace permissions on Windows.
        self.workspace = WORKSPACE / ".tools/mass_results_catalog_tests" / uuid.uuid4().hex
        self.workspace.mkdir(parents=True)

    def write(self, relative, data):
        path = self.workspace / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")
        return path

    def add_comparison(self, data=None):
        return self.write(".tools/mass_comparison_20260917/full/results.json", data or comparison_data())

    def test_comparison_keeps_values_pairing_and_provenance(self):
        self.add_comparison()
        run = load_run(self.workspace, "full-comparison")
        row = run["rows"][0]
        self.assertEqual(run["samples"], 100)
        self.assertEqual(row["mass"], 4.)
        self.assertEqual(row["paired"]["norm_ratio"]["squared_log"]["mean"], -0.1)
        self.assertEqual(row["paired"]["norm_ratio"]["squared_log"]["mcse"], 0.01)
        self.assertEqual(row["tolerance"], {"factors": [1., 2.], "methods": {method: [0., 0.8] for method in row["values"]}})
        self.assertEqual(row["diagnostics"]["flat_joint"]["bias"], 0.1)
        self.assertEqual(run["metadata"]["source_sha256"], {"engine.py": "abc123"})
        self.assertEqual(next(method for method in run["methods"] if method["id"] == "norm_ratio")["label"], "Recorded norm ratio")

    def test_population_risk_annotation_survives_finite_sample(self):
        self.add_comparison()
        row = load_run(self.workspace, "full-comparison")["rows"][0]
        value = row["values"]["positive_profile"]["squared_log"]
        self.assertEqual(value["mean"], 0.4)
        self.assertEqual(value["status"], "infinite_population_risk")
        self.assertEqual(row["paired"]["positive_profile"]["squared_log"]["status"], "infinite_population_risk")
        self.assertNotEqual(row["values"]["positive_profile"]["outside_factor_2"]["status"], "infinite_population_risk")
        self.assertNotEqual(row["paired"]["positive_profile"]["outside_factor_2"]["status"], "infinite_population_risk")

    def test_additional_law_keeps_paired_scores_and_oracle_label(self):
        data = comparison_data()
        source = data["rows"][0]
        source["uncertainty"]["calibrated_joint"] = source["uncertainty"]["flat_joint"]
        source["uncertainty"]["oracle_joint"] = source["uncertainty"]["flat_joint"]
        source["paired_uncertainty_flat_minus_other"] = {"calibrated_joint": {"coverage_95": {
            "mean_difference": -.04, "mcse": .006}}}
        self.add_comparison(data)
        run = load_run(self.workspace, "full-comparison")
        row = run["rows"][0]
        self.assertEqual(row["paired"]["calibrated_joint"]["coverage_95"]["mean"], -.04)
        self.assertEqual(row["paired"]["calibrated_joint"]["coverage_95"]["mcse"], .006)
        self.assertEqual(row["values"]["calibrated_joint"]["coverage_95"]["mean"], .9)
        oracle = next(method for method in run["methods"] if method["id"] == "oracle_joint")
        self.assertIn("oracle", oracle["label"])

    def test_repeated_reading_snr_and_calibration_provenance_survive_adapter(self):
        data = comparison_data()
        measurement = {"repeats": 25, "calibration": {"mode": "estimated", "samples": 20, "seed": 12}}
        data["experiment_config"] = {"name": "Calibration experiment", "measurement": measurement}
        data["measurement"] = dict(measurement, true_mean_force_sd=[.2, .2, .2],
                                   supplied_mean_force_sd_summary={"mean": [.19, .21, .2]})
        data["rows"][0].update(repeats=25, effective_force_snr=40., effective_acceleration_snr=10.,
                                effective_total_signal_snr=1700.**.5)
        self.add_comparison(data)
        run = load_run(self.workspace, "full-comparison")
        row = run["rows"][0]
        self.assertEqual((row["f"], row["a"]), (8., 2.))
        self.assertEqual((row["effective_force_snr"], row["effective_acceleration_snr"]), (40., 10.))
        self.assertEqual(row["repeats"], 25)
        self.assertEqual(row["effective_total_signal_snr"], 1700.**.5)
        self.assertEqual(run["config"]["measurement"], measurement)
        self.assertEqual(run["metadata"]["measurement"], data["measurement"])

    def test_heavy_tail_mcse_is_descriptive(self):
        self.add_comparison()
        row = load_run(self.workspace, "full-comparison")["rows"][0]
        for values in (row["values"]["norm_ratio"], row["paired"]["norm_ratio"]):
            self.assertEqual(values["squared_mass"]["status"], "descriptive_mcse")
            self.assertIn("infinite", values["squared_mass"]["note"])
        self.assertEqual(row["values"]["norm_ratio"]["heldout_known_acceleration_relative_squared_force_error"]["status"], "descriptive_mcse")

    def test_null_failure_is_never_zero(self):
        self.add_comparison()
        row = load_run(self.workspace, "full-comparison")["rows"][0]
        self.assertIsNone(row["values"]["forward_ols"]["squared_log"]["mean"])
        self.assertEqual(row["values"]["forward_ols"]["squared_log"]["status"], "invalid_estimates")
        self.assertEqual(row["diagnostics"]["forward_ols"]["invalid_count"], 2)

    def test_study_discovery_uses_only_completed_unchanged_attempts(self):
        base = ".tools/mass_estimator_studies/pilot/"
        relative = "experiments/n1/attempt-001"
        result = self.write(base + relative + "/results.json", comparison_data())
        record = {"id": "n1", "status": "completed", "attempt_directory": relative,
                  "output_sha256": {"results.json": hashlib.sha256(result.read_bytes()).hexdigest()}}
        checkpoint = {"schema_version": 1, "experiments": [record,
                      {**record, "id": "incomplete", "status": "failed"}]}
        self.write(base + "study_checkpoint.json", checkpoint)
        runs = list_runs(self.workspace)
        self.assertEqual([run["id"] for run in runs], ["study-pilot-n1"])
        self.assertEqual(load_run(self.workspace, "study-pilot-n1")["rows"][0]["mass"], 4.)
        result.write_text("{}", encoding="utf-8")
        self.assertEqual(list_runs(self.workspace), [])

    def test_null_inapplicability_is_not_relabelled_as_infinite_risk(self):
        data = comparison_data()
        metric = "heldout_known_acceleration_relative_squared_force_error"
        data["rows"][0]["methods"]["positive_profile"]["additional_tasks"] = {
            metric: {"mean": None, "empirical_mcse": None, "status": "not_applicable", "reason": "Zero true force"}}
        self.add_comparison(data)
        value = load_run(self.workspace, "full-comparison")["rows"][0]["values"]["positive_profile"][metric]
        self.assertEqual(value["status"], "not_applicable")
        self.assertEqual(value["note"], "Zero true force")

    def test_uncertainty_uses_target_and_context_semantics(self):
        self.add_comparison()
        run = load_run(self.workspace, "full-comparison")
        metrics = {metric["id"]: metric for metric in run["metrics"]}
        self.assertEqual(metrics["coverage_95"]["better"], "target")
        self.assertEqual(metrics["coverage_95"]["target"], 0.95)
        self.assertEqual(metrics["log_width_95"]["better"], "context")
        self.assertEqual(metrics["interval_score_95"]["better"], "lower")
        self.assertEqual(run["rows"][0]["values"]["flat_joint"]["log_density_score"]["mean"], -0.3)
        self.assertNotIn("coverage_95", run["rows"][0]["paired"].get("tube_joint", {}))

    def test_boundary_only_exposes_recorded_statistics(self):
        data = {"n": 1000, "seed": 15, "antithetic": True, "rows": [{
            "f": 8., "a": 2., "truth": 4., "ours": {"squared_log": {"mean": 0.4, "mcse": 0.01}, "bias": 0.1, "variance": 0.39, "success2": 0.8},
            "norm": {"squared_log": {"mean": 0.5, "mcse": 0.02}, "success2": 0.75},
            "delta": {"squared_log": {"mean": -0.1, "mcse": 0.005}, "root": {"mean": -0.02, "mcse": 0.001},
                      "factor_two_failure": {"mean": -0.05, "mcse": 0.01}, "squared_log_vs_floor": {"mean": -0.09, "mcse": 0.006}},
            "high_signal_gap_prediction": 0.001, "refinement": 1e-10}]}
        self.write(".tools/mass_boundary_20260918/high_antithetic.json", data)
        self.write(".tools/mass_boundary_20260918/known_force_limit.json", {"factor2_lower": 0.18})
        run = load_run(self.workspace, "boundary-high_antithetic")
        row = run["rows"][0]
        self.assertEqual(run["metadata"]["independent_units"], 500)
        self.assertEqual(row["paired"]["norm_ratio"]["reciprocal_root"]["mean"], -0.02)
        self.assertEqual(row["paired"]["norm_floor"]["squared_log"]["mean"], -0.09)
        self.assertNotIn("reciprocal_root", row["values"]["flat_joint"])
        self.assertNotIn("norm_floor", row["values"])
        self.assertAlmostEqual(row["values"]["flat_joint"]["outside_factor_2"]["mean"], 0.2)
        self.assertIsNone(row["values"]["flat_joint"]["outside_factor_2"]["mcse"])
        self.assertNotIn("invalid_count", row["diagnostics"]["flat_joint"])
        self.assertIn("known_force_limit_reference", run["metadata"])
        self.assertEqual(len(list_runs(self.workspace)), 1)

    def test_saved_config_name_status_and_history_id(self):
        self.write(".tools/mass_estimator_runs/run-123/results.json", comparison_data())
        self.write(".tools/mass_estimator_runs/run-123/config.json", {"name": "My experiment", "samples": 100})
        self.write(".tools/mass_estimator_runs/run-123/status.json", {"status": "completed"})
        runs = list_runs(self.workspace)
        self.assertEqual(runs[0]["id"], "run-123")
        self.assertEqual(runs[0]["title"], "My experiment")
        self.assertTrue(runs[0]["config_available"])
        run = load_run(self.workspace, "run-123")
        self.assertEqual(run["config"]["name"], "My experiment")
        self.assertEqual(run["metadata"]["run_status"]["status"], "completed")
        self.assertEqual(run["metadata"]["config_provenance"], "saved_config")

    def test_legacy_config_reconstruction_requires_complete_preset(self):
        data = comparison_data()
        self.add_comparison(data)
        self.assertIsNone(load_run(self.workspace, "full-comparison")["config"])
        data["rows"] *= 16
        self.add_comparison(data)
        run = load_run(self.workspace, "full-comparison")
        self.assertEqual(run["config"]["scenario"], {"type": "preset", "selection": "original"})
        self.assertEqual(run["config"]["samples"], 100)
        self.assertEqual(run["metadata"]["config_provenance"], "reconstructed_from_recorded_settings")

    def test_executed_config_wins_over_edited_sidecar(self):
        data = comparison_data()
        data["experiment_config"] = {"name": "Executed experiment", "samples": 100}
        self.write(".tools/mass_estimator_runs/run-executed/results.json", data)
        self.write(".tools/mass_estimator_runs/run-executed/config.json", {"name": "Later edit", "samples": 200})
        run = load_run(self.workspace, "run-executed")
        self.assertEqual(run["config"]["samples"], 100)
        self.assertEqual(run["title"], "Executed experiment")
        self.assertEqual(run["metadata"]["config_provenance"], "embedded_executed_config")
        self.assertEqual(list_runs(self.workspace)[0]["title"], "Executed experiment")

    def test_changed_file_invalidates_cache(self):
        data = comparison_data()
        self.add_comparison(data)
        first = load_run(self.workspace, "full-comparison")
        data["samples_per_scenario"] = 1024
        self.add_comparison(data)
        second = load_run(self.workspace, "full-comparison")
        self.assertEqual(first["samples"], 100)
        self.assertEqual(second["samples"], 1024)

    def test_unknown_paths_pending_and_broken_runs(self):
        self.add_comparison()
        self.write(".tools/mass_estimator_runs/pending/config.json", {"name": "Pending"})
        broken = self.write(".tools/mass_estimator_runs/broken/results.json", {})
        broken.write_text("{", encoding="utf-8")
        self.assertEqual([run["id"] for run in list_runs(self.workspace)], ["full-comparison"])
        for run_id in ("../full-comparison", "C:/Windows", "does-not-exist", "boundary-known_force_limit", "pending"):
            with self.assertRaises(KeyError):
                load_run(self.workspace, run_id)


if __name__ == "__main__":
    unittest.main()

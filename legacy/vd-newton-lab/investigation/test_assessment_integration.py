"""Scientific checks of true-null experiments and their saved evidence."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import unittest
from uuid import uuid4

import numpy as np

from experiment_config import run_config
from results_catalog import _comparison


class AssessmentIntegrationTests(unittest.TestCase):
    def run_experiment(self, **overrides):
        output = Path(__file__).resolve().parents[3] / ".tools" / "assessment_tests" / uuid4().hex
        config = {"samples": 8, "batch_size": 3, "save_trials": True,
                  "training_seed": 101, "heldout_seed": 202,
                  "scenario": {"type": "mass_excitation", "masses": [.25, 1., 4.],
                               "acceleration_snrs": [0.]}, **overrides}
        with redirect_stdout(io.StringIO()):
            report = run_config(config, output)
        arrays = []
        for row in report["rows"]:
            with np.load(output / "trials" / (row["id"] + ".npz"), allow_pickle=False) as saved:
                arrays.append({key: saved[key].copy() for key in saved.files})
        return report, arrays, output

    def test_true_null_observations_and_estimates_cannot_depend_on_mass(self):
        report, arrays, output = self.run_experiment()
        self.assertEqual([row["true_mass"] for row in report["rows"]], [.25, 1., 4.])
        for saved in arrays[1:]:
            for key in arrays[0]:
                if key.startswith(("observation__", "heldout__", "point__")):
                    np.testing.assert_array_equal(saved[key], arrays[0][key], err_msg=key)
        self.assertNotEqual(report["rows"][0]["methods"]["flat_joint"]["means"]["squared_log"],
                            report["rows"][1]["methods"]["flat_joint"]["means"]["squared_log"])
        self.assertTrue(all(row["refinement"]["passed"] for row in report["rows"]))
        self.assertEqual(json.loads((output/"trial_manifest.json").read_text())["schema_version"], 1)
        for filename in report["source_sha256"]:
            self.assertTrue((output/"source"/filename).is_file())

    def test_null_metrics_are_inapplicable_including_in_browser(self):
        report, _, _ = self.run_experiment()
        rows, _, _, _ = _comparison(report)
        for raw, adapted in zip(report["rows"], rows):
            self.assertEqual(len(raw["metric_applicability"]), 3)
            for method, summary in raw["methods"].items():
                for metric in raw["metric_applicability"]:
                    self.assertIsNone(summary["additional_tasks"][metric]["mean"])
                    value = adapted["values"][method][metric]
                    self.assertEqual(value["status"], "not_applicable")
                    self.assertTrue(value["note"])
                    if method != "flat_joint":
                        self.assertEqual(adapted["paired"][method][metric]["status"], "not_applicable")
            for law in raw["uncertainty"].values():
                for interval in law["intervals"].values():
                    probabilities = [interval[key]["mean"] for key in ("coverage", "below_interval", "above_interval")]
                    self.assertAlmostEqual(sum(probabilities), 1.)

    def test_more_readings_of_true_null_do_not_create_information(self):
        specification = {"type": "mass_excitation", "masses": [4.], "acceleration_snrs": [0.]}
        single, one, _ = self.run_experiment(scenario=specification)
        repeated, four, _ = self.run_experiment(scenario=specification, measurement={"repeats": 4})
        np.testing.assert_array_equal(four[0]["observation__force"]*2, one[0]["observation__force"])
        for key in one[0]:
            if key.startswith(("point__", "posterior__")):
                np.testing.assert_allclose(four[0][key], one[0][key], rtol=2e-12, atol=2e-12, equal_nan=True, err_msg=key)
        self.assertEqual(single["rows"][0]["uncertainty"]["flat_joint"]["intervals"]["95"]["coverage"],
                         repeated["rows"][0]["uncertainty"]["flat_joint"]["intervals"]["95"]["coverage"])

    def test_diagonal_null_has_no_mass_or_direction_information(self):
        settings = {"noise_model": "diagonal_gaussian_3d",
                    "noise": {"force_sd": [.5, 1., 2.], "acceleration_sd": [2., 1., .5]}}
        report, arrays, _ = self.run_experiment(**settings, direction=[1., 0., 0.])
        _, rotated, _ = self.run_experiment(**settings, direction=[0., 1., 1.])
        for current, other in zip(arrays, rotated):
            for key in current:
                if key.startswith(("observation__", "point__")):
                    np.testing.assert_array_equal(current[key], other[key], err_msg=key)
                    np.testing.assert_array_equal(current[key], arrays[0][key], err_msg=key)
        for row in report["rows"]:
            self.assertTrue(all(value is None for value in row["heldout_reference_values"].values()))

    def test_positive_scenario_matches_legacy_and_archives_heldout_readings(self):
        new, arrays, _ = self.run_experiment(scenario={"type": "mass_excitation", "masses": [4.], "acceleration_snrs": [2.]})
        old, old_arrays, _ = self.run_experiment(scenario={"type": "pairs", "pairs": [[8., 2.]]})
        self.assertEqual(new["rows"][0]["methods"], old["rows"][0]["methods"])
        self.assertEqual(new["rows"][0]["uncertainty"], old["rows"][0]["uncertainty"])
        for key in arrays[0]:
            np.testing.assert_array_equal(arrays[0][key], old_arrays[0][key], err_msg=key)
        saved = arrays[0]
        point = saved["point__flat_joint"]
        squared = np.sum((point[:, None]*saved["heldout__acceleration"]-saved["heldout__force"])**2, axis=1)
        denominator = np.sum(saved["heldout__true_force"]**2, axis=1)
        recorded = new["rows"][0]["methods"]["flat_joint"]["additional_tasks"]["heldout_noisy_acceleration_relative_squared_force_error"]["mean"]
        self.assertAlmostEqual(float(np.mean(squared/denominator)), recorded)


if __name__ == "__main__":
    unittest.main()

"""Independent calibration and repeated-reading protocol checks."""
from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import unittest
from unittest.mock import patch
from uuid import uuid4

import numpy as np

from compare_estimators import calibration_refinement_indices, measurement_noise
from experiment_config import run_config


class CalibrationGenerationTests(unittest.TestCase):
    def test_pooled_calibration_degrees_fresh_channels_and_shared_coordinates(self):
        count, k, repeats = 100000, 10, 4
        sf, sa = np.full(3, 2.), np.full(3, .5)
        config = {"measurement": {"repeats": repeats, "calibration": {"mode": "pooled_isotropic", "samples": k, "seed": 123}}}
        force, acceleration, metadata = measurement_noise(config, count, sf, sa)
        for values in (force, acceleration):
            np.testing.assert_array_equal(values[:, 0], values[:, 1])
            np.testing.assert_array_equal(values[:, 0], values[:, 2])
        ratios = np.column_stack((force[:, 0] / (sf[0]/math.sqrt(repeats)),
                                  acceleration[:, 0] / (sa[0]/math.sqrt(repeats))))**2
        np.testing.assert_allclose(np.mean(ratios, axis=0), 1., atol=.004)
        np.testing.assert_allclose(np.var(ratios, axis=0), 2/(3*(k-1)), atol=.002)
        self.assertLess(abs(np.corrcoef(ratios, rowvar=False)[0, 1]), .012)
        self.assertLess(abs(np.corrcoef(ratios[:-1, 0], ratios[1:, 0])[0, 1]), .012)
        self.assertEqual(metadata["calibration_degrees"], 3*(k-1))
        self.assertIn("no repeated-reading residual scatter", metadata["inference_evidence"])
        replay_force, replay_acceleration, _ = measurement_noise(config, count, sf, sa)
        np.testing.assert_array_equal(replay_force, force)
        np.testing.assert_array_equal(replay_acceleration, acceleration)
        with self.assertRaisesRegex(ValueError, "isotropic true noise"):
            measurement_noise(config, 3, np.array([1., 1., 2.]), sa)

    def test_refinement_includes_extreme_calibrations_outside_first_batch(self):
        force, acceleration = np.ones((40, 3)), np.ones((40, 3))
        sf, sa = np.ones((40, 3)), np.ones((40, 3))
        sf[33] = [1., 2., 8.]
        sa[34] = [1., 3., 16.]
        force[35] = [100., 100., 100.]
        acceleration[36] = [100., 100., 100.]
        sf[37] = [.01]*3
        sa[38] = [.001]*3
        indices = calibration_refinement_indices(force, acceleration, sf, sa)
        self.assertTrue(set(range(32)).issubset(indices))
        self.assertTrue({33, 34, 35, 37, 38}.issubset(indices))
        self.assertEqual(list(indices), sorted(set(indices)))

    def test_unbiased_variance_and_calibration_sampling_uncertainty(self):
        # For Gaussian calibration with estimated mean, sample variance has
        # expectation sigma^2 and relative variance 2/(k-1). This checks the
        # generating law without reproducing its individual random draws.
        count, calibration_samples = 100000, 30
        sf, sa = np.array([.5, 1., 2.]), np.array([2., 1., .5])
        config = {"measurement": {"repeats": 5, "calibration": {"mode": "estimated", "samples": calibration_samples, "seed": 123}}}
        force, acceleration, metadata = measurement_noise(config, count, sf, sa)
        ratios = np.concatenate((force/(sf/math.sqrt(5)), acceleration/(sa/math.sqrt(5))), axis=1)**2
        np.testing.assert_allclose(np.mean(ratios, axis=0), 1., atol=.004)
        np.testing.assert_allclose(np.var(ratios, axis=0), 2/(calibration_samples-1), atol=.002)
        off_diagonal = np.corrcoef(ratios, rowvar=False) - np.eye(6)
        self.assertLess(np.max(np.abs(off_diagonal)), .012)
        self.assertEqual(metadata["calibration"]["samples"], calibration_samples)
        self.assertIn("not integrated", metadata["calibration_scope"])

    def test_unsupported_calibration_fails_whole_run_without_clipping(self):
        class ExtremeCalibration:
            def chisquare(self, degrees, size):
                values = np.full(size, degrees, dtype=float)
                values[1, 2] *= 33**2
                return values
        config = {"measurement": {"repeats": 1, "calibration": {"mode": "estimated", "samples": 10, "seed": 123}}}
        with patch("compare_estimators.np.random.default_rng", return_value=ExtremeCalibration()):
            with self.assertRaisesRegex(FloatingPointError, "entire run.*without clipping or dropping"):
                measurement_noise(config, 3, np.ones(3), np.ones(3))


class MeasurementIntegrationTests(unittest.TestCase):
    def run_experiment(self, config):
        output = Path(__file__).resolve().parents[3] / ".tools" / "measurement_tests" / uuid4().hex
        with redirect_stdout(io.StringIO()):
            report = run_config({"samples": 8, "batch_size": 3, "save_trials": True,
                                 "training_seed": 101, "heldout_seed": 202,
                                 "scenario": {"type": "pairs", "pairs": [[8., 2.]]}, **config}, output)
        with np.load(output / "trials" / (report["rows"][0]["id"] + ".npz")) as saved:
            arrays = {key: saved[key].copy() for key in saved.files}
        return report, arrays

    def test_explicit_one_known_reading_preserves_legacy_results(self):
        old, old_arrays = self.run_experiment({})
        new, new_arrays = self.run_experiment({"measurement": {"repeats": 1, "calibration": {"mode": "known"}}})
        self.assertEqual(old["rows"][0]["methods"], new["rows"][0]["methods"])
        self.assertEqual(old["rows"][0]["uncertainty"], new["rows"][0]["uncertainty"])
        for key, value in old_arrays.items():
            np.testing.assert_array_equal(new_arrays[key], value)

    def test_pooled_plugin_integrated_and_oracle_share_mean_evidence_and_archive(self):
        from calibration_posterior import infer_batch as calibrated_infer
        from comparison_diagonal_posterior import infer_batch as diagonal_infer
        settings = {"repeats": 4, "calibration": {"mode": "pooled_isotropic", "samples": 30,
                                                 "seed": 303, "calibration_order": 24}}
        known, known_arrays = self.run_experiment({"measurement": {"repeats": 4}})
        with patch("calibration_posterior.infer_batch", wraps=calibrated_infer) as fitted:
            report, arrays = self.run_experiment({"measurement": settings})
        row = report["rows"][0]
        for channel in ("force", "acceleration"):
            np.testing.assert_array_equal(arrays["observation__" + channel], known_arrays["observation__" + channel])
        for suffix in ("joint", "geometric", "median"):
            np.testing.assert_array_equal(arrays["point__oracle_" + suffix], known_arrays["point__flat_" + suffix])
        # The first inference calls cover each batch exactly once. Calibration
        # degrees belong to the external experiment and do not include n=4.
        self.assertEqual([len(call.args[0]) for call in fitted.call_args_list[:3]], [3, 3, 2])
        for call in fitted.call_args_list:
            self.assertEqual(call.kwargs["degrees"], 87)
            self.assertEqual(np.ndim(call.kwargs["force_sd"]), 1)
        for index in (0, 3, 7):
            f, a = arrays["observation__force"][index:index+1], arrays["observation__acceleration"][index:index+1]
            sf = arrays["noise__supplied_mean_force_sd"][index, 0]
            sa = arrays["noise__supplied_mean_acceleration_sd"][index, 0]
            integrated = calibrated_infer(f, a, 4., force_sd=sf, acceleration_sd=sa, degrees=87)
            plugin = diagonal_infer(f, a, 4., force_sd=[sf]*3, acceleration_sd=[sa]*3)
            for name, values in integrated["points"].items():
                np.testing.assert_allclose(arrays["point__" + name][index], values[0], rtol=2e-12)
            np.testing.assert_allclose(arrays["point__flat_joint"][index], plugin["points"]["flat_joint"][0], rtol=2e-12)
        self.assertFalse(np.allclose(arrays["point__flat_joint"], arrays["point__calibrated_joint"]))
        self.assertEqual(set(row["uncertainty"]), {"flat_joint", "tube_joint", "calibrated_joint", "oracle_joint"})
        self.assertEqual(set(row["refinement"]["components"]),
                         {"plugin_angular", "integrated_angular", "integrated_calibration", "oracle_angular"})
        legacy = row["paired_uncertainty_flat_minus_tube"]
        self.assertEqual(legacy, {key: row["paired_uncertainty_flat_minus_other"]["tube_joint"][key] for key in legacy})
        self.assertIn("calibrated_joint", row["paired_uncertainty_flat_minus_other"])
        for law in ("calibrated_joint", "oracle_joint", "tube_joint"):
            pairs = row["paired_uncertainty_flat_minus_other"][law]
            for level in (50, 95):
                def covers(method):
                    return ((arrays[f"posterior__{method}__log_lower_{level}"] <= math.log(row["true_mass"])) &
                            (arrays[f"posterior__{method}__log_upper_{level}"] >= math.log(row["true_mass"]))).astype(float)
                difference = covers("flat_joint")-covers(law)
                self.assertEqual(pairs[f"coverage_{level}"]["mean_difference"], np.mean(difference))
                self.assertAlmostEqual(pairs[f"coverage_{level}"]["mcse"], np.std(difference, ddof=1)/np.sqrt(len(difference)))
                self.assertIn(f"log_width_{level}", pairs)
                self.assertNotIn(f"coverage_{level}", row["uncertainty"][law]["scores"])
        self.assertIn("posterior__calibrated_joint__log_lower_95", arrays)
        self.assertIn("posterior__oracle_joint__log_upper_50", arrays)
        self.assertIn("calibration_posterior.py", report["source_sha256"])
        self.assertIn("Only named oracle", report["assumptions"])

    def test_repeat_mean_matches_increased_snr_and_retains_heldout_task(self):
        repeats = 5
        repeated, arrays = self.run_experiment({"measurement": {"repeats": repeats}})
        single, single_arrays = self.run_experiment({})
        boosted, boosted_arrays = self.run_experiment({"scenario": {"type": "pairs", "pairs": [[8*math.sqrt(repeats), 2*math.sqrt(repeats)]]}})
        row = repeated["rows"][0]
        self.assertEqual(row["true_mass"], 4.)
        self.assertEqual(row["true_force_snr"], 8.)
        self.assertAlmostEqual(row["effective_force_snr"], 8*math.sqrt(repeats))
        self.assertAlmostEqual(row["effective_acceleration_snr"], 2*math.sqrt(repeats))
        self.assertEqual(row["heldout_reference_values"], single["rows"][0]["heldout_reference_values"])
        for channel, truth in (("force", 8.), ("acceleration", 2.)):
            key = "observation__" + channel
            latent = np.array([truth, 0., 0.])
            np.testing.assert_allclose(arrays[key]-latent, (single_arrays[key]-latent)/math.sqrt(repeats), atol=1e-15)
            np.testing.assert_allclose(arrays[key]*math.sqrt(repeats), boosted_arrays[key], atol=1e-14)
        for method in single["rows"][0]["methods"]:
            np.testing.assert_allclose(arrays["point__"+method], boosted_arrays["point__"+method], rtol=2e-12, atol=1e-13)
        self.assertIn("same latent", repeated["measurement"]["protocol"])

    def test_calibration_changes_fit_without_changing_observations_or_truth(self):
        known, known_arrays = self.run_experiment({"measurement": {"repeats": 2}})
        settings = {"repeats": 2, "calibration": {"mode": "estimated", "samples": 30, "seed": 303}}
        estimated, arrays = self.run_experiment({"measurement": settings})
        other, other_arrays = self.run_experiment({"measurement": {**settings, "calibration": {**settings["calibration"], "seed": 404}}})
        for key in ("observation__force", "observation__acceleration", "noise__true_force_sd", "noise__true_acceleration_sd"):
            np.testing.assert_array_equal(arrays[key], known_arrays[key])
            np.testing.assert_array_equal(arrays[key], other_arrays[key])
        for method in ("norm_ratio", "forward_ols", "reverse_ols"):
            np.testing.assert_array_equal(arrays["point__"+method], known_arrays["point__"+method])
        self.assertFalse(np.allclose(arrays["point__flat_joint"], known_arrays["point__flat_joint"]))
        self.assertFalse(np.allclose(arrays["point__flat_joint"], other_arrays["point__flat_joint"]))
        # Full calibration reaches every row of every inference batch; compare
        # individual calls to catch a reused first-batch calibration slice.
        from comparison_diagonal_posterior import infer_batch
        from comparison_diagonal_baselines import baseline_points
        for index in (0, 3, 7):
            f, a = arrays["observation__force"][index:index+1], arrays["observation__acceleration"][index:index+1]
            sf, sa = arrays["noise__supplied_mean_force_sd"][index], arrays["noise__supplied_mean_acceleration_sd"][index]
            expected = infer_batch(f, a, 4., force_sd=sf, acceleration_sd=sa)
            np.testing.assert_allclose(arrays["point__flat_joint"][index], expected["points"]["flat_joint"][0], rtol=2e-12)
            for method, values in baseline_points(f, a, sf, sa).items():
                np.testing.assert_allclose(arrays["point__"+method][index], values[0], rtol=2e-12, equal_nan=True)
        metadata = estimated["measurement"]
        np.testing.assert_allclose(metadata["supplied_force_sd_summary"]["mean"], arrays["noise__supplied_force_sd"].mean(axis=0))
        self.assertEqual(estimated["rows"][0]["heldout_reference_values"], known["rows"][0]["heldout_reference_values"])
        self.assertIn("supplied", estimated["baseline_metadata"]["flat_joint"]["assumptions"])
        self.assertIn("unknown direction", estimated["baseline_metadata"]["flat_joint"]["assumptions"])

    def test_scaled_calibration_is_applied_to_both_channels(self):
        config = {"noise_model": "diagonal_gaussian_3d",
                  "noise": {"force_sd": [.5, 1., 2.], "acceleration_sd": [2., 1., .5]},
                  "measurement": {"repeats": 10, "calibration": {"mode": "scaled", "force_scale": [.75, 1., 1.25], "acceleration_scale": [1.5]*3}}}
        report, arrays = self.run_experiment(config)
        for channel in ("force", "acceleration"):
            expected = np.array(config["noise"][channel+"_sd"])*config["measurement"]["calibration"][channel+"_scale"]
            np.testing.assert_allclose(arrays["noise__supplied_"+channel+"_sd"], np.broadcast_to(expected, (8, 3)))
            np.testing.assert_allclose(arrays["noise__supplied_mean_"+channel+"_sd"], np.broadcast_to(expected/math.sqrt(10), (8, 3)))
        self.assertEqual(report["measurement"]["calibration"]["mode"], "scaled")

    def test_rms_control_respects_batch_bound_and_per_trial_calibration(self):
        from comparison_posterior import infer_batch
        for calibration in ({"mode": "known"},
                            {"mode": "estimated", "samples": 30, "seed": 303}):
            with self.subTest(calibration=calibration):
                with patch("compare_estimators.infer_batch", wraps=infer_batch) as control:
                    report, arrays = self.run_experiment({"measurement": {
                        "repeats": 4, "calibration": calibration}})
                self.assertEqual([len(call.args[0]) for call in control.call_args_list], [3, 3, 2])
                sf = np.linalg.norm(arrays["noise__supplied_mean_force_sd"], axis=1)/np.sqrt(3)
                sa = np.linalg.norm(arrays["noise__supplied_mean_acceleration_sd"], axis=1)/np.sqrt(3)
                mass_scale = sf/sa
                expected = infer_batch(arrays["observation__force"]/sf[:, None],
                                       arrays["observation__acceleration"]/sa[:, None],
                                       report["rows"][0]["true_mass"]/mass_scale)
                np.testing.assert_allclose(arrays["point__flat_joint_rms"],
                                           expected["points"]["flat_joint"]*mass_scale,
                                           rtol=2e-12, atol=1e-13)


if __name__ == "__main__":
    unittest.main()

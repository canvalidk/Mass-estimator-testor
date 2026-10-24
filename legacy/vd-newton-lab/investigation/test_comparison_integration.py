"""Bounded integration checks for the external-truth comparison harness."""

from contextlib import nullcontext, redirect_stdout
from copy import deepcopy
import io
import itertools
import json
from pathlib import Path
import unittest
from uuid import uuid4

import numpy as np

from compare_estimators import (
    TRAIN_SEED, predictive_losses, refinement_check, run, scenarios,
)
from comparison_estimators import baseline_points
from comparison_posterior import infer_batch


POINT_METHODS = {
    "flat_joint", "tube_joint", "flat_geometric", "flat_median",
    "flat_reciprocal_root", "norm_ratio", "norm_floor", "positive_profile",
    "forward_ols", "reverse_ols", "noise_corrected_ols",
}


class ScenarioTests(unittest.TestCase):
    def test_scenario_counts_and_physical_identities(self):
        original = scenarios("original")
        fixed = scenarios("fixed_mass")
        combined = scenarios("all")
        self.assertEqual((len(original), len(fixed), len(combined)), (16, 30, 46))
        self.assertEqual(len({row["id"] for row in combined}), 46)
        for row in combined:
            force = row["true_force_snr"]
            acceleration = row["true_acceleration_snr"]
            self.assertAlmostEqual(force / acceleration, row["true_mass"])
            self.assertAlmostEqual(np.hypot(force, acceleration), row["total_signal_snr"])
        for mass in (0.25, 1, 4, 16, 32):
            group = [row for row in fixed if row["true_mass"] == mass]
            self.assertEqual(len(group), 6)
            self.assertEqual({row["total_signal_snr"] for row in group}, {0.5, 1, 2, 4, 8, 16})
            np.testing.assert_allclose(
                [row["true_force_snr"] for row in group],
                mass * np.array([row["true_acceleration_snr"] for row in group]),
            )
        equal_channel_originals = [row for row in original if row["true_force_snr"] == row["true_acceleration_snr"]]
        self.assertGreaterEqual(len(equal_channel_originals), 5)
        self.assertTrue(all(row["true_mass"] == 1 for row in equal_channel_originals))

    def test_unknown_selection_is_rejected(self):
        with self.assertRaises(ValueError):
            scenarios("missing")


class PredictiveLossTests(unittest.TestCase):
    def test_symmetric_noise_quadrature_matches_analytic_expectations(self):
        # Quadratic losses depend only on these moments. All 64 sign vectors
        # give exact zero means, identity covariance, and independent channels,
        # so this is an exact Gaussian-moment oracle rather than a random test.
        noise = np.array(list(itertools.product((-1.0, 1.0), repeat=6)))
        np.testing.assert_array_equal(noise.mean(axis=0), np.zeros(6))
        np.testing.assert_array_equal(noise.T @ noise / len(noise), np.eye(6))
        for truth, acceleration, relative_point in ((0.25, 2, 1), (3, 0.25, 0.4), (32, 0.25, 2)):
            with self.subTest(truth=truth, acceleration=acceleration, relative_point=relative_point):
                point = truth * relative_point
                force = truth * acceleration
                losses = predictive_losses(
                    np.full(len(noise), point), truth, acceleration,
                    noise[:, :3], noise[:, 3:],
                )
                known = (relative_point - 1) ** 2 + 3 / (4 * force ** 2)
                noisy = known + 3 * point ** 2 / (4 * force ** 2)
                self.assertAlmostEqual(
                    losses["heldout_known_acceleration_relative_squared_force_error"].mean(), known,
                )
                self.assertAlmostEqual(
                    losses["heldout_noisy_acceleration_relative_squared_force_error"].mean(), noisy,
                )

    def test_invalid_points_receive_infinite_prediction_loss(self):
        points = np.array([0, -1, np.nan, np.inf, -np.inf, 2.0])
        losses = predictive_losses(points, 2.0, 3.0, np.zeros((6, 3)), np.zeros((6, 3)))
        for values in losses.values():
            self.assertTrue(np.isposinf(values[:5]).all())
            self.assertEqual(values[5], 0)


class RefinementCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        force = np.array([[0, 0, 0], [8, 1, 0], [3, 0, 0]], float)
        acceleration = np.array([[0, 0, 0], [0.25, -0.3, 0.2], [-2, 0, 0]], float)
        cls.coarse = infer_batch(force, acceleration, truth=4, order=96)
        cls.fine = infer_batch(force, acceleration, truth=4, order=192)

    def test_actual_posterior_refinement_passes(self):
        result = refinement_check(self.coarse, self.fine)
        self.assertTrue(result["passed"])
        self.assertLessEqual(result["maximum_relative_point_difference"], result["relative_point_limit"])
        self.assertTrue(result["maximum_absolute_distribution_difference"])
        for key, difference in result["maximum_absolute_distribution_difference"].items():
            self.assertLessEqual(difference, result["absolute_distribution_limits"][key])

    def test_material_point_and_distribution_differences_fail(self):
        altered = deepcopy(self.fine)
        altered["points"]["flat_joint"][0] *= 1.01
        self.assertFalse(refinement_check(self.coarse, altered)["passed"])
        altered = deepcopy(self.fine)
        altered["distributions"]["tube_joint"]["log_lower_95"][0] += 0.01
        self.assertFalse(refinement_check(self.coarse, altered)["passed"])

    def test_nonfinite_point_or_distribution_is_an_error(self):
        for position in ("fine_point", "coarse_point", "fine_distribution", "coarse_distribution"):
            with self.subTest(position=position):
                coarse, fine = deepcopy(self.coarse), deepcopy(self.fine)
                target = fine if position.startswith("fine") else coarse
                if position.endswith("point"):
                    target["points"]["flat_joint"][0] = np.nan
                else:
                    target["distributions"]["flat_joint"]["log_crps"][0] = np.inf
                with self.assertRaises(FloatingPointError):
                    refinement_check(coarse, fine)


class EndToEndComparisonTests(unittest.TestCase):
    def test_small_run_accounts_for_every_trial_and_round_trips_outputs(self):
        scratch = Path(__file__).resolve().parents[3] / ".tools" / "mass_comparison_tests"
        scratch.mkdir(parents=True, exist_ok=True)
        # Keep generated test artifacts in the declared workspace scratch area;
        # do not require system-temp permissions or remove existing artifacts.
        with nullcontext(scratch / ("integration-" + uuid4().hex)) as temporary:
            output = Path(temporary).resolve()
            self.assertEqual(output.parent, scratch.resolve())
            with redirect_stdout(io.StringIO()):
                report = run(16, "original", 96, 7, output, save_trials=True, limit=1)

            def reject_nonfinite_json(value):
                raise AssertionError("Non-standard nonfinite JSON constant: " + value)

            saved_report = json.loads(
                (output / "results.json").read_text(encoding="utf-8"),
                parse_constant=reject_nonfinite_json,
            )
            self.assertEqual(saved_report, report)
            json.dumps(report, allow_nan=False)
            self.assertEqual(report["samples_per_scenario"], 16)
            self.assertEqual(len(report["rows"]), 1)
            factors = np.array(report["tolerance_factors"])
            self.assertEqual(len(factors), 401)
            self.assertEqual((factors[0], factors[-1]), (1, 20))
            row = report["rows"][0]
            self.assertTrue(row["refinement"]["passed"])
            self.assertEqual(set(row["methods"]), POINT_METHODS)
            self.assertEqual(set(row["paired_flat_minus_other"]), POINT_METHODS - {"flat_joint"})
            self.assertEqual(set(row["uncertainty"]), {"flat_joint", "tube_joint"})
            self.assertEqual(
                set(row["paired_uncertainty_flat_minus_tube"]),
                {"log_crps", "log_density_score", "interval_score_50", "interval_score_80", "interval_score_95"},
            )

            with np.load(output / "trials" / (row["id"] + ".npz"), allow_pickle=False) as trials:
                for name, summary in row["methods"].items():
                    with self.subTest(method=name):
                        values = trials["point__" + name]
                        self.assertEqual(values.shape, (16,))
                        valid = np.isfinite(values) & (values > 0)
                        self.assertEqual(summary["n"], 16)
                        self.assertEqual(summary["valid_count"], int(valid.sum()))
                        self.assertEqual(summary["invalid_count"], int((~valid).sum()))
                        self.assertEqual(summary["valid_count"] + summary["invalid_count"], 16)
                        curve = np.array(summary["tolerance_curve"])
                        self.assertEqual(curve.shape, (401,))
                        self.assertTrue(np.all(np.diff(curve) >= 0))
                        self.assertTrue(np.all((curve >= 0) & (curve <= valid.mean())))
                        if not valid.all():
                            self.assertIsNone(summary["means"]["squared_log"])
                        for factor, result in summary["tolerance_success"].items():
                            tolerance = float(factor)
                            success = valid & (values >= row["true_mass"] / tolerance) & (values <= row["true_mass"] * tolerance)
                            self.assertAlmostEqual(result["fraction"], success.mean())

                # Check that the binary artifact preserves actual invalid
                # outputs, including NaN, instead of substituting a mass.
                noise = np.random.default_rng(TRAIN_SEED).standard_normal((16, 6))
                expected_baselines = baseline_points(
                    noise[:, :3] + [row["true_force_snr"], 0, 0],
                    noise[:, 3:] + [row["true_acceleration_snr"], 0, 0],
                )
                for name, expected in expected_baselines.items():
                    np.testing.assert_array_equal(trials["point__" + name], expected)
                self.assertTrue(np.isnan(trials["point__noise_corrected_ols"]).any())

                for law, uncertainty in row["uncertainty"].items():
                    self.assertEqual(set(uncertainty["intervals"]), {"50", "80", "95"})
                    for level, interval in uncertainty["intervals"].items():
                        self.assertEqual(interval["nominal_coverage"], int(level) / 100)
                        coverage = interval["coverage"]["mean"]
                        self.assertGreaterEqual(coverage, 0)
                        self.assertLessEqual(coverage, 1)
                    for field in ("log_crps", "log_density_score", "cdf_truth", "log_lower_95", "log_upper_95"):
                        values = trials[f"posterior__{law}__{field}"]
                        self.assertEqual(values.shape, (16,))
                        self.assertTrue(np.isfinite(values).all())

                # Pairing must survive the estimator, scorer, and serializer
                # boundaries; bounded errors permit every original trial.
                flat = trials["point__flat_joint"]
                flat_loss = np.minimum(np.log(flat / row["true_mass"]) ** 2, np.log(2) ** 2)
                for name, comparisons in row["paired_flat_minus_other"].items():
                    values = trials["point__" + name]
                    good = np.isfinite(values) & (values > 0)
                    other_loss = np.full(16, np.log(2) ** 2)
                    other_loss[good] = np.minimum(np.log(values[good] / row["true_mass"]) ** 2, np.log(2) ** 2)
                    comparison = comparisons["capped_squared_log_factor_2"]
                    self.assertEqual(comparison["n"], 16)
                    self.assertEqual(comparison["finite_pair_count"], 16)
                    self.assertAlmostEqual(comparison["mean_difference"], (flat_loss - other_loss).mean())
                    self.assertAlmostEqual(comparison["mcse"], np.std(flat_loss - other_loss, ddof=1) / 4)


if __name__ == "__main__":
    unittest.main()

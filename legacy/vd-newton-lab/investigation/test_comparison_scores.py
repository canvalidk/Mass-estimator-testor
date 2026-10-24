"""Analytic and failure-accounting checks for external-truth scoring."""

import json
import unittest

import numpy as np

from comparison_scores import (
    CAP_FACTORS, UNBOUNDED_LOSSES, interval_score, paired_summary,
    point_losses, score_point,
)


class PointScoreTests(unittest.TestCase):
    def test_exact_estimates_have_zero_loss(self):
        losses = point_losses(np.full((2, 3), 7.0), 7.0)
        for values in losses.values():
            np.testing.assert_array_equal(values, np.zeros((2, 3)))
        summary = score_point(np.full(6, 7.0), 7.0)
        self.assertEqual(summary["invalid_count"], 0)
        self.assertTrue(all(x == 0 for x in summary["means"].values()))
        self.assertEqual(summary["log_diagnostics"]["signed_log_bias"], 0)
        self.assertEqual(summary["log_diagnostics"]["log_variance"], 0)

    def test_known_losses_and_inclusive_tolerances(self):
        values = np.array([5.0, 10.0, 20.0])
        losses = point_losses(values, 10.0)
        ln2 = np.log(2.0)
        np.testing.assert_allclose(losses["abs_log"], [ln2, 0, ln2])
        np.testing.assert_allclose(losses["squared_log"], [ln2 ** 2, 0, ln2 ** 2])
        root_error = np.sqrt(2) + 1 / np.sqrt(2) - 2
        np.testing.assert_allclose(losses["reciprocal_root"], [root_error, 0, root_error])
        np.testing.assert_allclose(losses["abs_relative"], [0.5, 0, 1])
        np.testing.assert_allclose(losses["squared_relative"], [0.25, 0, 1])
        np.testing.assert_allclose(losses["absolute_mass"], [5, 0, 10])
        np.testing.assert_allclose(losses["squared_mass"], [25, 0, 100])
        summary = score_point(values, 10)
        self.assertAlmostEqual(summary["means"]["absolute_mass"], 5)
        self.assertAlmostEqual(summary["means"]["squared_mass"], 125 / 3)
        self.assertEqual(summary["tolerance_success"]["2"]["fraction"], 1)
        self.assertAlmostEqual(summary["tolerance_success"]["1.5"]["fraction"], 1 / 3)
        diagnostics = summary["log_diagnostics"]
        self.assertAlmostEqual(diagnostics["signed_log_bias"], 0)
        self.assertAlmostEqual(diagnostics["log_variance"], 2 * ln2 ** 2 / 3)
        self.assertAlmostEqual(
            summary["means"]["squared_log"],
            diagnostics["signed_log_bias"] ** 2 + diagnostics["log_variance"],
        )

    def test_invalid_estimates_are_never_dropped(self):
        values = np.array([10.0, 0, -3, np.nan, np.inf, -np.inf])
        losses = point_losses(values, 10)
        for name in UNBOUNDED_LOSSES:
            self.assertEqual(losses[name][0], 0)
            self.assertTrue(np.isposinf(losses[name][1:]).all())
        for factor in CAP_FACTORS:
            name = "capped_squared_log_factor_" + format(factor, "g")
            np.testing.assert_allclose(losses[name][1:], np.log(factor) ** 2)
        summary = score_point(values, 10)
        self.assertEqual(summary["n"], 6)
        self.assertEqual(summary["invalid_count"], 5)
        for name in UNBOUNDED_LOSSES:
            self.assertIsNone(summary["means"][name])
            self.assertEqual(summary["mean_status"][name], "invalid_estimates")
        for result in summary["tolerance_success"].values():
            self.assertAlmostEqual(result["fraction"], 1 / 6)
        self.assertIsNone(summary["log_diagnostics"]["signed_log_bias"])
        self.assertIsNone(summary["absolute_log_quantiles"]["0.5"])
        json.dumps(summary, allow_nan=False)

    def test_capped_loss_has_declared_units_and_maximum(self):
        losses = point_losses([1, 2, 10, np.nan], 1)
        np.testing.assert_allclose(
            losses["capped_squared_log_factor_5"],
            [0, np.log(2) ** 2, np.log(5) ** 2, np.log(5) ** 2],
        )

    def test_quantiles_include_failures_without_nan_interpolation(self):
        summary = score_point([1, np.nan], 1)
        self.assertTrue(all(x is None for x in summary["absolute_log_quantiles"].values()))
        partial = score_point(np.r_[np.ones(9), np.nan], 1)
        self.assertEqual(partial["absolute_log_quantiles"]["0.5"], 0)
        self.assertIsNone(partial["absolute_log_quantiles"]["0.99"])
        finite = score_point(np.exp([0, 1, 2, 3, 4]), 1)
        self.assertAlmostEqual(finite["absolute_log_quantiles"]["0.5"], 2)
        self.assertAlmostEqual(finite["absolute_log_quantiles"]["0.9"], 3.6)

    def test_reciprocal_symmetry_for_symmetric_scores(self):
        values, truth = np.array([0.1, 0.4, 1, 1.7, 20]), 2.0
        forward, inverse = point_losses(values, truth), point_losses(1 / values, 1 / truth)
        symmetric = ["abs_log", "squared_log", "reciprocal_root"]
        symmetric.extend(name for name in forward if name.startswith(("capped_", "outside_")))
        for name in symmetric:
            np.testing.assert_allclose(forward[name], inverse[name], rtol=1e-13, atol=1e-14)

    def test_unit_changes_preserve_dimensionless_scores(self):
        values, truth, scale = np.array([0.2, 0.7, 1, 4, 8]), 1.7, 1000.0
        original, scaled = point_losses(values, truth), point_losses(values * scale, truth * scale)
        for name in original:
            multiplier = scale if name == "absolute_mass" else scale ** 2 if name == "squared_mass" else 1
            np.testing.assert_allclose(scaled[name], original[name] * multiplier, rtol=1e-13)

    def test_extreme_finite_values_do_not_overflow_log_scores(self):
        summary = score_point([1e-300, 1e300], 1e-300)
        self.assertEqual(summary["invalid_count"], 0)
        self.assertTrue(np.isfinite(summary["means"]["squared_log"]))
        self.assertIsNone(summary["means"]["squared_mass"])
        self.assertEqual(summary["mean_status"]["squared_mass"], "nonfinite_loss")
        json.dumps(summary, allow_nan=False)

    def test_truth_and_empty_validation(self):
        for truth in (0, -1, np.inf, np.nan, [1]):
            with self.subTest(truth=truth), self.assertRaises(ValueError):
                point_losses([1], truth)
        with self.assertRaises(ValueError):
            score_point([], 1)

    def test_single_trial_has_no_estimated_mcse(self):
        summary = score_point([2], 1)
        self.assertTrue(all(x is None for x in summary["mcse"].values()))


class PairedScoreTests(unittest.TestCase):
    def test_analytic_difference_and_mcse(self):
        result = paired_summary([1, 2, 3], [0, 2, 4])
        self.assertEqual(result["mean_difference"], 0)
        self.assertAlmostEqual(result["mcse"], 1 / np.sqrt(3))
        better = paired_summary([0, 0], [1, 1])
        self.assertEqual(better["mean_difference"], -1)
        self.assertEqual(better["mcse"], 0)

    def test_nonfinite_pairs_do_not_create_conditional_comparison(self):
        result = paired_summary([0, np.inf], [1, 2])
        self.assertEqual(result["finite_pair_count"], 1)
        self.assertIsNone(result["mean_difference"])
        self.assertIsNone(result["mcse"])
        self.assertEqual(result["status"], "nonfinite_pairs")

    def test_shape_validation_and_single_pair(self):
        for first, second in (([], []), ([1], [1, 2]), ([[1, 2]], [1, 2])):
            with self.assertRaises(ValueError):
                paired_summary(first, second)
        result = paired_summary([1], [2])
        self.assertEqual(result["mean_difference"], -1)
        self.assertIsNone(result["mcse"])

    def test_large_finite_losses_use_stable_sample_statistics(self):
        result = paired_summary([1e308, 1e308], [0, 0])
        self.assertEqual(result["mean_difference"], 1e308)
        self.assertEqual(result["mcse"], 0)


class IntervalScoreTests(unittest.TestCase):
    def test_analytic_inside_and_outside_scores(self):
        np.testing.assert_allclose(interval_score(0, 1, [-1, 0, 0.5, 1, 2], 0.2), [11, 1, 1, 1, 11])
        self.assertEqual(float(interval_score(0, 1, 2, 0.1)), 21)
        # A tiny valid alpha must not turn an exactly zero penalty into NaN.
        self.assertEqual(float(interval_score(0, 1, 0.5, np.nextafter(0.0, 1.0))), 1)

    def test_log_unit_shift_and_reciprocal_invariance(self):
        lower, upper, truth = np.array([-2, 1]), np.array([0, 2]), np.array([-3, 1.5])
        original = interval_score(lower, upper, truth, 0.1)
        np.testing.assert_allclose(interval_score(lower + 7, upper + 7, truth + 7, 0.1), original)
        np.testing.assert_allclose(interval_score(-upper, -lower, -truth, 0.1), original)

    def test_invalid_predicted_intervals_receive_infinite_loss(self):
        score = interval_score([1, -np.inf, np.nan], [0, np.inf, 1], 0, 0.1)
        self.assertTrue(np.isposinf(score).all())

    def test_alpha_and_truth_validation(self):
        for alpha in (0, 1, -0.1, np.nan, [0.1]):
            with self.subTest(alpha=alpha), self.assertRaises(ValueError):
                interval_score(0, 1, 0, alpha)
        with self.assertRaises(ValueError):
            interval_score(0, 1, np.nan, 0.1)


if __name__ == "__main__":
    unittest.main()

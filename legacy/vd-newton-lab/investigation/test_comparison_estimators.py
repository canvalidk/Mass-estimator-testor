"""Checks of comparator objectives, singular cases and coordinate invariance."""

import unittest

import numpy as np

from comparison_estimators import METHOD_METADATA, baseline_points


class ComparisonEstimatorTests(unittest.TestCase):
    def test_aligned_exact_readings(self):
        acceleration = np.array([[2., -1., 3.], [-3., 2., 1.], [4., 0., 0.]])
        masses = np.array([.2, 2., 20.])
        result = baseline_points(masses[:, None] * acceleration, acceleration)
        for name in ("norm_ratio", "forward_ols", "reverse_ols", "positive_profile"):
            np.testing.assert_allclose(result[name], masses, rtol=2e-15)
        # Exact readings do not cancel a declared nonzero measurement variance.
        self.assertTrue(np.all(result["noise_corrected_ols"] > masses))
        self.assertEqual(set(result), set(METHOD_METADATA))

    def test_singular_and_nonphysical_results_remain_visible(self):
        force = [[0., 0., 0.], [2., 0., 0.], [0., 0., 0.],
                 [2., 0., 0.], [1., 0., 0.], [-4., 0., 0.]]
        accel = [[0., 0., 0.], [0., 0., 0.], [2., 0., 0.],
                 [0., 1., 0.], [0., 1., 0.], [2., 0., 0.]]
        result = baseline_points(force, accel)
        np.testing.assert_equal(result["norm_ratio"], [np.nan, np.inf, 0., 2., 1., 2.])
        np.testing.assert_equal(result["forward_ols"], [np.nan, np.nan, 0., 0., 0., -2.])
        np.testing.assert_equal(result["reverse_ols"], [np.nan, np.inf, np.nan, np.inf, np.inf, -2.])
        np.testing.assert_equal(result["positive_profile"], [np.nan, np.inf, 0., np.inf, np.nan, np.inf])
        np.testing.assert_equal(result["noise_corrected_ols"], [np.nan, np.nan, 0., np.nan, np.nan, -8.])
        self.assertTrue(np.all(np.isfinite(result["norm_floor"])))
        self.assertTrue(np.all(result["norm_floor"] > 0))

    def test_noise_correction_rejects_nonpositive_energy(self):
        result = baseline_points([[1., 1., 1.], [2., 0., 0.], [0., 1., 0.]],
                                 [[1., 1., 1.], [1., 0., 0.], [2., 0., 0.]])
        np.testing.assert_equal(result["noise_corrected_ols"], [np.nan, np.nan, 0.])

    def test_profile_solves_independent_svd_problem(self):
        rng = np.random.default_rng(2026091701)
        accel = rng.normal(size=(60, 3))
        force = 1.7 * accel + .3 * rng.normal(size=(60, 3))
        result = baseline_points(force, accel)["positive_profile"]
        for f, a, point in zip(force, accel, result):
            # Rank-one fit to columns [a,F] is an independent TLS construction.
            _, _, vh = np.linalg.svd(np.column_stack((a, f)), full_matrices=False)
            svd_slope = vh[0, 1] / vh[0, 0]
            if svd_slope <= 0:
                expected = np.inf if np.dot(f, f) >= np.dot(a, a) else 0.
                self.assertEqual(point, expected)
                continue
            self.assertAlmostEqual(point, svd_slope, places=12)
            value = np.dot(f - point * a, f - point * a) / (1 + point**2)
            for other in (.7 * point, 1.3 * point):
                self.assertLessEqual(value, np.dot(f - other*a, f - other*a)/(1 + other**2))

    def test_ols_matches_least_squares_and_corrected_moments(self):
        rng = np.random.default_rng(2026091703)
        force, accel = rng.normal(size=(2, 12, 3))
        result = baseline_points(force, accel)
        for i, (f, a) in enumerate(zip(force, accel)):
            forward = np.linalg.lstsq(a[:, None], f, rcond=None)[0][0]
            inverse_mass = np.linalg.lstsq(f[:, None], a, rcond=None)[0][0]
            self.assertAlmostEqual(result["forward_ols"][i], forward, places=13)
            self.assertAlmostEqual(result["reverse_ols"][i], 1/inverse_mass, places=11)
        # A symmetric finite quadrature integrates moments of independent unit
        # errors exactly: this validates the correction, not unbiased ratios.
        errors = np.vstack((np.eye(3), -np.eye(3))) * np.sqrt(3.)
        true_accel = np.array([2., -.5, 1.])
        mass = 4.
        observed_accel = np.repeat(true_accel[None, :] + errors, 6, axis=0)
        observed_force = np.tile(mass*true_accel[None, :] + errors, (6, 1))
        dot = np.sum(observed_force * observed_accel, axis=1)
        corrected_energy = np.sum(observed_accel**2, axis=1) - 3.
        self.assertAlmostEqual(np.mean(corrected_energy), np.dot(true_accel, true_accel))
        self.assertAlmostEqual(np.mean(dot - mass*corrected_energy), 0., places=12)

    def test_rotation_and_channel_exchange(self):
        rng = np.random.default_rng(2026091702)
        force, accel = rng.normal(size=(2, 40, 3))
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        # Keep the exchange check away from boundary and degenerate cases.
        force = 2 * accel + .4 * force
        original = baseline_points(force, accel)
        rotated = baseline_points(force @ rotation, accel @ rotation)
        swapped = baseline_points(accel, force)
        for name in original:
            np.testing.assert_allclose(rotated[name], original[name], rtol=1e-12, atol=1e-12)
        for name in ("norm_ratio", "norm_floor", "positive_profile"):
            np.testing.assert_allclose(original[name] * swapped[name], 1., rtol=1e-12)
        np.testing.assert_allclose(original["forward_ols"] * swapped["reverse_ols"], 1., rtol=1e-12)

    def test_input_validation_and_single_pair(self):
        self.assertEqual(baseline_points([2, 0, 0], [1, 0, 0])["forward_ols"].shape, (1,))
        for force, accel in (([1, 2], [1, 2]), ([[1, 2, 3]], [[1, 2, 3], [4, 5, 6]]),
                             ([1, np.inf, 3], [1, 2, 3])):
            with self.assertRaises(ValueError):
                baseline_points(force, accel)
        self.assertEqual(baseline_points(np.empty((0, 3)), np.empty((0, 3)))["norm_ratio"].size, 0)


if __name__ == "__main__":
    unittest.main()

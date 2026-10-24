"""Physical-unit and global-objective checks for diagonal-noise comparators."""

import math
import unittest

import numpy as np

from comparison_estimators import baseline_points as isotropic_points
from comparison_diagonal_baselines import METHOD_METADATA, baseline_points, covariance_profile


def objective(mass, force, acceleration, sf, sa):
    """Independent direct evaluation in physical coordinates."""
    mass = np.asarray(mass)
    with np.errstate(invalid="ignore", divide="ignore"):
        values = np.sum((force - mass[..., None]*acceleration)**2 /
                        (sf**2 + mass[..., None]**2*sa**2), axis=-1)
    return np.where(np.isposinf(mass), np.sum((acceleration/sa)**2), values)


class DiagonalBaselineTests(unittest.TestCase):
    def test_unit_noise_reproduces_all_existing_baselines(self):
        rng = np.random.default_rng(63014)
        force, acceleration = rng.normal(size=(2, 60, 3))
        force = np.vstack((force, [0, 0, 0], [1, 0, 0], [1, 0, 0], [2, 0, 0]))
        acceleration = np.vstack((acceleration, [0, 0, 0], [0, 0, 0], [0, 1, 0], [-1, 0, 0]))
        old = isotropic_points(force, acceleration)
        new = baseline_points(force, acceleration, [1, 1, 1], [1, 1, 1])
        self.assertEqual(len(new), 10)
        self.assertEqual(set(new), set(METHOD_METADATA))
        for name, expected in old.items():
            np.testing.assert_array_equal(new[name], expected, err_msg=name)
        for added, previous in (("covariance_profile", "positive_profile"),
                                ("weighted_forward_ols", "forward_ols"),
                                ("weighted_reverse_ols", "reverse_ols"),
                                ("weighted_corrected_ols", "noise_corrected_ols")):
            np.testing.assert_array_equal(new[added], old[previous], err_msg=added)

    def test_isotropic_nonunit_noise_matches_scalar_standardization(self):
        rng = np.random.default_rng(309)
        force, acceleration = rng.normal(size=(2, 40, 3))
        expected = 2.5/.7 * isotropic_points(force/2.5, acceleration/.7)["positive_profile"]
        actual = covariance_profile(force, acceleration, [2.5]*3, [.7]*3)
        np.testing.assert_allclose(actual, expected, rtol=1e-14, equal_nan=True)

    def test_weighted_formulas_and_corrected_unavailability(self):
        force = np.array([[2., 4., 10.], [2., 4., 10.]])
        acceleration = np.array([[1., 2., 3.], [3., 4., 6.]])
        result = baseline_points(force, acceleration, [1, 2, 4], [2, 1, .5])
        self.assertAlmostEqual(result["weighted_forward_ols"][0], 5.875/2.5625)
        self.assertAlmostEqual(result["weighted_reverse_ols"][0], 417/128.5)
        self.assertTrue(math.isnan(result["weighted_corrected_ols"][0]))
        self.assertAlmostEqual(result["weighted_corrected_ols"][1], 13.75/(15.25-4.265625))
        self.assertAlmostEqual(result["noise_corrected_ols"][0], 40/(14-5.25))

    def test_physical_unit_changes_preserve_a_single_mass(self):
        rng = np.random.default_rng(701)
        force, acceleration = rng.normal(size=(2, 32, 3))
        sf, sa = np.array([.4, 1., 3.]), np.array([2., .3, 1.])
        original = baseline_points(force, acceleration, sf, sa)
        changed = baseline_points(force*7, acceleration*.2, sf*7, sa*.2)
        for name in original:
            np.testing.assert_allclose(changed[name], original[name]*35, rtol=3e-10, atol=1e-10,
                                       equal_nan=True, err_msg=name)

    def test_covariance_profile_coordinate_permutation_and_channel_swap(self):
        rng = np.random.default_rng(7319)
        force, acceleration = rng.normal(size=(2, 40, 3))
        sf, sa = np.array([.2, 1.5, 3.]), np.array([2., .4, 1.])
        direct = covariance_profile(force, acceleration, sf, sa)
        permutation = [2, 0, 1]
        permuted = covariance_profile(force[:, permutation], acceleration[:, permutation],
                                      sf[permutation], sa[permutation])
        reverse = covariance_profile(acceleration, force, sa, sf)
        with np.errstate(divide="ignore", invalid="ignore"):
            np.testing.assert_allclose(permuted, direct, rtol=2e-9, atol=1e-10, equal_nan=True)
            np.testing.assert_allclose(reverse, 1/direct, rtol=2e-9, atol=1e-10, equal_nan=True)

    def test_profile_compares_multiple_interior_minima_globally(self):
        # This profile has minima near 0.102 and 10.816, separated by a
        # maximum. A local optimization started at small mass gives the
        # wrong answer; the larger-mass minimum has the lower objective.
        force, acceleration = np.array([.1, 11., 0]), np.array([1., 1., 0])
        sf, sa = np.array([.1, 10., 1.]), np.ones(3)
        grid = np.geomspace(.001, 1000, 20001)
        values = objective(grid, force, acceleration, sf, sa)
        minima = np.flatnonzero((values[1:-1] < values[:-2]) & (values[1:-1] < values[2:])) + 1
        self.assertEqual(len(minima), 2)
        fitted = covariance_profile(force, acceleration, sf, sa)[0]
        self.assertAlmostEqual(fitted, 10.8162541016, places=7)
        self.assertLess(objective(fitted, force, acceleration, sf, sa), values[minima[0]] - .2)
        self.assertLessEqual(objective(fitted, force, acceleration, sf, sa), values.min() + 1e-12)

    def test_profile_beats_dense_grid_for_varied_covariances(self):
        rng = np.random.default_rng(91625)
        grid = np.r_[0., np.geomspace(1e-7, 1e7, 18001), np.inf]
        for _ in range(60):
            sf, sa = np.exp(rng.uniform(-3, 3, size=(2, 3)))
            force, acceleration = rng.normal(size=(2, 3))*np.exp(rng.uniform(-1, 2, size=(2, 3)))
            fitted = covariance_profile(force, acceleration, sf, sa)[0]
            actual = objective(fitted, force, acceleration, sf, sa)
            reference = np.min(objective(grid, force, acceleration, sf, sa))
            self.assertLessEqual(actual, reference + 2e-9*max(1., reference))

    def test_noise_free_codirectional_pair_recovers_physical_mass(self):
        acceleration = np.array([[1., -2., .4], [0., 3., -1.]])
        force = 3.75*acceleration
        profile = covariance_profile(force, acceleration, [.25, 2, 4], [3, .2, 1])
        np.testing.assert_allclose(profile, 3.75, rtol=1e-12)

    def test_explicit_boundaries_and_flatness(self):
        sf, sa = [.5, 2, 3], [1, .3, 2]
        force = np.array([[0, 0, 0], [1, 2, 3], [0, 0, 0]])
        acceleration = np.array([[1, 2, 3], [0, 0, 0], [0, 0, 0]])
        result = covariance_profile(force, acceleration, sf, sa)
        self.assertEqual(result[0], 0)
        self.assertTrue(np.isposinf(result[1]))
        self.assertTrue(np.isnan(result[2]))

    def test_bad_shapes_and_noise_scales_are_rejected(self):
        for sf, sa in (([1, 2], [1, 1, 1]), ([1, 0, 1], [1, 1, 1]),
                       ([1, 1, 1], [1, -1, 1]), ([np.inf, 1, 1], [1, 1, 1]),
                       ([1, 1, 1], [1, np.nan, 1])):
            with self.subTest(sf=sf, sa=sa):
                with self.assertRaises(ValueError):
                    baseline_points([1, 2, 3], [1, 2, 3], sf, sa)
        with self.assertRaises(ValueError):
            covariance_profile([1, 2], [1, 2], [1, 1, 1], [1, 1, 1])


if __name__ == "__main__":
    unittest.main()

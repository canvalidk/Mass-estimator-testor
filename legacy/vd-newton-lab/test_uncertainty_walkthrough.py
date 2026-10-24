"""Tests whose targets are interpretation and downstream propagation."""
import math
import unittest
import numpy as np
from estimator import estimate, isotropic_covariance, profile, isotropic_qmin
from uncertainty_walkthrough import null_draws, half_cauchy_quantile, positive_interval_product, detached_force_quantile, informative_summary


class UncertaintyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f, cls.a, cls.mass = null_draws(400000)

    def test_same_pair_force_reconstruction_preserves_every_draw(self):
        np.testing.assert_allclose(self.mass * self.a, self.f, rtol=3e-16)

    def test_same_pair_acceleration_reconstruction_preserves_every_draw(self):
        np.testing.assert_allclose(self.f / self.mass, self.a, rtol=3e-16)

    def test_joint_log_covariance_restores_finite_force_variance(self):
        lm, la, lf = np.log(self.mass), np.log(self.a), np.log(self.f)
        self.assertAlmostEqual(lm.var() + la.var() + 2 * np.cov(lm, la, ddof=0)[0, 1], lf.var(), places=12)
        self.assertAlmostEqual(np.corrcoef(lm, la)[0, 1], -1 / math.sqrt(2), delta=.006)
        self.assertAlmostEqual(lf.var(), math.pi ** 2 / 8, delta=.02)

    def test_detaching_mass_from_original_acceleration_changes_prediction(self):
        _, fresh_a, _ = null_draws(len(self.a), seed=20260908)
        detached = self.mass * fresh_a
        variance_ratio = np.var(np.log(detached)) / np.var(np.log(self.f))
        self.assertAlmostEqual(variance_ratio, 3, delta=.04)
        self.assertAlmostEqual(np.quantile(detached, .975) / detached_force_quantile(.975), 1, delta=.03)
        self.assertGreater(detached_force_quantile(.975), 9 * np.quantile(self.f, .975))

    def test_null_point_preserves_marginal_force_law_but_changes_pairing(self):
        point_force = self.a  # Unit scale mass point.
        # Similar marginals do not make the reconstructed joint law correct.
        np.testing.assert_allclose(np.quantile(point_force, [.1, .5, .9]), np.quantile(self.f, [.1, .5, .9]), rtol=.015)
        self.assertAlmostEqual(np.corrcoef(self.f, self.a)[0, 1], 0, delta=.006)
        self.assertAlmostEqual(np.corrcoef(point_force, self.a)[0, 1], 1, places=12)
        self.assertGreater(np.mean((point_force - self.f) ** 2), .7)

    def test_new_exact_positive_acceleration_propagates_mass_quantiles(self):
        samples = 3 * self.mass
        for p in (.025, .5, .975):
            self.assertAlmostEqual(np.quantile(samples, p) / (3 * half_cauchy_quantile(p)), 1, delta=.025)

    def test_new_exact_zero_acceleration_determines_zero_force_despite_mass_uncertainty(self):
        np.testing.assert_array_equal(self.mass * 0, np.zeros_like(self.mass))

    def test_new_exact_zero_force_determines_zero_acceleration_despite_mass_uncertainty(self):
        np.testing.assert_array_equal(0 / self.mass, np.zeros_like(self.mass))

    def test_new_unit_force_has_reciprocal_null_acceleration_law(self):
        acceleration = 1 / self.mass
        for p in (.025, .5, .975):
            self.assertAlmostEqual(np.quantile(acceleration, p) / half_cauchy_quantile(p), 1, delta=.025)

    def test_null_no_finite_mean_does_not_mean_uniform_probability(self):
        # Formula probabilities, not inference about a calibrated physical mass.
        probability_near_scale = 2 / math.pi * (math.atan(2) - math.atan(.5))
        probability_near_1000 = 2 / math.pi * (math.atan(2000) - math.atan(500))
        self.assertGreater(probability_near_scale, .4)
        self.assertLess(probability_near_1000, .001)
        self.assertEqual(profile([0, 0], [0, 0], np.eye(4), 1000), 0)

    def test_zero_true_acceleration_generates_identical_readings_for_different_masses(self):
        noise = np.random.default_rng(191).normal(size=(10, 2, 2))
        trials = [noise + np.array([[m * 0., 0], [0., 0.]]) for m in (.001, 1, 1000)]
        np.testing.assert_array_equal(trials[0], trials[2])
        # Feed the same readings through the actual estimator: no hidden true mass input.
        estimates = [estimate(t[0, 0], t[0, 1], np.eye(4)).mass for t in trials]
        self.assertEqual(estimates[0], estimates[2])

    def test_interval_endpoint_product_is_not_joint_central_interval(self):
        interval = positive_interval_product([half_cauchy_quantile(.025), half_cauchy_quantile(.975)], [.03133798, 2.24140273])
        self.assertLess(interval[0], .002)
        self.assertGreater(interval[1], 50)
        self.assertLess(np.quantile(self.mass * self.a, .975), 2.3)

    def test_ratio_summary_preserves_first_magnitude_moment_for_same_input(self):
        point = self.f.mean() / self.a.mean()
        self.assertAlmostEqual(np.mean(point * self.a), self.f.mean(), places=12)
        self.assertGreater(np.mean((point * self.a - self.f) ** 2), .7)

    def test_one_exact_mass_law_is_safe_to_collapse_for_same_pair(self):
        f = 2 * self.a
        point = f.mean() / self.a.mean()
        np.testing.assert_allclose(point * self.a, f, rtol=1e-15)

    def test_angle_changes_can_reject_without_changing_mass_point(self):
        for angle, reject in ((45, False), (55, True)):
            t = math.radians(angle)
            f, a = [5 * math.cos(t), 5 * math.sin(t)], [5, 0]
            point = estimate(f, a, np.eye(4), direction_order=256).mass
            self.assertAlmostEqual(point, 1, places=10)
            self.assertEqual(isotropic_qmin(f, a, 1, 1) > 9, reject)

    def test_rejected_fit_withholds_the_propagated_prediction(self):
        rejected = informative_summary([-5, 0], [5, 0])
        self.assertEqual(rejected["fit_status"], "rejected_at_9")
        self.assertIsNone(rejected["force_interval_after_global_gate"])
        self.assertIsNotNone(rejected["force_interval_at_exact_acceleration_3_if_model_assumed"])

    def test_moderate_angle_crosses_point_compatibility_without_global_rejection(self):
        p = math.sqrt(6)
        for angle, point_allowed in ((119, True), (121, False)):
            t = math.radians(angle)
            f, a = [p * math.cos(t), p * math.sin(t)], [p, 0]
            self.assertLess(isotropic_qmin(f, a, 1, 1), 9)
            self.assertEqual(profile(f, a, np.eye(4), 1) <= 9, point_allowed)


if __name__ == "__main__":
    unittest.main()

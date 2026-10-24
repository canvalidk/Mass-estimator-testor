"""Behavioral tests with analytic oracles, symmetry, independent integration,
simulation, and deliberate failure cases. Run: python -m unittest discover -v.
"""
import math
import unittest
import numpy as np
from estimator import (estimate, exact_recovery, isotropic_covariance, profile,
                       isotropic_qmin, compatibility_intervals, gate,
                       positive_normal_mean, radial_logs, gauss)


class BehaviorTests(unittest.TestCase):
    def assertRelative(self, got, expected, tolerance=1e-5):
        self.assertLessEqual(abs(got / expected - 1), tolerance, (got, expected))

    def test_exact_branches(self):
        cases = [([0, 0], [0, 0], ("not_identified", None)),
                 ([0, 0], [1, 0], ("contradiction", None)),
                 ([1, 0], [0, 0], ("contradiction", None)),
                 ([-2, 0], [1, 0], ("contradiction", None)),
                 ([0, 2], [1, 0], ("contradiction", None)),
                 ([6, 8], [3, 4], ("mass", 2.)),
                 ([-6, -8], [-3, -4], ("mass", 2.))]
        for f, a, wanted in cases:
            with self.subTest(f=f, a=a):
                self.assertEqual(exact_recovery(f, a), wanted)

    def test_zero_covariance_is_not_a_zero_reading(self):
        with self.assertRaises(ValueError):
            estimate([0, 0], [0, 0], np.zeros((4, 4)))
        self.assertRelative(estimate([0, 0], [0, 0], np.eye(4)).mass, 1)

    def test_null_in_one_two_three_dimensions(self):
        for d in (1, 2, 3):
            with self.subTest(d=d):
                e = estimate(np.zeros(d), np.zeros(d), isotropic_covariance(d, 3, 2), direction_order=16)
                self.assertRelative(e.mass, 1.5)
                self.assertRelative(e.mean_force, 3 * math.sqrt(2 / math.pi), 1e-6)
                for p in (.025, .1, .5, .9, .975):
                    self.assertRelative(e.quantile(p), 1.5 * math.tan(math.pi * p / 2), 1e-10)

    def test_null_log_density_and_finite_log_variance(self):
        e = estimate([0, 0], [0, 0], np.eye(4))
        ell, density = e.log_density()
        np.testing.assert_allclose(density, 1 / (math.pi * np.cosh(ell)), rtol=1e-12)
        x, w = gauss(400)
        x, w = 24 * x, 24 * w
        self.assertRelative(float(w @ (x * x / (math.pi * np.cosh(x)))), math.pi ** 2 / 4, 2e-8)

    def test_null_common_precision_scaling(self):
        for factor in (1e-6, .1, 1, 10, 1e6):
            e = estimate([0, 0], [0, 0], isotropic_covariance(2, 3 * factor, factor))
            self.assertRelative(e.mass, 3)
            self.assertRelative(e.quantile(.975), 3 * math.tan(.975 * math.pi / 2), 1e-10)

    def test_null_apparatus_ratio_changes_the_readout(self):
        self.assertRelative(estimate([0, 0], [0, 0], isotropic_covariance(2, 10, 1)).mass, 10)
        self.assertRelative(estimate([0, 0], [0, 0], isotropic_covariance(2, 1, 10)).mass, .1)

    def test_null_profile_all_masses_even_with_correlation(self):
        cov = np.array([[2, .1, .5, 0], [.1, 1, 0, .2], [.5, 0, 1, .1], [0, .2, .1, 3]])
        for mass in (1e-12, .01, 1, 100, 1e12):
            self.assertEqual(profile([0, 0], [0, 0], cov, mass), 0)
        self.assertEqual(compatibility_intervals([0, 0], [0, 0], 2, 1, 9), [[0., math.inf]])

    def test_known_direction_independent_truncated_normal_oracle(self):
        for f, a, sf, sa in [(0, 0, 2, 1), (0, 10, 1, 1), (10, 0, 1, 1),
                             (3, 5, 1, 2), (-5, 5, 1, 1), (-20, -30, 2, 3)]:
            with self.subTest(f=f, a=a):
                e = estimate([f], [a], isotropic_covariance(1, sf, sa), known_direction=[1], ratio_order=2048)
                self.assertRelative(e.mass, positive_normal_mean(f, sf) / positive_normal_mean(a, sa), 5e-5)

    def test_radial_integral_independent_positive_axis_quadrature(self):
        x, w = gauss(500)
        for t in (-50, -5, -4, 0, 5, 30):
            with self.subTest(t=t):
                upper = max(12., t + 12.)
                r, weights = (x + 1) * upper / 2, w * upper / 2
                kernel = np.exp(-.5 * r * r + t * r - max(t, 0) ** 2 / 2)
                logs = radial_logs(np.array([t]))
                for power in (1, 2):
                    oracle = math.log(float(weights @ (r ** power * kernel))) + max(t, 0) ** 2 / 2
                    self.assertAlmostEqual(logs[power - 1][0], oracle, delta=2e-9)

    def test_resolved_one_zero_limits(self):
        c = math.sqrt(2 / math.pi)
        for d in (1, 2, 3):
            f, a = np.zeros(d), np.zeros(d)
            a[0] = 20
            e = estimate(f, a, np.eye(2 * d), direction_order=48, ratio_order=2048)
            self.assertRelative(e.mass, c / 20, .005)
            swapped = estimate(a, f, np.eye(2 * d), direction_order=48, ratio_order=2048)
            self.assertRelative(swapped.mass, 20 / c, .005)

    def test_one_zero_compatibility_bounds_and_unresolved_exception(self):
        upper = 3 / math.sqrt(100 - 9)
        self.assertRelative(compatibility_intervals([0, 0], [10, 0], 1, 1, 9)[0][1], upper)
        self.assertRelative(compatibility_intervals([10, 0], [0, 0], 1, 1, 9)[0][0], 1 / upper)
        self.assertEqual(compatibility_intervals([0, 0], [2, 0], 1, 1, 9), [[0, math.inf]])

    def test_strong_aligned_limit(self):
        e = estimate([20, 0], [10, 0], np.eye(4), direction_order=256, ratio_order=2048)
        self.assertRelative(e.mass, 2, .002)
        self.assertLess(e.quantile(.975) - e.quantile(.025), 1)
        self.assertEqual(isotropic_qmin([20, 0], [10, 0], 1, 1), 0)

    def test_all_angles_same_point_under_exchange_symmetry(self):
        for degrees in (0, 15, 45, 90, 120, 170, 180):
            theta = math.radians(degrees)
            e = estimate([6 * math.cos(theta), 6 * math.sin(theta)], [3, 0], isotropic_covariance(2, 2, 1), direction_order=256)
            self.assertRelative(e.mass, 2, 1e-9)

    def test_anti_aligned_gate_dominates_conditional_point(self):
        e = estimate([-10, 0], [10, 0], isotropic_covariance(2, 2, 2))
        q = isotropic_qmin([-10, 0], [10, 0], 2, 2)
        self.assertEqual(q, 25)
        self.assertRelative(e.mass, 1)
        self.assertEqual(gate(e.mass, q, 9)["returned_mass"], None)

    def test_balanced_anti_alignment_can_pass_global_gate_but_fail_at_point(self):
        p = math.sqrt(6)
        f, a = [-p, 0], [p, 0]
        e = estimate(f, a, np.eye(4))
        self.assertLess(isotropic_qmin(f, a, 1, 1), 9)
        self.assertGreater(profile(f, a, np.eye(4), e.mass), 9)
        sets = compatibility_intervals(f, a, 1, 1, 9)
        self.assertEqual(len(sets), 2)
        self.assertRelative(sets[0][1], 2 - math.sqrt(3))
        self.assertRelative(sets[1][0], 2 + math.sqrt(3))
        self.assertLess(sets[0][1], e.mass)
        self.assertGreater(sets[1][0], e.mass)

    def test_balanced_perpendicular_profile_flat_but_nonzero(self):
        for mass in (.001, .1, 1, 10, 1000):
            self.assertAlmostEqual(profile([0, 2], [2, 0], np.eye(4), mass), 4)
        self.assertEqual(compatibility_intervals([0, 2], [2, 0], 1, 1, 9), [[0., math.inf]])
        self.assertEqual(compatibility_intervals([0, 4], [4, 0], 1, 1, 9), [])

    def test_boundary_infimum_not_stationary_maximum(self):
        self.assertEqual(isotropic_qmin([-10, 0], [10, 0], 2, 2), 25)
        self.assertEqual(profile([-10, 0], [10, 0], isotropic_covariance(2, 2, 2), 1), 50)
        self.assertEqual(isotropic_qmin([0, 0], [10, 0], 1, 1), 0)

    def test_quadratic_compatibility_matches_direct_profile(self):
        rng = np.random.default_rng(240906)
        for _ in range(150):
            f, a = rng.normal(size=(2, 2)) * 3
            sf, sa = np.exp(rng.normal(size=2))
            threshold = float(rng.uniform(.1, 15))
            intervals = compatibility_intervals(f, a, sf, sa, threshold)
            cov = isotropic_covariance(2, sf, sa)
            for mass in np.geomspace(.001, 1000, 25):
                self.assertEqual(any(lo <= mass <= hi for lo, hi in intervals), profile(f, a, cov, mass) <= threshold)

    def test_tangent_singleton_and_threshold_equality(self):
        self.assertEqual(compatibility_intervals([2, 0], [1, 0], 1, 1, 0), [[2., 2.]])
        self.assertEqual(gate(1, 9, 9)["returned_mass"], 1)

    def test_channel_swap_reciprocates_point_and_interval(self):
        cov = np.array([[4, .2, .3, .1], [.2, 1, -.1, .2], [.3, -.1, 1, .1], [.1, .2, .1, 2]])
        f, a = [4, 2], [1, -.5]
        e = estimate(f, a, cov)
        ix = [2, 3, 0, 1]
        swapped = estimate(a, f, cov[np.ix_(ix, ix)])
        self.assertRelative(e.mass * swapped.mass, 1, 1e-10)
        self.assertRelative(e.quantile(.025) * swapped.quantile(.975), 1, 1e-8)

    def test_units_rescale_readout_not_information(self):
        f, a, cov = np.array([3., 2]), np.array([2., -.5]), np.eye(4)
        e = estimate(f, a, cov)
        sf, sa = 1000., .01
        scales = np.diag([sf, sf, sa, sa])
        rescaled = estimate(f * sf, a * sa, scales @ cov @ scales)
        self.assertRelative(rescaled.mass, sf / sa * e.mass, 1e-10)
        self.assertRelative(rescaled.quantile(.975), sf / sa * e.quantile(.975), 1e-9)
        self.assertRelative(profile(f * sf, a * sa, scales @ cov @ scales, e.mass * sf / sa), profile(f, a, cov, e.mass), 1e-10)

    def test_rotation_with_covariance_is_invariant(self):
        angle = .431
        r = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
        big_r = np.kron(np.eye(2), r)
        f, a = np.array([4., 1]), np.array([1., 2])
        cov = np.diag([.5, 3, 2, 1])
        e = estimate(f, a, cov, direction_order=512)
        rotated = estimate(r @ f, r @ a, big_r @ cov @ big_r.T, direction_order=512)
        self.assertRelative(e.mass, rotated.mass, 1e-9)
        self.assertRelative(profile(f, a, cov, 2), profile(r @ f, r @ a, big_r @ cov @ big_r.T, 2), 1e-12)

    def test_rotating_only_readings_changes_anisotropic_problem(self):
        cov = np.diag([.2, 5, 1, 1])
        x = estimate([3, 0], [1, 0], cov)
        y = estimate([0, 3], [0, 1], cov)
        self.assertGreater(abs(x.mass - y.mass), .1)

    def test_correlated_fixed_direction_matches_independent_monte_carlo(self):
        rng = np.random.default_rng(12421)
        cov = np.array([[1., .7], [.7, 1.]])
        data = rng.multivariate_normal([.5, 1.5], cov, size=500000)
        data = data[(data > 0).all(axis=1)]
        e = estimate([.5], [1.5], cov, known_direction=[1])
        self.assertRelative(e.mass, float(data[:, 0].mean() / data[:, 1].mean()), .006)
        for p in (.025, .5, .975):
            self.assertRelative(e.quantile(p), float(np.quantile(data[:, 0] / data[:, 1], p)), .04)

    def test_null_correlations_change_spread_without_changing_symmetric_point(self):
        independent = estimate([0], [0], np.eye(2), known_direction=[1])
        correlated = estimate([0], [0], [[1, .8], [.8, 1]], known_direction=[1])
        self.assertRelative(correlated.mass, 1)
        self.assertLess(correlated.quantile(.975), independent.quantile(.975))

    def test_known_direction_and_unknown_direction_are_different_inputs(self):
        cov = np.eye(4)
        known = estimate([0, 4], [2, 0], cov, known_direction=[1, 0])
        unknown = estimate([0, 4], [2, 0], cov)
        self.assertGreater(abs(known.mass - unknown.mass), .5)

    def test_numerical_refinement_extreme_and_ordinary_cases(self):
        for f, a in [([3, 1], [2, 0]), ([-5, 0], [5, 0]), ([0, 0], [20, 0]), ([20, 0], [10, 0])]:
            with self.subTest(f=f, a=a):
                coarse = estimate(f, a, np.eye(4), direction_order=256, ratio_order=1024)
                fine = estimate(f, a, np.eye(4), direction_order=512, ratio_order=2048)
                self.assertRelative(coarse.mass, fine.mass, 2e-4)
                self.assertRelative(coarse.quantile(.975), fine.quantile(.975), .004)

    def test_three_dimensional_rotation_and_refinement(self):
        f, a = [4, 1, 2], [1, 2, 1]
        e = estimate(f, a, np.eye(6), direction_order=24)
        rotated = estimate([2, 4, 1], [1, 1, 2], np.eye(6), direction_order=40)
        self.assertRelative(e.mass, rotated.mass, 1e-6)

    def test_invalid_inputs_fail_explicitly(self):
        for f, a, cov in [([1], [1, 0], np.eye(2)), ([math.nan], [1], np.eye(2)),
                          ([1], [1], [[1, 2], [2, 1]]), ([1], [1], [[1, .5], [0, 1]])]:
            with self.subTest(f=f, a=a, cov=cov), self.assertRaises(ValueError):
                estimate(f, a, cov)
        for s in (0, -1, math.inf):
            with self.assertRaises(ValueError):
                isotropic_covariance(2, s, 1)

    def test_ratio_of_means_is_not_mean_of_ratios(self):
        f, a = np.array([1., 100]), np.array([1., 10])
        self.assertRelative(f.mean() / a.mean(), 101 / 11)
        self.assertEqual(float((f / a).mean()), 5.5)
        self.assertNotAlmostEqual(float(f.mean() / a.mean()), float((f / a).mean()))

    def test_composition_and_same_law_regrouping(self):
        rng = np.random.default_rng(121)
        for _ in range(100):
            f1, f2, a1, a2 = np.exp(rng.normal(size=(4, 20)))
            p = rng.dirichlet(np.ones(20))
            r = lambda f, a: float(p @ f / (p @ a))
            self.assertAlmostEqual(r(f1 + f2, a1), r(f1, a1) + r(f2, a1))
            self.assertAlmostEqual(1 / r(f1, a1 + a2), 1 / r(f1, a1) + 1 / r(f1, a2))
            grouped = sum(p[i:i + 5] @ f1[i:i + 5] for i in range(0, 20, 5)) / sum(p[i:i + 5] @ a1[i:i + 5] for i in range(0, 20, 5))
            self.assertAlmostEqual(r(f1, a1), grouped)

    def test_same_null_repetitions_do_not_create_precision(self):
        for n in (1, 10, 10000):
            e = estimate([0, 0], [0, 0], isotropic_covariance(2, 2 / math.sqrt(n), 1 / math.sqrt(n)))
            self.assertRelative(e.quantile(.975) / e.quantile(.025), math.tan(.975 * math.pi / 2) / math.tan(.025 * math.pi / 2), 1e-10)

    def test_common_mass_detects_disagreeing_trials(self):
        cov = isotropic_covariance(2, .1, .1)
        self.assertEqual(isotropic_qmin([2, 0], [1, 0], .1, .1), 0)
        self.assertEqual(isotropic_qmin([8, 0], [1, 0], .1, .1), 0)
        scores = [profile([2, 0], [1, 0], cov, m) + profile([8, 0], [1, 0], cov, m) for m in np.geomspace(.01, 1000, 2000)]
        self.assertGreater(min(scores), 40)

    def test_null_ratio_tail_does_not_have_a_finite_mean(self):
        # Truncated integrals grow with the truncation. This is a regression
        # witness; the unbounded conclusion follows from the source's tail fact.
        truncated_mean = lambda bound: math.log1p(bound * bound) / math.pi
        self.assertGreater(truncated_mean(1e9), truncated_mean(1e6) + 4)
        self.assertGreater(truncated_mean(1e6), truncated_mean(1e3) + 4)

    def test_informative_readings_still_have_nonzero_boundary_density(self):
        e = estimate([8], [8], np.eye(2), known_direction=[1], ratio_order=1024)
        self.assertGreater(e.angle_density[-1], 0)
        self.assertLess(e.angle_density[-1], 1e-10)
        self.assertLess(e.quantile(.975) / e.quantile(.025), 3)

    def test_large_fixed_mass_avoids_overflow(self):
        self.assertAlmostEqual(profile([2, 0], [1, 0], np.eye(4), 1e300), 1)

    def test_full_vector_against_direct_magnitude_direction_integral(self):
        # Independent coordinates and integration: f, alpha, direction instead
        # of the implementation's ratio angle and analytic radial integral.
        f, a = np.array([3., 1]), np.array([1., 2])
        cov = np.array([[1, .2, .3, 0], [.2, 2, 0, .1], [.3, 0, 1, .2], [0, .1, .2, 2]])
        precision = np.linalg.inv(cov)
        x, w = gauss(72)
        x, w = 8 * (x + 1), 8 * w
        ff, aa = np.meshgrid(x, x, indexing="ij")
        quadrature = np.outer(w, w)
        totals = np.zeros(3)
        for angle in np.arange(160) * (2 * math.pi / 160):
            u = np.array([math.cos(angle), math.sin(angle)])
            residual = np.r_[f, a] - np.concatenate((ff[..., None] * u, aa[..., None] * u), axis=2)
            probability = np.exp(-.5 * np.einsum("...i,ij,...j->...", residual, precision, residual)) * quadrature
            totals += [probability.sum(), (probability * ff).sum(), (probability * aa).sum()]
        e = estimate(f, a, cov, direction_order=256, ratio_order=1024)
        self.assertRelative(e.mass, totals[1] / totals[2], 1e-7)
        self.assertRelative(e.mean_force, totals[1] / totals[0], 1e-7)

    def test_profile_against_direct_generalized_least_squares(self):
        cov = np.array([[1, .2, .3, 0], [.2, 2, 0, .1], [.3, 0, 1, .2], [0, .1, .2, 2]])
        precision, y = np.linalg.inv(cov), np.array([3., 1, 1, 2])
        for mass in (.01, .5, 2, 100):
            design = np.vstack((mass * np.eye(2), np.eye(2)))
            latent = np.linalg.solve(design.T @ precision @ design, design.T @ precision @ y)
            residual = y - design @ latent
            self.assertRelative(profile(y[:2], y[2:], cov, mass), residual @ precision @ residual, 1e-10)

    def test_fixed_mass_can_reject_one_zero_that_free_mass_accepts(self):
        self.assertEqual(isotropic_qmin([0, 0], [10, 0], 1, 1), 0)
        self.assertEqual(profile([0, 0], [10, 0], np.eye(4), 1), 50)

    def test_approach_to_origin_depends_on_information_path(self):
        # Fixed uncertainties lose excitation and approach apparatus scale 1.
        e = estimate([.0002, 0], [.0001, 0], np.eye(4))
        self.assertRelative(e.mass, 1, 1e-6)
        # Shrink readings AND uncertainties: keep SNR and the complete ratio law.
        original = estimate([20, 0], [10, 0], np.eye(4), ratio_order=1024, direction_order=256)
        shrunk = estimate([.002, 0], [.001, 0], isotropic_covariance(2, .0001, .0001), ratio_order=1024, direction_order=256)
        self.assertRelative(original.mass, shrunk.mass, 1e-10)
        self.assertRelative(original.quantile(.975), shrunk.quantile(.975), 1e-9)


if __name__ == "__main__":
    unittest.main()

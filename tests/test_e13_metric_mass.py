"""Tests for E13 Push 6 (Insights 11–13): GL(d)-covariant Metric Mass g in Sym+_d.

Pure stdlib unittest (also pytest-discoverable) so it runs without external packages.
"""
import math
import random
import unittest

from explore.e13_general_clr import swap_reading
from explore.e13_metric_mass import (
    MetricReading,
    _make_metric_world,
    compose_metric_common_force,
    compose_metric_lump,
    eye,
    mat_inv_spd,
    mat_vec,
    metric_clr_pvalue,
    nested_surprise_ladder,
    t50_metric_score_test,
    test_scalar_isotropy_l22,
    transform_gl,
    transform_metric_gl,
)
from explore.e13_unnormalized_evidence import Reading, isotropic


class TestE13MetricMass(unittest.TestCase):
    def test_covector_residual_matches_qr_and_is_gl_d_and_swap_invariant(self):
        # Correlated 2D reading with cross-channel covariance
        F2 = [3.1, -1.4]
        a2 = [1.2, 0.7]
        cov2 = [
            [0.25, 0.04, 0.02, -0.01],
            [0.04, 0.18, -0.015, 0.02],
            [0.02, -0.015, 0.11, 0.01],
            [-0.01, 0.02, 0.01, 0.14],
        ]
        r2 = Reading(F2, a2, cov2)
        mr2 = MetricReading(r2)
        g2 = [[2.2, 0.4], [0.4, 3.5]]
        k_fast = mr2.residual_cost(g2)
        _, _, _, _, h_qr, k_qr, _ = mr2.split_at_matrix(g2)
        self.assertAlmostEqual(k_fast, k_qr, places=11)

        # 1. Exact GL(2) invariance (non-orthogonal shear + unequal axis scaling; NO Euclidean delta!)
        A2 = [[2.3, 0.7], [-0.4, 1.5]]
        r2_gl = transform_gl(r2, A2)
        g2_gl = transform_metric_gl(g2, A2)
        mr2_gl = MetricReading(r2_gl)
        self.assertAlmostEqual(mr2_gl.residual_cost(g2_gl), k_fast, places=11)
        _, _, _, _, h_gl, k_gl, _ = mr2_gl.split_at_matrix(g2_gl)
        self.assertAlmostEqual(h_gl, h_qr, places=11)
        self.assertAlmostEqual(k_gl, k_qr, places=11)
        # T50 score test is also GL(2)-invariant on the symmetric tangent space
        t_base = t50_metric_score_test([mr2], g2)["T2"]
        t_gl = t50_metric_score_test([mr2_gl], g2_gl)["T2"]
        self.assertAlmostEqual(t_gl, t_base, places=10)

        # 2. Exact force-acceleration swap duality: g <-> g^{-1}
        mr2_sw = MetricReading(swap_reading(r2))
        g2_inv = mat_inv_spd(g2)
        self.assertAlmostEqual(mr2_sw.residual_cost(g2_inv), k_fast, places=11)
        self.assertAlmostEqual(t50_metric_score_test([mr2_sw], g2_inv)["T2"], t_base, places=10)

        # 3. Correlated 3D reading + GL(3) invariance (sheared 3x3 basis change!)
        cov3 = [[0.0] * 6 for _ in range(6)]
        for i in range(6):
            cov3[i][i] = 0.2 + 0.03 * i
            if i < 3:
                cov3[i][i + 3] = cov3[i + 3][i] = 0.04
        cov3[0][1] = cov3[1][0] = 0.02
        r3 = Reading([2.5, -1.0, 4.2], [1.1, 0.4, -0.8], cov3)
        g3 = [[2.0, 0.3, -0.1], [0.3, 4.0, 0.2], [-0.1, 0.2, 6.0]]
        mr3 = MetricReading(r3)
        k3 = mr3.residual_cost(g3)
        self.assertAlmostEqual(k3, mr3.split_at_matrix(g3)[5], places=11)
        A3 = [[1.5, 0.4, -0.2], [0.1, 2.0, 0.5], [-0.3, 0.2, 0.8]]
        mr3_gl = MetricReading(transform_gl(r3, A3))
        g3_gl = transform_metric_gl(g3, A3)
        self.assertAlmostEqual(mr3_gl.residual_cost(g3_gl), k3, places=10)

    def test_partial_axis_sigma_infinity_and_additive_combination(self):
        # 2D x-only push: K_x(diag(2, g_yy)) == 0 for ALL g_yy > 0 (sigma = inf in unpushed direction!)
        mr_x = MetricReading(isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.1))
        mr_y = MetricReading(isotropic([0.0, 5.0], [0.0, 1.0], 0.1, 0.1))
        for gyy in (1e-3, 0.5, 5.0, 1e4):
            self.assertAlmostEqual(mr_x.residual_cost([[2.0, 0.0], [0.0, gyy]]), 0.0, places=12)
        # Adding K_x + K_y pins down diag(2, 5) with zero prior
        lad = nested_surprise_ladder([mr_x, mr_y])
        self.assertAlmostEqual(lad["g_hat"][0][0], 2.0, places=2)
        self.assertAlmostEqual(lad["g_hat"][1][1], 5.0, places=2)
        self.assertAlmostEqual(lad["g_hat"][0][1], 0.0, places=2)
        self.assertLess(lad["K_metric"], 1e-4)
        self.assertGreater(lad["K_scalar"], 20.0)

        # 3D x-only push stays identically 0 across arbitrary g_yy, g_zz
        mr3_x = MetricReading(isotropic([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 0.1))
        self.assertAlmostEqual(mr3_x.residual_cost([[2.0, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 500.0]]), 0.0, places=12)

    def test_three_rung_surprise_ladder_and_t50_full_rank_hilltop_ghost(self):
        def make_clean(M):
            return [
                MetricReading(isotropic(mat_vec(M, [math.cos(a), math.sin(a)]), [math.cos(a), math.sin(a)], 0.1, 0.1))
                for a in (0.0, 0.6, 1.2, 2.1)
            ]

        # 1. Isotropic -> all 3 rungs ~ 0
        l_iso = nested_surprise_ladder(make_clean([[2.0, 0.0], [0.0, 2.0]]))
        self.assertLess(l_iso["K_scalar"], 1e-4)
        self.assertLess(l_iso["K_metric"], 1e-4)
        self.assertLess(l_iso["K_linear"], 1e-4)

        # 2. Anisotropic metric -> scalar fails, metric & linear succeed
        l_ani = nested_surprise_ladder(make_clean([[2.0, 0.0], [0.0, 5.0]]))
        self.assertGreater(l_ani["delta_iso"], 40.0)
        self.assertLess(l_ani["K_metric"], 1e-3)

        # 3. Gyroscopic non-reciprocal M = [[2, -1], [1, 2]] (E7) -> metric fails, linear succeeds!
        l_gyr = nested_surprise_ladder(make_clean([[2.0, -1.0], [1.0, 2.0]]))
        self.assertGreater(l_gyr["delta_recip"], 40.0)
        self.assertLess(l_gyr["K_linear"], 1e-3)

        # 4. Full-rank (rank=3) 2D anti-parallel hilltop ghost at g0 = I_2:
        # T50 accepts with T2=0, rank=3, p=1; Metric-GCLR & Z-safe Metric-Z-GCLR both reject!
        mrs_anti = [
            MetricReading(isotropic([-2.0, 0.0], [1.0, 0.0], 0.1, 0.1)),
            MetricReading(isotropic([-1.0, 0.0], [2.0, 0.0], 0.1, 0.1)),
            MetricReading(isotropic([0.0, -2.0], [0.0, 1.0], 0.1, 0.1)),
            MetricReading(isotropic([0.0, -1.0], [0.0, 2.0], 0.1, 0.1)),
            MetricReading(isotropic([-2.0, -2.0], [1.0, 1.0], 0.1, 0.1)),
            MetricReading(isotropic([-1.0, -1.0], [2.0, 2.0], 0.1, 0.1)),
        ]
        t50 = t50_metric_score_test(mrs_anti, eye(2))
        self.assertEqual(t50["rank"], 3)
        self.assertAlmostEqual(t50["T2"], 0.0, places=10)
        self.assertAlmostEqual(t50["p_value"], 1.0, places=10)
        p_gclr = metric_clr_pvalue(mrs_anti, eye(2), family="sym", weight_mode="unweighted", n_draws=49, opt_iter=20, seed=5)["p_value"]
        p_zgclr = metric_clr_pvalue(mrs_anti, eye(2), family="sym", weight_mode="z_smooth", n_draws=49, opt_iter=20, seed=5)["p_value"]
        self.assertLessEqual(p_gclr, 0.04)
        self.assertLessEqual(p_zgclr, 0.04)

    def test_z_invariance_partial_rank_and_heavy_tail_coverage(self):
        base = [
            MetricReading(isotropic([2.1, 0.2], [1.0, 0.1], 0.5, 0.5)),
            MetricReading(isotropic([0.1, 4.8], [0.0, 1.0], 0.5, 0.5)),
        ]
        with_zeros = base + [MetricReading(isotropic([0.0, 0.0], [0.0, 0.0], 1.0, 1.0))] * 5
        g0 = [[2.0, 0.0], [0.0, 5.0]]

        # Bit-for-bit Z invariance for T50, Metric-Z-GCLR+, and Metric-SCLR+
        self.assertEqual(t50_metric_score_test(base, g0)["T2"], t50_metric_score_test(with_zeros, g0)["T2"])
        self.assertEqual(
            metric_clr_pvalue(base, g0, weight_mode="z_smooth", null_mode="gaussian", n_draws=39, opt_iter=15, seed=9)["p_value"],
            metric_clr_pvalue(with_zeros, g0, weight_mode="z_smooth", null_mode="gaussian", n_draws=39, opt_iter=15, seed=9)["p_value"],
        )
        self.assertEqual(
            metric_clr_pvalue(base, g0, weight_mode="unweighted", null_mode="spherical", n_draws=39, opt_iter=15, seed=9)["p_value"],
            metric_clr_pvalue(with_zeros, g0, weight_mode="unweighted", null_mode="spherical", n_draws=39, opt_iter=15, seed=9)["p_value"],
        )

        # Heavy-tail Cauchy-like nu=2 on 2D anisotropic g=diag(2,5): T50 & Gaussian collapse, SCLR+ holds ~95%!
        rng = random.Random(2029)
        reps = 80
        c_t50 = c_g = c_s = 0
        for rep in range(reps):
            rds = _make_metric_world(rng, g0, 8, 3.0, sigma=1.0, nu_df=2)
            mrs = [MetricReading(r) for r in rds]
            if t50_metric_score_test(mrs, g0)["p_value"] > 0.05:
                c_t50 += 1
            if metric_clr_pvalue(mrs, g0, family="metric", weight_mode="z_smooth", null_mode="gaussian", n_draws=39, opt_iter=15, seed=rep * 13 + 1)["p_value"] > 0.05:
                c_g += 1
            if metric_clr_pvalue(mrs, g0, family="metric", weight_mode="unweighted", null_mode="spherical", n_draws=39, opt_iter=15, seed=rep * 13 + 2)["p_value"] > 0.05:
                c_s += 1
        self.assertLess(c_t50 / reps, 0.70)
        self.assertLess(c_g / reps, 0.70)
        self.assertGreater(c_s / reps, 0.87)

    def test_l22_scalar_isotropy_and_simultaneous_l5_l14_composition(self):
        # 1. L22 scalar isotropy test: accepts isotropic & x-only; rejects anisotropic diag(2, 5)
        rng = random.Random(77)
        mrs_iso = [MetricReading(r) for r in _make_metric_world(rng, eye(2, 2.0), 8, 15.0, sigma=0.3)]
        mrs_ani = [MetricReading(r) for r in _make_metric_world(rng, [[2.0, 0.0], [0.0, 5.0]], 8, 15.0, sigma=0.3)]
        mrs_x = [MetricReading(isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.1))] * 3
        self.assertGreater(test_scalar_isotropy_l22(mrs_iso, n_draws=49, seed=3)["p_value"], 0.05)
        self.assertLessEqual(test_scalar_isotropy_l22(mrs_ani, n_draws=49, seed=3)["p_value"], 0.05)
        self.assertEqual(test_scalar_isotropy_l22(mrs_x, n_draws=49, seed=3)["p_value"], 1.0)

        # 2. Simultaneous exact L5 lumping and L14 common force on anisotropic metrics (resolving O17!)
        rds1 = [isotropic([2.0, 0.0], [1.0, 0.0], 0.1, 0.01), isotropic([0.0, 5.0], [0.0, 1.0], 0.1, 0.01)]
        rds2 = [isotropic([3.0, 0.0], [1.0, 0.0], 0.1, 0.01), isotropic([0.0, 1.0], [0.0, 1.0], 0.1, 0.01)]
        mrs1, mrs2 = [MetricReading(r) for r in rds1], [MetricReading(r) for r in rds2]

        self.assertLess(compose_metric_lump(mrs1, mrs2, [[5.0, 0.0], [0.0, 6.0]])["K_lump"], 1e-4)
        self.assertGreater(compose_metric_lump(mrs1, mrs2, [[5.4, 0.0], [0.0, 6.0]])["K_lump"], 4.0)
        self.assertLess(compose_metric_common_force(rds1, rds2, [[1.2, 0.0], [0.0, 5.0 / 6.0]])["K_red"], 1e-4)
        self.assertGreater(compose_metric_common_force(rds1, rds2, [[1.5, 0.0], [0.0, 5.0 / 6.0]])["K_red"], 4.0)

        # 3. Partial-axis gluing: diag(2, 5) + (g2_xx=3, g2_yy unmeasured):
        # rules out G_yy < 5 while leaving G_yy > 5 flat!
        mrs2_x = [MetricReading(isotropic([3.0, 0.0], [1.0, 0.0], 0.1, 0.01))]
        k_low = compose_metric_lump(mrs1, mrs2_x, [[5.0, 0.0], [0.0, 4.5]])["K_lump"]
        k_hi = compose_metric_lump(mrs1, mrs2_x, [[5.0, 0.0], [0.0, 15.0]])["K_lump"]
        self.assertGreater(k_low, 10.0)
        self.assertLess(k_hi, 1e-4)


if __name__ == "__main__":
    unittest.main()

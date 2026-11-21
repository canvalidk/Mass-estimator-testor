"""Tests for E13: unnormalized evidence ('bell with sigma = inf') prototype.

Pure stdlib unittest (also pytest-discoverable) so it runs without external packages.
"""
import math
import random
import unittest

from explore.e13_unnormalized_evidence import (
    CHI2_95_1,
    Reading,
    analyze_case,
    chi2_sf,
    compose_common_force,
    compose_lump,
    dot,
    evaluate_series,
    isotropic,
    moreira_clr_pvalue,
)


class TestE13UnnormalizedEvidence(unittest.TestCase):
    def test_chi2_sf_exact_values(self):
        self.assertAlmostEqual(chi2_sf(CHI2_95_1, 1), 0.05, places=12)
        self.assertAlmostEqual(chi2_sf(4.0, 2), math.exp(-2.0), places=12)
        self.assertAlmostEqual(chi2_sf(6.0, 4), 4.0 * math.exp(-3.0), places=12)
        self.assertAlmostEqual(chi2_sf(10.0, 6), 18.5 * math.exp(-5.0), places=11)

    def test_moreira_clr_pvalue_limits_and_monte_carlo(self):
        # Strong-signal limit -> chi^2_1; zero-H limit -> chi^2_M
        self.assertAlmostEqual(moreira_clr_pvalue(CHI2_95_1, 1e8, 2), 0.05, places=5)
        self.assertAlmostEqual(moreira_clr_pvalue(5.0, 1e-8, 6), chi2_sf(5.0, 6), places=4)
        # Finite-H check against direct Monte Carlo (df=6, H=4.0, LR=3.5)
        rng = random.Random(42)
        h_test, df_test, c_test, n_mc = 4.0, 6, 3.5, 25000
        exceed = 0
        for _ in range(n_mc):
            t = rng.gauss(0.0, 1.0)
            w = sum(rng.gauss(0.0, 1.0) ** 2 for _ in range(df_test - 1))
            k_val = t * t + w
            lr = 0.5 * (k_val - h_test + math.hypot(k_val - h_test, 2.0 * math.sqrt(h_test) * t))
            if lr >= c_test:
                exceed += 1
        self.assertAlmostEqual(exceed / n_mc, moreira_clr_pvalue(c_test, h_test, df_test), delta=0.01)

    def test_all_zero_is_flat_sigma_infinity_and_combines_as_exact_identity(self):
        r_zero = [isotropic([0, 0, 0], [0, 0, 0], 2.0, 1.0)]
        res = analyze_case(r_zero, common_ratio=2.0, scale_hint=2.0, n_grid=400)
        self.assertIsNone(res["m_hat"])
        self.assertIsNone(res["door_b_hat"])
        self.assertEqual(res["K_min"], 0.0)
        self.assertEqual(res["surprise_p"], 1.0)
        self.assertEqual(res["t45_set"], "(0, inf)")
        self.assertEqual(res["clr_set"], "(0, inf)")
        self.assertEqual(res["door_b_set"], "(0, inf)")

        # Combining stationary 0/0 with density 5.0 +- 0.5 kg leaves density untouched
        K_density = lambda m: ((m - 5.0) / 0.5) ** 2
        for m in (3.0, 4.02, 5.0, 5.98, 7.0):
            k_zero = evaluate_series(r_zero, math.atan(m))["K"]
            self.assertEqual(k_zero, 0.0)
            self.assertEqual(K_density(m) + k_zero, K_density(m))

    def test_units_swap_and_rotation_invariance_with_correlated_covariance(self):
        # 2D reading with cross-channel correlation and nontrivial subspace rotation
        F = [3.2, -1.1]
        a = [1.4, 0.5]
        cov = [
            [0.25, 0.03, 0.02, -0.01],
            [0.03, 0.16, -0.01, 0.015],
            [0.02, -0.01, 0.09, 0.005],
            [-0.01, 0.015, 0.005, 0.12],
        ]
        r = Reading(F, a, cov)
        m = 2.3
        e_base = evaluate_series([r], math.atan(m))

        # 1. Units scaling: F -> lam*F, a -> mu*a => same K, T^2 at m' = (lam/mu)*m
        lam, mu = 3.5, 0.7
        scales = [lam, lam, mu, mu]
        cov_u = [[cov[i][j] * scales[i] * scales[j] for j in range(4)] for i in range(4)]
        r_u = Reading([lam * x for x in F], [mu * x for x in a], cov_u)
        e_u = evaluate_series([r_u], math.atan(m * lam / mu))
        self.assertAlmostEqual(e_u["K"], e_base["K"], places=11)
        self.assertAlmostEqual(e_u["T2"], e_base["T2"], places=11)

        # 2. Force-acceleration swap => same K, T^2 at 1/m
        idx = [2, 3, 0, 1]
        cov_s = [[cov[idx[i]][idx[j]] for j in range(4)] for i in range(4)]
        r_s = Reading(a, F, cov_s)
        e_s = evaluate_series([r_s], math.atan(1.0 / m))
        self.assertAlmostEqual(e_s["K"], e_base["K"], places=11)
        self.assertAlmostEqual(e_s["T2"], e_base["T2"], places=11)

        # 3. Spatial rotation R on both channels simultaneously
        ang = 0.73
        c, s = math.cos(ang), math.sin(ang)
        rot2 = [[c, -s], [s, c]]
        rot4 = [
            [c, -s, 0.0, 0.0],
            [s, c, 0.0, 0.0],
            [0.0, 0.0, c, -s],
            [0.0, 0.0, s, c],
        ]
        F_r = [rot2[0][0] * F[0] + rot2[0][1] * F[1], rot2[1][0] * F[0] + rot2[1][1] * F[1]]
        a_r = [rot2[0][0] * a[0] + rot2[0][1] * a[1], rot2[1][0] * a[0] + rot2[1][1] * a[1]]
        cov_r = [
            [sum(rot4[i][p] * cov[p][q] * rot4[j][q] for p in range(4) for q in range(4)) for j in range(4)]
            for i in range(4)
        ]
        r_rot = Reading(F_r, a_r, cov_r)
        e_rot = evaluate_series([r_rot], math.atan(m))
        self.assertAlmostEqual(e_rot["K"], e_base["K"], places=11)
        self.assertAlmostEqual(e_rot["T2"], e_base["T2"], places=11)

    def test_cans_scenes_door_a_vs_door_b_including_50_sigma(self):
        # Scene A (near): (5.5, 5.6) +- 0.1
        scene_a = analyze_case(
            [isotropic(5.5, 1.0, 0.1, 0.002), isotropic(5.6, 1.0, 0.1, 0.002)],
            common_ratio=50.0, scale_hint=5.5, n_grid=800,
        )
        self.assertAlmostEqual(scene_a["m_hat"], 5.55, places=2)
        self.assertGreater(scene_a["surprise_p"], 0.5)

        # Scene B (10 sigma): (5.0, 6.0) +- 0.1
        scene_b = analyze_case(
            [isotropic(5.0, 1.0, 0.1, 0.002), isotropic(6.0, 1.0, 0.1, 0.002)],
            common_ratio=50.0, scale_hint=5.5, n_grid=800,
        )
        self.assertAlmostEqual(scene_b["m_hat"], 5.50, places=2)
        self.assertLess(scene_b["surprise_p"], 1e-10)
        # Door B widens cleanly to cover [5, 6] without any 79-decade blowup
        self.assertIn("[4.9", scene_b["door_b_set"])
        self.assertIn("6.0", scene_b["door_b_set"])

        # Scene B extreme (50 sigma): does not crash; Door B splits around 5 and 6
        scene_50 = analyze_case(
            [isotropic(5.0, 1.0, 0.02, 0.0004), isotropic(6.0, 1.0, 0.02, 0.0004)],
            common_ratio=50.0, scale_hint=5.5, n_grid=800,
        )
        self.assertLess(scene_50["surprise_p"], 1e-200)
        self.assertIn(" U ", scene_50["door_b_set"])

    def test_right_angle_scene_c_and_antiparallel_ghost(self):
        # Right-angle Scene C (G29): both T45 and CLR stay (0, inf) with huge surprise
        scene_c = analyze_case(
            [isotropic(1.0, 0.0, 0.05, 0.05), isotropic(0.0, 1.0, 0.05, 0.05)],
            common_ratio=1.0, scale_hint=1.0, n_grid=400,
        )
        self.assertIsNone(scene_c["m_hat"])
        self.assertAlmostEqual(scene_c["K_min"], 400.0, places=6)
        self.assertEqual(scene_c["t45_set"], "(0, inf)")
        self.assertEqual(scene_c["clr_set"], "(0, inf)")

        # Anti-parallel pair: T45 has hilltop ghost near 1; Moreira CLR returns EMPTY
        anti = analyze_case(
            [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)],
            common_ratio=1.0, scale_hint=1.0, n_grid=600,
        )
        self.assertAlmostEqual(anti["K_min"], 500.0, places=5)
        self.assertNotEqual(anti["t45_set"], "EMPTY")
        self.assertEqual(anti["clr_set"], "EMPTY")

    def test_infimal_convolution_lumping_and_common_force(self):
        # Exact Gaussian variance addition: N(2, 0.3^2) + N(3, 0.4^2) = N(5, 0.5^2)
        K1 = lambda m: ((m - 2.0) / 0.3) ** 2
        K2 = lambda m: ((m - 3.0) / 0.4) ** 2
        for M, val in zip((4.0, 5.0, 6.0), compose_lump(K1, K2, (4.0, 5.0, 6.0))):
            self.assertAlmostEqual(val, ((M - 5.0) / 0.5) ** 2, places=6)

        # Swap duality: common force (L14) on reciprocal curves matches lumping (L5)
        Ms = (0.5, 1.2, 2.0)
        lump_recip = compose_lump(lambda w: K1(1.0 / w), lambda w: K2(1.0 / w), [1.0 / M for M in Ms])
        cf_vals = compose_common_force(K1, K2, Ms)
        for a, b in zip(lump_recip, cf_vals):
            self.assertAlmostEqual(a, b, places=10)

        # Gluing 5.0 +- 0.5 kg onto stationary 0/0 (K_zero == 0):
        # zero cost for M > 5, exact Gaussian cost ((M - 5)/0.5)^2 for M <= 5
        K_dens = lambda m: ((m - 5.0) / 0.5) ** 2
        K_zero = lambda m: 0.0
        vals = compose_lump(K_dens, K_zero, (3.5, 4.02, 5.0, 6.5, 20.0))
        self.assertAlmostEqual(vals[0], ((3.5 - 5.0) / 0.5) ** 2, places=3)
        self.assertAlmostEqual(vals[1], ((4.02 - 5.0) / 0.5) ** 2, places=3)
        self.assertLess(vals[2], 1e-4)
        self.assertLess(vals[3], 1e-6)
        self.assertLess(vals[4], 1e-6)

    def test_monte_carlo_coverage_in_hard_regimes(self):
        rng = random.Random(2026)

        def rand_unit_vec(d):
            v = [rng.gauss(0.0, 1.0) for _ in range(d)]
            nrm = math.sqrt(dot(v, v))
            return [x / nrm for x in v]

        def make_world(d, N, r_len, m_true):
            th = math.atan(m_true)
            s, co = math.sin(th), math.cos(th)
            return [
                isotropic(
                    [r_len * s * x + rng.gauss(0.0, 1.0) for x in u],
                    [r_len * co * x + rng.gauss(0.0, 1.0) for x in u],
                    1.0, 1.0,
                )
                for u in (rand_unit_vec(d) for _ in range(N))
            ]

        def check_world(gen, m_true, reps=250):
            phi0 = math.atan(m_true)
            c_t45 = c_clr = 0
            for _ in range(reps):
                rds = gen()
                df_tot = sum(r.k for r in rds)
                e0 = evaluate_series(rds, phi0)
                if e0["T2"] <= CHI2_95_1:
                    c_t45 += 1
                s_ff = sum(dot(r.force, r.force) for r in rds)
                s_aa = sum(dot(r.acceleration, r.acceleration) for r in rds)
                s_fa = sum(dot(r.force, r.acceleration) for r in rds)
                k_w = max(0.0, 0.5 * ((s_ff + s_aa) - math.hypot(s_ff - s_aa, 2.0 * s_fa)))
                if moreira_clr_pvalue(max(0.0, e0["K"] - k_w), e0["H"], df_tot) >= 0.05:
                    c_clr += 1
            return c_t45 / reps, c_clr / reps

        # 1. Pure noise, far-off mass m = 20 (Gleser-Hwang corner)
        t45_n, clr_n = check_world(lambda: make_world(3, 20, 0.0, 20.0), 20.0)
        self.assertGreater(t45_n, 0.90)
        self.assertGreater(clr_n, 0.90)

        # 2. Weak signal 3D N=20 (killed ours27 and ours31)
        t45_w, clr_w = check_world(lambda: make_world(3, 20, 1.0, 2.0), 2.0)
        self.assertGreater(t45_w, 0.90)
        self.assertGreater(clr_w, 0.90)

        # 3. T43 leverage (1 huge reading + 99 nulls)
        def make_t43():
            th = math.atan(2.0)
            s, co = math.sin(th), math.cos(th)
            return [isotropic(20.0 * s + rng.gauss(0.0, 1.0), 20.0 * co + rng.gauss(0.0, 1.0), 1.0, 1.0)] + [
                isotropic(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0), 1.0, 1.0) for _ in range(99)
            ]

        t45_l, clr_l = check_world(make_t43, 2.0)
        self.assertGreater(t45_l, 0.90)
        self.assertGreater(clr_l, 0.90)


if __name__ == "__main__":
    unittest.main()

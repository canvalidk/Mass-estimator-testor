"""Tests for E13 Push 4 & 5: Z-safe weight Pareto frontier (E2/E3) and Piecewise Exact Door B (O14/O18).

Pure stdlib unittest (also pytest-discoverable) so it runs without external packages.
"""
import math
import random
import unittest

from explore.e13_general_clr import (
    SpanReading,
    _make_correlated_world,
    _make_scale_mixture_world,
    _weight_fn,
    gclr_confidence_set,
    swap_reading,
)
from explore.e13_push4_5 import (
    door_b_piecewise_confidence_set,
    door_b_piecewise_exact_pvalue,
    evaluate_weighted_t45,
    exact_t_ar_quantile,
    stein_pareto_metrics,
    student_t_series_cost,
)
from explore.e13_unnormalized_evidence import CHI2_95_1, Reading, isotropic


class TestE13Push4And5(unittest.TestCase):
    def test_weighted_t45_invariance_z_safety_and_correlated_coverage(self):
        F = [3.2, -1.1]
        a = [1.4, 0.5]
        cov = [
            [0.25, 0.03, 0.02, -0.01],
            [0.03, 0.16, -0.01, 0.015],
            [0.02, -0.01, 0.09, 0.005],
            [-0.01, 0.015, 0.005, 0.12],
        ]
        r1 = Reading(F, a, cov)
        r2 = Reading([1.5, 0.8], [0.7, 0.3], cov)
        rds = [r1, r2]
        m = 2.1

        for wm in ("z_smooth", "z_linear", "z_cap10"):
            base = evaluate_weighted_t45(rds, math.atan(m), weight_mode=wm)["T2"]

            # Exact bit-for-bit Z invariance (appending exact zeros)
            rds_z = rds + [isotropic([0.0, 0.0], [0.0, 0.0], 1.0, 1.0)] * 5
            self.assertEqual(evaluate_weighted_t45(rds_z, math.atan(m), weight_mode=wm)["T2"], base)

            # Units scaling
            lam, mu = 2.5, 0.4
            scales = [lam, lam, mu, mu]
            cov_u = [[cov[i][j] * scales[i] * scales[j] for j in range(4)] for i in range(4)]
            rds_u = [
                Reading([lam * x for x in r.force], [mu * x for x in r.acceleration], cov_u)
                for r in rds
            ]
            self.assertAlmostEqual(
                evaluate_weighted_t45(rds_u, math.atan(m * lam / mu), weight_mode=wm)["T2"],
                base,
                places=11,
            )

            # Force-acceleration swap
            rds_s = [swap_reading(r) for r in rds]
            self.assertAlmostEqual(
                evaluate_weighted_t45(rds_s, math.atan(1.0 / m), weight_mode=wm)["T2"],
                base,
                places=11,
            )

        # Full-type correlated coverage (d=3, k=2, rho=0.6, mixed ratios)
        rng = random.Random(2026)
        phi0 = math.atan(2.5)
        for r_len in (0.0, 1.0):
            reps = 200
            c_lin = c_cap = 0
            for _ in range(reps):
                world = _make_correlated_world(rng, 10, 2.5, r_len, d=3, k=2, rho=0.6, ratio_spread=True)
                if evaluate_weighted_t45(world, phi0, weight_mode="z_linear")["T2"] <= CHI2_95_1:
                    c_lin += 1
                if evaluate_weighted_t45(world, phi0, weight_mode="z_cap10")["T2"] <= CHI2_95_1:
                    c_cap += 1
            self.assertGreater(c_lin / reps, 0.90)
            self.assertGreater(c_cap / reps, 0.90)

    def test_stein_pareto_closed_form_identities(self):
        # Exact analytical identities from Stein's lemma + Cauchy-Schwarz:
        #   w = 1:   ARE_1 = 1,    ARE_3 = 1,    N_drown(k=1, r=20) = 400
        #   w = H:   ARE_1 = 3/5,  ARE_3 = 5/7,  N_drown(k=1, r=20) = 400^3 / 15 = 4,266,666.667
        #   w = H^2: ARE_1 = 5/21,               N_drown(k=1, r=20) = 400^5 / 945
        m1_u = stein_pareto_metrics("unweighted", k=1, r_strong=20.0)
        m1_l = stein_pareto_metrics("z_linear", k=1, r_strong=20.0)
        m3_l = stein_pareto_metrics("z_linear", k=3, r_strong=20.0)
        m1_p2 = stein_pareto_metrics("z_pow2", k=1, r_strong=20.0)

        self.assertAlmostEqual(m1_u["are"], 1.0, places=10)
        self.assertAlmostEqual(m1_u["n_drown"], 400.0, places=8)
        self.assertAlmostEqual(m1_l["are"], 0.60, places=10)
        self.assertAlmostEqual(m3_l["are"], 5.0 / 7.0, places=10)
        self.assertAlmostEqual(m1_l["n_drown"], (400.0 ** 3) / 15.0, places=3)
        self.assertAlmostEqual(m1_p2["are"], 5.0 / 21.0, places=10)
        self.assertAlmostEqual(m1_p2["e_hw2"], 945.0, places=7)

    def test_three_regime_pareto_tradeoff_and_antiparallel_ghost_elimination(self):
        rng = random.Random(2026)
        th = math.atan(2.0)
        s, co = math.sin(th), math.cos(th)
        phi0 = math.atan(2.0)
        phi_4 = math.atan(4.0)
        phi_235 = math.atan(2.35)

        def score_1d(rds, ph, wm):
            sp, cp = math.sin(ph), math.cos(ph)
            num = den = 0.0
            for f, a in rds:
                u = sp * f + cp * a
                v = cp * f - sp * a
                w = _weight_fn(u * u, wm)
                num += w * u * v
                den += (w * w) * (u * u)
            return 0.0 if den <= 1e-30 else (num * num) / den

        # 1. Buried sparse (1x r=20 + 2000 random noisy zeros): w=H crushes T45 while keeping ~95% cov
        reps = 120
        cov_l = pwr_t45 = pwr_lin = 0
        for _ in range(reps):
            rds = [(20.0 * s + rng.gauss(0, 1), 20.0 * co + rng.gauss(0, 1))] + [
                (rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(2000)
            ]
            if score_1d(rds, phi0, "z_linear") <= CHI2_95_1:
                cov_l += 1
            if score_1d(rds, phi_4, "unweighted") > CHI2_95_1:
                pwr_t45 += 1
            if score_1d(rds, phi_4, "z_linear") > CHI2_95_1:
                pwr_lin += 1
        self.assertGreater(cov_l / reps, 0.89)
        self.assertGreater(pwr_lin / reps, 0.90)
        self.assertLess(pwr_t45 / reps, 0.55)

        # 2. Mixed strong (1x r=20 + 25x r=5): T45 beats unbounded w=H (genuine Pareto trade-off!)
        reps2 = 200
        m_t45 = m_lin = 0
        for _ in range(reps2):
            rds = [(20.0 * s + rng.gauss(0, 1), 20.0 * co + rng.gauss(0, 1))] + [
                (5.0 * s + rng.gauss(0, 1), 5.0 * co + rng.gauss(0, 1)) for _ in range(25)
            ]
            if score_1d(rds, phi_235, "unweighted") > CHI2_95_1:
                m_t45 += 1
            if score_1d(rds, phi_235, "z_linear") > CHI2_95_1:
                m_lin += 1
        self.assertGreater(m_t45, m_lin + 15)

        # 3. Ghost-free Z-GCLR(RP^1) with w=H still eliminates the anti-parallel hilltop ghost
        anti = [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)]
        self.assertEqual(
            gclr_confidence_set(anti, domain="whole", weight_mode="z_linear", scale_hint=1.0, n_grid=61, n_draws=99),
            "EMPTY",
        )

    def test_exact_t_ar_quantile_and_piecewise_door_b_closing_1d_n2_gap(self):
        # 1. Exact closed-form check at N=1, k=1, nu=10: F_{1,10; 0.95} = 4.9646027 -> 11*log(1.49646) = 4.4342
        q1 = exact_t_ar_quantile([1], nu=10.0, n_mc=50000, seed=7)
        self.assertAlmostEqual(q1, 11.0 * math.log1p(4.9646027 / 10.0), delta=0.06)

        # 2. Honest Z-failure of AR-bounded t-LR: appending a reading raises q95_t_AR strictly
        q2 = exact_t_ar_quantile([1, 1], nu=10.0, n_mc=40000, seed=7)
        q3 = exact_t_ar_quantile([1, 1, 1], nu=10.0, n_mc=40000, seed=7)
        self.assertGreater(q2, q1 + 1.5)
        self.assertGreater(q3, q2 + 1.5)

        # 3. 1D N=2 Scene B 10-sigma & 50-sigma: closes the (0, inf) SCLR+ gap AND avoids t-AR EMPTY
        s_0 = door_b_piecewise_confidence_set([isotropic(0.0, 0.0, 1.0, 1.0)] * 2, nu=10.0, scale_hint=1.0)
        self.assertEqual(s_0, "(0, inf)")

        s_10 = door_b_piecewise_confidence_set(
            [isotropic(5.0, 1.0, 0.1, 0.002), isotropic(6.0, 1.0, 0.1, 0.002)],
            nu=10.0, scale_hint=5.5,
        )
        self.assertNotEqual(s_10, "(0, inf)")
        self.assertNotEqual(s_10, "EMPTY")
        self.assertIn("[4.8", s_10)
        self.assertIn("6.1", s_10)

        s_50 = door_b_piecewise_confidence_set(
            [isotropic(5.0, 1.0, 0.02, 0.0004), isotropic(6.0, 1.0, 0.02, 0.0004)],
            nu=10.0, scale_hint=5.5,
        )
        self.assertIn(" U ", s_50)

        # 4. Proved >= 95% coverage on 1D N=2 (ar_bounded_t_lr branch) AND 3D N=15 weak (sclr_t branch)
        rng = random.Random(2028)
        phi0 = math.atan(2.0)
        reps_1d = 250
        cov_1d = 0
        for rep in range(reps_1d):
            rds = _make_scale_mixture_world(rng, 1, 2, 0.0 if rep % 2 == 0 else 3.0, 2.0, nu_df=10)
            res = door_b_piecewise_exact_pvalue(
                [SpanReading(r) for r in rds], phi0, nu=10.0, n_grid=41, scale_hint=2.0, q_t_ar=q2
            )
            self.assertEqual(res["branch"], "ar_bounded_t_lr")
            if res["accepted"]:
                cov_1d += 1
        self.assertGreater(cov_1d / reps_1d, 0.93)

        reps_3d = 100
        cov_3d = 0
        for rep in range(reps_3d):
            rds = _make_scale_mixture_world(rng, 3, 15, 1.0, 2.0, nu_df=10)
            res = door_b_piecewise_exact_pvalue(
                [SpanReading(r) for r in rds], phi0, nu=10.0, n_grid=25, n_draws=99, seed=rep * 17 + 1, scale_hint=2.0
            )
            self.assertEqual(res["branch"], "sclr_t")
            if res["accepted"]:
                cov_3d += 1
        self.assertGreater(cov_3d / reps_3d, 0.89)


if __name__ == "__main__":
    unittest.main()

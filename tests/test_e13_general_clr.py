"""Tests for E13 Push 1–3: General-type GCLR/Z-GCLR, Exact Two-Body Composition, and SCLR+/O7 Floor.

Pure stdlib unittest (also pytest-discoverable) so it runs without external packages.
"""
import math
import random
import unittest

from explore.e13_general_clr import (
    SpanReading,
    _make_correlated_world,
    _make_scale_mixture_world,
    _rand_unit_sphere,
    gclr_confidence_set,
    general_clr_pvalue,
    precompute_grid,
    swap_reading,
    two_body_common_force_conditional_pvalue,
    two_body_lump_conditional_pvalue,
)
from explore.e13_unnormalized_evidence import (
    CHI2_95_1,
    Reading,
    compose_lump,
    dot,
    evaluate_series,
    isotropic,
)


class TestE13GeneralCLR(unittest.TestCase):
    def test_span_reading_matches_evaluate_series_on_correlated_subspace(self):
        rng = random.Random(101)
        rds = _make_correlated_world(rng, 4, 2.5, 2.0, d=3, k=2, rho=0.6, ratio_spread=True)
        srs = [SpanReading(r) for r in rds]
        for m in (0.3, 1.0, 2.5, 7.0):
            phi = math.atan(m)
            ev = evaluate_series(rds, phi)
            h_sum = 0.0
            k_sum = 0.0
            for sr in srs:
                _, _, _, _, h_i, k_i = sr.uv_at_phi(phi)
                h_sum += h_i
                k_sum += k_i
            self.assertAlmostEqual(h_sum, ev["H"], places=11)
            self.assertAlmostEqual(k_sum, ev["K"], places=11)

    def test_units_swap_and_rotation_invariance_of_gclr_zgclr_sclr(self):
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
        # Symmetric grid under scaling/swap: evaluate lr_obs via fine dense grid on RP^1
        sr = SpanReading(r)
        _, _, _, _, h0, k0 = sr.uv_at_phi(math.atan(m))

        # 1. Units scaling (F, a) -> (lam*F, mu*a) at m' = (lam/mu)*m
        lam, mu = 3.5, 0.7
        scales = [lam, lam, mu, mu]
        cov_u = [[cov[i][j] * scales[i] * scales[j] for j in range(4)] for i in range(4)]
        sr_u = SpanReading(Reading([lam * x for x in F], [mu * x for x in a], cov_u))
        _, _, _, _, hu, ku = sr_u.uv_at_phi(math.atan(m * lam / mu))
        self.assertAlmostEqual(hu, h0, places=11)
        self.assertAlmostEqual(ku, k0, places=11)

        # 2. Force-acceleration swap at 1/m
        sr_s = SpanReading(swap_reading(r))
        _, _, _, _, hs, ks = sr_s.uv_at_phi(math.atan(1.0 / m))
        self.assertAlmostEqual(hs, h0, places=11)
        self.assertAlmostEqual(ks, k0, places=11)

        # 3. Spatial rotation R on both channels simultaneously (exact grid match at same phi!)
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
        sr_rot = SpanReading(Reading(F_r, a_r, cov_r))
        for kw in (
            dict(weight_mode="unweighted", null_mode="gaussian"),
            dict(weight_mode="z_smooth", null_mode="gaussian"),
            dict(weight_mode="unweighted", null_mode="spherical"),
        ):
            r0 = general_clr_pvalue([sr], math.atan(m), n_grid=41, n_draws=99, seed=7, **kw)
            rr = general_clr_pvalue([sr_rot], math.atan(m), n_grid=41, n_draws=99, seed=7, **kw)
            self.assertAlmostEqual(rr["lr_obs"], r0["lr_obs"], places=11)
            self.assertAlmostEqual(rr["k_obs"], r0["k_obs"], places=11)
            self.assertAlmostEqual(rr["weights"][0], r0["weights"][0], places=11)

    def test_exact_z_invariance_and_sigma_infinity_limit(self):
        base = [isotropic(4.5, 2.0, 0.5, 0.5), isotropic(3.5, 1.5, 0.5, 0.5)]
        with_zeros = base + [isotropic(0.0, 0.0, 1.0, 1.0)] * 25
        with_big_sig = base + [isotropic(1.0, 2.0, 1e6, 1e6)] * 25
        sr_b = [SpanReading(r) for r in base]
        sr_z = [SpanReading(r) for r in with_zeros]
        sr_s = [SpanReading(r) for r in with_big_sig]
        phi_t = math.atan(3.2)

        # Z-GCLR+ and SCLR+ are BIT-FOR-BIT invariant to appended exact zeros
        for kw in (
            dict(weight_mode="z_smooth", null_mode="gaussian"),
            dict(weight_mode="unweighted", null_mode="spherical"),
        ):
            pb = general_clr_pvalue(sr_b, phi_t, n_grid=41, n_draws=199, seed=42, scale_hint=2.2, **kw)["p_value"]
            pz = general_clr_pvalue(sr_z, phi_t, n_grid=41, n_draws=199, seed=42, scale_hint=2.2, **kw)["p_value"]
            self.assertEqual(pb, pz)

        # Z-GCLR+ also stays invariant as sigma -> infinity (within MC error)
        pb_hi = general_clr_pvalue(sr_b, phi_t, weight_mode="z_smooth", null_mode="gaussian", n_grid=41, n_draws=999, seed=11, scale_hint=2.2)["p_value"]
        ps_hi = general_clr_pvalue(sr_s, phi_t, weight_mode="z_smooth", null_mode="gaussian", n_grid=41, n_draws=999, seed=11, scale_hint=2.2)["p_value"]
        self.assertAlmostEqual(pb_hi, ps_hi, delta=0.03)

        # Unweighted GCLR+ shifts on appended zeros AND on sigma=1e6 (same whitened z=0)
        p_u_base = general_clr_pvalue(sr_b, phi_t, weight_mode="unweighted", null_mode="gaussian", n_grid=41, n_draws=499, seed=11, scale_hint=2.2)["p_value"]
        p_u_zero = general_clr_pvalue(sr_z, phi_t, weight_mode="unweighted", null_mode="gaussian", n_grid=41, n_draws=499, seed=11, scale_hint=2.2)["p_value"]
        self.assertNotEqual(p_u_base, p_u_zero)

    def test_scene_c_antiparallel_ghost_and_rank1_pinhole(self):
        # Right-angle Scene C (G29): all four variants return (0, inf)
        scene_c = [isotropic(1.0, 0.0, 0.05, 0.05), isotropic(0.0, 1.0, 0.05, 0.05)]
        for dom in ("positive", "whole"):
            for wm in ("unweighted", "z_smooth"):
                self.assertEqual(
                    gclr_confidence_set(scene_c, domain=dom, weight_mode=wm, scale_hint=1.0, n_grid=61, n_draws=99),
                    "(0, inf)",
                )

        # Multi-direction anti-parallel pair: T45 accepts hilltop at m=1; GCLR & Z-GCLR kill it
        anti_pair = [isotropic(2.0, -1.0, 0.1, 0.1), isotropic(1.0, -2.0, 0.1, 0.1)]
        self.assertLess(evaluate_series(anti_pair, math.atan(1.0))["T2"], CHI2_95_1)
        self.assertEqual(
            gclr_confidence_set(anti_pair, domain="whole", weight_mode="unweighted", scale_hint=1.0, n_grid=61, n_draws=99),
            "EMPTY",
        )
        self.assertEqual(
            gclr_confidence_set(anti_pair, domain="whole", weight_mode="z_smooth", scale_hint=1.0, n_grid=61, n_draws=99),
            "EMPTY",
        )

        # Rank-1 single anti-parallel reading (2, -1): orthogonal mass m=0.5 has U(0.5)=0
        # Unweighted GCLR still rejects m=0.5; Z-GCLR has the theorem-forced pinhole ghost (p=1)
        sr_single = [SpanReading(isotropic(2.0, -1.0, 0.1, 0.1))]
        p_u = general_clr_pvalue(sr_single, math.atan(0.5), domain="whole", weight_mode="unweighted", n_grid=61, n_draws=199, scale_hint=0.5)["p_value"]
        p_z = general_clr_pvalue(sr_single, math.atan(0.5), domain="whole", weight_mode="z_smooth", n_grid=61, n_draws=199, scale_hint=0.5)["p_value"]
        self.assertLessEqual(p_u, 0.01)
        self.assertEqual(p_z, 1.0)

    def test_full_type_correlated_subspace_coverage(self):
        rng = random.Random(2026)
        phi0 = math.atan(2.5)
        for r_len in (0.0, 1.0):
            reps = 120
            c_gp = c_zgp = c_wilks = 0
            for rep in range(reps):
                rds = _make_correlated_world(rng, 10, 2.5, r_len, d=3, k=2, rho=0.6, ratio_spread=True)
                srs = [SpanReading(r) for r in rds]
                pre = precompute_grid(srs, domain="positive", n_grid=25, scale_hint=2.5)
                e0 = evaluate_series(rds, phi0)
                k_min_g = min(
                    sum(
                        sum(sum(a * b for a, b in zip(row, srs[i].z)) ** 2 for row in pre[1][i][g])
                        for i in range(len(srs))
                    )
                    for g in range(len(pre[0]))
                )
                if (e0["K"] - min(k_min_g, e0["K"])) <= CHI2_95_1:
                    c_wilks += 1
                if general_clr_pvalue(srs, phi0, domain="positive", weight_mode="unweighted", n_draws=99, seed=rep * 19 + 1, precomputed_grid=pre)["p_value"] > 0.05:
                    c_gp += 1
                if general_clr_pvalue(srs, phi0, domain="positive", weight_mode="z_smooth", n_draws=99, seed=rep * 19 + 2, precomputed_grid=pre)["p_value"] > 0.05:
                    c_zgp += 1
            self.assertGreater(c_gp / reps, 0.89)
            self.assertGreater(c_zgp / reps, 0.89)
            self.assertLess(c_wilks / reps, 0.80)

    def test_push2_two_body_lumping_common_force_and_zero_gluing(self):
        # 1. Exact L5 <-> L14 swap duality
        r1 = [isotropic(4.0, 2.0, 0.3, 0.2), isotropic(3.8, 2.1, 0.3, 0.2)]
        r2 = [isotropic(6.0, 2.0, 0.3, 0.2), isotropic(6.3, 1.9, 0.3, 0.2)]
        p_cf = two_body_common_force_conditional_pvalue(r1, r2, 1.2, 2.0, weight_mode="z_smooth", n_draws=99, n_grid=25, n_split=25, seed=77)
        p_lump_swapped = two_body_lump_conditional_pvalue(
            [SpanReading(swap_reading(r)) for r in r1],
            [SpanReading(swap_reading(r)) for r in r2],
            1.0 / 1.2,
            1.0 / 2.0,
            weight_mode="z_smooth",
            n_draws=99,
            n_grid=25,
            n_split=25,
            seed=77,
        )
        self.assertEqual(p_cf, p_lump_swapped)

        # 2. Gluing 5.0 +- 0.5 kg onto stationary 3D 0/0 with z_smooth recovers the exact
        # one-sided Chernoff 95% boundary [5.0 - 1.645 * 0.5, inf) = [4.178, inf):
        # M0 = 4.00 is 2 sigma below 5.0 -> exact one-sided p = Phi(-2) ~ 0.023 (rejected);
        # M0 = 4.35 is 1.3 sigma below 5.0 -> exact one-sided p = Phi(-1.3) ~ 0.097 (accepted);
        # M0 = 5.0 and 50.0 -> p = 1.0.
        sr_dens = [SpanReading(isotropic(5.0, 1.0, 0.5, 0.001))]
        sr_zero = [SpanReading(isotropic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1.0))]
        p_40 = max(two_body_lump_conditional_pvalue(sr_dens, sr_zero, 4.0, c, weight_mode="z_smooth", n_draws=399, seed=13) for c in (3.99, 3.8, 2.0))
        p_435 = two_body_lump_conditional_pvalue(sr_dens, sr_zero, 4.35, 4.34, weight_mode="z_smooth", n_draws=399, seed=13)
        p_50 = two_body_lump_conditional_pvalue(sr_dens, sr_zero, 50.0, 5.0, weight_mode="z_smooth", n_draws=199, seed=13)
        self.assertLessEqual(p_40, 0.05)
        self.assertGreater(p_435, 0.05)
        self.assertEqual(p_50, 1.0)

        # 3. Weak-signal two-body lumping coverage (M = 2 + 3 = 5, 3D N=5 each, r=1)
        rng = random.Random(2027)
        reps = 100
        c_cond = 0
        th1, th2 = math.atan(2.0), math.atan(3.0)
        for rep in range(reps):
            b1 = [
                isotropic(
                    [math.sin(th1) * x + rng.gauss(0, 1) for x in u],
                    [math.cos(th1) * x + rng.gauss(0, 1) for x in u],
                    1.0, 1.0,
                )
                for u in (_rand_unit_sphere(3, rng) for _ in range(5))
            ]
            b2 = [
                isotropic(
                    [math.sin(th2) * x + rng.gauss(0, 1) for x in u],
                    [math.cos(th2) * x + rng.gauss(0, 1) for x in u],
                    1.0, 1.0,
                )
                for u in (_rand_unit_sphere(3, rng) for _ in range(5))
            ]
            if two_body_lump_conditional_pvalue(
                [SpanReading(r) for r in b1],
                [SpanReading(r) for r in b2],
                5.0, 2.0, weight_mode="z_smooth", n_draws=99, n_grid=21, n_split=21, seed=rep * 23 + 3
            ) > 0.05:
                c_cond += 1
        self.assertGreater(c_cond / reps, 0.89)

    def test_push3_sclr_heavy_tails_1d_small_n_floor_and_o7_floor_union(self):
        # 1. Exact 1D N <= 4 small-N floor for SCLR+ (1/2^N >= 1/16 = 0.0625 > 0.05)
        for N in (2, 3, 4):
            rds = [isotropic(2.0 + 0.005 * ((i % 3) - 1), 1.0, 0.01, 0.01) for i in range(N)]
            self.assertEqual(
                gclr_confidence_set(rds, domain="positive", null_mode="spherical", scale_hint=2.0, n_grid=61, n_draws=199),
                "(0, inf)",
            )
        # At N=5, exact 1/32 = 0.03125 < 0.05 rejects wildly wrong m=10
        rds5 = [isotropic(2.0 + 0.005 * ((i % 3) - 1), 1.0, 0.01, 0.01) for i in range(5)]
        p5_far = general_clr_pvalue([SpanReading(r) for r in rds5], math.atan(10.0), domain="positive", null_mode="spherical", n_grid=61, scale_hint=2.0)["p_value"]
        self.assertAlmostEqual(p5_far, 1.0 / 32.0, places=12)

        # 2. Scale-blindness on over-agreeing data vs O7 floor union
        rds_tight10 = [isotropic(5.5 + 0.002 * (i - 4.5), 1.0, 0.1, 0.002) for i in range(10)]
        srs_t10 = [SpanReading(r) for r in rds_tight10]
        # m = 5.46 is well inside Gaussian 95% (~[5.44, 5.56]) but rejected by scale-blind SCLR+
        p_g = general_clr_pvalue(srs_t10, math.atan(5.46), domain="positive", weight_mode="z_smooth", null_mode="gaussian", n_grid=91, n_draws=199, scale_hint=5.5)["p_value"]
        p_s = general_clr_pvalue(srs_t10, math.atan(5.46), domain="positive", weight_mode="unweighted", null_mode="spherical", n_grid=91, scale_hint=5.5)["p_value"]
        self.assertGreater(p_g, 0.05)
        self.assertLessEqual(p_s, 0.05)
        self.assertGreater(max(p_g, p_s), 0.05)  # O7 floor union respects declared sigma!

        # 3. 2D N=2 discordant pair ((5,0) & (6,0)): SCLR+ and Floor Union widen past Gaussian
        rds_2d_far2 = [isotropic([5.0, 0.0], [1.0, 0.0], 0.1, 0.002), isotropic([6.0, 0.0], [1.0, 0.0], 0.1, 0.002)]
        srs_2d = [SpanReading(r) for r in rds_2d_far2]
        # At m = 5.8 (between 5.5 and 6.0): Gaussian rejects, Spherical & Floor Union accept!
        p_2d_g = general_clr_pvalue(srs_2d, math.atan(5.8), domain="positive", weight_mode="z_smooth", null_mode="gaussian", n_grid=91, n_draws=199, scale_hint=5.5)["p_value"]
        p_2d_s = general_clr_pvalue(srs_2d, math.atan(5.8), domain="positive", weight_mode="unweighted", null_mode="spherical", n_grid=91, n_draws=199, scale_hint=5.5)["p_value"]
        self.assertLessEqual(p_2d_g, 0.05)
        self.assertGreater(p_2d_s, 0.05)

        # 4. Heavy-tail Cauchy-like nu=2 world (3D, N=15): Gaussian GCLR+ collapses, SCLR+ holds ~95%
        rng = random.Random(2028)
        reps = 100
        c_g = c_s = 0
        phi0 = math.atan(2.0)
        for rep in range(reps):
            rds = _make_scale_mixture_world(rng, 3, 15, 3.0, 2.0, nu_df=2)
            srs = [SpanReading(r) for r in rds]
            pre = precompute_grid(srs, domain="positive", n_grid=25, scale_hint=2.0)
            if general_clr_pvalue(srs, phi0, domain="positive", null_mode="gaussian", n_draws=99, seed=rep * 13 + 1, precomputed_grid=pre)["p_value"] > 0.05:
                c_g += 1
            if general_clr_pvalue(srs, phi0, domain="positive", null_mode="spherical", n_draws=99, seed=rep * 13 + 2, precomputed_grid=pre)["p_value"] > 0.05:
                c_s += 1
        self.assertLess(c_g / reps, 0.65)
        self.assertGreater(c_s / reps, 0.89)


if __name__ == "__main__":
    unittest.main()

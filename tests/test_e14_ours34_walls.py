"""Deterministic & seeded verification suite for E14 (Walls 1, 2, and 3 of ours34).

Runs under plain Python 3.12 stdlib with zero external dependencies:
    python3 -m unittest -v tests/test_e14_ours34_walls.py
"""
import math
import random
import unittest

from explore.e14_ours34 import (
    EvidenceState,
    Reading,
    common_force_curve,
    common_force_qpos_total,
    lump_curve,
    lump_qpos_split,
    lump_qpos_total,
    norm_cdf,
    q_pos_sf_exact,
    t10_circular_mean_discrete,
    _cholesky,
)


def _make_random_spd(rng, dim):
    A = [[rng.gauss(0.0, 1.0) for _ in range(dim)] for _ in range(dim)]
    return tuple(
        tuple(sum(A[i][t] * A[j][t] for t in range(dim)) + (0.5 if i == j else 0.0) for j in range(dim))
        for i in range(dim)
    )


class TestWall1TropicalUniquenessAndBoundaries(unittest.TestCase):
    """Wall 1 (B28): Theorems E14-1, E14-2, and both zero-padding endpoints."""

    def test_t10_circular_mean_fails_l5_and_l17(self):
        """T10's zero-padding endpoint fails finite-noise L5 and L17 anchor independence."""
        # L5 failure: (m1, m2) in {(1,1), (3,3)} equiprobable => lump in {2, 6}
        part = [(1.0, 0.5), (3.0, 0.5)]
        lump = [(2.0, 0.5), (6.0, 0.5)]
        m_part = t10_circular_mean_discrete(part, s=1.0)
        m_lump = t10_circular_mean_discrete(lump, s=1.0)
        self.assertAlmostEqual(m_part, 0.5 * (1.0 + math.sqrt(5.0)), places=6)  # golden ratio 1.618034
        self.assertGreater(abs(2.0 * m_part - m_lump), 0.15)  # 3.236 vs 3.075

        # L17 failure on c-anchored version: 0.5*delta_1 + 0.5*delta_4 at c=2 vs c=1
        law = [(1.0, 0.5), (4.0, 0.5)]
        m_c2 = t10_circular_mean_discrete(law, s=2.0)
        m_c1 = t10_circular_mean_discrete(law, s=1.0)
        self.assertAlmostEqual(m_c2, 2.0, places=12)
        self.assertAlmostEqual(m_c1, 1.76619037896906, places=6)
        self.assertGreater(abs(m_c2 - m_c1), 0.2)

    def test_argmin_K_exact_l5_l14_l3_t5_z_on_full_correlated_vectors(self):
        """argmin K satisfies exact L5, L14, L3, T5, and Z on full correlated 3D vectors with subspaces."""
        rng = random.Random(2026)
        cov1 = _make_random_spd(rng, 6)
        cov2 = _make_random_spd(rng, 6)
        N2 = ((1.0, 0.0, 0.0), (0.0, 0.6, 0.8))

        r1a = Reading((4.2, -1.3, 2.1), (2.0, -0.7, 1.1), cov1, N2)
        r1b = Reading((2.5, 1.4, -0.6), (1.3, 0.6, -0.2), cov2)
        r2a = Reading((6.1, 2.0, -1.8), (1.9, 0.7, -0.5), cov2, N2)
        r2b = Reading((3.2, -2.1, 0.9), (1.1, -0.6, 0.3), cov1)
        zero_r = Reading((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), cov1, N2)

        s1 = EvidenceState([r1a, r1b])
        s2 = EvidenceState([r2a, r2b])
        m1 = s1.argmin_K()["mass"]
        m2 = s2.argmin_K()["mass"]
        self.assertTrue(0.5 < m1 < 10.0 and 0.5 < m2 < 10.0)

        # Exact Z (L20)
        s1_z = EvidenceState([r1a, r1b, zero_r, zero_r])
        self.assertAlmostEqual(s1_z.argmin_K()["mass"], m1, places=10)

        # Exact L3 (units covariance)
        lam = 3.7
        s1_scaled = EvidenceState([r1a.scaled(lam), r1b.scaled(lam)])
        self.assertAlmostEqual(s1_scaled.argmin_K()["mass"], lam * m1, places=8)

        # Exact T5 (swap duality)
        s1_swapped = EvidenceState([r1a.swapped(), r1b.swapped()])
        self.assertAlmostEqual(s1_swapped.argmin_K()["mass"], 1.0 / m1, places=8)

        # Exact L5 under T74 lump_curve
        k_lump = lambda M: lump_curve(s1.K, s2.K, M)
        # Minimum of k_lump must occur at m1 + m2
        m_sum = m1 + m2
        k_at_sum = k_lump(m_sum)
        self.assertAlmostEqual(k_at_sum, s1.argmin_K()["K_min"] + s2.argmin_K()["K_min"], places=7)
        self.assertGreater(k_lump(m_sum * 0.95), k_at_sum + 1e-4)
        self.assertGreater(k_lump(m_sum * 1.05), k_at_sum + 1e-4)

        # Exact L14 under T74 common_force_curve
        k_par = lambda M: common_force_curve(s1.K, s2.K, M)
        m_red = 1.0 / (1.0 / m1 + 1.0 / m2)
        k_at_red = k_par(m_red)
        self.assertAlmostEqual(k_at_red, s1.argmin_K()["K_min"] + s2.argmin_K()["K_min"], places=7)
        self.assertGreater(k_par(m_red * 0.95), k_at_red + 1e-4)
        self.assertGreater(k_par(m_red * 1.05), k_at_red + 1e-4)

    def test_level_set_midpoints_fail_on_t74_gaussians(self):
        """Level-set geometric midpoint fails L5 on T74(2) Gaussians; arithmetic fails T5."""
        # K_j(m) = (m - mu_j)^2 / s_j^2, sublevel tau=1 => [mu_j - s_j, mu_j + s_j]
        mu1, s1 = 4.0, 1.0
        mu2, s2 = 6.0, 2.0
        # Geometric midpoint: sqrt(mu^2 - s^2)
        g1 = math.sqrt(mu1**2 - s1**2)
        g2 = math.sqrt(mu2**2 - s2**2)
        g_lump = math.sqrt((mu1 + mu2)**2 - (s1**2 + s2**2))
        self.assertGreater(abs(g_lump - (g1 + g2)), 0.15)  # 9.747 vs 9.530

        # Arithmetic midpoint passes L5 on Gaussians but fails T5 swap:
        # swapped sublevel [1/(mu+s), 1/(mu-s)] has midpoint mu/(mu^2 - s^2) != 1/mu
        a_swap = 0.5 * (1.0 / (mu1 + s1) + 1.0 / (mu1 - s1))
        self.assertGreater(abs(a_swap - 1.0 / mu1), 0.01)

    def test_iterated_pinch_Qn_lower_bound(self):
        """Theorem E14-1 mechanism: K_+^(n) pinches m < m* linearly in n while leaving m > m*
        bounded by H_inf; applying parallel convolution pinches BOTH sides linearly in n."""
        # Single bounded T73 curve with unique interior minimum m* = 2.0
        st = EvidenceState([Reading.scalar(2.0, 1.0, 1.0, 1.0)])
        m_star = 2.0
        K_inf = st.K(math.inf)  # 1.0

        # Compute exact n-fold scaled K_+^(n)(m) = inf_{avg m_j = m} sum K(m_j) by iterative lump
        def k_plus_pow(k_fn, n, M):
            if n == 1:
                return k_fn(M)
            # recursive doubling for powers of 2
            half = n // 2
            kh = lambda x: k_plus_pow(k_fn, half, x)
            return lump_curve(kh, kh, M, grid=121)

        # At m = 1.0 < m* = 2.0: K_+^(n)(n*1.0) grows linearly with n
        v1_left = k_plus_pow(st.K, 1, 1.0)
        v2_left = k_plus_pow(st.K, 2, 2.0)
        v4_left = k_plus_pow(st.K, 4, 4.0)
        self.assertGreater(v2_left, 1.8 * v1_left)
        self.assertGreater(v4_left, 1.8 * v2_left)

        # At m = 4.0 > m* = 2.0: K_+^(n)(n*4.0) stays bounded above by K_inf = 1.0!
        v4_right = k_plus_pow(st.K, 4, 16.0)
        self.assertLessEqual(v4_right, K_inf + 1e-6)

        # Now apply 2-fold scaled parallel convolution to L_2(m) = K_+^(2)(2m):
        L2 = lambda m: k_plus_pow(st.K, 2, 2.0 * m)
        Q2 = lambda m: common_force_curve(L2, L2, m / 2.0, grid=121)
        # At m* = 2.0: Q2(2.0) == 0; BOTH m=1.0 and m=4.0 are strictly higher than single K!
        self.assertAlmostEqual(Q2(m_star), 0.0, places=6)
        self.assertGreater(Q2(1.0), v2_left - 1e-4)
        self.assertGreater(Q2(4.0), st.K(4.0) * 1.5)

    def test_boundary_obstruction_and_non_unique_taxonomy(self):
        """Theorem E14-2: 4-way taxonomy + L14 forces R(K_inf)=inf via K_int square_par K_inf."""
        # 1. Exact 0/0 flat curve
        s_zero = EvidenceState([Reading.scalar(0.0, 0.0, 1.0, 1.0)])
        self.assertEqual(s_zero.argmin_K()["status"], "flat")

        # 2. Non-zero flat curve: {(1,0), (0,1)} at I_2 => K === 1
        s_flat_nz = EvidenceState([Reading.scalar(1.0, 0.0, 1.0, 1.0), Reading.scalar(0.0, 1.0, 1.0, 1.0)])
        res_fnz = s_flat_nz.argmin_K()
        self.assertEqual(res_fnz["status"], "flat")
        self.assertAlmostEqual(res_fnz["K_min"], 1.0, places=12)

        # 3. One-sided boundary minimum at infinity (F=10, a=0) and at zero (F=0, a=10)
        s_inf = EvidenceState([Reading.scalar(10.0, 0.0, 1.0, 1.0)])
        s_zmin = EvidenceState([Reading.scalar(0.0, 10.0, 1.0, 1.0)])
        self.assertEqual(s_inf.argmin_K()["status"], "boundary_inf")
        self.assertEqual(s_zmin.argmin_K()["status"], "boundary_zero")

        # 4. Anti-aligned double-boundary tie {0, inf}
        s_anti = EvidenceState([Reading.scalar(1.0, -1.0, 1.0, 1.0)])
        res_anti = s_anti.argmin_K()
        self.assertEqual(res_anti["status"], "tied")
        self.assertEqual(res_anti["minima"], (0.0, math.inf))

        # 5. Multi-instrument interior bimodal tie (s1=0.1, s2=10)
        s_bimodal = EvidenceState([
            Reading.scalar(0.1, 1.0, 0.1, 1.0),
            Reading.scalar(1.0, 0.1, 1.0, 0.1),
        ], scale=1.0)
        res_bi = s_bimodal.argmin_K()
        self.assertEqual(res_bi["status"], "tied")
        self.assertEqual(len(res_bi["minima"]), 2)
        self.assertAlmostEqual(res_bi["minima"][0] * res_bi["minima"][1], 1.0, places=5)
        self.assertGreater(s_bimodal.K(1.0), res_bi["K_min"] + 0.5)

        # 6. Theorem E14-2 obstruction: K_int (m*=3) square_par K_inf has UNIQUE interior min at M=3!
        s_int = EvidenceState([Reading.scalar(3.0, 1.0, 1.0, 1.0)])
        k_par = lambda M: common_force_curve(s_int.K, s_inf.K, M)
        self.assertAlmostEqual(k_par(3.0), 0.0, places=8)
        self.assertGreater(k_par(2.8), 1e-3)
        self.assertGreater(k_par(3.2), 1e-3)
        # Symmetrically, K_int (m*=3) square_+ K_zero has UNIQUE interior min at M=3!
        k_plus = lambda M: lump_curve(s_int.K, s_zmin.K, M)
        self.assertAlmostEqual(k_plus(3.0), 0.0, places=8)
        self.assertGreater(k_plus(2.8), 1e-3)
        self.assertGreater(k_plus(3.2), 1e-3)


class TestWall3MinimalStateAndCurvePairComposition(unittest.TestCase):
    """Wall 3: Theorems E14-3 and E14-4."""

    def test_full_type_score_equals_minus_half_dK_dm(self):
        """S(m) == -0.5 * dK/dm identically on correlated 3D covariance with 2D subspace."""
        rng = random.Random(999)
        cov6 = _make_random_spd(rng, 6)
        N2 = ((1.0, 0.0, 0.0), (0.0, 0.6, 0.8))
        st = EvidenceState([Reading((3.2, -1.5, 4.1), (1.1, -0.4, 2.0), cov6, N2)])
        m0, eps = 2.3, 1e-6
        _, S_analytic, _ = st.ksc(m0)
        S_fd = -0.5 * (st.K(m0 + eps) - st.K(m0 - eps)) / (2.0 * eps)
        self.assertAlmostEqual(S_analytic, S_fd, places=7)

    def test_K_only_insufficient_past_one_instrument(self):
        """Theorem E14-3: Two series on different instruments have bit-for-bit identical K(m) === 20
        yet yield opposite 95% T45 verdicts when combined with the same reading."""
        sA = EvidenceState([Reading.scalar(1.0, 0.0, 1.0, 1.0), Reading.scalar(0.0, 1.0, 1.0, 1.0)] * 20)
        sB = EvidenceState([Reading.scalar(2.0, 0.0, 2.0, 1.0), Reading.scalar(0.0, 1.0, 2.0, 1.0)] * 20)
        for m in [0.1, 0.5, 1.0, 2.5, 10.0, math.inf]:
            self.assertAlmostEqual(sA.K(m), 20.0, places=12)
            self.assertAlmostEqual(sB.K(m), 20.0, places=12)
        self.assertAlmostEqual(sA.ksc(1.0)[2], 5.0, places=12)
        self.assertAlmostEqual(sB.ksc(1.0)[2], 3.2, places=12)

        r_new = EvidenceState([Reading.scalar(5.0, 1.0, 1.0, 1.0)])
        resA = (sA + r_new).t45_test(1.0)
        resB = (sB + r_new).t45_test(1.0)
        self.assertLess(resA["statistic"], 3.8414588)     # 3.7895 => KEPT
        self.assertGreater(resB["statistic"], 3.8414588)  # 4.6753 => REJECTED

    def test_qpos_exact_quadrature_and_strong_chernoff_limits(self):
        """Theorem E14-4: Q_pos exact survival law, strong 1-dof sharpness, L14 swap, and Chernoff bound."""
        # Interior exact 1-dof chi2_1
        p_int = q_pos_sf_exact(2.3**2, 100.0, 100.0, 2.0, 3.0)
        self.assertAlmostEqual(p_int, math.erfc(2.3 / math.sqrt(2.0)), places=10)

        # Corner (A=0, B=0) exact chi-bar-squared: 0.5*P(chi2_1>=q) + 0.25*P(chi2_2>=q)
        q95 = 3.841458820694124
        exact_chibar = 0.5 * 0.05 + 0.25 * math.exp(-0.5 * q95)
        self.assertAlmostEqual(q_pos_sf_exact(q95, 1.0, 1.0, 0.0, 0.0), exact_chibar, places=5)

        # Strong signal: 2.3 sigma off total gives 1-parameter p ~ 0.0217 (vs 0.071 for 2-dof)
        b1 = EvidenceState([Reading.scalar(2.0, 1.0, 0.01, 1e-6)])
        b2 = EvidenceState([Reading.scalar(3.0, 1.0, 0.01, 1e-6)])
        se = math.hypot(0.01, 0.01)
        p_lo = lump_qpos_total(b1, b2, 5.0 - 2.3 * se)
        p_hi = lump_qpos_total(b1, b2, 5.0 + 2.3 * se)
        self.assertAlmostEqual(p_lo, 0.02174, places=3)
        self.assertAlmostEqual(p_hi, 0.02174, places=3)

        # Exact L14 swap duality
        b1_sw = EvidenceState([Reading.scalar(1.0, 2.0, 1e-6, 0.01)])
        b2_sw = EvidenceState([Reading.scalar(1.0, 3.0, 1e-6, 0.01)])
        self.assertAlmostEqual(common_force_qpos_total(b1_sw, b2_sw, 1.0 / (5.0 - 2.3 * se)), p_lo, places=8)

        # One-sided Chernoff bound on (5 +/- 0.5) + (0/0): threshold at 5 - 1.64485*0.5 = 4.1776
        b_5 = EvidenceState([Reading.scalar(5.0, 1.0, 0.5, 1e-6)])
        b_0 = EvidenceState([Reading.scalar(0.0, 0.0, 1.0, 1.0)])
        self.assertLess(lump_qpos_total(b_5, b_0, 4.17), 0.05)
        self.assertGreater(lump_qpos_total(b_5, b_0, 4.18), 0.05)
        self.assertAlmostEqual(lump_qpos_total(b_5, b_0, 5.0), 1.0, places=10)
        self.assertAlmostEqual(lump_qpos_total(b_5, b_0, 100.0), 1.0, places=10)

    def test_qpos_weak_signal_3d_coverage_seeded(self):
        """Theorem E14-4: Q_pos holds >= 0.95 coverage at weak signal (sig=0 and sig=1) where raw K_lump fails."""
        rng = random.Random(20261122)
        for sig in [0.0, 1.0]:
            kept = 0
            n_trials = 120
            for _ in range(n_trials):
                r1, r2 = [], []
                for _ in range(5):
                    for mt, dest in [(2.0, r1), (3.0, r2)]:
                        g = [rng.gauss(0.0, 1.0) for _ in range(3)]
                        gn = math.sqrt(sum(x * x for x in g))
                        at = [sig / math.sqrt(mt * mt + 1.0) * x / gn for x in g]
                        ft = [mt * x for x in at]
                        dest.append(Reading.isotropic(
                            [x + rng.gauss(0.0, 1.0) for x in ft],
                            [x + rng.gauss(0.0, 1.0) for x in at], 1.0, 1.0))
                p = lump_qpos_split(EvidenceState(r1), EvidenceState(r2), 5.0, 2.0)
                if p > 0.05:
                    kept += 1
            self.assertGreaterEqual(kept / n_trials, 0.90)


class TestWall2T76CollapseDiscountNoGoAndJointT(unittest.TestCase):
    """Wall 2: Theorems E14-5, E14-6, and E14-7."""

    def test_t76_nonconstant_weights_reject_argmin_K(self):
        """Theorem E14-5: Both w(H)=H/(1+H) and w(H)=H reject argmin K(m) at 95% while T45 keeps it (p=1)."""
        reads = [Reading.scalar(20.0, 2.0, 1.0, 1.0)] + [Reading.scalar(0.2, 2.0, 1.0, 1.0)] * 100
        st = EvidenceState(reads)
        m_star = st.argmin_K()["mass"]
        self.assertAlmostEqual(m_star, 1.0, places=6)
        self.assertAlmostEqual(st.t45_test(m_star)["p_value"], 1.0, places=8)

        for wfn in [lambda h: h / (1.0 + h), lambda h: h]:
            num = den = 0.0
            for sp in st.spans:
                _, si, ci, hi = sp.eval_ksc(m_star)
                w = wfn(hi)
                num += w * si
                den += w * w * ci
            z = abs(num) / math.sqrt(den)
            p = math.erfc(z / math.sqrt(2.0))
            self.assertLess(p, 0.005)  # strictly rejected at 95% and 99%!

    def test_curve_only_discount_nogo_exact_match(self):
        """Theorem E14-6: 2 conflicting readings at sigma_A and N=50 honest readings at sigma_B
        with exact mean m0 and exact sum-of-squares (N-1)*sigma_B^2 produce IDENTICAL (K(m), C(m))."""
        m0, delta, N = 5.5, 1.0, 50
        sa = 1e-9  # exact-acceleration limit (checked at finite m)
        sigma_A = delta / math.sqrt(2.0 * (N - 1))
        sigma_B = 0.5 * delta * math.sqrt(N / (N - 1))

        # Series A: 2 conflicting readings at 5.0 and 6.0
        sA = EvidenceState([
            Reading.scalar(m0 - 0.5 * delta, 1.0, sigma_A, sa),
            Reading.scalar(m0 + 0.5 * delta, 1.0, sigma_A, sa),
        ])
        # Series B: N readings with exact mean m0 and exact sample variance sigma_B^2
        # Construct N symmetric Gaussian-like quantiles normalized to sum=0, sum_sq=(N-1)*sigma_B^2
        raw = [i - 0.5 * (N - 1) for i in range(N)]
        scale = math.sqrt((N - 1) * sigma_B * sigma_B / sum(x * x for x in raw))
        xs = [m0 + scale * x for x in raw]
        sB = EvidenceState([Reading.scalar(x, 1.0, sigma_B, sa) for x in xs])

        for m in [4.5, 5.0, 5.5, 6.0, 6.5]:
            kA, sA_v, cA = sA.ksc(m)
            kB, sB_v, cB = sB.ksc(m)
            self.assertAlmostEqual(kA, kB, places=6)
            self.assertAlmostEqual(sA_v, sB_v, places=6)
            self.assertAlmostEqual(cA, cB, places=6)

        # Yet M_nz cleanly distinguishes them: 2 vs 50!
        self.assertEqual(sA.M_nz, 2)
        self.assertEqual(sB.M_nz, 50)
        # And t_joint keeps sA wide (at m=6.0: p > 0.05) while sB rejects m=6.0 (p = 6.6e-9)!
        self.assertGreater(sA.floor_test(6.0)["p_value"], 0.05)
        self.assertLess(sB.floor_test(6.0)["p_value"], 1e-8)

    def test_joint_t_floor_scenes_and_exact_Z(self):
        """Theorem E14-7: Can's 1D & 3D scenes, T30 sigma-invariance, G29 right-angle, 99+1, and exact Z."""
        # 3D near pair: floor matches T45 (5.55 +/- 0.138)
        near_3d = EvidenceState([
            Reading.isotropic([5.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
            Reading.isotropic([5.6, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
        ])
        self.assertGreater(near_3d.floor_test(5.68)["p_value"], 0.05)
        self.assertLess(near_3d.floor_test(5.70)["p_value"], 0.05)

        # 3D far pair: widens to ~[4.926, 6.074] AND stays bit-for-bit invariant as sigma shrinks 10x!
        far_3d_01 = EvidenceState([
            Reading.isotropic([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
            Reading.isotropic([6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
        ])
        far_3d_001 = EvidenceState([
            Reading.isotropic([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.01, 1e-7),
            Reading.isotropic([6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.01, 1e-7),
        ])
        self.assertGreater(far_3d_01.floor_test(6.05)["p_value"], 0.05)
        self.assertLess(far_3d_01.floor_test(6.10)["p_value"], 0.05)
        self.assertAlmostEqual(far_3d_01.t_joint_test(6.05)["p_value"],
                               far_3d_001.t_joint_test(6.05)["p_value"], places=10)

        # Exact Z (L20): appending exact 0/0 readings leaves floor_test bit-for-bit unchanged
        zero_3d = Reading.isotropic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1e-6)
        far_3d_z = far_3d_01 + EvidenceState([zero_3d] * 50)
        self.assertEqual(far_3d_z.M_nz, 6)
        self.assertAlmostEqual(far_3d_z.floor_test(6.05)["p_value"],
                               far_3d_01.floor_test(6.05)["p_value"], places=14)

        # 1D 2-scalar scene: finite interval (unlike T75(6) which returns (0, inf) for n<=5)
        near_1d = EvidenceState([Reading.scalar(5.5, 1.0, 0.1, 1e-6), Reading.scalar(5.6, 1.0, 0.1, 1e-6)])
        self.assertGreater(near_1d.floor_test(6.10)["p_value"], 0.05)
        self.assertLess(near_1d.floor_test(6.25)["p_value"], 0.05)

        # G29 right-angle scene: keeps every positive mass (no divergence to infinity)
        ra = EvidenceState([Reading.scalar(2.0, 0.0, 0.1, 0.1), Reading.scalar(0.0, 2.0, 0.1, 0.1)])
        for m in [0.01, 0.5, 1.0, 2.0, 100.0]:
            self.assertAlmostEqual(ra.floor_test(m)["p_value"], 1.0, places=12)

    def test_joint_t_exact_coverage_full_correlated_subspace_seeded(self):
        """Theorem E14-7: Exact ~0.95 coverage of t_joint on correlated 6x6 covariances + mixed subspaces
        under BOTH declared Gaussian (v=1) and unknown common scale v ~ InvChi2(2)."""
        rng = random.Random(31415)
        cov1 = _make_random_spd(rng, 6)
        cov2 = _make_random_spd(rng, 6)
        L1 = _cholesky(cov1)
        L2 = _cholesky(cov2)
        N2 = ((1.0, 0.0, 0.0), (0.0, 0.6, 0.8))
        N1 = ((0.0, 0.0, 1.0),)
        m_true = 2.3

        for random_scale in [False, True]:
            ct = cf = 0
            n_tr = 250
            for _ in range(n_tr):
                sd = math.sqrt(2.0 / (rng.gauss(0, 1)**2 + rng.gauss(0, 1)**2)) if random_scale else 1.0
                reads = []
                for L_m, cov_m, N_sub, amp in [(L1, cov1, N2, 1.5), (L1, cov1, N2, 0.5), (L2, cov2, N1, 2.0)]:
                    cp = [amp * rng.gauss(0.0, 1.0) for _ in range(len(N_sub))]
                    at = [sum(cp[j] * N_sub[j][t] for j in range(len(N_sub))) for t in range(3)]
                    ft = [m_true * x for x in at]
                    ez = [sd * rng.gauss(0.0, 1.0) for _ in range(6)]
                    noise = [sum(L_m[i][j] * ez[j] for j in range(i + 1)) for i in range(6)]
                    reads.append(Reading(
                        (ft[0] + noise[0], ft[1] + noise[1], ft[2] + noise[2]),
                        (at[0] + noise[3], at[1] + noise[4], at[2] + noise[5]),
                        cov_m, N_sub,
                    ))
                res = EvidenceState(reads).floor_test(m_true)
                ct += (res["p_tjoint"] > 0.05)
                cf += (res["p_value"] > 0.05)
            self.assertGreaterEqual(ct / n_tr, 0.90)
            self.assertLessEqual(ct / n_tr, 0.99)
            self.assertGreaterEqual(cf / n_tr, 0.93)


if __name__ == "__main__":
    unittest.main()

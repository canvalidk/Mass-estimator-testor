"""Executable pure-stdlib checks for E17:
  - Part I (L25): Aczel curve-addition uniqueness, smooth L16+L20+L17 quadratic collapse without L13,
    and the sharp 1D-discontinuous / 2D-continuous split of T75 (4)'s spherical null at 0/0.
  - Part II (G25): Global non-Gaussian full-type vector consistency of T81's argmin K (E17-G1),
    global-ray counterexample to B27's local weighted rate + global chord bound (E17-G2),
    and heavy-tail t_2.5 (2 < nu < 4) breakdown of T81 vs Door B robust log-curve (E17-G3).

Runnable with `python3 -m unittest -v tests/test_e17_l25_and_g25.py`.
"""
import math
import random
import unittest

import explore.e17_l25_g25 as e17
from met.ours34 import Reading, State


class TestE17L25AndG25(unittest.TestCase):

    def test_e17_l1_aczel_unnormalized_monoid_vs_normalized_posterior(self):
        """Theorem E17-L1: Unnormalized K_A + K_B == K_{A U B} with 0/0 identity, whereas
        normalized posteriors p_A, p_B with half-Cauchy prior pi fail pointwise combination."""
        # Two scalar readings on unit instrument
        rA = Reading((6.0,), (2.0,), ((1.0, 0.0), (0.0, 1.0)))
        rB = Reading((9.0,), (3.0,), ((1.0, 0.0), (0.0, 1.0)))
        r0 = Reading((0.0,), (0.0,), ((1.0, 0.0), (0.0, 1.0)))
        sA, sB, s0, sAB = State([rA]), State([rB]), State([r0]), State([rA, rB])
        for m in (0.5, 1.5, 3.0, 8.0):
            self.assertAlmostEqual(sA.curve(m) + sB.curve(m), sAB.curve(m), places=13)
            self.assertAlmostEqual(sA.curve(m) + s0.curve(m), sA.curve(m), places=14)

        # Normalized posterior with half-Cauchy prior pi(m) = (2/pi)/(1+m^2): multiplying p_A(m)*p_B(m)
        # double-counts pi(m) and gives a DIFFERENT normalized law than p_{A U B}(m) (breaks L25(b)).
        grid = [0.05 * (200.0 ** (i / 120.0)) for i in range(121)]
        dm = [grid[i + 1] - grid[i] for i in range(len(grid) - 1)]
        mids = [0.5 * (grid[i] + grid[i + 1]) for i in range(len(grid) - 1)]
        pi = [(2.0 / math.pi) / (1.0 + m * m) for m in mids]

        def norm_post(st):
            raw = [p * math.exp(-0.5 * st.curve(m)) for p, m in zip(pi, mids)]
            Z = math.fsum(r * d for r, d in zip(raw, dm))
            return [r / Z for r in raw]

        pA, pB, pAB = norm_post(sA), norm_post(sB), norm_post(sAB)
        prod_raw = [a * b for a, b in zip(pA, pB)]
        Z_prod = math.fsum(r * d for r, d in zip(prod_raw, dm))
        p_prod = [r / Z_prod for r in prod_raw]
        tv_err = 0.5 * math.fsum(abs(a - b) * d for a, b, d in zip(p_prod, pAB, dm))
        self.assertGreater(tv_err, 0.03)  # normalized product fails to recover pooled posterior

    def test_e17_l2a_smooth_l16_l20_splitting_forces_quadratic_gram_without_l13(self):
        """Theorem E17-L2: Any C^2 sign-symmetric (L16) additive map Phi(y) with Phi(0)=0 (L20)
        satisfies n * Phi(y / sqrt(n)) -> Q(y), a pure quadratic function of the 3 Gram entries
        (||F||^2, F.a, ||a||^2) — recovering the Gaussian curve K(m) WITHOUT assuming L13."""
        # Test three very different non-Gaussian C^2 L16-symmetric losses (Student-t, tanh, rational)
        def k_i(f, a, m):
            num = math.fsum((fi - m * ai) ** 2 for fi, ai in zip(f, a))
            return num / (1.0 + m * m)

        non_gauss_losses = [
            ("student_nu4", lambda f, a, m: 2.5 * math.log(1.0 + k_i(f, a, m) / 4.0), 2.5 / 4.0),
            ("tanh", lambda f, a, m: math.tanh(k_i(f, a, m)), 1.0),
            ("rational", lambda f, a, m: k_i(f, a, m) / (1.0 + 0.3 * k_i(f, a, m)), 1.0),
        ]
        f = (3.2, -1.7, 0.9)
        a = (1.1, -0.4, 0.6)
        n = 10000
        s = math.sqrt(n)
        f_split = tuple(x / s for x in f)
        a_split = tuple(x / s for x in a)
        for name, phi, slope0 in non_gauss_losses:
            for m in (0.4, 2.0, 5.0):
                split_sum = n * phi(f_split, a_split, m)
                quad_limit = slope0 * k_i(f, a, m)
                self.assertAlmostEqual(split_sum, quad_limit, delta=1e-3 * max(1.0, quad_limit), msg=name)

    def test_e17_l2b_spherical_t75_1d_discontinuous_vs_2d_continuous_at_zero(self):
        """Theorem E17-L2(4): Door B spherical null T75 (4) is exact Z at 0/0, DISCONTINUOUS at 0/0
        in 1D scalars (one tiny reading breaks discrete sign ties: p = 0.0625 -> 0.03125), yet
        CONTINUOUS at 0/0 in 2D vectors when the base spherical null has no atom at T_obs."""
        # 1D scalar discontinuity
        base5 = [(5.0, 1.0)] * 5
        p5, _ = e17.spherical_p_1d_exact(base5, m0=1.0)
        p5_z, _ = e17.spherical_p_1d_exact(base5 + [(0.0, 0.0)], m0=1.0)
        p5_tiny, _ = e17.spherical_p_1d_exact(base5 + [(1e-6, 0.0)], m0=1.0)
        self.assertAlmostEqual(p5, 2.0 / 32.0, places=12)
        self.assertAlmostEqual(p5_z, 2.0 / 32.0, places=12)
        self.assertAlmostEqual(p5_tiny, 2.0 / 64.0, places=12)

        # 2D vector continuity under 20 tiny (1e-6, 1e-6) readings (continuous S^1 null)
        base_2d = [((5.0, 0.2), (1.0, 0.1)), ((4.8, -0.3), (1.0, -0.05)), ((5.2, 0.1), (1.0, 0.02))]
        p2_base, _ = e17.spherical_p_2d_mc(base_2d, m0=2.5, draws=3000, seed=17)
        p2_tiny, _ = e17.spherical_p_2d_mc(
            base_2d + [((1e-6, -1e-6), (1e-6, 1e-6))] * 20, m0=2.5, draws=3000, seed=17
        )
        # Note: RNG stream consumes extra draws for the 20 tiny readings, so compare within MC SE (~0.02)
        self.assertAlmostEqual(p2_base, p2_tiny, delta=0.035)
        self.assertLess(p2_base, 0.15)

    def test_e17_g1_global_nongaussian_vector_consistency_full_type(self):
        """Theorem E17-G1: Exact zero bias E[K(m) - K(m0)] == G(m) and global convergence of T81's
        argmin_{[0,inf]} K on correlated full-type (3D k=2 + 2D + 1D) mixed series with strongly
        skewed Exp(1)-1, asymmetric discrete, and heavy-tailed t_5 non-Gaussian world errors."""
        # 1. Verify exact zero-bias identity E[K(m) - K(m0)] == G(m) empirically under skewed Exp(1)-1
        cov2 = ((1.2, 0.3, 0.2, -0.1), (0.3, 0.9, -0.1, 0.25), (0.2, -0.1, 0.8, 0.15), (-0.1, 0.25, 0.15, 1.1))
        m0, m_alt = 2.0, 0.7
        a_true = (1.5, -0.8)
        f_true = (m0 * a_true[0], m0 * a_true[1])
        G_exact = State([Reading(f_true, a_true, cov2)]).curve(m_alt)
        rng = random.Random(1710)
        diffs = []
        for _ in range(1500):
            err = e17.sample_correlated_error(rng, cov2, kind="exp")
            st = State([Reading((f_true[0] + err[0], f_true[1] + err[1]),
                                (a_true[0] + err[2], a_true[1] + err[3]), cov2)])
            diffs.append(st.curve(m_alt) - st.curve(m0))
        mean_diff = math.fsum(diffs) / len(diffs)
        self.assertAlmostEqual(mean_diff, G_exact, delta=0.25 * max(1.0, G_exact))

        # 2. Cross-check fast scatter-matrix curve against met.ours34.State on n_blocks=3 to 1e-12
        chk = e17.simulate_full_type_t81_error(3, 1.0, m0=2.0, seed=77, return_readings=True)
        st_ref = State(chk["readings"])
        for m in (0.0, 0.4, 2.0, 6.0, math.inf):
            self.assertAlmostEqual(chk["curve_fn"](m), st_ref.curve(m), delta=1e-11 * max(1.0, st_ref.curve(m)))

        # 3. Global convergence of T81 on correlated 3D+2D+1D mixed non-Gaussian errors as D grows
        #    Schedule: amp_scale = 1.6 * n_blocks^(-0.15) => Delta ~ n_blocks^0.70, Delta/sqrt(D) ~ n_blocks^0.20 -> inf
        med_errs = []
        for n_blocks in (4, 32, 256):
            amp = 1.6 * (n_blocks ** (-0.15))
            errs = [
                e17.simulate_full_type_t81_error(n_blocks, amp, m0=2.0, seed=1000 * n_blocks + j)["ang_err"]
                for j in range(30)
            ]
            errs.sort()
            med_errs.append(errs[len(errs) // 2])
        self.assertLess(med_errs[1], med_errs[0])
        self.assertLess(med_errs[2], med_errs[1])
        self.assertLess(med_errs[2], 0.05)

    def test_e17_g2_global_ray_chord_bound_and_counterexample_to_local_b27(self):
        """Theorem E17-G2: (a) Exact global chord inequality sin^2(theta_i(m) - theta_i(m0)) >= w_i^2 * rho(m, m0)^2
        for all m in [0, inf] and all s_i > 0; (b) Global ray counterexample where B27's local weighted
        rate J^2/(W+J) -> inf (compact accuracy -> 1) yet global argmin_{[0,inf]} K escapes to >= N/2."""
        # (a) Exact global chord lower bound across wide (m, m0, s) grid
        for m0 in (0.1, 1.0, 7.5):
            for s in (1e-4, 0.1, 1.0, 10.0, 1e4):
                w = m0 * s / (m0 * m0 + s * s)
                th0 = math.atan(m0 / s)
                for m in (0.0, 1e-4, 0.05, 0.5, 1.0, 2.5, 20.0, 1e4, math.inf):
                    th = 0.5 * math.pi if math.isinf(m) else math.atan(m / s)
                    lhs = math.sin(th - th0) ** 2
                    rho = 1.0 if math.isinf(m) else abs(m - m0) / max(m, m0)
                    rhs = (w * rho) ** 2
                    self.assertGreaterEqual(lhs + 1e-15, rhs)

        # (b) Global ray counterexample: compact accuracy climbs while global escape climbs toward 1/2
        rng = random.Random(1720)
        rates = {}
        for N in (100, 1000):
            r1_sq = 8.0 * math.log(N)
            esc = cmp_ok = 0
            trials = 250
            for _ in range(trials):
                res = e17.exact_2ratio_counterexample_trial(rng, N, r1_sq, m0=1.0)
                if res["escaped"]:
                    esc += 1
                if abs(res["m_compact"] - 1.0) < 0.5:
                    cmp_ok += 1
            rates[N] = (cmp_ok / trials, esc / trials)
        self.assertGreater(rates[1000][0], rates[100][0])  # compact interval converges
        self.assertGreater(rates[1000][0], 0.88)
        self.assertGreater(rates[1000][1], rates[100][1])  # global argmin K escapes MORE as N grows!
        self.assertGreater(rates[1000][1], 0.35)

    def test_e17_g3_heavy_tail_t25_breaks_t81_at_gaussian_rate_and_recovers_above_nu_threshold(self):
        """Theorem E17-G3: Under standardized t_2.5 errors (exact Cov = I_2, infinite 4th moment,
        2/nu = 0.8): (a) at Delta_N = 8 * N^0.65 (so Delta/sqrt(N) -> inf yet Delta/N^0.8 -> 0),
        T81's argmin K fails consistency while Door B's log-curve converges to 0 error;
        (b) above the E17-G3 threshold (Delta_N = 6 * N^0.95 >> N^0.8), T81 recovers consistency."""
        m0 = 2.0
        u0_f, u0_a = m0 / math.hypot(m0, 1.0), 1.0 / math.hypot(m0, 1.0)
        grid_m = [0.2 * (50.0 ** (j / 30.0)) for j in range(31)]
        rng = random.Random(1730)

        # (a) Sub-threshold regime Delta_N = 8 * N^0.65: T81 fails while Door B log-curve converges
        t81_fail = {}
        rob_fail = {}
        for N in (100, 1000):
            amp = math.sqrt(8.0 * (N ** 0.65) / N)
            ft = fr = 0
            trials = 160
            for _ in range(trials):
                s_ff = s_fa = s_aa = 0.0
                pairs = []
                for _ in range(N):
                    ef, ea = e17.sample_unit_t_vec(rng, 2, 2.5)
                    f, a = amp * u0_f + ef, amp * u0_a + ea
                    s_ff += f * f
                    s_fa += f * a
                    s_aa += a * a
                    pairs.append((f, a))
                m_hat = e17.argmin_k_1instr_gram(s_ff, s_fa, s_aa)
                if not (math.isfinite(m_hat) and abs(m_hat - m0) < 0.5):
                    ft += 1
                best_m = min(
                    grid_m,
                    key=lambda m: math.fsum(math.log(1.0 + (f - m * a) ** 2 / (4.0 * (1.0 + m * m))) for f, a in pairs),
                )
                if abs(best_m - m0) >= 0.5:
                    fr += 1
            t81_fail[N] = ft / trials
            rob_fail[N] = fr / trials
        self.assertGreater(t81_fail[1000], 0.30)       # T81 does NOT converge
        self.assertLess(rob_fail[1000], rob_fail[100]) # Door B robust log-curve converges
        self.assertLess(rob_fail[1000], 0.05)

        # (b) Super-threshold regime Delta_N = 6 * N^0.95 (Delta / N^0.8 = 6 * N^0.15 -> inf): T81 recovers!
        ft_hi = {}
        for N in (50, 800):
            amp = math.sqrt(6.0 * (N ** 0.95) / N)
            ft = 0
            trials = 160
            for _ in range(trials):
                s_ff = s_fa = s_aa = 0.0
                for _ in range(N):
                    ef, ea = e17.sample_unit_t_vec(rng, 2, 2.5)
                    f, a = amp * u0_f + ef, amp * u0_a + ea
                    s_ff += f * f
                    s_fa += f * a
                    s_aa += a * a
                m_hat = e17.argmin_k_1instr_gram(s_ff, s_fa, s_aa)
                if not (math.isfinite(m_hat) and abs(m_hat - m0) < 0.5):
                    ft += 1
            ft_hi[N] = ft / trials
        self.assertLess(ft_hi[800], ft_hi[50])
        self.assertLess(ft_hi[800], 0.12)


if __name__ == "__main__":
    unittest.main()

"""Executable pure-stdlib checks for E18 (Walking all the way through Door B):
  - E18-1: Exact-covariance symmetric non-Gaussian (W_{ell,Sigma}) breaks Door A T45 (~74-78% cov);
    universal small-scalar impossibility on W_{sph} (N <= 4, or N <= 5 symmetric); Chebyshev-min /
    Gauss-unimodal two-island behavior on Can's near/far scalar scenes.
  - E18-2: Shattering of F2 uniqueness inside Door B spherical null T75 (4): Gaussian-LMP,
    Cauchy-LMP, Laplace/spatial-sign-LMP, and spherical LR disagree sharply on 10 noisy agreeing
    2D readings + 1 wild outlier (Cauchy-LMP rejects T81's point argmin K!).
  - E18-3: (a) Spherical LR T75 (4) CANNOT kill anti-aligned hilltops — even on domain='line'
    (allowing rejection of Newton!), p_line(1.0) == 1.000 on (10, -10); (b) Pure spherical T75 (4)
    collapses to < 0.02 sigma on over-agreeing readings, violating believed sigma (L18);
    (c) EIV leverage U_i breaks Bakirov-Szekely / Ibragimov-Muller for t_joint under per-reading
    heavy tails even while equal-weight Student-t stays >= 0.95!

Runnable with `python3 -m unittest -v tests/test_e18_door_b_spherical.py`.
"""
import math
import random
import unittest

import explore.e18_door_b as e18


class TestE18DoorBSpherical(unittest.TestCase):

    def test_e18_1a_exact_covariance_symmetric_nongaussian_breaks_door_a_t45(self):
        """Theorem E18-1(1): Dropping L13 Gaussianity while keeping exact believed Var(e)=1
        (W_{ell,Sigma}) breaks Door A T45 coverage down to ~74-78% for all N = 1, 2, 10, 50."""
        for N in (1, 2, 10, 50):
            res = e18.simulate_t45_vs_cheb_min_exact_cov(N, trials=8000, seed=1800 + N)
            self.assertLess(res["t45_cov"], 0.80, f"expected T45 break at N={N}: {res}")
            self.assertGreaterEqual(res["cheb_min_cov"], 0.95)

    def test_e18_1b_chebyshev_min_and_gauss_unimodal_on_cans_2scalar_scenes(self):
        """Theorem E18-1(3): Exact-variance W_{ell,Sigma} bounds work on 2 scalar readings
        (unlike W_{sph}, which is universally vacuous for N <= 5 symmetric): tight single interval
        on Can's near pair (5.5, 5.6) and two disconnected islands rejecting 5.5 on far (5.0, 6.0)."""
        near_c = e18.cheb_min_2scalar_intervals(5.5, 5.6, 0.1, unimodal=False)
        far_c = e18.cheb_min_2scalar_intervals(5.0, 6.0, 0.1, unimodal=False)
        self.assertEqual(len(near_c), 1)
        self.assertAlmostEqual(near_c[0][0], 5.5 - 0.211474, places=4)
        self.assertAlmostEqual(near_c[0][1], 5.6 + 0.211474, places=4)
        self.assertEqual(len(far_c), 2)  # disconnected islands around 5.0 and 6.0; rejects 5.5!
        self.assertLess(far_c[0][1], 5.5)
        self.assertGreater(far_c[1][0], 5.5)

        near_g = e18.cheb_min_2scalar_intervals(5.5, 5.6, 0.1, unimodal=True)
        far_g = e18.cheb_min_2scalar_intervals(5.0, 6.0, 0.1, unimodal=True)
        self.assertEqual(len(near_g), 1)
        self.assertLess(near_g[0][1] - near_g[0][0], 0.40)  # 5.55 +- 0.19
        self.assertEqual(len(far_g), 2)

    def test_e18_2_door_b_spherical_f2_uniqueness_shatters_on_10_agree_plus_1_outlier(self):
        """Theorem E18-2: Four exact admissible spherical-null tests (Gaussian-LMP, Cauchy-LMP,
        Laplace/spatial-sign-LMP, and spherical LR) give conflicting real-data verdicts on
        10 noisy agreeing 2D readings at m=5 + 1 outlier at m=50; Cauchy-LMP rejects T81's point."""
        rng = random.Random(182)
        agree = []
        for i in range(10):
            ax, ay = 2.0, 0.5 * math.cos(i)
            agree.append(((5.0 * ax + rng.gauss(0.0, 0.3), 5.0 * ay + rng.gauss(0.0, 0.3)),
                          (ax + rng.gauss(0.0, 0.1), ay + rng.gauss(0.0, 0.1))))
        data = agree + [((50.0, 0.0), (1.0, 0.0))]
        m_t81 = e18.argmin_k_2d_unit(data)
        self.assertGreater(m_t81, 12.0)  # outlier pulls T81 far away from true m=5

        p_true = e18.sph_2d_family_pvalues(data, 5.0, draws=2500, seed=18)
        p_mid = e18.sph_2d_family_pvalues(data, 6.5, draws=2500, seed=18)
        p_t81 = e18.sph_2d_family_pvalues(data, m_t81, draws=2500, seed=18)
        p_hi = e18.sph_2d_family_pvalues(data, 25.0, draws=2500, seed=18)

        # Cauchy-LMP centers at true m=5, REJECTS m=6.5 AND REJECTS T81's point, keeps m=25!
        self.assertGreater(p_true["cauchy_score"], 0.40)
        self.assertLess(p_mid["cauchy_score"], 0.01)
        self.assertLess(p_t81["cauchy_score"], 0.05)
        self.assertGreater(p_hi["cauchy_score"], 0.20)

        # Unweighted Gaussian-LMP & Spherical LR do the exact opposite: keep m=6.5 & T81, reject m=25!
        self.assertGreater(p_mid["sph_lr_pos"], 0.20)
        self.assertGreater(p_t81["sph_lr_pos"], 0.99)
        self.assertLess(p_hi["sph_lr_pos"], 0.05)

        # Spatial-sign Laplace-LMP keeps all four!
        for ps in (p_true, p_mid, p_t81, p_hi):
            self.assertGreater(ps["sign_score"], 0.15)

    def test_e18_3a_spherical_lr_keeps_anti_aligned_hilltops_even_on_whole_line(self):
        """Theorem E18-3(1): Conditioning on ||V_i|| freezes Q_VV in every spherical draw, so
        T75 (4) spherical LR CANNOT kill anti-aligned hilltops — even on domain='line' (allowing
        rejection of Newton!), p_line(1.0) == 1.0 on (10,-10), and p_pos(0.305) ~ 0.50 on (10,-3)."""
        anti10 = [((10.0, 0.0), (-10.0, 0.0)), ((10.0, 0.5), (-10.0, -0.5)), ((10.0, -0.4), (-10.0, 0.4))]
        p_hill = e18.sph_2d_family_pvalues(anti10, 1.0, draws=2000, seed=18)
        p_side = e18.sph_2d_family_pvalues(anti10, 2.0, draws=2000, seed=18)
        self.assertGreater(p_hill["sph_lr_line"], 0.99)  # hilltop survives even on domain='line'!
        self.assertLess(p_side["sph_lr_line"], 0.01)     # while non-hilltop masses are rejected!

        anti3 = [((10.0, 0.0), (-3.0, 0.0)), ((10.0, 0.5), (-3.1, -0.15))]
        p_305 = e18.sph_2d_family_pvalues(anti3, 0.305, draws=2000, seed=18)
        p_1 = e18.sph_2d_family_pvalues(anti3, 1.0, draws=2000, seed=18)
        p_inf = e18.sph_2d_family_pvalues(anti3, math.inf, draws=2000, seed=18)
        self.assertGreater(p_305["sph_lr_pos"], 0.35)    # hilltop ghost at global max K=218.6 kept!
        self.assertLess(p_1["sph_lr_pos"], 0.01)         # valley between hilltop & inf rejected
        self.assertGreater(p_inf["sph_lr_pos"], 0.35)    # T81 boundary end kept by LR
        self.assertLess(p_inf["unw_score"], 0.01)        # yet rejected by unweighted spherical score

    def test_e18_3b_pure_spherical_collapses_below_believed_sigma_on_over_agreement(self):
        """Theorem E18-3(2): Pure spherical T75 (4) discards absolute sigma magnitude (L18) and
        collapses to < 0.02 Gaussian sigma on over-agreeing readings (forcing an unforced union
        'floor' with Gaussian to restore L18)."""
        over_agree = [
            ((55.00, 0.01), (10.0, 0.00)),
            ((55.02, -0.01), (10.0, 0.00)),
            ((54.98, 0.01), (10.0, 0.00)),
            ((55.01, -0.02), (10.0, 0.00)),
        ]
        # Believed Gaussian 1-sigma is ~0.28; m0 = 5.505 is only 0.02 Gaussian sigma from 5.50!
        p_center = e18.sph_2d_family_pvalues(over_agree, 5.50, draws=2000, seed=18)
        p_tiny_off = e18.sph_2d_family_pvalues(over_agree, 5.505, draws=2000, seed=18)
        self.assertGreater(p_center["sph_lr_pos"], 0.50)
        self.assertLess(p_tiny_off["sph_lr_pos"], 0.01)

    def test_e18_3c_eiv_leverage_breaks_bakirov_szekely_for_tjoint(self):
        """Theorem E18-3(3): Under independent per-reading Cauchy scales (nu=1), equal-weight
        Student-t_{N-1} keeps >= 95% coverage (Bakirov-Szekely 2005 / Ibragimov-Muller 2010),
        whereas EIV t_joint under-covers sharply (~0.77) because outlier leverage U_i scales with sqrt(v_i)."""
        res = e18.compare_bakirov_equal_weight_vs_eiv_tjoint(N=5, nu=1.0, signal=0.5, trials=4000, seed=183)
        self.assertGreaterEqual(res["direct_equal_weight_cov"], 0.945)
        self.assertLess(res["eiv_tjoint_cov"], 0.85)


if __name__ == "__main__":
    unittest.main()

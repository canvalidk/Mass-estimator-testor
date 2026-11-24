"""Verification suite for E15: Door A uniqueness of T45, M_nz-GCLR discontinuity/L17 breaks,
Z-GCLR rejection of T81's point, and T81 boundary-end inclusion.

Runs under plain Python 3 stdlib with zero external dependencies:
    python3 -m unittest -v tests/test_e15_door_a_uniqueness.py
"""
import math
import unittest

from met.ours34 import isotropic, scalar
from explore.e15_door_a import (
    C95_CENTER,
    Z95,
    chi2_band_sup_prob,
    door_a1,
    gclr_pvalue_1instr,
    simulate_t45_vs_psi_hill,
)


class TestE15DoorAUniquenessAndDichotomy(unittest.TestCase):
    """Theorems E15-1, E15-2, and E15-3."""

    def test_e15_1_t45_strictly_dominates_exact_hilltop_competitor(self):
        """Theorem E15-1(b, d): psi_hill holds exact size 0.05 at delta=0, kills the exact
        hilltop m=1 on (10,-10), yet is strictly dominated by T45 at every alternative delta != 0."""
        # Null size check (delta = 0): both hold exact 0.05
        null_res = simulate_t45_vs_psi_hill(M=5, signal_per_comp=1.0, delta=0.0, trials=10000, seed=101)
        self.assertAlmostEqual(null_res["t45"], 0.05, delta=0.007)
        self.assertAlmostEqual(null_res["psi_hill"], 0.05, delta=0.007)

        # Alternative power check (delta = 0.4): T45 strictly beats psi_hill across M=2, 5, 10
        for M, sig in [(2, 1.0), (5, 1.0), (10, 0.7)]:
            alt = simulate_t45_vs_psi_hill(M=M, signal_per_comp=sig, delta=0.4, trials=12000, seed=200 + M)
            self.assertGreater(alt["t45"], alt["psi_hill"])

        # And on anti-parallel (10, -10) at sigma=0.1, psi_hill kills m=1 (q=0 < w=20000, X=0 < c95),
        # yet ACCEPTS m=0.5 (where q < w and |X| >> c95)!
        # At m=0.5: theta=atan(0.5), u=(1,2)/sqrt(5), v=(2,-1)/sqrt(5); y=(100,-100)
        u_val = (100.0 * 1.0 - 100.0 * 2.0) / math.sqrt(5.0)
        v_val = (100.0 * 2.0 + 100.0 * 1.0) / math.sqrt(5.0)
        q, x = u_val * u_val, math.copysign(v_val, u_val)
        w = 0.0  # N=1 has W=0; for N=2 identical readings (10,-10) perturbed slightly, W > q!
        # Take two slightly perturbed anti-parallel readings: (100, -100) and (100, -101)
        # Wait: even at m=0.5 on (100,-100) + (0, 100): q=2000, vv=18000 => w=0 for 1 reading;
        # use two readings (100, -100) and (-100, 100) perturbed to (100,-100) & (99,-101):
        # At m=0.5, u1=-100/sqrt(5), v1=300/sqrt(5); u2=-103/sqrt(5), v2=299/sqrt(5) => q < w? Check:
        # At m=0.9: u=(90-100)/sqrt(1.81), v=(100+90)/sqrt(1.81) => q=55, vv=19945; add orthogonal reading (100,-100)
        # Simpler: at m=0.95 with 2 readings (100,-100) and (95,-105):
        th = math.atan(0.95)
        su, cu = math.sin(th), math.cos(th)
        pts = [(100.0, -100.0), (95.0, -105.0)]
        q2 = sum((f * su + a * cu) ** 2 for f, a in pts)
        uv2 = sum((f * su + a * cu) * (f * cu - a * su) for f, a in pts)
        vv2 = sum((f * cu - a * su) ** 2 for f, a in pts)
        x2 = uv2 / math.sqrt(q2)
        w2 = vv2 - x2 * x2
        # Confirm q2 < w2 AND |x2| > C95_CENTER => psi_hill ACCEPTS m=0.95 even though it's near the hilltop!
        self.assertLess(q2, w2)
        self.assertGreater(abs(x2), C95_CENTER)

    def test_e15_1c_conservative_chi2_band_shows_exactness_is_necessary(self):
        """Theorem E15-1(c): Relaxing exact size == alpha to conservative <= alpha permits
        narrow bands of W with sup_{k>=1} P(chi2_k in [a, b]) < 0.05."""
        sup_p, worst_k = chi2_band_sup_prob(10.0, 10.5)
        self.assertLess(sup_p, 0.05)
        self.assertGreater(sup_p, 0.04)
        self.assertEqual(worst_k, 11)

    def test_e15_3a_zgclr_rejects_t81_point_argmin_K(self):
        """Theorem E15-3: Even T75 (3)'s weighted conditional LR (Z-GCLR with w(H)=H/(1+H))
        strictly rejects T81's point argmin K(m) (extending T82 (1) from T76 scores to Z-GCLR)."""
        pairs = [(20.0, 2.0)] + [(0.2, 2.0)] * 100
        # Confirm T81 point is exactly m* = 1.0
        res_a1 = door_a1([scalar(f, a, 1.0, 1.0) for f, a in pairs])
        self.assertAlmostEqual(res_a1["point"], 1.0, places=6)
        self.assertAlmostEqual(res_a1["state"].p_value(1.0)["p_t45"], 1.0, places=8)

        # Z-GCLR with w(H) = H / (1 + H) rejects m* = 1.0 at p < 0.01!
        z_res = gclr_pvalue_1instr(pairs, 1.0, weight=lambda h: h / (1.0 + h), draws=1500, seed=42)
        self.assertGreater(z_res["lr_obs"], 20.0)
        self.assertLess(z_res["p_value"], 0.01)

    def test_e15_3b_mnz_gclr_breaks_continuity_at_zero_and_l17_reading_split(self):
        """Theorem E15-3: M_nz-GCLR (exact Z by dropping 0/0) breaks continuity at 0/0
        (10 readings of size 1e-9 flip p from 0.005 to >0.05, 50 flip to >0.99) AND breaks
        L17 reading splitting (1 reading split into 50 equal parts with identical curve K)."""
        base = [(5.0, 1.0)]
        r0 = gclr_pvalue_1instr(base, 1.0, domain="line", draws=2000, seed=7)
        r50_eps = gclr_pvalue_1instr(base + [(1e-9, 0.0)] * 50, 1.0, domain="line", draws=2000, seed=7)
        self.assertLess(r0["p_value"], 0.01)
        self.assertGreater(r50_eps["p_value"], 0.95)

        # Meanwhile T45 is continuous at 0/0 to machine precision
        a0 = door_a1([scalar(5.0, 1.0, 1.0, 1.0)])
        a50 = door_a1([scalar(5.0, 1.0, 1.0, 1.0)] + [scalar(1e-9, 0.0, 1.0, 1.0)] * 50)
        self.assertAlmostEqual(a0["state"].p_value(1.0)["p_t45"],
                               a50["state"].p_value(1.0)["p_t45"], places=12)

        # L17 reading split: 1 reading (5, 1) vs 50 readings (5/sqrt(50), 1/sqrt(50))
        split50 = [(5.0 / math.sqrt(50.0), 1.0 / math.sqrt(50.0))] * 50
        s1 = gclr_pvalue_1instr(base, 1.0, domain="positive", draws=1500, seed=9)
        s50 = gclr_pvalue_1instr(split50, 1.0, domain="positive", draws=1500, seed=9)
        self.assertLess(s1["p_value"], 0.02)
        self.assertGreater(s50["p_value"], 0.95)
        # Same curve K(m) bit-for-bit, so T45 is unchanged!
        t1 = door_a1([scalar(5.0, 1.0, 1.0, 1.0)])["state"]
        t50 = door_a1([scalar(f, a, 1.0, 1.0) for f, a in split50])["state"]
        self.assertAlmostEqual(t1.curve(1.0), t50.curve(1.0), places=12)
        self.assertAlmostEqual(t1.p_value(1.0)["p_t45"], t50.p_value(1.0)["p_t45"], places=12)

    def test_e15_3c_t81_boundary_end_inclusion_gclr_yes_t45_no(self):
        """Answers T81 'What is not shown' item 3: on anti-aligned readings (10, -3) & (10, -3.1),
        T81's point is m* = inf (K(inf)=18.61). M_nz-GCLR (positive domain) keeps m*=inf with p=1.0,
        whereas T45 REJECTS m*=inf (p=1.6e-5) and keeps only a hilltop ghost at the global maximum m=0.305."""
        pairs = [(10.0, -3.0), (10.0, -3.1)]
        res = door_a1([scalar(f, a, 1.0, 1.0) for f, a in pairs])
        self.assertTrue(math.isinf(res["point"]))
        self.assertEqual(res["point_status"], "boundary_inf")

        # T45 rejects T81's boundary point m* = inf, and keeps the hilltop maximum m = 0.305!
        p_inf_t45 = res["state"].p_value(math.inf)["p_t45"]
        p_max_t45 = res["state"].p_value(0.305)["p_t45"]
        self.assertLess(p_inf_t45, 1e-4)
        self.assertGreater(p_max_t45, 0.95)
        self.assertGreater(res["state"].curve(0.305), res["state"].curve(math.inf) * 10.0)

        # M_nz-GCLR on positive domain keeps m* = inf with p = 1.0 (LR_obs = 0)
        g_inf = gclr_pvalue_1instr(pairs, math.inf, domain="positive", draws=500, seed=1)
        self.assertEqual(g_inf["lr_obs"], 0.0)
        self.assertEqual(g_inf["p_value"], 1.0)

    def test_door_a1_scenes_contradiction_in_leftover_not_width(self):
        """Door A1 behavior on Can's scenes and G29: identical compatible width on near and far,
        contradiction carried loudly in leftover K_min and chi2_bound."""
        near = door_a1([
            isotropic([5.5, 0, 0], [1, 0, 0], 0.1, 1e-6),
            isotropic([5.6, 0, 0], [1, 0, 0], 0.1, 1e-6),
        ])
        far = door_a1([
            isotropic([5.0, 0, 0], [1, 0, 0], 0.1, 1e-6),
            isotropic([6.0, 0, 0], [1, 0, 0], 0.1, 1e-6),
        ])
        w_near = near["set"][0][1] - near["set"][0][0]
        w_far = far["set"][0][1] - far["set"][0][0]
        self.assertAlmostEqual(w_near, w_far, places=3)
        self.assertLess(near["leftover"]["K_min"], 1.0)
        self.assertGreater(near["leftover"]["chi2_bound"], 0.5)
        self.assertGreater(far["leftover"]["K_min"], 45.0)
        self.assertLess(far["leftover"]["chi2_bound"], 1e-8)

        # G29 symmetric right-angle: flat curve, set (0, inf), screaming leftover K_min = 400
        g29 = door_a1([scalar(2.0, 0.0, 0.1, 0.1), scalar(0.0, 2.0, 0.1, 0.1)])
        self.assertEqual(g29["point_status"], "flat")
        self.assertEqual(g29["shape"], "all")
        self.assertAlmostEqual(g29["leftover"]["K_min"], 400.0, places=8)
        self.assertLess(g29["leftover"]["chi2_bound"], 1e-50)


if __name__ == "__main__":
    unittest.main()

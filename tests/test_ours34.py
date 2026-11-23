"""ours34, the probe (met/ours34.py; explore/e14-ours34-walls): construction, premise report,
and explicit adversarial breaks (`GOAL.md` standing rule 1).

Runs under plain Python 3 stdlib with zero external dependencies:
    python3 -m unittest -v tests/test_ours34.py
"""
import math
import random
import unittest

from met.ours34 import (
    Reading,
    State,
    common_force,
    isotropic,
    lump,
    ours34,
    scalar,
    _cholesky,
)


def _random_spd(rng, dim):
    A = [[rng.gauss(0.0, 1.0) for _ in range(dim)] for _ in range(dim)]
    return tuple(
        tuple(sum(A[i][t] * A[j][t] for t in range(dim)) + (0.5 if i == j else 0.0) for j in range(dim))
        for i in range(dim)
    )


class _DuckFullReading:
    """Mimics `met.full.Reading` (numpy-style `.tolist()`, d x k subspace, optional window)."""

    def __init__(self, force, acceleration, cov, subspace=None, window=None):
        self.force = force
        self.acceleration = acceleration
        self.cov = cov
        self.subspace = subspace
        self.window = window


class TestOurs34ConstructionAndSymmetries(unittest.TestCase):
    """Exact algebraic properties of the zero-QR Gram engine and additive state."""

    def test_duck_typed_full_reading_and_additive_state(self):
        rng = random.Random(42)
        cov1 = _random_spd(rng, 6)
        cov2 = _random_spd(rng, 6)
        # 3 x 2 orthonormal subspace (columns [1,0,0] and [0, 0.6, 0.8])
        N_dxk = ((1.0, 0.0), (0.0, 0.6), (0.0, 0.8))
        r1 = _DuckFullReading([3.2, -1.5, 4.1], [1.1, -0.4, 2.0], cov1, N_dxk, window=(0.0, 1.0))
        r2 = _DuckFullReading([2.4, 0.8, -1.1], [0.9, 0.3, -0.4], cov2)

        s1 = State([r1])
        s2 = State([r2])
        s12 = s1 + s2
        s_direct = State([r1, r2])
        self.assertEqual(s12.M_nz, s_direct.M_nz)
        self.assertEqual(s12.M_nz, 5)
        for m in [0.0, 0.5, 2.3, 10.0, math.inf]:
            self.assertAlmostEqual(s12.curve(m), s_direct.curve(m), places=12)
            self.assertAlmostEqual(s12.p_value(m)["p_value"], s_direct.p_value(m)["p_value"], places=12)

        # Verify S_m == -0.5 * dK/dm on full correlated 6x6 + 2D subspace
        m0, eps = 2.3, 1e-6
        _, s_analytic, _ = s1.ksc_m(m0)
        s_fd = -0.5 * (s1.curve(m0 + eps) - s1.curve(m0 - eps)) / (2.0 * eps)
        self.assertAlmostEqual(s_analytic, s_fd, places=7)

    def test_l3_units_and_t5_swap_exact(self):
        rd = [
            isotropic([4.0, 0.3, -0.2], [2.0, 0.1, -0.1], 0.5, 0.2),
            isotropic([6.1, -0.4, 0.5], [2.9, -0.2, 0.2], 0.4, 0.3),
        ]
        base = ours34(rd)
        lam = 3.5
        scaled = ours34([r.scaled(lam) for r in rd])
        self.assertAlmostEqual(scaled["point"], lam * base["point"], places=8)
        self.assertAlmostEqual(scaled["set"][0][0], lam * base["set"][0][0], places=5)
        self.assertAlmostEqual(scaled["set"][0][1], lam * base["set"][0][1], places=5)
        self.assertAlmostEqual(scaled["leftover"]["K_min"], base["leftover"]["K_min"], places=10)

        swapped = ours34([r.swapped() for r in rd])
        self.assertAlmostEqual(swapped["point"], 1.0 / base["point"], places=8)
        self.assertAlmostEqual(swapped["set"][0][0], 1.0 / base["set"][0][1], places=5)
        self.assertAlmostEqual(swapped["set"][0][1], 1.0 / base["set"][0][0], places=5)
        self.assertAlmostEqual(swapped["leftover"]["K_min"], base["leftover"]["K_min"], places=10)

    def test_validation_rejects_non_spd_and_non_orthonormal(self):
        with self.assertRaises(ValueError):
            Reading((1.0,), (1.0,), ((1.0, 2.0), (2.0, 1.0)))  # indefinite cov
        with self.assertRaises(ValueError):
            Reading((1.0, 0.0), (1.0, 0.0), ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
                    subspace=((2.0,), (0.0,)))  # non-unit column


class TestOurs34PremisePasses(unittest.TestCase):
    """Premise rows and historical probe breakers that ours34 passes."""

    def test_all_zero_and_exact_z_invariance(self):
        z_res = ours34([scalar(0.0, 0.0, 1.0, 1.0)] * 3)
        self.assertEqual(z_res["point_status"], "flat")
        self.assertTrue(z_res["point_is_convention"])
        self.assertEqual(z_res["shape"], "all")
        self.assertEqual(z_res["set"], ((0.0, math.inf),))
        self.assertEqual(z_res["leftover"]["K_min"], 0.0)
        self.assertEqual(z_res["leftover"]["M_nz"], 0)

        base_rd = [
            isotropic([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
            isotropic([6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
        ]
        z_3d = isotropic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1e-6)
        r_base = ours34(base_rd)
        r_pad = ours34(base_rd + [z_3d] * 25)
        self.assertAlmostEqual(r_pad["point"], r_base["point"], places=12)
        self.assertAlmostEqual(r_pad["set"][0][0], r_base["set"][0][0], places=10)
        self.assertAlmostEqual(r_pad["set"][0][1], r_base["set"][0][1], places=10)
        self.assertEqual(r_pad["leftover"]["M_nz"], r_base["leftover"]["M_nz"])

    def test_t80_unbounded_one_sided_and_interval_shapes(self):
        # T73 (1): F = 10 +/- 1, a = 1 +/- 1 => [3.28, inf)
        r = ours34([scalar(10.0, 1.0, 1.0, 1.0)])
        self.assertAlmostEqual(r["point"], 10.0, places=8)
        self.assertEqual(r["shape"], "at_least")
        self.assertAlmostEqual(r["set"][0][0], 3.2797, places=3)
        self.assertTrue(math.isinf(r["set"][0][1]))

    def test_cans_3d_scenes_t30_sigma_invariance_and_t78_conflict(self):
        near = ours34([
            isotropic([5.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
            isotropic([5.6, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
        ])
        self.assertEqual(near["shape"], "interval")
        self.assertAlmostEqual(near["set"][0][0], near["raw_set"][0][0], places=4)
        self.assertAlmostEqual(near["set"][0][1], near["raw_set"][0][1], places=4)

        far_01 = ours34([
            isotropic([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
            isotropic([6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-6),
        ])
        far_001 = ours34([
            isotropic([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.001, 1e-8),
            isotropic([6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.001, 1e-8),
        ])
        self.assertAlmostEqual(far_01["set"][0][0], 4.925, places=2)
        self.assertAlmostEqual(far_01["set"][0][1], 6.075, places=2)
        # Bit-for-bit sigma-invariance across 100x sigma reduction (T30)
        self.assertAlmostEqual(far_01["set"][0][0], far_001["set"][0][0], places=4)
        self.assertAlmostEqual(far_01["set"][0][1], far_001["set"][0][1], places=4)

        # T78 (2) 3D conflict {0.2, 1.0, 1.4}: sigma-invariant interval [0.601, 1.664]
        c_01 = ours34([
            isotropic([0.2, 0, 0], [1, 0, 0], 0.01, 0.01),
            isotropic([1.0, 0, 0], [1, 0, 0], 0.01, 0.01),
            isotropic([1.4, 0, 0], [1, 0, 0], 0.01, 0.01),
        ])
        c_001 = ours34([
            isotropic([0.2, 0, 0], [1, 0, 0], 0.001, 0.001),
            isotropic([1.0, 0, 0], [1, 0, 0], 0.001, 0.001),
            isotropic([1.4, 0, 0], [1, 0, 0], 0.001, 0.001),
        ])
        self.assertAlmostEqual(c_01["set"][0][0], c_001["set"][0][0], places=3)
        self.assertAlmostEqual(c_01["set"][0][1], c_001["set"][0][1], places=3)

    def test_passes_t39_and_g29_ours33_breakers(self):
        # T39 opposed pair (where ours33 blew up to >60 decades!): finite interval covering 2.0
        th0 = math.atan(2.0)
        u0, u1 = math.sin(th0), math.cos(th0)
        v0, v1 = math.cos(th0), -math.sin(th0)
        t39 = ours34([
            scalar(1000.0 * u0 + 3.0 * v0, 1000.0 * u1 + 3.0 * v1, 1.0, 1.0),
            scalar(1000.0 * u0 - 3.0 * v0, 1000.0 * u1 - 3.0 * v1, 1.0, 1.0),
        ])
        self.assertAlmostEqual(t39["point"], 2.0, places=3)
        self.assertEqual(t39["shape"], "interval")
        lo, hi = t39["set"][0]
        self.assertTrue(0.0 < lo < 2.0 < hi < 5.0)
        self.assertLess(math.log10(hi / lo), 2.0)  # vs >60 decades for ours33!

        # G29 right-angle pair (where ours33 diverged to infinity!): returns (0, inf)
        g29 = ours34([scalar(2.0, 0.0, 0.1, 0.1), scalar(0.0, 2.0, 0.1, 0.1)])
        self.assertEqual(g29["point_status"], "flat")
        self.assertEqual(g29["shape"], "all")

    def test_l5_and_l14_composition_discordant_asymmetric_and_boundary(self):
        # Discordant Body 1 + Clean Body 2: Q_joint preserves widening without dilution
        b1_disc = ours34([
            isotropic([5.0, 0, 0], [1, 0, 0], 0.01, 1e-7),
            isotropic([6.0, 0, 0], [1, 0, 0], 0.01, 1e-7),
        ])
        b2_clean = ours34([isotropic([10.0, 0, 0], [1, 0, 0], 0.01, 1e-7)] * 5)
        comp = lump(b1_disc, b2_clean, method="joint")
        self.assertAlmostEqual(comp["point"], 15.5, places=6)
        self.assertEqual(comp["shape"], "interval")
        self.assertLess(comp["set"][0][0], 14.9)
        self.assertGreater(comp["set"][0][1], 16.1)

        # Asymmetric 100 + 10: Q_joint is completely ghost-free (around 110)
        b1_asym = ours34([scalar(100.0, 1.0, 0.01, 0.01)])
        b2_asym = ours34([scalar(10.0, 1.0, 0.01, 0.01)])
        comp_asym = lump(b1_asym, b2_asym, method="joint")
        self.assertAlmostEqual(comp_asym["point"], 110.0, places=5)
        self.assertEqual(comp_asym["shape"], "interval")
        self.assertGreater(comp_asym["set"][0][0], 105.0)
        self.assertLess(comp_asym["set"][0][1], 115.0)

        # Boundary chaining: (F=10, a=0) [pt=inf] with (3 +/- 0.1)
        b_inf = ours34([scalar(10.0, 0.0, 1.0, 1.0)])
        b_3 = ours34([scalar(3.0, 1.0, 0.1, 1e-6)])
        c_sum = lump(b_inf, b_3, method="joint")
        self.assertTrue(math.isinf(c_sum["point"]))
        self.assertEqual(c_sum["shape"], "at_least")
        self.assertGreater(c_sum["set"][0][0], 6.5)

        c_par = common_force(b_inf, b_3, method="joint")
        self.assertAlmostEqual(c_par["point"], 3.0, places=6)
        # Main finite interval sits right around [1.70, 3.25] (plus the [0, 0] pinhole!)
        finite_seg = [seg for seg in c_par["set"] if seg[1] > 1.0][0]
        self.assertAlmostEqual(finite_seg[0], 1.703, places=2)
        self.assertAlmostEqual(finite_seg[1], 3.245, places=2)

        # (5 +/- 0.5) + (0/0)
        b_5 = ours34([scalar(5.0, 1.0, 0.5, 1e-6)])
        b_0 = ours34([scalar(0.0, 0.0, 1.0, 1.0)])
        c_null_j = lump(b_5, b_0, method="joint")
        self.assertEqual(c_null_j["shape"], "at_least")
        self.assertAlmostEqual(c_null_j["set"][0][0], 4.02, places=2)
        c_null_q = lump(b_5, b_0, method="qpos")
        self.assertEqual(c_null_q["shape"], "at_least")
        self.assertAlmostEqual(c_null_q["set"][0][0], 4.178, places=2)


class TestOurs34DocumentedBreaks(unittest.TestCase):
    """Explicit adversarial breaks of ours34 — each test asserts the break exists so
    every flaw is recorded and never silently forgotten (`GOAL.md` standing rule 1)."""

    def test_break_1_hilltop_and_right_angle_pinhole_ghosts(self):
        """Break 1 (score-test limitation): stationary hilltops (K'(m)=0) and right-angle
        pinholes (U(m)=0) give p=1 even where K(m) is maximized."""
        # Anti-parallel (10, -10) at sig=0.1: keeps a razor-thin ghost at the MAXIMUM m=1!
        anti = ours34([isotropic([10.0, 0, 0], [-10.0, 0, 0], 0.1, 0.1)])
        self.assertAlmostEqual(anti["state"].p_value(1.0)["p_value"], 1.0, places=12)
        self.assertGreater(anti["state"].curve(1.0), anti["state"].curve(0.0) * 1.9)

        # Bimodal 2-instrument conflict: raw_set has 3 islands including a hilltop ghost at m=1!
        bimodal = ours34([scalar(0.1, 1.0, 0.01, 0.1), scalar(1.0, 0.1, 0.1, 0.01)])
        self.assertEqual(len(bimodal["raw_set"]), 3)
        self.assertTrue(bimodal["raw_set"][1][0] < 1.0 < bimodal["raw_set"][1][1])

        # One-sided F=10, a=0: keeps isolated right-angle pinhole at m=0 alongside [5.003, inf)
        one_sided = ours34([scalar(10.0, 0.0, 1.0, 1.0)])
        self.assertEqual(one_sided["shape"], "disconnected")
        self.assertEqual(one_sided["set"][0][0], 0.0)
        self.assertLess(one_sided["set"][0][1], 1e-10)
        self.assertAlmostEqual(one_sided["set"][1][0], 5.003, places=2)

    def test_break_2_discontinuity_at_zero_via_M_nz(self):
        """Break 2 (Theorem E14-6): exact 0/0 leaves M_nz unchanged, whereas 20 readings of
        size 1e-9 jump M_nz from 6 to 66 and collapse the far-pair widening back to T45!"""
        base = [
            isotropic([5.0, 0, 0], [1.0, 0, 0], 0.1, 1e-6),
            isotropic([6.0, 0, 0], [1.0, 0, 0], 0.1, 1e-6),
        ]
        r_exact = ours34(base + [isotropic([0, 0, 0], [0, 0, 0], 0.1, 1e-6)] * 20)
        r_eps = ours34(base + [isotropic([1e-9, 0, 0], [0, 0, 0], 0.1, 1e-6)] * 20)
        w_exact = r_exact["set"][0][1] - r_exact["set"][0][0]
        w_eps = r_eps["set"][0][1] - r_eps["set"][0][0]
        self.assertGreater(w_exact, 1.1)     # widened (~1.15)
        self.assertLess(w_eps, 0.30)         # collapsed to T45 (~0.278)!

    def test_break_3_qpos_fake_split_and_discordant_collapse_vs_qjoint_2dof_tax(self):
        """Break 3 (Wall 3 composition trade-off):
        - 1-dof linearized Q_pos keeps fake split M=30 on asymmetric 100+10 (p=1.0) and
          collapses on discordant+clean bodies;
        - 2-dof Q_joint kills both breaks, at the cost of +25% width (2.448 vs 1.960 sigma)."""
        b1_asym = ours34([scalar(100.0, 1.0, 0.01, 0.01)])
        b2_asym = ours34([scalar(10.0, 1.0, 0.01, 0.01)])
        # Q_pos breaks at M = 30 (70 kg away from true 110!):
        c_qpos = lump(b1_asym, b2_asym, method="qpos")
        self.assertGreater(c_qpos["p_value"](30.0), 0.99)
        self.assertLess(c_qpos["set"][0][0], 20.0)

        # Q_joint pays the exact 2-dof tax at strong signal: sqrt(chi2_2(0.95)) / z_0.975 = 2.4477 / 1.9600
        b1_g = ours34([scalar(2.0, 1.0, 0.01, 1e-6)])
        b2_g = ours34([scalar(3.0, 1.0, 0.01, 1e-6)])
        c_j = lump(b1_g, b2_g, method="joint")
        hw = 0.5 * (c_j["set"][0][1] - c_j["set"][0][0])
        se = math.hypot(0.01, 0.01)
        self.assertAlmostEqual(hw / se, math.sqrt(-2.0 * math.log(0.05)), places=2)  # 2.448 sigma

    def test_break_4_per_reading_heavy_tails_undercoverage(self):
        """Break 4 (Wall 2 common-scale vs per-reading scale): under independent per-reading
        v_i ~ InvChi2(2), ours34's common-scale floor under-covers (below 0.90), while
        maintaining >= 0.95 under common scale v ~ InvChi2(2)."""
        rng = random.Random(2026)
        n_tr = 160
        m_true = 2.5
        cov_common = cov_per_reading = 0
        for _ in range(n_tr):
            v_c = 2.0 / (rng.gauss(0, 1)**2 + rng.gauss(0, 1)**2)
            rd_c, rd_p = [], []
            for idx in range(5):
                sf, sa = 1.0 + 0.4 * idx, 1.0
                v_i = 2.0 / (rng.gauss(0, 1)**2 + rng.gauss(0, 1)**2)
                g = [rng.gauss(0, 1) for _ in range(3)]
                gn = math.sqrt(sum(x * x for x in g))
                asc = 3.0 / math.sqrt((m_true / sf)**2 + (1.0 / sa)**2)
                at = [asc * x / gn for x in g]
                ft = [m_true * x for x in at]
                ef = [rng.gauss(0, sf) for _ in range(3)]
                ea = [rng.gauss(0, sa) for _ in range(3)]
                sc, sp = math.sqrt(v_c), math.sqrt(v_i)
                rd_c.append(isotropic([ft[j] + sc * ef[j] for j in range(3)],
                                      [at[j] + sc * ea[j] for j in range(3)], sf, sa))
                rd_p.append(isotropic([ft[j] + sp * ef[j] for j in range(3)],
                                      [at[j] + sp * ea[j] for j in range(3)], sf, sa))
            cov_common += (State(rd_c).p_value(m_true)["p_value"] > 0.05)
            cov_per_reading += (State(rd_p).p_value(m_true)["p_value"] > 0.05)
        self.assertGreaterEqual(cov_common / n_tr, 0.93)
        self.assertLess(cov_per_reading / n_tr, 0.91)


if __name__ == "__main__":
    unittest.main()

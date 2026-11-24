"""Executable premise, adversarial, coverage, and known-cost suite for `met/ours35.py` (Door A1).

Exploration branch `explore/e16-ours35-door-a` (24 Nov 2026).
Pure Python 3 stdlib — runnable with `python3 -m unittest -v tests/test_ours35.py` or `pytest`.
"""
import math
import random
import unittest

import met.ours34 as o34
import met.ours35 as o35


def _rand_spd(rng, n):
    A = [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(n)]
    return [
        [math.fsum(A[i][t] * A[j][t] for t in range(n)) / n + (0.25 if i == j else 0.0) for j in range(n)]
        for i in range(n)
    ]


def _rand_ortho_cols(rng, d, k):
    cols = []
    for _ in range(k):
        v = [rng.gauss(0.0, 1.0) for _ in range(d)]
        for q in cols:
            r = math.fsum(a * b for a, b in zip(v, q))
            v = [v[i] - r * q[i] for i in range(d)]
        nrm = math.sqrt(math.fsum(x * x for x in v))
        cols.append([x / nrm for x in v])
    return [[cols[j][i] for j in range(k)] for i in range(d)]


def _chol_sample(rng, cov):
    L = o35._cholesky(cov)
    z = [rng.gauss(0.0, 1.0) for _ in range(len(cov))]
    return [math.fsum(L[i][j] * z[j] for j in range(i + 1)) for i in range(len(cov))]


def _in_set(intervals, m, tol=1e-12):
    for lo, hi in intervals:
        if lo - tol <= m <= hi + tol:
            return True
    return False


class TestOurs35DoorA(unittest.TestCase):

    def test_01_cross_check_against_ours34_full_type(self):
        """Curve K, score S_phi, conditional variance C_phi, Z_T45, and argmin_K match met.ours34
        (which is cross-checked against main's met.evidence and met.conditional_score) to ~1e-13
        across 20 random mixed-dimension correlated-Sigma subspace-N series."""
        rng = random.Random(3501)
        for _ in range(20):
            rs35, rs34 = [], []
            for _ in range(rng.randint(1, 5)):
                d = rng.randint(1, 3)
                k = rng.randint(1, d)
                cov = _rand_spd(rng, 2 * d)
                sub = _rand_ortho_cols(rng, d, k)
                a = [rng.gauss(0.0, 3.0) for _ in range(d)]
                f = [2.3 * x + rng.gauss(0.0, 0.5) for x in a]
                rs35.append(o35.Reading(f, a, cov, sub))
                rs34.append(o34.Reading(f, a, cov, sub))
            st35, st34 = o35.State(rs35), o34.State(rs34)
            for m in (0.0, 0.25, 1.0, 2.3, 7.5, math.inf):
                k35, s35, c35 = st35.ksc_phi(m)
                k34, s34, c34 = st34.ksc_phi(m)
                self.assertAlmostEqual(k35, k34, delta=1e-12 * max(1.0, abs(k34)))
                if c34 > 1e-20:
                    self.assertAlmostEqual(s35 / math.sqrt(c35), s34 / math.sqrt(c34), delta=1e-11)
            r35 = o35.ours35(rs35)
            r34 = o34.ours34(rs34)
            if math.isfinite(r35["point"]) and math.isfinite(r34["point"]):
                self.assertAlmostEqual(r35["point"], r34["point"], delta=1e-10 * max(1.0, r34["point"]))
            self.assertAlmostEqual(r35["leftover"], r34["leftover"]["K_min"], delta=1e-11 * max(1.0, r34["leftover"]["K_min"]))

    def test_02_exact_zero_noise_branch_t1(self):
        """Exact branch T1 for sigma = 0 handles one mass, every mass, and contradiction."""
        self.assertEqual(o35.exact_mass((0.0, 0.0), (0.0, 0.0)), ("every", None))
        tag, m = o35.exact_mass((6.0, -9.0), (2.0, -3.0))
        self.assertEqual(tag, "one")
        self.assertAlmostEqual(m, 3.0, places=14)
        self.assertEqual(o35.exact_mass((-6.0,), (2.0,))[0], "none")
        self.assertEqual(o35.exact_mass((1.0, 0.0), (0.0, 1.0))[0], "none")
        self.assertEqual(o35.exact_mass((1.0,), (0.0,))[0], "none")
        self.assertEqual(o35.exact_mass((0.0,), (1.0,))[0], "none")
        self.assertEqual(o35.exact_series([((0.0,), (0.0,)), ((6.0,), (2.0,)), ((9.0,), (3.0,))]), ("one", 3.0))
        self.assertEqual(o35.exact_series([((6.0,), (2.0,)), ((8.0,), (2.0,))])[0], "none")
        self.assertEqual(o35.exact_series([((0.0,), (0.0,))]), ("every", None))

    def test_03_invariances_l3_t5_l4_order_window_and_extreme_scale(self):
        """Exact L3 unit scaling (including 1e+-200 ratio!), T5 channel swap, L4 3D rotation,
        L16 reading-order invariance, and L18 window indifference."""
        cov3 = _rand_spd(random.Random(42), 6)
        sub = ((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))
        r1 = o35.Reading((4.0, -2.0, 1.0), (1.5, -0.8, 0.3), cov3, subspace=sub, window=(0.0, 1.0))
        r2 = o35.scalar(7.0, 2.5, 1.2, 0.4, corr=0.35, window=(1.0, 2.0))
        base = o35.ours35([r1, r2])

        # L3 unit scaling by lam = 13.7
        lam = 13.7
        sc = o35.ours35([r1.scaled(lam), r2.scaled(lam)])
        self.assertAlmostEqual(sc["point"], lam * base["point"], delta=1e-11 * lam * base["point"])
        self.assertAlmostEqual(sc["leftover"], base["leftover"], delta=1e-12 * max(1.0, base["leftover"]))
        self.assertEqual(len(sc["set"]), len(base["set"]))
        for (a, b), (sa, sb) in zip(base["set"], sc["set"]):
            self.assertAlmostEqual(sa, lam * a, delta=1e-11 * max(1.0, lam * a))
            self.assertAlmostEqual(sb, lam * b, delta=1e-11 * max(1.0, lam * b))

        # Extreme L3 dynamic range: s = 1e+200 and s = 1e-200 (10^400 span, zero overflow/underflow)
        ref = o35.ours35([o35.scalar(15.0, 5.0, 1.0, 1.0)])
        hi_s = o35.ours35([o35.scalar(15.0e100, 5.0e-100, 1e100, 1e-100)])
        lo_s = o35.ours35([o35.scalar(15.0e-100, 5.0e100, 1e-100, 1e100)])
        self.assertAlmostEqual(hi_s["point"] / 1e200, ref["point"], places=12)
        self.assertAlmostEqual(lo_s["point"] * 1e200, ref["point"], places=12)
        self.assertAlmostEqual(hi_s["set"][0][0] / 1e200, ref["set"][0][0], places=12)
        self.assertAlmostEqual(hi_s["set"][0][1] / 1e200, ref["set"][0][1], places=12)
        self.assertAlmostEqual(lo_s["set"][0][0] * 1e200, ref["set"][0][0], places=12)
        self.assertAlmostEqual(lo_s["set"][0][1] * 1e200, ref["set"][0][1], places=12)

        # T5 channel swap: m <-> 1/m
        sw = o35.ours35([r1.swapped(), r2.swapped()])
        self.assertAlmostEqual(sw["point"], 1.0 / base["point"], delta=1e-11 / base["point"])
        self.assertAlmostEqual(sw["leftover"], base["leftover"], delta=1e-12 * max(1.0, base["leftover"]))
        inv_set = tuple(sorted((1.0 / hi, 1.0 / lo) for lo, hi in base["set"]))
        for (a, b), (sa, sb) in zip(inv_set, sw["set"]):
            self.assertAlmostEqual(sa, a, delta=1e-11 * max(1.0, a))
            self.assertAlmostEqual(sb, b, delta=1e-11 * max(1.0, b))

        # L4 3D rotation invariance on correlated 3D reading + subspace
        R = _rand_ortho_cols(random.Random(99), 3, 3)
        rot = o35.ours35([r1.rotated(R), r2])
        self.assertAlmostEqual(rot["point"], base["point"], delta=1e-11 * base["point"])
        self.assertAlmostEqual(rot["leftover"], base["leftover"], delta=1e-12 * max(1.0, base["leftover"]))
        for (a, b), (ra, rb) in zip(base["set"], rot["set"]):
            self.assertAlmostEqual(ra, a, delta=1e-11 * max(1.0, a))
            self.assertAlmostEqual(rb, b, delta=1e-11 * max(1.0, b))

        # L16 order invariance and L18 window indifference
        r1_win = o35.Reading(r1.force, r1.acceleration, r1.cov, r1.subspace, window=(10.0, 99.0))
        rev = o35.ours35([r2, r1_win])
        self.assertAlmostEqual(rev["point"], base["point"], delta=1e-12)
        self.assertAlmostEqual(rev["leftover"], base["leftover"], delta=1e-12)
        self.assertAlmostEqual(rev["set"][0][0], base["set"][0][0], delta=1e-12)
        self.assertAlmostEqual(rev["set"][0][1], base["set"][0][1], delta=1e-12)

    def test_04_exact_z_continuity_at_zero_and_l17_reading_splitting(self):
        """Every output field is exact Z (L20), continuous at 0/0 (immune to 50x 1e-9 blowup
        that breaks M_nz-GCLR and ours34's t_joint), and exact under L17 reading splitting."""
        r = o35.scalar(5.0, 1.0, 1.0, 1.0)
        base = o35.ours35([r])

        # Exact Z: adding 50 exact 0/0 readings with arbitrary sigma changes NOTHING
        zeros = [o35.scalar(0.0, 0.0, 100.0, 0.01) for _ in range(50)]
        with_z = o35.ours35([r] + zeros)
        self.assertEqual(with_z["point"], base["point"])
        self.assertEqual(with_z["set"], base["set"])
        self.assertEqual(with_z["leftover"], base["leftover"])
        self.assertEqual(with_z["state"].p_value(1.0), base["state"].p_value(1.0))

        # Continuity at 0/0: 50 readings of (1e-9, 1e-9) barely nudge anything (< 1e-12)
        tiny = [o35.scalar(1e-9, 1e-9, 1.0, 1.0) for _ in range(50)]
        with_tiny = o35.ours35([r] + tiny)
        self.assertAlmostEqual(with_tiny["state"].p_value(1.0), base["state"].p_value(1.0), delta=1e-12)
        self.assertAlmostEqual(with_tiny["point"], base["point"], delta=1e-12)
        self.assertAlmostEqual(with_tiny["set"][0][0], base["set"][0][0], delta=1e-11)
        self.assertAlmostEqual(with_tiny["leftover"], base["leftover"], delta=1e-12)

        # Exact L17 reading splitting: 1 reading == 50 equal split parts (same K, same C)
        for n in (10, 50):
            s = math.sqrt(n)
            split_rs = [o35.scalar(5.0 / s, 1.0 / s, 1.0, 1.0) for _ in range(n)]
            res_split = o35.ours35(split_rs)
            self.assertAlmostEqual(res_split["point"], base["point"], delta=1e-12)
            self.assertAlmostEqual(res_split["leftover"], base["leftover"], delta=1e-12)
            self.assertAlmostEqual(res_split["state"].p_value(1.0), base["state"].p_value(1.0), delta=1e-12)
            self.assertAlmostEqual(res_split["set"][0][0], base["set"][0][0], delta=1e-11)
            self.assertEqual(math.isinf(res_split["set"][0][1]), math.isinf(base["set"][0][1]))

    def test_05_pure_additive_monoid_state_no_hidden_counts(self):
        """State = (K, C) is an exact associative additive monoid with 0/0 identity and NO M/M_nz."""
        b1 = [o35.scalar(6.0, 2.0, 1.0, 0.5)]
        b2 = [o35.isotropic((4.0, -2.0), (1.4, -0.6), 0.8, 0.3)]
        b3 = [o35.scalar(0.0, 0.0, 10.0, 10.0)]
        s1, s2, s3 = o35.State(b1), o35.State(b2), o35.State(b3)
        self.assertFalse(hasattr(s1, "M_nz"))
        self.assertFalse(hasattr(s1, "M_total"))
        self.assertTrue(s3.is_zero)

        comb_left = (s1 + s2) + s3
        comb_right = s1 + (s2 + s3)
        direct = o35.State(b1 + b2 + b3)
        for m in (0.0, 0.5, 2.0, 3.0, 10.0, math.inf):
            for st in (comb_left, comb_right):
                k1, sc1, c1 = st.ksc_phi(m)
                kd, scd, cd = direct.ksc_phi(m)
                self.assertAlmostEqual(k1, kd, delta=1e-13 * max(1.0, kd))
                self.assertAlmostEqual(sc1, scd, delta=1e-12 * max(1.0, abs(scd)))
                self.assertAlmostEqual(c1, cd, delta=1e-12 * max(1.0, cd))

    def test_06_extreme_snr_tiny_signal_and_flat_zero(self):
        """Resolves ultra-narrow SNR=1e10 intervals analytically, locates argmin_K at tiny
        signal 1e-10 without falsely flagging flat, and handles empty / 0/0 cleanly."""
        # Ultra-strong SNR = 1e10
        res_hi = o35.ours35([o35.scalar(3.0e10, 1.0e10, 1.0, 1.0)])
        self.assertAlmostEqual(res_hi["point"], 3.0, places=12)
        self.assertEqual(res_hi["shape"], "interval")
        lo, hi = res_hi["set"][0]
        # Exact theoretical 95% T45 width at F=3e10, a=1e10, sf=sa=1: 2 * z95 * sqrt(1 + m^2) / a
        z95 = 1.959963984540054
        expected_width = 2.0 * z95 * math.sqrt(10.0) / 1.0e10
        self.assertAlmostEqual(hi - lo, expected_width, delta=1e-6 * expected_width)

        # Ultra-tiny signal (3e-10, 1e-10): NOT flat!
        res_lo = o35.ours35([o35.scalar(3.0e-10, 1.0e-10, 1.0, 1.0)])
        self.assertEqual(res_lo["point_status"], "interior")
        self.assertFalse(res_lo["point_is_convention"])
        self.assertAlmostEqual(res_lo["point"], 3.0, places=12)
        self.assertEqual(res_lo["shape"], "all")

        # Empty and 0/0
        for rs in ([], [o35.scalar(0.0, 0.0, 2.0, 0.5)]):
            rz = o35.ours35(rs)
            self.assertTrue(rz["point_is_convention"])
            self.assertEqual(rz["point_status"], "flat")
            self.assertEqual(rz["set"], ((0.0, math.inf),))
            self.assertEqual(rz["shape"], "all")
            self.assertEqual(rz["leftover"], 0.0)

    def test_07_off_grid_pinhole_equivariance_and_narrow_hilltop(self):
        """Off-grid single-reading pinhole (10, -7.13) is resolved and exactly L3/T5 equivariant,
        and ultra-narrow high-SNR hilltop (1000,-713)&(1000,-714) is never falsely dropped."""
        r_pin = o35.scalar(10.0, -7.13, 1.0, 1.0)
        res = o35.ours35([r_pin])
        self.assertEqual(res["point_status"], "boundary_inf")
        self.assertEqual(len(res["set"]), 1)
        self.assertAlmostEqual(res["set"][0][0], 0.713, places=11)
        self.assertAlmostEqual(res["set"][0][1], 0.713, places=11)

        res_sc = o35.ours35([r_pin.scaled(13.7)])
        self.assertAlmostEqual(res_sc["set"][0][0], 0.713 * 13.7, places=10)
        res_sw = o35.ours35([r_pin.swapped()])
        self.assertAlmostEqual(res_sw["set"][0][0], 1.0 / 0.713, places=10)

        # Ultra-narrow hilltop width ~0.0011 on anti-aligned high-SNR pair
        res_hill = o35.ours35([o35.scalar(1000.0, -713.0, 1.0, 1.0), o35.scalar(1000.0, -714.0, 1.0, 1.0)])
        self.assertEqual(res_hill["point_status"], "boundary_inf")
        self.assertEqual(res_hill["shape"], "interval")
        lo, hi = res_hill["set"][0]
        self.assertTrue(0.713 < lo < 0.7135 < hi < 0.714)
        self.assertGreater(res_hill["leftover"], 1e6)

    def test_08_exact_frequentist_coverage_full_type_correlated_mixed(self):
        """Uniform 95% exact coverage across null, weak, and strong signal on correlated 3D (k=2)
        + correlated 1D mixed-instrument series, with and without 0/0 contamination."""
        cov3 = (
            (1.0, 0.3, -0.2, 0.25, 0.0, 0.1),
            (0.3, 1.4, 0.1, 0.0, -0.2, 0.05),
            (-0.2, 0.1, 0.8, 0.1, 0.05, 0.15),
            (0.25, 0.0, 0.1, 0.5, 0.1, -0.05),
            (0.0, -0.2, 0.05, 0.1, 0.7, 0.08),
            (0.1, 0.05, 0.15, -0.05, 0.08, 0.6),
        )
        cov1 = ((4.0, 0.6), (0.6, 0.25))
        sub = ((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))
        z_read = o35.scalar(0.0, 0.0, 50.0, 0.2)
        rng = random.Random(3508)
        m_true = 2.5
        trials = 800
        for amp in (0.0, 0.7, 5.0):
            hits = 0
            for _ in range(trials):
                e3 = _chol_sample(rng, cov3)
                a3 = (amp, -0.5 * amp, 0.0)
                r3 = o35.Reading(
                    tuple(m_true * a3[i] + e3[i] for i in range(3)),
                    tuple(a3[i] + e3[3 + i] for i in range(3)),
                    cov3,
                    subspace=sub,
                )
                e1 = _chol_sample(rng, cov1)
                a1 = 0.8 * amp
                r1 = o35.Reading((m_true * a1 + e1[0],), (a1 + e1[1],), cov1)
                if o35.State([r3, r1, z_read]).p_value(m_true) > 0.05:
                    hits += 1
            rate = hits / trials
            self.assertGreater(rate, 0.925, f"under-coverage at amp={amp}: {rate}")
            self.assertLess(rate, 0.975, f"over-coverage at amp={amp}: {rate}")

    def test_09_can_scenes_g29_and_door_a_disagreement_cost(self):
        """Can's near (5.5, 5.6) vs far (5, 6), G29 right-angle, and wildly disagreeing (100,10)&(10,100):
        Door A1 keeps the exact Gaussian set width and puts the entire contradiction in `leftover`."""
        near = o35.ours35([o35.scalar(5.5, 1.0, 0.1, 1e-6), o35.scalar(5.6, 1.0, 0.1, 1e-6)])
        far = o35.ours35([o35.scalar(5.0, 1.0, 0.1, 1e-6), o35.scalar(6.0, 1.0, 0.1, 1e-6)])
        self.assertAlmostEqual(near["point"], 5.55, places=6)
        self.assertAlmostEqual(far["point"], 5.50, places=6)
        w_near = near["set"][0][1] - near["set"][0][0]
        w_far = far["set"][0][1] - far["set"][0][0]
        self.assertAlmostEqual(w_near, w_far, places=6)
        self.assertAlmostEqual(near["leftover"], 0.5, places=4)
        self.assertAlmostEqual(far["leftover"], 50.0, places=4)

        # G29 right-angle: flat K == 100, set == (0, inf), screaming leftover == 100
        g29 = o35.ours35([o35.scalar(10.0, 0.0, 1.0, 1.0), o35.scalar(0.0, 10.0, 1.0, 1.0)])
        self.assertTrue(g29["point_is_convention"])
        self.assertEqual(g29["shape"], "all")
        self.assertAlmostEqual(g29["leftover"], 100.0, places=10)

        # Wild disagreement (100, 10) & (10, 100): narrow set around 1.0 + huge leftover 8100
        wild = o35.ours35([o35.scalar(100.0, 10.0, 1.0, 1.0), o35.scalar(10.0, 100.0, 1.0, 1.0)])
        self.assertAlmostEqual(wild["point"], 1.0, places=10)
        self.assertTrue(0.85 < wild["set"][0][0] < 0.92 and 1.08 < wild["set"][0][1] < 1.15)
        self.assertAlmostEqual(wild["leftover"], 8100.0, places=8)

    def test_10_known_door_a1_cost_anti_aligned_boundary_end_vs_aligned_end(self):
        """Documents Theorem E15-3's sharp Door A1 cost: on anti-aligned (10,-3)&(10,-3.1),
        T81's point is boundary_inf (K_min = 18.61), yet T45 rejects inf (p = 1.6e-5) and keeps
        only the hilltop around global maximum m = 0.305 (K = 218.61). On aligned (10, 0),
        boundary_inf IS kept (p = 1.0, shape = 'at_least')."""
        anti = o35.ours35([o35.scalar(10.0, -3.0, 1.0, 1.0), o35.scalar(10.0, -3.1, 1.0, 1.0)])
        self.assertEqual(anti["point_status"], "boundary_inf")
        self.assertTrue(math.isinf(anti["point"]))
        self.assertAlmostEqual(anti["leftover"], 18.61, places=2)
        self.assertLess(anti["state"].p_value(math.inf), 1e-4)
        self.assertEqual(anti["shape"], "interval")
        self.assertFalse(_in_set(anti["set"], math.inf))
        self.assertTrue(_in_set(anti["set"], 0.305))

        aligned = o35.ours35([o35.scalar(10.0, 0.0, 1.0, 1.0)])
        self.assertEqual(aligned["point_status"], "boundary_inf")
        self.assertEqual(aligned["shape"], "at_least")
        self.assertAlmostEqual(aligned["state"].p_value(math.inf), 1.0, places=12)

        aligned_zero = o35.ours35([o35.scalar(0.0, 10.0, 1.0, 1.0)])
        self.assertEqual(aligned_zero["point_status"], "boundary_zero")
        self.assertEqual(aligned_zero["shape"], "at_most")
        self.assertAlmostEqual(aligned_zero["state"].p_value(0.0), 1.0, places=12)

        # Two aligned boundary bodies also lump cleanly to at_least without a fake m=0 pinhole sliver
        lp_al = o35.lump(aligned, aligned)
        self.assertEqual(lp_al["shape"], "at_least")
        self.assertGreater(lp_al["set"][0][0], 5.0)

    def test_11_two_body_lump_common_force_known_mass_and_coverage(self):
        """L5 lump, L14 common_force, exact 1-dof lump_known/common_force_known, exact 0/0
        one-sided reduction, ghost-free asymmetric 100+10, ultra-narrow SNR=1e6 precision,
        and two-body frequentist coverage >= 95%."""
        b1 = o35.ours35([o35.scalar(50.0, 10.0, 1.0, 0.2)])
        b2 = o35.ours35([o35.scalar(100.0, 10.0, 1.0, 0.2)])
        b0 = o35.ours35([o35.scalar(0.0, 0.0, 1.0, 1.0)])

        # L5 lump 5 + 10
        lp = o35.lump(b1, b2)
        self.assertAlmostEqual(lp["point"], 15.0, places=10)
        self.assertAlmostEqual(lp["leftover"], 0.0, places=10)
        self.assertAlmostEqual(lp["curve"](15.0), 0.0, places=10)
        self.assertGreater(lp["curve"](20.0), 10.0)

        # L5 lump 5 + 0/0 collapses cleanly to exact 1-dof lower bound (lo_1, inf)
        lp_z = o35.lump(b1, b0)
        self.assertTrue(lp_z["point_is_convention"])
        self.assertEqual(lp_z["shape"], "at_least")
        self.assertAlmostEqual(lp_z["set"][0][0], b1["set"][0][0], places=8)

        # L14 common_force 5 || 10 and 5 || 0/0
        cf = o35.common_force(b1, b2)
        self.assertAlmostEqual(cf["point"], 10.0 / 3.0, places=10)
        cf_z = o35.common_force(b1, b0)
        self.assertEqual(cf_z["shape"], "at_most")
        self.assertAlmostEqual(cf_z["set"][0][1], b1["set"][0][1], places=8)

        # Exact 1-dof known-mass shifts (no 2-dof tax)
        lk = o35.lump_known(b1, 10.0)
        self.assertAlmostEqual(lk["point"], 15.0, places=12)
        self.assertAlmostEqual(lk["set"][0][0], b1["set"][0][0] + 10.0, places=12)
        self.assertAlmostEqual(lk["set"][0][1], b1["set"][0][1] + 10.0, places=12)
        cfk = o35.common_force_known(b1, 10.0)
        self.assertAlmostEqual(cfk["point"], 10.0 / 3.0, places=12)
        self.assertAlmostEqual(cfk["set"][0][0], 1.0 / (1.0 / b1["set"][0][0] + 0.1), places=12)
        self.assertAlmostEqual(cfk["set"][0][1], 1.0 / (1.0 / b1["set"][0][1] + 0.1), places=12)

        # Asymmetric nonlinear 100 + 10: ghost-free (unlike linearized 1-dof Q_pos)
        a100 = o35.ours35([o35.scalar(100.0, 1.0, 1.0, 0.1)])
        a10 = o35.ours35([o35.scalar(10.0, 1.0, 1.0, 0.1)])
        lp_asym = o35.lump(a100, a10)
        self.assertAlmostEqual(lp_asym["point"], 110.0, places=6)
        self.assertLess(lp_asym["p_value"](15.0), 1e-10)

        # Ultra-narrow SNR=1e6 two-body lumping: analytic 2-dof half-width sqrt(2*ln(20))*sqrt(v1+v2)
        n1 = o35.ours35([o35.scalar(5.0e6, 1.0e6, 1.0, 1.0)])
        n2 = o35.ours35([o35.scalar(10.0e6, 1.0e6, 1.0, 1.0)])
        lp_n = o35.lump(n1, n2)
        # Var(m1) = (1 + 25)/1e12, Var(m2) = (1 + 100)/1e12 => sd = sqrt(127)/1e6
        expected_half = math.sqrt(2.0 * math.log(20.0)) * math.sqrt(127.0) / 1.0e6
        self.assertAlmostEqual(lp_n["set"][0][1] - lp_n["set"][0][0], 2.0 * expected_half, delta=1e-3 * expected_half)

        # Two-body coverage >= 0.95 (conservative 2-dof Q_joint guarantee)
        rng = random.Random(3511)
        m1_t, m2_t = 3.0, 7.0
        M_t = m1_t + m2_t
        hits = 0
        trials = 400
        for _ in range(trials):
            s1 = o35.State([o35.scalar(m1_t * 1.5 + rng.gauss(0, 1), 1.5 + rng.gauss(0, 1), 1.0, 1.0)])
            s2 = o35.State([o35.scalar(m2_t * 2.0 + rng.gauss(0, 1), 2.0 + rng.gauss(0, 1), 1.0, 1.0)])
            # Direct evaluation of exact Q_joint p-value at true split bounds sup_{m1+m2=M_t} from below
            if o35._lump_p_split(s1, s2, m1_t, m2_t) > 0.05:
                hits += 1
        self.assertGreater(hits / trials, 0.925)

    def test_12_invalid_inputs_rejected_l2(self):
        """Non-SPD cov, asymmetric cov, bad subspace, NaN, and invalid level raise ValueError."""
        with self.assertRaises(ValueError):
            o35.Reading((1.0,), (1.0,), ((1.0, 1.0), (1.0, 1.0)))  # singular cov (L2)
        with self.assertRaises(ValueError):
            o35.Reading((1.0,), (1.0,), ((1.0, 0.5), (0.0, 1.0)))  # asymmetric cov
        with self.assertRaises(ValueError):
            o35.Reading((1.0, 2.0), (1.0, 2.0), ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
                        subspace=((1.0,), (1.0,)))  # non-unit subspace column
        with self.assertRaises(ValueError):
            o35.Reading((float("nan"),), (1.0,), ((1.0, 0.0), (0.0, 1.0)))
        with self.assertRaises(ValueError):
            o35.ours35([o35.scalar(1.0, 1.0, 1.0, 1.0)], level=1.0)
        with self.assertRaises(ValueError):
            o35.State([], scale=0.0)


if __name__ == "__main__":
    unittest.main()

"""Useful limits and uncertainty expressions, checked against independent oracles."""
import math
import unittest
import numpy as np
from estimator import estimate, isotropic_covariance
from practical_formulas import (known_direction_mass, aligned_local_uncertainty,
                                log_mass_moments, summarize)


class PracticalTests(unittest.TestCase):
    def test_alignment_alone_is_not_the_zero_noise_limit(self):
        e = estimate([2, 0], [1, 0], np.eye(4), ratio_order=2048)
        self.assertAlmostEqual(e.mass, 1.520356039, places=7)
        self.assertGreater(abs(e.mass-2), .4)
        near = estimate([2*math.cos(.001), 2*math.sin(.001)], [1, 0], np.eye(4))
        self.assertLess(abs(near.mass-e.mass), 1e-6)

    def test_full_vector_aligned_small_noise_convergence(self):
        errors = []
        for sigma in (1, .5, .25):
            e = estimate([2, 0], [1, 0], isotropic_covariance(2, sigma, sigma),
                         direction_order=256, ratio_order=2048)
            errors.append(abs(e.mass-2))
        self.assertTrue(errors[0] > errors[1] > errors[2])
        self.assertLess(errors[-1], .001)

    def test_known_direction_exact_formula_and_limit(self):
        for f, a, sf, sa in [(2, 1, 1, 1), (-2, 1, 1, 2), (0, 0, 3, 2)]:
            e = estimate([f], [a], isotropic_covariance(1, sf, sa),
                         known_direction=[1], ratio_order=2048)
            self.assertAlmostEqual(known_direction_mass(f, a, sf, sa), e.mass, places=7)
        # Explicit exponential relative-error bound for both positive z >= k.
        for sigma in (1, .5, .25):
            k = 1/sigma
            bound = 2*math.exp(-k*k/2)/(math.sqrt(2*math.pi)*k)
            self.assertLess(abs(known_direction_mass(2, 1, sigma, sigma)/2-1), bound)

    def test_equal_standardized_aligned_signals_exact_ratio(self):
        for strength in (.1, 1, 4):
            e = estimate([3*strength, 0], [2*strength, 0], isotropic_covariance(2, 3, 2))
            self.assertAlmostEqual(e.mass, 1.5, places=12)

    def test_null_log_spread_and_interval_analytic(self):
        e = estimate([0, 0], [0, 0], isotropic_covariance(2, 3, 2))
        mean, spread = log_mass_moments(e)
        self.assertAlmostEqual(mean, 0, places=11)
        self.assertAlmostEqual(spread, math.pi/2, places=9)
        r = summarize(e)
        self.assertAlmostEqual(r.lower, 1.5*math.tan(.025*math.pi/2), places=10)
        self.assertAlmostEqual(r.upper, 1.5*math.tan(.975*math.pi/2), places=9)

    def test_reciprocal_uncertainty_and_units(self):
        f, a = [2, 1], [1, -.5]
        e = estimate(f, a, isotropic_covariance(2, .8, 1.2), ratio_order=2048)
        inv = estimate(a, f, isotropic_covariance(2, 1.2, .8), ratio_order=2048)
        scaled = estimate(np.array(f)*1000, a, isotropic_covariance(2, 800, 1.2), ratio_order=2048)
        r, ri, rs = summarize(e), summarize(inv), summarize(scaled)
        self.assertAlmostEqual(ri.lower, 1/r.upper, places=9)
        self.assertAlmostEqual(ri.upper, 1/r.lower, places=9)
        self.assertAlmostEqual(r.log_standard_deviation, ri.log_standard_deviation, places=10)
        self.assertAlmostEqual(r.log_standard_deviation, rs.log_standard_deviation, places=10)
        self.assertAlmostEqual(rs.mass/r.mass, 1000, places=8)

    def test_familiar_local_uncertainty_including_correlation(self):
        f, a, sf, sa = 20., 10., .6, .2
        for rho in (-.7, 0, .7):
            cross = rho*sf*sa
            cov = np.block([[sf**2*np.eye(2), cross*np.eye(2)],
                            [cross*np.eye(2), sa**2*np.eye(2)]])
            got = aligned_local_uncertainty([f, 0], [a, 0], cov)
            oracle = (f/a)*math.sqrt((sf/f)**2 + (sa/a)**2 - 2*cross/(f*a))
            self.assertAlmostEqual(got, oracle, places=12)

    def test_general_covariance_local_oracle_and_rotation(self):
        # At F=2,a=1 the tangent constraint is dFy-2*dAy=0. Gaussian
        # conditioning gives Sigma - Sigma*n*n'*Sigma/(n'*Sigma*n).
        b = np.array([[2, 0, 0, 0], [.6, 1, 0, 0], [.4, -.2, 1, 0], [.3, .5, .2, 1]])
        cov = b @ b.T * .001
        normal, gradient = np.array([0, 1, 0, -2.]), np.array([1, 0, -2, 0.])
        v = cov @ normal
        conditional = cov - np.outer(v, v)/(normal @ v)
        oracle = math.sqrt(gradient @ conditional @ gradient)
        got = aligned_local_uncertainty([2, 0], [1, 0], cov)
        self.assertAlmostEqual(got, oracle, places=12)
        rotation = np.array([[.6, -.8], [.8, .6]])
        transform = np.kron(np.eye(2), rotation)
        rotated = aligned_local_uncertainty(rotation @ [2, 0], rotation @ [1, 0], transform @ cov @ transform.T)
        self.assertAlmostEqual(rotated, oracle, places=12)

    def test_high_signal_log_spread_matches_local_approximation(self):
        e = estimate([20, 0], [10, 0], isotropic_covariance(2, .6, .2),
                     direction_order=512, ratio_order=4096)
        refined = estimate([20, 0], [10, 0], isotropic_covariance(2, .6, .2),
                           direction_order=768, ratio_order=8192)
        spread, fine = summarize(e).log_standard_deviation, summarize(refined).log_standard_deviation
        local = math.sqrt((.6/20)**2 + (.2/10)**2)
        self.assertLess(abs(fine/spread-1), .0001)
        self.assertLess(abs(fine/local-1), .003)

    def test_bad_local_formula_inputs_and_interval_probability(self):
        for f, a in [([0, 0], [1, 0]), ([1, 0], [0, 0]), ([-2, 0], [1, 0]), ([2, 1], [1, 0])]:
            with self.assertRaises(ValueError):
                aligned_local_uncertainty(f, a, np.eye(4))
        e = estimate([0], [0], np.eye(2))
        for probability in (0, 1, float('nan')):
            with self.assertRaises(ValueError):
                summarize(e, probability=probability)


if __name__ == '__main__':
    unittest.main()

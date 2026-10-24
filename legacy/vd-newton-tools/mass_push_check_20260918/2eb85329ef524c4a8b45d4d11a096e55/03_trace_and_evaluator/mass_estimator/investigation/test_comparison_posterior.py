"""Independent oracles and refinement for continuous posterior scoring."""

import math
import unittest

import numpy as np

from comparison_posterior import infer_batch
from crossover_vector_checks import infer


class ComparisonPosteriorTests(unittest.TestCase):
    def test_null_posterior_analytic(self):
        for truth in (.03, 1., 17.):
            result = infer_batch(np.zeros(3), np.zeros(3), truth)
            for point in result['points'].values():
                np.testing.assert_allclose(point, 1, atol=2e-8)
            y = math.log(truth)
            for summary in result['distributions'].values():
                self.assertAlmostEqual(summary['cdf_truth'][0], 2*math.atan(truth)/math.pi, places=9)
                self.assertAlmostEqual(summary['log_density_score'][0], math.log(math.pi*math.cosh(y)), places=9)
                self.assertAlmostEqual(summary['mean_log_mass'][0], 0, places=10)
                self.assertAlmostEqual(summary['log_sd'][0], math.pi/2, places=9)
                for level in (50, 80, 95):
                    tail = (1-level/100)/2
                    self.assertAlmostEqual(summary[f'log_lower_{level}'][0], math.log(math.tan(math.pi*tail/2)), places=7)
                    self.assertAlmostEqual(summary[f'log_upper_{level}'][0], -math.log(math.tan(math.pi*tail/2)), places=7)
                # Independent integral of exact CDF; not the implementation's
                # density-weighted pair-expectation formula.
                x, w = np.polynomial.legendre.leggauss(240)
                left = -32+(y+32)*(x+1)/2
                right = y+(32-y)*(x+1)/2
                fl = 2*np.arctan(np.exp(left))/math.pi
                survival = 2*np.arctan(np.exp(-right))/math.pi
                expected = (y+32)/2*np.dot(w, fl**2)+(32-y)/2*np.dot(w, survival**2)
                self.assertAlmostEqual(summary['log_crps'][0], expected, places=9)

    def test_existing_angular_integration(self):
        force = np.array([[0, 0, 0], [8, 0, 0], [3, 0, 0], [2, 0, 0], [3, 0, 0], [16, 1, 0]], float)
        acceleration = np.array([[0, 0, 0], [0, 0, 0], [1, .5, 0], [0, 1, 0], [-2, 0, 0], [.5, -1, 1]], float)
        result = infer_batch(force, acceleration, truth=3.7)
        reference, cdf = infer(force, acceleration, ratio=3.7, order=256)
        for column, name in enumerate(('flat_joint', 'tube_joint')):
            np.testing.assert_allclose(result['points'][name], reference[:, column], rtol=1e-7)
            np.testing.assert_allclose(result['distributions'][name]['cdf_truth'], cdf[:, column], atol=1e-7)

    def test_reciprocal_symmetry_and_truth_independence(self):
        rng = np.random.default_rng(90411)
        force, acceleration = rng.normal(size=(2, 12, 3))*4
        forward = infer_batch(force, acceleration, truth=2.7)
        reverse = infer_batch(acceleration, force, truth=1/2.7)
        alternative = infer_batch(force, acceleration, truth=.1)
        for name, value in forward['points'].items():
            np.testing.assert_allclose(value*reverse['points'][name], 1, atol=2e-8)
            np.testing.assert_array_equal(value, alternative['points'][name])
        for name, summary in forward['distributions'].items():
            other = reverse['distributions'][name]
            np.testing.assert_allclose(summary['cdf_truth']+other['cdf_truth'], 1, atol=1e-9)
            for key in ('log_crps', 'log_density_score', 'log_sd'):
                np.testing.assert_allclose(summary[key], other[key], atol=1e-9)
            for level in (50, 80, 95):
                np.testing.assert_allclose(summary[f'log_lower_{level}'], -other[f'log_upper_{level}'], atol=2e-8)

    def test_nonnull_crps_and_quantiles_from_angular_cdf(self):
        # Evaluate the old, independent theta integral at many log masses.
        # Integrating CDF squared directly checks the continuous CRPS formula
        # without sharing its density, spectral antiderivative or pair moment.
        f, a, truth = [4, 1, 0], [1, -.3, .5], 2.3
        result = infer_batch(f, a, truth)
        y = math.log(truth)
        nodes, weights = np.polynomial.legendre.leggauss(180)
        expected = np.zeros(2)
        for lo, hi, is_right in ((-32., y, False), (y, 32., True)):
            z = (lo+hi)/2+(hi-lo)*nodes/2
            cdfs = np.array([infer(f, a, ratio=math.exp(v), order=160)[1][0] for v in z])
            integrand = (1-cdfs)**2 if is_right else cdfs**2
            expected += (hi-lo)/2 * (weights @ integrand)
        for column, name in enumerate(('flat_joint', 'tube_joint')):
            summary = result['distributions'][name]
            self.assertAlmostEqual(summary['log_crps'][0], expected[column], places=7)
            for level in (50, 80, 95):
                for side, target in (('lower', (1-level/100)/2), ('upper', (1+level/100)/2)):
                    ratio = math.exp(summary[f'log_{side}_{level}'][0])
                    _, cdf = infer(f, a, ratio=ratio, order=256)
                    self.assertAlmostEqual(cdf[0, column], target, places=7)

    def test_refinement_observed_range(self):
        rng = np.random.default_rng(2026091703)
        force, acceleration = rng.normal(size=(2, 160, 3))
        force[:, 0] += np.tile([.25, 1, 4, 8, 16], 32)
        acceleration[:, 0] += np.tile([.1, 1, 8, -8], 40)
        # Adverse signs with equal strong magnitudes have two tail modes.
        force = np.vstack((force, [21, 0, 0], [21, 0, 0], [21, 0, 0]))
        acceleration = np.vstack((acceleration, [-21, 0, 0], [21, 0, 0], [0, 0, 0]))
        coarse = infer_batch(force, acceleration, truth=8., order=96)
        fine = infer_batch(force, acceleration, truth=8., order=192)
        for name, value in coarse['points'].items():
            np.testing.assert_allclose(value, fine['points'][name], rtol=1e-6, atol=1e-8)
        for name, summary in coarse['distributions'].items():
            for key, value in summary.items():
                np.testing.assert_allclose(value, fine['distributions'][name][key], atol=1e-4, rtol=1e-6,
                                           err_msg=f'{name}: {key}')
            self.assertTrue(np.all(summary['log_crps'] >= 0))

    def test_validation(self):
        for force, acceleration, truth in (([1, 2], [1, 2], 1), ([1, 2, 3], [1, 2, 3], 0),
                                           ([math.nan, 0, 0], [1, 2, 3], 1)):
            with self.assertRaises(ValueError):
                infer_batch(force, acceleration, truth)


if __name__ == '__main__':
    unittest.main()

"""Independent radial oracles, full-covariance checks and diagonal invariances."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np

from compare_estimators import refinement_check
from comparison_diagonal_posterior import _infer_chunk, gaussian_radius_moments, infer_batch
from comparison_posterior import infer_batch as isotropic_infer

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from estimator import estimate as full_covariance_estimate


class GaussianRadiusTests(unittest.TestCase):
    def test_noncentral_isotropic_analytic_oracles(self):
        radii = np.r_[0., np.geomspace(1e-9, 100, 80)]
        for sigma in (.15, 1., 7.):
            mean = np.column_stack([radii, 0*radii, 0*radii]) * sigma
            variance = np.full(3, sigma*sigma)
            inverse, first, second = gaussian_radius_moments(mean, variance, 48)
            erf = np.array([math.erf(float(radius)/math.sqrt(2)) for radius in radii])
            reference_inverse = np.divide(erf, radii, out=np.full_like(radii, math.sqrt(2/math.pi)), where=radii != 0) / sigma
            reference_first = sigma*((radii*radii+1)*reference_inverse*sigma + math.sqrt(2/math.pi)*np.exp(-radii*radii/2))
            np.testing.assert_allclose(inverse, reference_inverse, rtol=3e-13)
            np.testing.assert_allclose(first, reference_first, rtol=3e-13)
            np.testing.assert_allclose(second, sigma*sigma*(radii*radii+3), rtol=3e-15)

    def test_anisotropic_radius_refinement_at_ratio_sixteen(self):
        rng = np.random.default_rng(2026091901)
        directions = np.vstack([np.eye(3), np.ones(3)/math.sqrt(3), rng.normal(size=(20,3))])
        directions /= np.linalg.norm(directions, axis=1)[:,None]
        means = np.concatenate([directions*strength for strength in (0., .1, 1., 3., 8., 32., 64.)])
        for sd in ([.25,.25,4.], [.25,1.,4.], [.25,4.,4.]):
            sd = np.asarray(sd)/np.sqrt(np.mean(np.asarray(sd)**2))
            ordinary = gaussian_radius_moments(means, sd*sd, 48)
            refined = gaussian_radius_moments(means, sd*sd, 128)
            for first, second in zip(ordinary, refined):
                np.testing.assert_allclose(first, second, rtol=2e-9, atol=1e-12)

    def test_radius_rejects_nonfinite_gaussians(self):
        for mean, variance in (([math.nan,0.,0.], [1.,1.,1.]),
                               ([0.,0.,0.], [1.,math.inf,1.]),
                               ([0.,0.,0.], [1.,0.,1.])):
            with self.assertRaises(ValueError):
                gaussian_radius_moments(mean, variance)


class DiagonalPosteriorTests(unittest.TestCase):
    def assert_summaries_equal(self, first, second, atol=3e-8):
        for method, point in first['points'].items():
            np.testing.assert_allclose(point, second['points'][method], rtol=2e-8, atol=atol,
                                       err_msg=method)
        for method, distribution in first['distributions'].items():
            for key, value in distribution.items():
                np.testing.assert_allclose(value, second['distributions'][method][key], rtol=2e-8, atol=atol,
                                           err_msg=f'{method}: {key}')

    def test_raw_diagonal_path_matches_isotropic_oracle(self):
        # Call _infer_chunk explicitly: public isotropic specialization must
        # not turn this into a comparison of the isotropic routine with itself.
        force = np.array([[0,0,0], [8,0,0], [2,-1,3], [3,0,0], [16,1,0], [8,0,0], [32,0,0]], float)
        acceleration = np.array([[0,0,0], [0,0,0], [.1,.7,-.4], [-2,0,0], [.5,-1,1], [-8,0,0], [-32,0,0]], float)
        raw = _infer_chunk(force, acceleration, np.ones(3), np.ones(3), 2.3, 96)
        # Refine the independent oracle too: its degree-16 interval endpoint
        # interpolation can differ by ~6e-8 for the zero-acceleration case.
        oracle = isotropic_infer(force, acceleration, 2.3, 192)
        self.assert_summaries_equal(raw, oracle)

    def test_zero_data_with_anisotropy_has_analytic_law(self):
        # Equal covariance shape in both RMS-scaled channels leaves the
        # zero-data theta law uniform, although every coordinate SD differs.
        result = infer_batch(np.zeros(3), np.zeros(3), truth=2.,
                             force_sd=[.75,1.5,3.], acceleration_sd=[.25,.5,1.])
        shift = math.log(3)
        for point in result['points'].values():
            np.testing.assert_allclose(point, 3., rtol=1e-8)
        for distribution in result['distributions'].values():
            self.assertAlmostEqual(distribution['cdf_truth'][0], 2*math.atan(2/3)/math.pi, places=9)
            self.assertAlmostEqual(distribution['mean_log_mass'][0], shift, places=10)
            self.assertAlmostEqual(distribution['log_sd'][0], math.pi/2, places=9)
            self.assertAlmostEqual(distribution['log_density_score'][0], math.log(math.pi*math.cosh(math.log(2/3))), places=9)
            for level in (50,80,95):
                tail = (1-level/100)/2
                expected = math.log(math.tan(math.pi*tail/2))
                self.assertAlmostEqual(distribution[f'log_lower_{level}'][0], shift+expected, places=7)
                self.assertAlmostEqual(distribution[f'log_upper_{level}'][0], shift-expected, places=7)

    def test_strong_adverse_symmetric_readings_keep_exact_median(self):
        # A nearly flat CDF between two well-separated tail modes makes
        # numerical inversion unstable. Channel symmetry fixes the median.
        for sign in (-1., 1.):
            result = infer_batch([32.,0.,0.], [sign*32.,0.,0.], truth=1.,
                                 force_sd=[.5,1.,2.], acceleration_sd=[.5,1.,2.])
            for point in result['points'].values():
                np.testing.assert_allclose(point, 1., rtol=1e-8)
            for summary in result['distributions'].values():
                np.testing.assert_allclose(summary['log_median'], 0., atol=1e-10)

    def test_full_covariance_reference_converges_to_flat_point(self):
        cases = [([2.,-1.,.5], [.8,.1,.2], [.7,1.,1.4], [1.2,.9,.8]),
                 ([0.,0.,0.], [0.,0.,0.], [.5,1.,2.], [2.,1.,.5])]
        for force, acceleration, sf, sa in cases:
            covariance = np.diag(np.r_[np.asarray(sf)**2, np.asarray(sa)**2])
            references = [full_covariance_estimate(force, acceleration, covariance,
                                                   direction_order=order, ratio_order=256).mass
                          for order in (24,48)]
            self.assertLess(abs(references[0]/references[1]-1), 2e-6)
            result = infer_batch(force, acceleration, force_sd=sf, acceleration_sd=sa)
            np.testing.assert_allclose(result['points']['flat_joint'], references[1], rtol=1e-8)

    def test_coordinate_permutation_and_common_sign_flips(self):
        force = np.array([[3.,-1.,.5], [.1,2.,-1.]])
        acceleration = np.array([[1.,.2,.6], [3.,-1.,.2]])
        sf, sa = np.array([.5,1.,2.]), np.array([2.,1.,.5])
        ordinary = infer_batch(force, acceleration, truth=2.3, force_sd=sf, acceleration_sd=sa)
        permutation, signs = [2,0,1], np.array([-1.,1.,-1.])
        transformed = infer_batch(force[:,permutation]*signs, acceleration[:,permutation]*signs,
                                  truth=2.3, force_sd=sf[permutation], acceleration_sd=sa[permutation])
        self.assert_summaries_equal(ordinary, transformed)

    def test_channel_swap_and_truth_independence(self):
        force = np.array([[3.,-1.,.5], [.1,2.,-1.]])
        acceleration = np.array([[1.,.2,.6], [3.,-1.,.2]])
        sf, sa = [.5,1.,2.], [2.,1.,.5]
        forward = infer_batch(force, acceleration, truth=2.7, force_sd=sf, acceleration_sd=sa)
        backward = infer_batch(acceleration, force, truth=1/2.7, force_sd=sa, acceleration_sd=sf)
        changed_truth = infer_batch(force, acceleration, truth=.03, force_sd=sf, acceleration_sd=sa)
        for method, point in forward['points'].items():
            np.testing.assert_allclose(point*backward['points'][method], 1., rtol=2e-8)
            np.testing.assert_array_equal(point, changed_truth['points'][method])
        for method, summary in forward['distributions'].items():
            reverse = backward['distributions'][method]
            np.testing.assert_allclose(summary['mean_log_mass'], -reverse['mean_log_mass'], atol=1e-9)
            np.testing.assert_allclose(summary['cdf_truth']+reverse['cdf_truth'], 1., atol=1e-9)
            for key in ('log_crps', 'log_density_score', 'log_sd'):
                np.testing.assert_allclose(summary[key], reverse[key], atol=1e-9)
            for level in (50,80,95):
                np.testing.assert_allclose(summary[f'log_lower_{level}'], -reverse[f'log_upper_{level}'], atol=2e-8)
            for key in summary:
                if key not in ('cdf_truth', 'log_crps', 'log_density_score'):
                    np.testing.assert_array_equal(summary[key], changed_truth['distributions'][method][key])

    def test_channel_units_rescale_mass_and_log_summaries(self):
        force, acceleration = np.array([3.,-1.,.5]), np.array([1.,.2,.6])
        sf, sa = np.array([.5,1.,2.]), np.array([2.,1.,.5])
        base = infer_batch(force, acceleration, truth=2.7, force_sd=sf, acceleration_sd=sa)
        scaled = infer_batch(7*force, .4*acceleration, truth=2.7*7/.4, force_sd=7*sf, acceleration_sd=.4*sa)
        mass_scale, shift = 7/.4, math.log(7/.4)
        for method, point in base['points'].items():
            np.testing.assert_allclose(scaled['points'][method], point*mass_scale, rtol=2e-8)
        for method, summary in base['distributions'].items():
            for key, values in summary.items():
                expected = values+shift if key in ('mean_log_mass','log_median') or key.startswith(('log_lower_','log_upper_')) else values
                np.testing.assert_allclose(scaled['distributions'][method][key], expected, atol=2e-8)

    def test_isotropic_channels_with_different_noise_scales(self):
        force, acceleration = np.array([4.,-1.,.5]), np.array([1.,.2,.6])
        reference = isotropic_infer(force, acceleration, truth=2.3)
        result = infer_batch(force*2., acceleration*.5, truth=2.3*4.,
                             force_sd=[2.]*3, acceleration_sd=[.5]*3)
        for method, point in reference['points'].items():
            np.testing.assert_allclose(result['points'][method], 4*point, rtol=1e-12)
        for method, summary in reference['distributions'].items():
            for key, values in summary.items():
                expected = values+math.log(4) if key in ('mean_log_mass','log_median') or key.startswith(('log_lower_','log_upper_')) else values
                np.testing.assert_allclose(result['distributions'][method][key], expected, atol=1e-12)

    def test_refinement_for_weak_strong_and_axis_dependent_readings(self):
        rng = np.random.default_rng(20260919)
        directions = np.vstack([np.eye(3), np.ones(3)/math.sqrt(3)])
        for sf, sa in (([.5,1,2],[2,1,.5]), ([.25,.25,4],[4,.25,.25]), ([.25,1,4],[1,4,.25])):
            sf, sa = np.asarray(sf), np.asarray(sa)
            force = np.concatenate([directions*strength for strength in (.25,2.,8.,16.,32.)])
            acceleration = np.concatenate([directions*8. for _ in range(5)])
            force += rng.normal(size=force.shape)*sf
            acceleration += rng.normal(size=acceleration.shape)*sa
            coarse = infer_batch(force, acceleration, force_sd=sf, acceleration_sd=sa, order=96)
            fine = infer_batch(force, acceleration, force_sd=sf, acceleration_sd=sa, order=192)
            check = refinement_check(coarse, fine)
            self.assertTrue(check['passed'], check)
            for summary in coarse['distributions'].values():
                self.assertTrue(np.all(summary['log_crps'] >= 0))


if __name__ == '__main__':
    unittest.main()

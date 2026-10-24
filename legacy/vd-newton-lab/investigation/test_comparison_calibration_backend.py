"""Per-trial noise calibration preserves each independent inference problem."""

import math
import unittest

import numpy as np

from compare_estimators import refinement_check
from comparison_diagonal_baselines import baseline_points, covariance_profile
from comparison_diagonal_posterior import infer_batch
from comparison_posterior import infer_batch as isotropic_infer


class CalibrationBackendTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(2026091911)
        # More than one internal posterior chunk, with varying channel ratios.
        self.force = rng.normal(size=(19, 3)) + [4., 1., .2]
        self.acceleration = rng.normal(size=(19, 3)) + [1., .2, .05]
        self.sf = np.exp(rng.uniform(-.7, .7, size=(19, 3)))
        self.sa = np.exp(rng.uniform(-.7, .7, size=(19, 3)))

    def assert_posteriors_close(self, first, second, atol=3e-8):
        for method, value in first['points'].items():
            np.testing.assert_allclose(value, second['points'][method], rtol=2e-8,
                                       atol=atol, err_msg=method)
        for method, distribution in first['distributions'].items():
            for key, value in distribution.items():
                np.testing.assert_allclose(value, second['distributions'][method][key],
                                           rtol=2e-8, atol=atol, err_msg=f'{method}: {key}')

    @staticmethod
    def combine(results):
        return {'points': {name: np.concatenate([r['points'][name] for r in results])
                           for name in results[0]['points']},
                'distributions': {name: {key: np.concatenate([r['distributions'][name][key] for r in results])
                                          for key in distribution}
                                  for name, distribution in results[0]['distributions'].items()}}

    def test_trial_specific_noise_matches_separate_inferences(self):
        combined = infer_batch(self.force, self.acceleration, 2.7,
                               force_sd=self.sf, acceleration_sd=self.sa)
        separate = self.combine([infer_batch(f, a, 2.7, force_sd=sf, acceleration_sd=sa)
                                 for f, a, sf, sa in zip(self.force, self.acceleration, self.sf, self.sa)])
        self.assert_posteriors_close(combined, separate, atol=1e-10)
        actual_baselines = baseline_points(self.force, self.acceleration, self.sf, self.sa)
        reference_baselines = [baseline_points(f, a, sf, sa)
                               for f, a, sf, sa in zip(self.force, self.acceleration, self.sf, self.sa)]
        for method, actual in actual_baselines.items():
            reference = np.concatenate([row[method] for row in reference_baselines])
            np.testing.assert_allclose(actual, reference, rtol=2e-13, equal_nan=True, err_msg=method)

    def test_fixed_channel_broadcast_and_identical_rows(self):
        for sf, sa in ((self.sf, [.5, 1., 2.]), ([.5, 1., 2.], self.sa)):
            direct = infer_batch(self.force, self.acceleration, 1.3, force_sd=sf, acceleration_sd=sa)
            repeated = infer_batch(self.force, self.acceleration, 1.3,
                                   force_sd=np.broadcast_to(sf, self.force.shape),
                                   acceleration_sd=np.broadcast_to(sa, self.force.shape))
            self.assert_posteriors_close(direct, repeated, atol=1e-12)
        direct = infer_batch(self.force, self.acceleration, 1.3,
                             force_sd=[.5, 1., 2.], acceleration_sd=[2., .5, 1.])
        repeated = infer_batch(self.force, self.acceleration, 1.3,
                               force_sd=np.tile([.5, 1., 2.], (19, 1)),
                               acceleration_sd=np.tile([2., .5, 1.], (19, 1)))
        self.assert_posteriors_close(direct, repeated, atol=1e-10)

    def test_variable_isotropic_scales_match_individual_oracle(self):
        sf = np.repeat(self.sf[:, :1], 3, axis=1)
        sa = np.repeat(self.sa[:, :1], 3, axis=1)
        combined = infer_batch(self.force, self.acceleration, 2.7, force_sd=sf, acceleration_sd=sa)
        separate = self.combine([infer_batch(f, a, 2.7, force_sd=f_sd, acceleration_sd=a_sd)
                                 for f, a, f_sd, a_sd in zip(self.force, self.acceleration, sf, sa)])
        self.assert_posteriors_close(combined, separate, atol=1e-12)

    def test_variable_noise_channel_units_and_truth_independence(self):
        base = infer_batch(self.force, self.acceleration, 2.7, force_sd=self.sf, acceleration_sd=self.sa)
        changed = infer_batch(self.force, self.acceleration, .03, force_sd=self.sf, acceleration_sd=self.sa)
        scaled = infer_batch(self.force*7, self.acceleration*.4, 2.7*7/.4,
                             force_sd=self.sf*7, acceleration_sd=self.sa*.4)
        for method, values in base['points'].items():
            np.testing.assert_array_equal(changed['points'][method], values)
            np.testing.assert_allclose(scaled['points'][method], values*7/.4, rtol=2e-8)
        for method, distribution in base['distributions'].items():
            for key, values in distribution.items():
                shift = math.log(7/.4) if key in ('mean_log_mass', 'log_median') or key.startswith(('log_lower_', 'log_upper_')) else 0
                np.testing.assert_allclose(scaled['distributions'][method][key], values+shift, atol=3e-8)
                if key not in ('cdf_truth', 'log_crps', 'log_density_score'):
                    np.testing.assert_array_equal(changed['distributions'][method][key], values)

    def test_isotropic_vector_truth_matches_scalar_scoring_including_tails(self):
        truth = np.geomspace(math.exp(-40), math.exp(40), len(self.force))
        combined = isotropic_infer(self.force, self.acceleration, truth)
        separate = self.combine([isotropic_infer(f, a, t)
                                 for f, a, t in zip(self.force, self.acceleration, truth)])
        self.assert_posteriors_close(combined, separate, atol=1e-11)
        changed = isotropic_infer(self.force, self.acceleration, 1.)
        for method, values in combined['points'].items():
            np.testing.assert_array_equal(values, changed['points'][method])

    def test_invalid_per_trial_scales_and_truth_shapes(self):
        for invalid in (np.ones((1, 3)), np.ones((19, 1)), np.ones((18, 3)),
                        np.full((19, 3), np.nan), np.full((19, 3), np.inf),
                        np.zeros((19, 3)), -np.ones((19, 3))):
            for sf, sa in ((invalid, self.sa), (self.sf, invalid)):
                with self.assertRaises(ValueError):
                    infer_batch(self.force, self.acceleration, force_sd=sf, acceleration_sd=sa)
                with self.assertRaises(ValueError):
                    baseline_points(self.force, self.acceleration, sf, sa)
                with self.assertRaises(ValueError):
                    covariance_profile(self.force, self.acceleration, sf, sa)
        for truth in ([1.], np.ones((19, 1)), np.full(19, np.nan), np.zeros(19)):
            with self.assertRaises(ValueError):
                isotropic_infer(self.force, self.acceleration, truth)

    def test_fifty_replicates_high_signal_and_weak_acceleration_refine(self):
        rng = np.random.default_rng(2026091912)
        for sd in (np.ones(3), np.array([.5, 2., 2.])/math.sqrt(2.75)):
            sd = sd/math.sqrt(50)
            for acceleration_strength in (.1, 8.):
                force = np.array([8., 0., 0.]) + rng.normal(size=(12, 3))*sd
                acceleration = np.array([acceleration_strength, 0., 0.]) + rng.normal(size=(12, 3))*sd
                coarse = infer_batch(force, acceleration, 8/acceleration_strength,
                                     force_sd=sd, acceleration_sd=sd, order=96)
                fine = infer_batch(force, acceleration, 8/acceleration_strength,
                                   force_sd=sd, acceleration_sd=sd, order=192)
                checked = refinement_check(coarse, fine)
                self.assertTrue(checked['passed'], checked)

    def test_per_trial_contrasts_up_to_thirty_two_refine(self):
        rng = np.random.default_rng(2026091913)
        shapes = np.array([[.125, 1., 4.], [.125, .125, 4.], [4., .125, .125]])
        shapes /= np.sqrt(np.mean(shapes**2, axis=1))[:, None]
        sf = np.repeat(shapes, 9, axis=0)
        sa = sf[:, ::-1].copy()
        directions = np.tile(np.eye(3), (9, 1))
        strength = np.tile(np.repeat([.1, 8., 8*math.sqrt(50)], 3), 3)
        force = strength[:, None]*directions + rng.normal(size=sf.shape)*sf
        acceleration = 8*directions + rng.normal(size=sa.shape)*sa
        coarse = infer_batch(force, acceleration, force_sd=sf, acceleration_sd=sa, order=96)
        fine = infer_batch(force, acceleration, force_sd=sf, acceleration_sd=sa, order=192)
        checked = refinement_check(coarse, fine)
        self.assertTrue(checked['passed'], checked)

    def test_planned_calibration_extremes_refine_without_clipping(self):
        # These calibration draws are independent of the measurement errors.
        # Retain the extreme tail draw rather than selecting a gentler seed.
        scale = np.array([.5, 2., 2.])/math.sqrt(2.75)
        factors = np.sqrt(np.random.default_rng(2026091923).chisquare(9, size=(1024, 6))/9)
        sf, sa = factors[:, :3]*scale, factors[:, 3:]*scale
        mass_scale = np.sqrt(np.sum(sf*sf, axis=1)/np.sum(sa*sa, axis=1))
        indices = np.unique([np.argmax(np.max(sf, axis=1)/np.min(sf, axis=1)),
                             np.argmax(np.max(sa, axis=1)/np.min(sa, axis=1)),
                             np.argmin(mass_scale), np.argmax(mass_scale)])
        self.assertGreater(max(np.max(sf/np.min(sf, axis=1)[:, None]),
                               np.max(sa/np.min(sa, axis=1)[:, None])), 16.)
        errors = np.random.default_rng(2026091921).standard_normal((1024, 6))[indices]
        for repeats in (1, 50):
            effective_scale = scale/math.sqrt(repeats)
            for a in (.1, .25, 2., 8.):
                force = np.array([8., 0., 0.]) + errors[:, :3]*effective_scale
                acceleration = np.array([a, 0., 0.]) + errors[:, 3:]*effective_scale
                noise = {'force_sd': sf[indices]/math.sqrt(repeats),
                         'acceleration_sd': sa[indices]/math.sqrt(repeats)}
                coarse = infer_batch(force, acceleration, 8/a, **noise, order=96)
                fine = infer_batch(force, acceleration, 8/a, **noise, order=192)
                checked = refinement_check(coarse, fine)
                self.assertTrue(checked['passed'], checked)


if __name__ == '__main__':
    unittest.main()

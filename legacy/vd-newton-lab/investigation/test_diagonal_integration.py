"""Trial generation, physical units and scoring across the diagonal workflow."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import unittest
import uuid

import numpy as np

from experiment_config import run_config
from results_catalog import load_run


WORKSPACE = Path(__file__).resolve().parents[3]


class DiagonalIntegrationTests(unittest.TestCase):
    def run_experiment(self, config):
        workspace = WORKSPACE / '.tools' / 'diagonal_integration_tests' / uuid.uuid4().hex
        output = workspace / '.tools' / 'mass_estimator_runs' / 'test-run'
        with redirect_stdout(io.StringIO()):
            report = run_config(config, output)
        return report, output, workspace

    def test_known_noise_generates_correct_vectors_and_physical_mass(self):
        config = {'samples': 16, 'save_trials': True, 'noise_model': 'diagonal_gaussian_3d',
                  'noise': {'force_sd': [1., 2., 3.], 'acceleration_sd': [.5, 1., 2.]},
                  'direction': [1., 2., -1.],
                  'scenario': {'type': 'fixed_mass', 'masses': [4.], 'signal_snrs': [4.]}}
        report, output, workspace = self.run_experiment(config)
        row = report['rows'][0]
        saved = report['experiment_config']
        self.assertAlmostEqual(row['true_mass'], 4.)
        self.assertTrue(row['refinement']['passed'])
        self.assertEqual(len(row['methods']), 16)
        self.assertIn('covariance_profile', row['methods'])
        self.assertIn('flat_joint_rms', row['methods'])
        self.assertEqual(set(row['paired_flat_minus_other']), set(row['methods']) - {'flat_joint'})
        errors = np.random.default_rng(saved['training_seed']).standard_normal((16, 6))
        with np.load(output / 'trials' / (row['id'] + '.npz')) as trials:
            expected_f = errors[:, :3] * saved['noise']['force_sd'] + row['true_force_magnitude'] * np.array(saved['direction'])
            expected_a = errors[:, 3:] * saved['noise']['acceleration_sd'] + row['true_acceleration_magnitude'] * np.array(saved['direction'])
            np.testing.assert_array_equal(trials['observation__force'], expected_f)
            np.testing.assert_array_equal(trials['observation__acceleration'], expected_a)
        browser = load_run(workspace, 'test-run')
        self.assertEqual(browser['metadata']['noise_model'], 'diagonal_gaussian_3d')
        self.assertEqual(browser['rows'][0]['direction'], saved['direction'])
        self.assertAlmostEqual(browser['rows'][0]['mass'], 4.)
        mass_metric = next(metric for metric in browser['metrics'] if metric['id'] == 'absolute_mass')
        self.assertIn('physical', mass_metric['description'])
        self.assertIn('independent coordinate noise', next(method for method in browser['methods'] if method['id'] == 'flat_joint')['assumptions'])
        json.dumps(report, allow_nan=False)

    def test_equal_noise_control_preserves_existing_results(self):
        common = {'samples': 16, 'scenario': {'type': 'force_path', 'force_snr': 8., 'acceleration_snrs': [2.]}}
        isotropic, _, _ = self.run_experiment(common)
        diagonal, _, _ = self.run_experiment(dict(common, noise_model='diagonal_gaussian_3d'))
        old, new = isotropic['rows'][0], diagonal['rows'][0]
        for method, values in old['methods'].items():
            self.assertEqual(values['means'], new['methods'][method]['means'], method)
        self.assertEqual(old['uncertainty'], new['uncertainty'])
        self.assertEqual(new['methods']['flat_joint']['means'], new['methods']['flat_joint_rms']['means'])

    def test_heldout_noise_and_magnitude_references_use_physical_units(self):
        config = {'samples': 8, 'noise_model': 'diagonal_gaussian_3d',
                  'noise': {'force_sd': [2., 2., 2.], 'acceleration_sd': [3., 3., 3.]},
                  'direction': [0., 0., 1.],
                  'scenario': {'type': 'pairs', 'pairs': [[8., 2.]]}}
        report, _, _ = self.run_experiment(config)
        row = report['rows'][0]
        self.assertAlmostEqual(row['true_mass'], 8/3)
        references = row['heldout_reference_values']
        self.assertAlmostEqual(references['known_acceleration_expected_loss_at_true_mass'], 12/(4*16**2))
        self.assertAlmostEqual(references['noisy_acceleration_expected_loss_at_true_mass'], (12+(8/3)**2*27)/(4*16**2))
        self.assertAlmostEqual(references['noisy_acceleration_optimal_fixed_coefficient'], (8/3)*144/(144+27))
        point = row['methods']['flat_joint']
        self.assertAlmostEqual(point['additional_tasks']['joint_root_loss']['mean'], np.sqrt(16*6)*point['means']['reciprocal_root'])


if __name__ == '__main__':
    unittest.main()

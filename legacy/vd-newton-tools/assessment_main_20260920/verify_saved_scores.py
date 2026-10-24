"""Independent arithmetic check of saved primary scores in the two new studies."""
import json
import math
from pathlib import Path
import sys
import numpy as np

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / '03_trace_and_evaluator/mass_estimator/investigation'))
from trial_archive import load_trial_manifest, load_trial_archive

conditions = 0
experiments = 0
simulations = 0
maximum_relative_point_refinement = 0.
maximum_absolute_interval_refinement = 0.
for study in ('assessment-main-20260920', 'assessment-calibration-20260920'):
    folder = root / '.tools/mass_estimator_studies' / study
    checkpoint = json.loads((folder / 'study_checkpoint.json').read_text())
    assert checkpoint['status'] == 'completed'
    for entry in checkpoint['experiments']:
        assert len(entry['attempts']) == 1
        path = folder / entry['attempt_directory']
        report = json.loads((path / 'results.json').read_text())
        manifest = load_trial_manifest(path)
        for row, archived in zip(report['rows'], manifest['entries'], strict=True):
            assert row['id'] == archived['id']
            arrays = load_trial_archive(path, archived)
            truth = row['true_mass']
            samples = report['samples_per_scenario']
            for name, method in row['methods'].items():
                points = arrays['point__' + name]
                valid = np.isfinite(points) & (points > 0)
                failure = ~valid | (points < truth/2) | (points > truth*2)
                log_errors = np.full(samples, np.inf)
                log_errors[valid] = np.log(points[valid])-math.log(truth)
                losses = np.minimum(log_errors**2, math.log(2)**2)
                assert method['n'] == samples
                assert method['invalid_count'] == int(np.sum(~valid))
                for metric, values in [('outside_factor_2',failure),('capped_squared_log_factor_2',losses)]:
                    np.testing.assert_allclose(method['means'][metric],np.mean(values),rtol=1e-12,atol=1e-15)
                    np.testing.assert_allclose(method['mcse'][metric],np.std(values,ddof=1)/math.sqrt(samples),rtol=1e-12,atol=1e-15)
            for law in ('flat_joint','tube_joint'):
                for level in (50,95):
                    lower = arrays[f'posterior__{law}__log_lower_{level}']
                    upper = arrays[f'posterior__{law}__log_upper_{level}']
                    covered = (lower <= math.log(truth)) & (math.log(truth) <= upper)
                    summary = row['uncertainty'][law]['intervals'][str(level)]['coverage']
                    np.testing.assert_allclose(summary['mean'],np.mean(covered),rtol=0,atol=1e-15)
            assert row['refinement']['passed']
            maximum_relative_point_refinement = max(maximum_relative_point_refinement,row['refinement']['maximum_relative_point_difference'])
            maximum_absolute_interval_refinement = max(maximum_absolute_interval_refinement,
                *(value for key,value in row['refinement']['maximum_absolute_distribution_difference'].items()
                  if key.startswith(('log_lower_', 'log_upper_'))))
            conditions += 1
            simulations += samples
        experiments += 1
print(json.dumps({'experiments':experiments,'conditions':conditions,'simulated_experiments':simulations,
                  'all_primary_point_scores_and_interval_coverages_reproduced':True,
                  'all_refinement_checks_passed':True,
                  'maximum_relative_point_refinement':maximum_relative_point_refinement,
                  'maximum_absolute_log_endpoint_refinement':maximum_absolute_interval_refinement},indent=2))

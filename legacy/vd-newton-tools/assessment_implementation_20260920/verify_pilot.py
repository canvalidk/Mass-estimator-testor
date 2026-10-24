import hashlib
import json
from pathlib import Path
import sys

workspace = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(workspace / '03_trace_and_evaluator/mass_estimator/investigation'))
from results_catalog import list_runs, load_run

study = workspace / '.tools/mass_estimator_studies/assessment-pilot-20260920'
checkpoint = json.loads((study / 'study_checkpoint.json').read_text())
assert checkpoint['status'] == 'completed'
hashed_files = 0
for experiment in checkpoint['experiments']:
    assert experiment['status'] == 'completed'
    assert len(experiment['attempts']) == 1, 'Resume must not start another attempt'
    run = study / experiment['attempt_directory']
    for relative, expected in experiment['output_sha256'].items():
        assert hashlib.sha256((run / relative).read_bytes()).hexdigest() == expected
        hashed_files += 1

original = json.loads((study / 'experiments/stable_n1/attempt-001/results.json').read_text())
rescored = json.loads((study / 'rescored-factor-3.json').read_text())
compared_methods = 0
compared_pairs = 0
for a, b in zip(original['rows'], rescored['rows'], strict=True):
    assert a['id'] == b['id']
    for method, result in b['methods'].items():
        for key, value in result['original_point_scores'].items():
            assert a['methods'][method][key] == value, (a['id'], method, key)
        compared_methods += 1
    for method, result in b['paired_reference_minus_other'].items():
        for metric, summary in result['original_point_losses'].items():
            for key in ['n', 'mean_difference', 'mcse', 'status']:
                assert a['paired_flat_minus_other'][method][metric][key] == summary[key]
            compared_pairs += 1

pilot_runs = [run for run in list_runs(workspace)
              if run['id'].startswith('study-assessment-pilot-20260920-')]
assert len(pilot_runs) == 2
for run in pilot_runs:
    normalized = load_run(workspace, run['id'])
    assert len(normalized['rows']) == 12
print(json.dumps({'resume_attempts_per_experiment': 1,
                  'verified_output_files': hashed_files,
                  'identical_rescored_method_summaries': compared_methods,
                  'identical_paired_summaries': compared_pairs,
                  'catalog_runs': [run['id'] for run in pilot_runs]}, indent=2))

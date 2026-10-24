"""Read saved study results and summarize declared scores, without refitting."""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
from comparison_scores import point_losses, paired_summary

plan = json.loads((HERE/'plan.json').read_text())
rows = []
for entry in plan['runs']:
    folder = ROOT/'.tools/mass_estimator_runs'/entry['id']
    if not (folder/'results.json').exists():
        continue
    report = json.loads((folder/'results.json').read_text())
    for row in report['rows']:
        methods = row['methods']
        ours = methods['flat_joint']
        coverage = row['uncertainty']['flat_joint']['intervals']['95']
        item = {
            'run': entry['id'].removeprefix('robust-20260919-'),
            'a': row['true_acceleration_snr'], 'n': row.get('repeats', 1),
            'ours_log': ours['means']['squared_log'],
            'norm_log': methods['norm_ratio']['means']['squared_log'],
            'ours_cap': ours['means']['capped_squared_log_factor_2'],
            'norm_cap': methods['norm_ratio']['means']['capped_squared_log_factor_2'],
            'profile_cap': methods['covariance_profile']['means']['capped_squared_log_factor_2'],
            'geometric_cap': methods['flat_geometric']['means']['capped_squared_log_factor_2'],
            'ours_success2': 1-ours['means']['outside_factor_2'],
            'norm_success2': 1-methods['norm_ratio']['means']['outside_factor_2'],
            'profile_success2': 1-methods['covariance_profile']['means']['outside_factor_2'],
            'coverage95': coverage['coverage']['mean'],
            'coverage95_mcse': coverage['coverage']['empirical_mcse'],
            'log_width95': coverage['mean_log_width']['mean'],
            'paired_vs_norm_cap': row['paired_flat_minus_other']['norm_ratio']['capped_squared_log_factor_2'],
            'paired_vs_profile_cap': row['paired_flat_minus_other']['covariance_profile']['capped_squared_log_factor_2'],
        }
        # Same training noise and calibration seeds couple runs by trial, so
        # calibration effects on our own point can be summarized as paired.
        geometry = next(g for g in ('precise-axis', 'crossed', 'equal') if item['run'].startswith(g+'-'))
        control = ROOT/'.tools/mass_estimator_runs'/f'robust-20260919-{geometry}-known'/'trials'/f"{row['id']}.npz"
        actual = folder/'trials'/f"{row['id']}.npz"
        if control.exists() and actual.exists():
            with np.load(control) as c, np.load(actual) as a:
                current_loss = point_losses(a['point__flat_joint'], row['true_mass'])
                control_loss = point_losses(c['point__flat_joint'], row['true_mass'])
                item['ours_vs_known_one_reading_cap'] = paired_summary(
                    current_loss['capped_squared_log_factor_2'], control_loss['capped_squared_log_factor_2'])
        rows.append(item)

(HERE/'analysis.json').write_text(json.dumps(rows, indent=2)+'\n')
print(json.dumps({'completed_runs': len(set(r['run'] for r in rows)), 'cells': len(rows)}))
for item in rows:
    if len(sys.argv) > 1 and sys.argv[1] not in item['run']:
        continue
    print(json.dumps(item))

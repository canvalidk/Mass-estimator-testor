"""Compact, cell-specific evidence for the conversation; no fitting or writes."""
from pathlib import Path
import json

workspace = Path(__file__).resolve().parents[2]
base = workspace / '.tools/mass_estimator_studies'
for name in ('assessment-main-20260920', 'assessment-calibration-20260920'):
    root = base / name
    checkpoint = json.loads((root / 'study_checkpoint.json').read_text())
    print(name, checkpoint['status'])
    for record in checkpoint['experiments']:
        if record['status'] != 'completed':
            print(record['id'], record['status'])
            continue
        report = json.loads((root / record['attempt_directory'] / 'results.json').read_text())
        print(record['id'], 'R', report['samples_per_scenario'])
        for row in report['rows']:
            if (row['true_mass'], row['true_acceleration_snr']) not in ((.25,2.),(1.,0.),(1.,2.),(1.,8.),(4.,8.)):
                continue
            point = row['methods']['flat_joint']
            uncertainty = row['uncertainty']['flat_joint']['intervals']
            comparisons = row['paired_flat_minus_other']
            print(json.dumps({
                'mass':row['true_mass'], 'a_snr':row['true_acceleration_snr'],
                'factor2_failure':[point['means']['outside_factor_2'], point['mcse']['outside_factor_2']],
                'capped_log_loss':[point['means']['capped_squared_log_factor_2'], point['mcse']['capped_squared_log_factor_2']],
                'log_diagnostics':point['log_diagnostics'],
                'coverage':{level:summary['coverage'] for level,summary in uncertainty.items()},
                'interval95':uncertainty['95'],
                'paired_loss':{method:comparisons[method]['capped_squared_log_factor_2'] for method in ('flat_geometric','flat_median','norm_ratio','positive_profile')},
            }))
    precision = root / 'precision_summary.json'
    if precision.exists():
        precision = json.loads(precision.read_text())
        counts = {}
        for exp in precision['experiments'].values():
            for row in exp['entries']:
                key = row['scope']+'/'+row['metric']
                count = counts.setdefault(key, [0,0])
                count[0] += row['target_met'] is True
                count[1] += 1
        print('precision_met',json.dumps(counts))

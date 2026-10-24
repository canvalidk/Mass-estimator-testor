from pathlib import Path
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
folder = Path(__file__).parent
new = json.loads((folder/'full/results.json').read_text(encoding='utf-8'))
old = json.loads((ROOT/'.tools/mass_crossover_20260916/vector_results.json').read_text())
source = ROOT/'03_trace_and_evaluator/mass_estimator/investigation'
for name, digest in new['source_sha256'].items():
    assert hashlib.sha256((source/name).read_bytes()).hexdigest() == digest, name
assert new['samples_per_scenario'] == old['sampling']['n_per_scenario'] == 32768
assert new['training_seed'] == old['sampling']['seed']
worst = 0.
for previous in old['sampling']['scenarios']:
    current = next(r for r in new['rows'] if r['family'] == 'original'
                   and r['true_force_snr'] == previous['true_force_snr']
                   and r['true_acceleration_snr'] == previous['true_acceleration_snr'])
    for name, before in previous['methods'].items():
        after = current['methods'][name]
        assert abs(before['finite_positive_fraction']-after['finite_positive_fraction']) < 1e-14
        assert abs(before['within_factor_two']-after['tolerance_success']['2']['fraction']) < 1e-14
        if before['log_mse'] is not None:
            difference = abs(before['log_mse']-after['means']['squared_log'])
            worst = max(worst, difference)
            assert difference < 1e-8, (name, difference)
for row in new['rows']:
    assert row['refinement']['passed']
    for name, summary in row['methods'].items():
        assert summary['n'] == new['samples_per_scenario']
        assert summary['valid_count'] + summary['invalid_count'] == summary['n']
        curve = summary['tolerance_curve']
        assert all(0 <= x <= 1 for x in curve)
        assert all(a <= b for a,b in zip(curve,curve[1:]))
        for key, value in summary['means'].items():
            assert value is None or math.isfinite(value)
print(json.dumps({'original_scenarios_reproduced':len(old['sampling']['scenarios']),
                  'maximum_absolute_log_mse_change':worst,
                  'source_hashes_match':True,
                  'all_trial_counts_and_tolerance_curves_valid':True},indent=2))
for row in new['rows'][:7]:
    print(json.dumps({'a':row['true_acceleration_snr'],
                      'log_mse':{name:round(v['means']['squared_log'],5)
                                 for name,v in row['methods'].items()
                                 if not v.get('population_unbounded_risk')},
                      'paired_abslog_norm':row['paired_flat_minus_other']['norm_ratio']['abs_log']['mean_difference'],
                      'paired_root_norm':row['paired_flat_minus_other']['norm_ratio']['reciprocal_root']['mean_difference']},indent=None))

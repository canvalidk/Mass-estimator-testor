"""Independent saved-array audit; imports readers only, never inference."""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / '03_trace_and_evaluator/mass_estimator/investigation'
sys.path.insert(0, str(ENGINE))
from assessment_summary import load_completed_study
from trial_archive import load_trial_manifest, load_trial_archive

study_name = sys.argv[1] if len(sys.argv) > 1 else 'assessment-calibration-integrated-20260920'
study = ROOT / '.tools/mass_estimator_studies' / study_name
loaded = load_completed_study(study)
checked = 0
missing_paired_coverage = []
rows = []
maximum_refinement = {}
current_source_differences = set()

def stats(values):
    values = np.asarray(values, dtype=float)
    return {'mean': float(np.mean(values)), 'mcse': float(np.std(values, ddof=1)/math.sqrt(len(values)))}

def same(actual, expected, label):
    global checked
    assert math.isclose(actual, expected, rel_tol=2e-11, abs_tol=2e-12), (label, actual, expected)
    checked += 1

for experiment in loaded['reports']:
    report = experiment['report']
    attempt = study / Path(experiment['relative_results']).parent
    manifest = load_trial_manifest(attempt)
    entries = {entry['id']: entry for entry in manifest['entries']}
    assert set(entries) == {row['id'] for row in report['rows']}
    for name, digest in report['source_sha256'].items():
        assert hashlib.sha256((attempt / 'source' / name).read_bytes()).hexdigest() == digest
        if hashlib.sha256((ENGINE / name).read_bytes()).hexdigest() != digest:
            current_source_differences.add(name)
    assert manifest['provenance']['source_sha256'] == report['source_sha256']
    for source in report['rows']:
        assert source['refinement']['passed']
        for name, component in source['refinement']['components'].items():
            assert component['passed']
            maximum_refinement[name] = max(maximum_refinement.get(name, 0.),
                component['maximum_relative_point_difference'],
                *component['maximum_absolute_distribution_difference'].values())
        arrays = load_trial_archive(attempt, entries[source['id']])
        n, truth = report['samples_per_scenario'], source['true_mass']
        assert n == 512 and np.array_equal(arrays['replicate_id'], np.arange(n))
        result = {'experiment': experiment['id'], 'id': source['id'], 'mass': truth,
                  'a_snr_per_reading': source['true_acceleration_snr'], 'n': source['repeats'],
                  'R': n, 'methods': {}, 'paired_flat_minus_other': {}}
        values_by_law = {}
        for name in ('flat_joint', 'calibrated_joint', 'oracle_joint'):
            points = arrays['point__' + name]
            valid = np.isfinite(points) & (points > 0)
            log_error = np.full(n, np.inf)
            log_error[valid] = np.log(points[valid]) - math.log(truth)
            losses = {
                'capped_squared_log_factor_2': np.minimum(log_error**2, math.log(2)**2),
                'outside_factor_2': (~(valid & (points >= truth/2) & (points <= truth*2))).astype(float)}
            record = {'invalid': int(np.count_nonzero(~valid)), 'point': {}, 'intervals': {}}
            assert record['invalid'] == source['methods'][name]['invalid_count']
            for metric, values in losses.items():
                record['point'][metric] = stats(values)
                for statistic, actual in record['point'][metric].items():
                    expected = source['methods'][name]['means' if statistic == 'mean' else 'mcse'][metric]
                    same(actual, expected, f'{name}/{metric}/{statistic}')
            values_by_law[name] = dict(losses)
            for metric in ('log_crps', 'log_density_score'):
                values_by_law[name][metric] = arrays[f'posterior__{name}__{metric}']
            for level in (50, 80, 95):
                lo, hi = (arrays[f'posterior__{name}__log_{endpoint}_{level}'] for endpoint in ('lower', 'upper'))
                assert np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)) and np.all(lo <= hi)
                y, alpha = math.log(truth), 1-level/100
                covered = ((lo <= y) & (y <= hi)).astype(float)
                below, above = (y < lo).astype(float), (y > hi).astype(float)
                np.testing.assert_array_equal(covered + below + above, 1.)
                interval_values = {'coverage': covered, 'below_interval': below, 'above_interval': above,
                    'mean_log_width': hi-lo, 'interval_score': hi-lo + 2/alpha*(lo-y)*below + 2/alpha*(y-hi)*above}
                values_by_law[name][f'interval_score_{level}'] = interval_values['interval_score']
                if level not in (50, 95):
                    continue
                record['intervals'][str(level)] = {}
                for metric, values in interval_values.items():
                    summary = stats(values)
                    record['intervals'][str(level)][metric] = summary
                    values_by_law[name][f'{metric}_{level}'] = values
                    for statistic, actual in summary.items():
                        same(actual, source['uncertainty'][name]['intervals'][str(level)][metric]['empirical_mcse' if statistic == 'mcse' else statistic],
                             f'{name}/{level}/{metric}/{statistic}')
            result['methods'][name] = record
        for name in ('calibrated_joint', 'oracle_joint'):
            paired = {}
            for metric, values in values_by_law['flat_joint'].items():
                summary = stats(values-values_by_law[name][metric])
                paired[metric] = {'mean_difference': summary['mean'], 'mcse': summary['mcse']}
                point_report = source['paired_flat_minus_other'][name]
                law_report = source['paired_uncertainty_flat_minus_other'][name]
                reported = point_report.get(metric) or law_report.get(metric)
                if reported:
                    same(summary['mean'], reported['mean_difference'], f'paired/{name}/{metric}/mean')
                    same(summary['mcse'], reported['mcse'], f'paired/{name}/{metric}/mcse')
                elif metric.startswith('coverage_'):
                    missing_paired_coverage.append([experiment['id'], source['id'], name, metric])
            result['paired_flat_minus_other'][name] = paired
        rows.append(result)
assert len(rows) == 24
assert 'calibration_posterior' not in sys.modules and 'compare_estimators' not in sys.modules
audit = {'status': 'passed', 'study': str(study), 'checkpoint_sha256': loaded['checkpoint_sha256'],
         'verified_artifacts': sum(item['verified_artifacts'] for item in loaded['reports']),
         'numeric_checks': checked, 'cells': len(rows), 'maximum_refinement_by_component': maximum_refinement,
         'current_source_differences': sorted(current_source_differences),
         'reporting_gap': ('Paired coverage/tail/width differences were not saved in original results; recomputed here from identical replicate arrays.'
                           if missing_paired_coverage else None),
         'missing_paired_coverage_summaries': len(missing_paired_coverage), 'rows': rows}
target = Path(__file__).with_name('independent_verified_pilot_audit.json' if study_name.endswith('-verified') else 'independent_pilot_audit.json')
with target.open('x', encoding='utf-8') as handle:
    json.dump(audit, handle, indent=2, allow_nan=False)
    handle.write('\n')
print(json.dumps({key: value for key, value in audit.items() if key != 'rows'}))
for row in rows:
    if row['mass'] == 1 and row['a_snr_per_reading'] == 8:
        print(json.dumps({'experiment': row['experiment'], 'mass': row['mass'],
                         'a_snr': row['a_snr_per_reading'],
                         'paired_coverage_95': row['paired_flat_minus_other']['calibrated_joint']['coverage_95']}))

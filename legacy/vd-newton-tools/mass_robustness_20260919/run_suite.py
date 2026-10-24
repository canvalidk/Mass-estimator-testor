"""Predeclared calibration/repeated-reading study; uses the shared engine."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / '03_trace_and_evaluator/mass_estimator/investigation'
sys.path.insert(0, str(ENGINE))
from experiment_config import normalize_config, run_config

HERE = Path(__file__).resolve().parent
STORE = ROOT / '.tools/mass_estimator_runs'


def configurations():
    common = np.array([.5, 2., 2.]) / np.sqrt(2.75)
    crossed = np.array([.5, 1., 2.]) / np.sqrt(1.75)
    geometries = {
        'equal': ([1., 1., 1.], [1., 1., 1.]),
        'precise-axis': (common.tolist(), common.tolist()),
        'crossed': (crossed.tolist(), crossed[::-1].tolist()),
    }
    settings = [('known', 1, {'mode': 'known'})]
    for factor in (.5, .75, .9, 1.1, 1.25, 1.5):
        settings.append((f'sd-{factor:g}', 1, {
            'mode': 'scaled', 'force_scale': [factor]*3, 'acceleration_scale': [factor]*3}))
    for factor in (.75, 1.25):
        settings.append((f'force-sd-{factor:g}', 1, {
            'mode': 'scaled', 'force_scale': [factor]*3, 'acceleration_scale': [1.]*3}))
    for count in (10, 30, 100):
        settings.append((f'cal-{count}', 1, {'mode': 'estimated', 'samples': count, 'seed': 2026091923}))
    for repeats in (2, 5, 10, 50):
        settings.append((f'repeats-{repeats}', repeats, {'mode': 'known'}))
    settings.append(('repeats-50-cal-10', 50, {'mode': 'estimated', 'samples': 10, 'seed': 2026091923}))
    for geometry, (sf, sa) in geometries.items():
        for tag, repeats, calibration in settings:
            run_id = f'robust-20260919-{geometry}-{tag}'.replace('.', 'p')
            label = ('Repeated readings' if repeats > 1 else 'Calibration')
            yield run_id, normalize_config({
                'name': f'{label} · {geometry} · {tag}', 'samples': 1024,
                'order': 96, 'batch_size': 128,
                'training_seed': 2026091921, 'heldout_seed': 2026091922,
                'save_trials': True, 'noise_model': 'diagonal_gaussian_3d',
                'noise': {'force_sd': sf, 'acceleration_sd': sa}, 'direction': [1., 0., 0.],
                'scenario': {'type': 'force_path', 'force_snr': 8., 'acceleration_snrs': [.1, .25, 2., 8.]},
                'measurement': {'repeats': repeats, 'calibration': deepcopy(calibration)},
            })


def main():
    suite = list(configurations())
    plan = {
        'aims': ['Sensitivity to supplied noise SDs, separately from true noise',
                 'Concentration with independent repeated readings of a fixed latent pair',
                 'Frequentist calibration of plug-in probability intervals'],
        'primary_loss': 'capped_squared_log_factor_2',
        'primary_comparators': ['norm_ratio', 'covariance_profile'],
        'secondary': ['squared_log (where finite)', 'outside_factor_2',
                      'coverage_95', 'mean_log_width_95', 'flat_geometric', 'all existing scores retained'],
        'interpretation': 'Report effects and paired Monte Carlo uncertainty per condition; no global winner or post-hoc score selection.',
        'repeated_protocol': 'Same fixed latent force/acceleration vectors; Gaussian sample mean is sufficient. Different excitation and shared sensor offsets are outside this study.',
        'estimated_calibration': 'Independent calibration per Monte Carlo trial; coordinate sample variances with k-1 degrees of freedom, plug-in rather than variance-integrated inference.',
        'coupling': 'Standardized measurement noise reused across conditions and run settings, independent calibration seed; cells are not independent replicates.',
        'runs': [{'id': key, 'config': config} for key, config in suite],
    }
    plan_path = HERE/'plan.json'
    if plan_path.exists():
        if json.loads(plan_path.read_text(encoding='utf-8')) != plan:
            raise RuntimeError('Existing predeclared plan differs; preserve it and use a new study.')
    else:
        plan_path.write_text(json.dumps(plan, indent=2)+'\n', encoding='utf-8')
    if '--plan-only' in sys.argv:
        print(json.dumps({'runs': len(suite), 'cells': len(suite)*4, 'trials_per_cell': 1024}), flush=True)
        return
    started = time.perf_counter()
    for index, (run_id, config) in enumerate(suite):
        output = STORE/run_id
        if (output/'results.json').exists():
            previous = json.loads((output/'results.json').read_text(encoding='utf-8'))
            if previous.get('experiment_config') != config:
                raise RuntimeError(f'Conflicting completed configuration: {run_id}')
            print(json.dumps({'skipped_completed': run_id}), flush=True)
            continue
        print(json.dumps({'starting': run_id, 'study': index+1, 'of': len(suite)}), flush=True)
        report = run_config(config, output)
        print(json.dumps({'finished': run_id, 'seconds': report['elapsed_seconds'],
                          'suite_seconds': time.perf_counter()-started}), flush=True)


if __name__ == '__main__':
    main()

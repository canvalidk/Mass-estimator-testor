"""Verify reporting-only rerun preserved every saved numeric trial array."""
from pathlib import Path
import json
import sys
import numpy as np

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / '03_trace_and_evaluator/mass_estimator/investigation'))
from assessment_summary import load_completed_study
from trial_archive import load_trial_manifest, load_trial_archive

old = root / '.tools/mass_estimator_studies/assessment-calibration-integrated-20260920'
new = old.with_name(old.name + '-verified')
first, second = load_completed_study(old), load_completed_study(new)
assert first['checkpoint']['specification'] == second['checkpoint']['specification']
count, elements = 0, 0
for experiment1, experiment2 in zip(first['reports'], second['reports']):
    assert experiment1['id'] == experiment2['id']
    run1, run2 = old / Path(experiment1['relative_results']).parent, new / Path(experiment2['relative_results']).parent
    manifests = load_trial_manifest(run1), load_trial_manifest(run2)
    assert [entry['id'] for entry in manifests[0]['entries']] == [entry['id'] for entry in manifests[1]['entries']]
    for entry1, entry2 in zip(manifests[0]['entries'], manifests[1]['entries']):
        a, b = load_trial_archive(run1, entry1), load_trial_archive(run2, entry2)
        assert a.keys() == b.keys()
        for name in a:
            np.testing.assert_array_equal(a[name], b[name])
            count += 1
            elements += a[name].size
    for row1, row2 in zip(experiment1['report']['rows'], experiment2['report']['rows']):
        for field in ('id', 'methods', 'uncertainty', 'paired_flat_minus_other', 'refinement'):
            assert row1[field] == row2[field], (experiment1['id'], row1['id'], field)
result = {'status': 'identical', 'array_count': count, 'numeric_elements': elements,
          'scope': 'Every saved trial array exactly equal including nonfinite values; point/law summaries, paired point summaries and refinement diagnostics exactly equal.',
          'original_checkpoint': first['checkpoint_sha256'], 'verified_checkpoint': second['checkpoint_sha256']}
target = Path(__file__).with_name('reporting_rerun_identity.json')
with target.open('x', encoding='utf-8') as handle:
    json.dump(result, handle, indent=2)
    handle.write('\n')
print(json.dumps(result))

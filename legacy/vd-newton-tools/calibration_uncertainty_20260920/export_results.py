from pathlib import Path
import json
import sys

root = Path(__file__).resolve().parents[2]
investigation = root / '03_trace_and_evaluator/mass_estimator/investigation'
sys.path.insert(0, str(investigation))
sys.path.append(str(root / '.tools/assessment_plot_dependencies'))
from assessment_summary import export_assessment
from null_truth_sweep import build_null_sweep, plot_null_sweep

study = root / '.tools/mass_estimator_studies/assessment-calibration-integrated-20260920-verified'
output = Path(__file__).resolve().parent / 'analysis'
summary = export_assessment(study, output, comparators=('calibrated_joint', 'oracle_joint', 'flat_geometric', 'flat_median'))
null = build_null_sweep(study)
(output / 'null_coverage.json').write_text(json.dumps(null, indent=2, allow_nan=False)+'\n', encoding='utf-8')
plot_null_sweep(null, output / 'null_coverage.png')
print(json.dumps({'output': str(output), 'figures': summary['figures'], 'null_status': null['status']}))

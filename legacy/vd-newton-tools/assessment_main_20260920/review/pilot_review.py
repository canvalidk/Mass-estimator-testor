import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / '.tools/mass_estimator_studies/assessment-pilot-20260920'
CAP = math.log(2)**2
METHODS = ['flat_joint','flat_geometric','flat_median','norm_ratio','positive_profile']
summary = []
checks = []
for label in ['stable_n1','stable_n4']:
    folder = BASE / 'experiments' / label / 'attempt-001'
    report = json.loads((folder/'results.json').read_text())
    print(label, 'm a coverage below above width logbias cap factor2; cap differences flat-other (MCSE), invalid profile')
    null_points = []
    null_intervals = []
    for row in report['rows']:
        m,a = row['true_mass'], row['true_acceleration_snr']
        y = math.log(m)
        data = np.load(folder/'trials'/(row['id']+'.npz'), allow_pickle=False)
        lower=data['posterior__flat_joint__log_lower_95']
        upper=data['posterior__flat_joint__log_upper_95']
        cover=(lower<=y)&(y<=upper)
        below=y<lower
        above=y>upper
        interval=row['uncertainty']['flat_joint']['intervals']['95']
        for key, values in [('coverage',cover),('below_interval',below),('above_interval',above),('mean_log_width',upper-lower)]:
            np.testing.assert_allclose(np.mean(values),interval[key]['mean'],rtol=1e-13,atol=1e-15)
            np.testing.assert_allclose(np.std(values,ddof=1)/np.sqrt(len(values)),interval[key]['empirical_mcse'],rtol=1e-13,atol=1e-15)
        losses={}
        independent={}
        for method in METHODS:
            values=data['point__'+method]
            valid=np.isfinite(values)&(values>0)
            errors=np.full(len(values),np.inf)
            errors[valid]=np.log(values[valid])-y
            losses[method]=np.minimum(errors**2,CAP)
            fail=(~valid)|(values<m/2)|(values>m*2)
            means=row['methods'][method]['means']
            np.testing.assert_allclose(np.mean(losses[method]),means['capped_squared_log_factor_2'],rtol=1e-13,atol=1e-15)
            np.testing.assert_allclose(np.mean(fail),means['outside_factor_2'],rtol=1e-13,atol=1e-15)
            independent[method]={'bounded_loss':float(losses[method].mean()),'invalid_fraction':float((~valid).mean()),'factor2_failure':float(fail.mean()),'median_estimate':float(np.median(values)),'mean_estimate':float(np.mean(values)) if valid.all() else None}
        deltas={}
        for method in METHODS[1:]:
            delta=losses['flat_joint']-losses[method]
            original=row['paired_flat_minus_other'][method]['capped_squared_log_factor_2']
            np.testing.assert_allclose(np.mean(delta),original['mean_difference'],rtol=1e-12,atol=1e-15)
            np.testing.assert_allclose(np.std(delta,ddof=1)/np.sqrt(len(delta)),original['mcse'],rtol=1e-12,atol=1e-15)
            deltas[method]={'mean':float(delta.mean()),'mcse':float(delta.std(ddof=1)/np.sqrt(len(delta)))}
        flat=row['methods']['flat_joint']
        width=float(np.mean(upper-lower))
        result={'experiment':label,'m':m,'a':a,'coverage':float(cover.mean()),'coverage_mcse':interval['coverage']['empirical_mcse'],'below':float(below.mean()),'above':float(above.mean()),'mean_log_width':width,'geometric_mean_endpoint_ratio':math.exp(width),'flat_log_bias':flat['log_diagnostics']['signed_log_bias'],'flat_log_variance':flat['log_diagnostics']['log_variance'],'methods':independent,'paired_capped_loss':deltas,'tube_coverage':row['uncertainty']['tube_joint']['intervals']['95']['coverage']['mean']}
        summary.append(result)
        display=f"{m:4g} {a:3g} {cover.mean():.4f} {below.mean():.4f} {above.mean():.4f} {width:.3f} {result['flat_log_bias']:+.3f} {losses['flat_joint'].mean():.4f} {independent['flat_joint']['factor2_failure']:.4f}"
        display+=' | '+ ' '.join(f"{k}={v['mean']:+.5f}({v['mcse']:.5f})" for k,v in deltas.items())
        print(display, 'TLS invalid',independent['positive_profile']['invalid_fraction'])
        if a == 0:
            null_points.append(data['point__flat_joint'])
            null_intervals.append(np.stack([lower,upper]))
    for values in (null_points,null_intervals):
        for another in values[1:]:
            np.testing.assert_array_equal(values[0],another)
    checks.append(label+': point summaries, paired MCSE, coverage/tails/width independently reproduced from NPZ; null point/interval arrays identical across masses')
precision=json.loads((BASE/'precision_summary.json').read_text())
planning={}
for name, entry in precision['experiments'].items():
    planning[name]={}
    for metric in ['capped_squared_log_factor_2','outside_factor_2','coverage_95']:
        es=[e for e in entry['entries'] if e['metric']==metric]
        planning[name][metric]={'max_mcse':max(e['mcse'] for e in es if e['mcse'] is not None),'max_pilot_R':max(e['pilot_variance_suggested_replicates'] or 0 for e in es),'conservative_R':max(e['conservative_bounded_variance_replicates'] for e in es),'num_met':sum(e['target_met'] is True for e in es),'num_entries':len(es)}
print('CHECKS',json.dumps(checks))
print('PLANNING',json.dumps(planning))
(Path(__file__).parent/'pilot_review.json').write_text(json.dumps({'cells':summary,'checks':checks,'planning':planning},indent=2,allow_nan=False)+'\n')

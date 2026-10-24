"""Rescore the existing IID Gaussian 3D trials without changing any estimator.

Every trial is retained. Infinite/zero profile outputs fail every finite
factor tolerance and receive the cap under explicitly capped log loss.
Scratch investigation only; production and managed sources are unchanged.
"""
from pathlib import Path
import sys
import json
import hashlib
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / '03_trace_and_evaluator/mass_estimator/investigation/crossover_vector_checks.py'
sys.path.insert(0, str(SOURCE.parent))
import crossover_vector_checks as cross

OUT = Path(__file__).parent
N = 32768
SEED = 2026091601
x, w = np.polynomial.legendre.leggauss(96)
x, w = (x+1)/2, w/2

def point(F, A):
    p, q = np.sum(F*F, axis=1), np.sum(A*A, axis=1)
    d = np.sum(F*A, axis=1)
    theta = np.where(d > 0, .5*np.arctan2(2*d, q-p),
                     np.where(p >= q, math.pi/2, 0.))
    t = np.concatenate((theta[:, None]*x,
                        theta[:, None]+(math.pi/2-theta[:, None])*x), axis=1)
    weights = np.concatenate((theta[:, None]*w,
                              (math.pi/2-theta[:, None])*w), axis=1)
    sn, cs = np.sin(t), np.cos(t)
    h2 = p[:, None]*sn*sn+q[:, None]*cs*cs+2*d[:, None]*sn*cs
    peak = np.where(d >= 0, (p+q+np.hypot(p-q,2*d))/2, np.maximum(p,q))
    density = np.exp((h2-peak[:, None])/2)*weights
    return np.sum(sn*density,axis=1)/np.sum(cs*density,axis=1)

def mean_se(v):
    return [float(np.mean(v)), float(np.std(v,ddof=1)/math.sqrt(len(v)))]

errors = np.random.default_rng(SEED).standard_normal((N,6))
factors = np.array([1.1,1.25,1.5,2,3,5,10,20])
curve_factors = np.geomspace(1,20,401)
old = json.loads((ROOT/'.tools/mass_crossover_20260916/vector_results.json').read_text())
rows, arrays = [], {}
for a in [.25,.5,1.,2.,3.,5.,8.]:
    f, truth = 8., 8./a
    F, A = errors[:,:3]+[f,0,0], errors[:,3:]+[a,0,0]
    ours = np.concatenate([point(F[i:i+512],A[i:i+512]) for i in range(0,N,512)])
    checked, _ = cross.infer(F[:128],A[:128],truth,order=128)
    discrepancy = float(np.max(abs(ours[:128]/checked[:,0]-1)))
    assert discrepancy < 1e-9, discrepancy
    methods = {'equation':ours, 'norm':np.linalg.norm(F,axis=1)/np.linalg.norm(A,axis=1),
               'profile':cross.profile_point(F,A)}
    row = {'force_snr':f,'acceleration_snr':a,'truth':truth,'methods':{},'paired':{},
           'quadrature_check':discrepancy}
    losses = {}
    reference = next(r for r in old['sampling']['scenarios']
                     if r['true_force_snr']==f and r['true_acceleration_snr']==a)
    for name, m in methods.items():
        valid = (m>0)&np.isfinite(m)
        logerr = np.full(N,np.inf)
        logerr[valid] = np.log(m[valid]/truth)
        ae = abs(logerr)
        success = [(ae<=np.log(k)).astype(float) for k in factors]
        stats = {'invalid_fraction':float(np.mean(~valid)),
                 'factor_success':{str(k):mean_se(s) for k,s in zip(factors,success)},
                 'capped_log_mse':{str(k):mean_se(np.minimum(ae**2,np.log(k)**2)) for k in [2,5,10,20]},
                 'log_mse':mean_se(ae**2) if valid.all() and name!='profile' else None}
        if name=='profile':
            stats['population_log_risk']='infinite: boundary outputs have positive probability even if absent from this finite sample'
        if valid.all():
            stats.update(log_bias=float(logerr.mean()),log_variance=float(logerr.var()),
                         median_ratio=float(np.median(m/truth)))
        old_name = {'equation':'flat_joint','norm':'norm_ratio','profile':'positive_profile'}[name]
        assert abs(stats['factor_success']['2.0'][0]-reference['methods'][old_name]['within_factor_two'])<1e-14
        if name!='profile':
            assert abs(stats['log_mse'][0]-reference['methods'][old_name]['log_mse'])<1e-10
        # Dimensionless squared distance of the TRUE standardized pair from
        # the estimated positive ray. Boundaries are retained as limiting rays.
        angle = np.arctan(m)
        ray_error = (f*np.cos(angle)-a*np.sin(angle))**2
        stats['true_pair_ray_loss']=mean_se(ray_error)
        row['methods'][name]=stats
        losses[name]=ae
        arrays[f'{a}_{name}']=m
        arrays[f'{a}_{name}_success_curve']=np.array([np.mean(ae<=np.log(k)) for k in curve_factors])
    for name in ['norm','profile']:
        row['paired'][name] = {'factor_gaps':{str(k):mean_se((losses['equation']<=np.log(k)).astype(float)-(losses[name]<=np.log(k)).astype(float)) for k in factors},
                              'capped_log_mse_gaps':{str(k):mean_se(np.minimum(losses['equation']**2,np.log(k)**2)-np.minimum(losses[name]**2,np.log(k)**2)) for k in [2,5,10,20]}}
    rows.append(row)
    print(json.dumps({'a':a,'methods':{name:{'success':[round(100*s[0],2) for s in stats['factor_success'].values()],
                                         'log_mse':stats['log_mse'],'bias':stats.get('log_bias'),'variance':stats.get('log_variance'),
                                         'capped_10':stats['capped_log_mse']['10'],'ray_loss':stats['true_pair_ray_loss']}
                                     for name,stats in row['methods'].items()}}),flush=True)

report = {'n':N,'seed':SEED,'factors':factors.tolist(),'managed_commit':'8782a293297330f7a354d87ebe62792018fa8866',
          'assumptions':'True codirectional 3D vectors, independent isotropic Gaussian component errors of unit SD; same draws as original study; each trial is one unpooled pair.',
          'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows}
(OUT/'scores.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
np.savez_compressed(OUT/'trials.npz',curve_factors=curve_factors,**arrays)

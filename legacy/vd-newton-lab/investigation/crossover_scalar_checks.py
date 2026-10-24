"""Fair crossover comparison for the known-direction scalar mass model.

Research contribution, 2026-09-16. Python + NumPy, no SciPy.
All methods see exactly the same independent Gaussian errors. Raw-ratio invalid
outputs remain failures in unconditional success and capped-log-loss metrics.
No raw-ratio sampling MSE is computed or claimed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LAB))
from estimator import positive_normal_mean

GRID=(.25,.5,1.,2.,3.,5.,8.)
C=math.sqrt(2/math.pi)
LOG2=math.log(2.)
LOG10=math.log(10.)
METHODS=('ratio_of_means','hard_floor','smooth_algebraic','posterior_geometric','raw_signed_ratio')


def cdf(x):
    x=np.asarray(x,dtype=float)
    return np.fromiter((.5*math.erfc(-float(t)/math.sqrt(2)) for t in x.flat),
                       dtype=float,count=x.size).reshape(x.shape)


def h(t):
    t=np.asarray(t,dtype=float)
    return t+np.exp(-t*t/2)/math.sqrt(2*math.pi)/cdf(t)


def expected_log(t,order=256):
    """E(log A|t,A>0), flat-positive prior and unit Gaussian noise.

    Integrate in u=log(A), removing the logarithmic endpoint singularity.
    Truncation u>=-40 omits under 1e-13 of the log integral on tested range;
    A<=max(t,0)+12 leaves an exponentially small Gaussian tail.
    """
    t=np.asarray(t,dtype=float)
    result=np.empty_like(t)
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for start in range(0,len(t),128):
        x=t[start:start+128]
        upper=np.log(np.maximum(x,0)+12.)
        width=upper+40.
        u=-40+(nodes[:,None]+1)*width[None,:]/2
        a=np.exp(u)
        # Divide by Phi(x) via log; positive-normal normalizer.
        density=np.exp(u-(a-x[None,:])**2/2-.5*math.log(2*math.pi)-np.log(cdf(x))[None,:])
        result[start:start+128]=(weights[:,None]*u*density).sum(0)*width/2
    return result


def metric(estimate,truth,guaranteed_positive):
    relative=estimate/truth
    valid=np.isfinite(relative)&(relative>0)
    signed_log=np.full(len(relative),np.nan)
    signed_log[valid]=np.log(relative[valid])
    absolute_log=np.where(valid,np.abs(signed_log),math.inf)
    loss=np.minimum(absolute_log,LOG10)**2
    success=valid&(absolute_log<=LOG2)
    # Infinity in median means >= half of outputs are unusable for log error.
    median_abs=float(np.median(absolute_log))
    result={
        'positive_finite_rate':float(valid.mean()),
        'nonpositive_or_nonfinite_rate':float((~valid).mean()),
        'factor2_probability':float(success.mean()),
        'factor10_failure_probability':float((absolute_log>LOG10).mean()),
        'capped_log_squared_error_mean':float(loss.mean()),
        'median_absolute_log_error_invalid_as_infinity':median_abs if math.isfinite(median_abs) else None,
        'relative_estimate_quantiles_05_50_95':np.quantile(relative,[.05,.5,.95]).tolist(),
        'mean_log_error_positive_outputs_only':float(signed_log[valid].mean()),
        'mean_log_squared_error_positive_outputs_only':float((signed_log[valid]**2).mean()),
    }
    # For methods positive on every input, these are unconditional log risks.
    result['mean_log_squared_error_unconditional']=float((signed_log**2).mean()) if guaranteed_positive else None
    return result,success,loss


def run(n,seed):
    rng=np.random.default_rng(seed)
    e=rng.standard_normal((n,2))
    observed={theta:theta+e for theta in GRID}
    assert min(float(x.min()) for x in observed.values()) > -10.
    assert max(float(x.max()) for x in observed.values()) < 16.
    # Smooth posterior decision comparator: deterministic interpolation of its
    # exact one-dimensional expectation, independent of simulation outcomes.
    table_x=np.linspace(-10.,16.,13001)
    table_log=expected_log(table_x,256)
    validation_x=np.r_[np.linspace(-10,16,111),0.]
    validation_reference=expected_log(validation_x,768)
    validation_interp=np.interp(validation_x,table_x,table_log)
    validation_error=float(np.max(np.abs(validation_interp-validation_reference)))
    halfnormal_expected_log=-(.5772156649015329+math.log(2))/2
    h_error=max(abs(float(h(x))-positive_normal_mean(float(x),1.)) for x in validation_x)
    assert h_error<2e-11,h_error
    assert validation_error<3e-5,validation_error
    assert abs(float(expected_log(np.array([0.]),768)[0])-halfnormal_expected_log)<2e-12
    transforms={}
    for theta in GRID:
        x=observed[theta]
        transforms[theta]={
            'ratio_of_means':h(x),
            'hard_floor':np.maximum(x,C),
            # Stable algebraic evaluation for large negative x.
            'smooth_algebraic':np.where(x>=0,(x+np.sqrt(x*x+4*C*C))/2,
                                      2*C*C/(np.sqrt(x*x+4*C*C)-x)),
            'posterior_geometric':np.exp(np.interp(x,table_x,table_log)),
            'raw_signed_ratio':x,
        }
    rows=[]
    for tf in GRID:
        for ta in GRID:
            truth=tf/ta
            result={};successes={};losses={}
            for method in METHODS:
                estimates=transforms[tf][method][:,0]/transforms[ta][method][:,1]
                result[method],successes[method],losses[method]=metric(estimates,truth,method!='raw_signed_ratio')
            paired={}
            for other in METHODS[1:]:
                ds=successes['ratio_of_means'].astype(float)-successes[other].astype(float)
                dl=losses['ratio_of_means']-losses[other]
                paired[other]={
                    'factor2_advantage_of_ratio_of_means':float(ds.mean()),
                    'factor2_paired_standard_error':float(ds.std(ddof=1)/math.sqrt(n)),
                    'capped_log_risk_excess_of_ratio_of_means':float(dl.mean()),
                    'capped_log_risk_paired_standard_error':float(dl.std(ddof=1)/math.sqrt(n)),
                }
            x,y=observed[tf][:,0],observed[ta][:,1]
            rows.append({'theta_force':tf,'theta_acceleration':ta,
                         'true_mass_over_instrument_ratio':truth,
                         'either_reading_nonpositive_rate':float(((x<=0)|(y<=0)).mean()),
                         'both_readings_nonpositive_rate':float(((x<=0)&(y<=0)).mean()),
                         'raw_nonpositive_probability_exact':float(cdf(-tf)*cdf(ta)+cdf(tf)*cdf(-ta)),
                         'metrics':result,'paired_comparisons':paired})
    deterministic=[]
    for t in (0.,.25,.5,1.,2.,3.,5.,8.):
        deterministic.append({'reading_over_sigma':t,'h':float(h(t)),
            'hard_floor':max(t,C),'smooth_algebraic':(t+math.sqrt(t*t+4*C*C))/2,
            'posterior_geometric_channel':math.exp(float(expected_log(np.array([t]),768)[0]))})
    return {'scope':'Known fixed direction; independent Gaussian scalar errors; known positive sigmas; positive true f,a; flat df da posterior for probability-based methods.',
            'managed_commit':'8782a293297330f7a354d87ebe62792018fa8866',
            'numpy_version':np.__version__,
            'estimator_sha256':hashlib.sha256((LAB/'estimator.py').read_bytes()).hexdigest(),
            'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'seed':seed,'trials_per_grid_cell':n,'signal_to_noise_grid':GRID,
            'error_draws':'IID Gaussian pairs reused across methods and signal settings; each individual cell consists of IID trials.',
            'scale':'Set sf=sa=1 without loss for scale-equivariant estimators; truth=tf/ta.',
            'hard_floor_constant':C,'capped_log_loss_cap':LOG10**2,
            'raw_invalid_policy':'Positive finite output required for log error; otherwise factor2 failure and maximum capped log loss. Relative-estimate quantiles retain negative raw outputs. Positive-subset log metrics are explicitly conditional; no unconditional raw log or ordinary MSE is claimed.',
            'posterior_geometric_definition':'exp(E(log f | X)-E(log a | Y)), the posterior action minimizing expected squared log error under flat df da. It is a different readout, not the composition theorem readout.',
            'validation':{'expected_log_table_max_error_against_order768':validation_error,
                          'h_max_error_against_existing_scalar_implementation':h_error,
                          'quadrature_halfnormal_log_moment_exact':halfnormal_expected_log},
            'deterministic_crossover':deterministic,'cells':rows}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--trials',type=int,default=500000)
    parser.add_argument('--seed',type=int,default=20260916)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run(args.trials,args.seed)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps({'output':str(args.output),'validation':result['validation'],
                      'cells':len(result['cells']),'trials_per_cell':args.trials},indent=2))

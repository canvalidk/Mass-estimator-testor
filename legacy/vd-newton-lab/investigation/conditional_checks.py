"""Conditional likelihood-ratio confidence sets and excitation-strength checks.

Research implementation for independent isotropic Gaussian channels with a
common instrument uncertainty ratio. Candidates are positive masses; the LR
reference alternative includes all signed slopes and its projective endpoints.
The construction adapts the known-covariance conditional LR method of Moreira
(2003). It does not define a posterior or filter the ratio-of-means readout.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from joint_checks import (SOURCE_COMMIT, chi_quantile, confidence_set,
                          coverage_flags, gamma_p, json_safe, nodes,
                          positive_roots, statistics)


def gamma_array(a,x):
    """Vectorized regularized incomplete gamma; fixed positive shape a."""
    x=np.asarray(x,dtype=float)
    result=np.zeros_like(x)
    positive=x>0
    series=positive&(x<a+1)
    for mask,is_series in ((series,True),(positive&~series,False)):
        xx=x[mask]
        if not xx.size: continue
        factor=np.exp(a*np.log(xx)-xx-math.lgamma(a))
        if is_series:
            term=np.full_like(xx,1/a)
            total=term.copy()
            for j in range(1,10000):
                term*=xx/(a+j)
                total+=term
                if np.all(abs(term)<abs(total)*3e-15): break
            else: raise ArithmeticError('Gamma series did not converge')
            result[mask]=factor*total
        else:
            tiny=1e-300
            b=xx+1-a
            c=np.full_like(xx,1/tiny)
            d=1/b
            h=d.copy()
            for j in range(1,10000):
                an=-j*(j-a)
                b+=2
                d=b+an*d
                c=b+an/c
                d=np.where(abs(d)<tiny,tiny,d)
                c=np.where(abs(c)<tiny,tiny,c)
                d=1/d
                change=d*c
                h*=change
                if np.all(abs(change-1)<3e-15): break
            else: raise ArithmeticError('Gamma fraction did not converge')
            result[mask]=1-factor*h
    return np.clip(result,0.,1.)


def cdf_weights(c,total,k,order):
    """Resolve the chi-square transition within |normal|<=8 integration.

    The discarded tail is at most erfc(8/sqrt(2))<1.25e-15, uniformly in c,s,K.
    The integration endpoint must depend on c, not c+s: using the latter can
    yield false refinement agreement when c+s is much larger than c.
    """
    c,total=np.broadcast_arrays(np.asarray(c,float),np.asarray(total,float))
    df=k-1
    spread=12*math.sqrt(2*df)
    qlow,qhigh=max(0.,df-spread),df+spread
    denominator=np.maximum(total,1e-300)
    endpoint=np.arcsin(np.minimum(1.,8/np.sqrt(np.maximum(c,1e-300))))
    begin=np.minimum(endpoint,np.arccos(np.sqrt(np.minimum(1.,qhigh/denominator))))
    end=np.minimum(endpoint,np.arccos(np.sqrt(np.minimum(1.,qlow/denominator))))
    bounds=np.stack([np.zeros_like(total),begin,end,endpoint],axis=-1)
    t,w=nodes(order)
    width=bounds[...,1:]-bounds[...,:-1]
    beta=bounds[...,:-1,None]+width[...,None]*t
    weights=width[...,None]*w
    cosine=np.cos(beta)
    w_cdf=gamma_array(df/2,total[...,None,None]*cosine**2/2)
    weighted=math.sqrt(2/math.pi)*weights*cosine*w_cdf
    shape=total.shape+(-1,)
    return (np.sin(beta)**2).reshape(shape),weighted.reshape(shape)


def conditional_cdf(c,s,k,order=64):
    """P{U/c+W/(c+s)<=1}, U~chi1 and W~chi(k-1), independent."""
    c,s=np.broadcast_arrays(np.asarray(c,float),np.asarray(s,float))
    if np.any(c<0) or np.any(s<0): raise ValueError('Nonnegative c and s required')
    if k==1: return gamma_array(.5,c/2)
    sine2,weighted=cdf_weights(c,c+s,k,order)
    integral=np.sum(weighted*np.exp(-c[...,None]*sine2/2),axis=-1)
    return np.clip(np.sqrt(c)*integral,0.,1.)


def eigenvalues(a,b,c):
    a,b,c=np.asarray(a),np.asarray(b),np.asarray(c)
    gap=np.hypot(a-b,2*c)
    maximum=(a+b+gap)/2
    # Determinant form avoids subtracting nearly equal eigenvalues.
    minimum=np.zeros_like(np.asarray(maximum,float))
    np.divide(np.maximum(0.,a*b-c*c),maximum,out=minimum,where=maximum>0)
    return minimum,maximum


def lr_statistic(t,s,z_squared):
    difference=t-s
    discriminant=np.hypot(difference,2*np.sqrt(np.maximum(0.,s*z_squared)))
    value=(difference+discriminant)/2
    denominator=discriminant-difference
    stable=np.zeros_like(np.asarray(value,float))
    np.divide(2*s*z_squared,denominator,out=stable,where=denominator>0)
    return np.where(difference<0,stable,value)


def strength_gain(q):
    """Twice per-coordinate profile gain of a learned nonnegative spike."""
    delta=np.maximum(0.,np.asarray(q)-1)
    return delta-np.log1p(delta)


def fitted_threshold(maximum,k,level=.95,order=64):
    """CLR accepts iff |S_r|²>=returned threshold. Exact integral, numeric root.

    For L=lambda_max<=chiK(level), threshold zero accepts every mass. Otherwise
    solve F(c;L-c)=level and return L-c. The analytic bounds q_chi1<=c<=q_chiK
    keep calibration in its valid range even for enormous observed eigenvalues.
    """
    maximum=np.asarray(maximum,float)
    result=np.zeros_like(maximum)
    mask=maximum>chi_quantile(k,level)
    if not np.any(mask): return result
    large=maximum[mask]
    if k==1:
        result[mask]=large-chi_quantile(1,level)
        return result
    lower=np.full_like(large,chi_quantile(1,level))
    upper=np.full_like(large,chi_quantile(k,level))
    for _ in range(52):
        c=(lower+upper)/2
        value=conditional_cdf(c,large-c,k,order)
        lower=np.where(value<level,c,lower)
        upper=np.where(value>=level,c,upper)
    result[mask]=large-(lower+upper)/2
    return result


def quadratic_set(a,b,c,threshold):
    """Positive r with (a-threshold)r²+2cr+b-threshold>=0.

    Prototype floating-point inversion; finite endpoints included, zero and
    infinity open. Degenerate fixtures and generic grids checked below.
    """
    if threshold<=0: return [[0.,math.inf]]
    roots=positive_roots([b-threshold,2*c,a-threshold])
    mapping={math.atan(r):r for r in roots}
    angles=sorted(set([0.,math.pi/2]+list(mapping)))
    def value(theta):
        sn,cs=math.sin(theta),math.cos(theta)
        return a*sn*sn+2*c*sn*cs+b*cs*cs-threshold
    accepted=[]
    for lo,hi in zip(angles[:-1],angles[1:]):
        if value((lo+hi)/2)>=0: accepted.append([lo,hi])
    tolerance=1e-9*max(1.,a+b)
    for theta in angles[1:-1]:
        if abs(value(theta))<=tolerance and not any(lo<=theta<=hi for lo,hi in accepted):
            accepted.append([theta,theta])
    merged=[]
    for lo,hi in sorted(accepted):
        if merged and lo<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],hi)
        else: merged.append([lo,hi])
    return [[0. if lo==0 else mapping[lo],math.inf if hi==math.pi/2 else mapping[hi]]
            for lo,hi in merged]


def conditional_set(a,b,c,k,level=.95,order=64):
    _,maximum=eigenvalues(a,b,c)
    threshold=float(fitted_threshold(maximum,k,level,order))
    return quadratic_set(a,b,c,threshold)


def checks():
    assert abs(float(conditional_cdf(1e8,1e12,3))-1)<2e-14
    # Threshold subtraction at L=1e12 is limited by the spacing of float64.
    extreme=1e12-float(fitted_threshold(1e12,3))
    assert abs(extreme-chi_quantile(1,.95))<np.spacing(1e12)
    for a in (.5,1.,1.5,11.5,95.5,767.5,6143.5):
        x=np.r_[0.,np.geomspace(.01,max(10.,3*a),61)]
        expected=np.array([gamma_p(a,float(v)) for v in x])
        assert np.max(abs(gamma_array(a,x)-expected))<3e-12
    for k in (1,2,3,24,192,1536):
        c=chi_quantile(k,.95)
        assert abs(float(conditional_cdf(c,0.,k))-.95)<2e-10
        c=np.array([.1,chi_quantile(1,.95),c])
        s=np.array([0.,float(k),10.*k])
        low=conditional_cdf(c,s,k,128)
        high=conditional_cdf(c,s,k,256)
        assert np.max(abs(low-high))<2e-9,(k,low,high)
        maximum=np.array([0.,.9*chi_quantile(k,.95),2.*k+10,10.*k+10])
        threshold=fitted_threshold(maximum,k)
        active=threshold>0
        assert np.max(abs(conditional_cdf(maximum[active]-threshold[active],threshold[active],k,256)-.95))<2e-9
    rng=np.random.default_rng(20260917)
    max_refinement=0.
    for k in (1,3,24,192):
        for _ in range(12):
            x,y=rng.normal(size=(2,k))
            a,b,c=float(x@x),float(y@y),float(x@y)
            intervals=conditional_set(a,b,c,k)
            minimum,maximum=eigenvalues(a,b,c)
            for r in np.geomspace(1e-4,1e4,33):
                t,z,w,s=statistics(a,b,c,r)
                lr=lr_statistic(t,s,z)
                assert abs(lr-(t-minimum))<1e-10*max(1.,a+b)
                f=float(conditional_cdf(lr,s,k))
                member=any(lo<=r<=hi for lo,hi in intervals)
                assert member==(f<=.95),(k,a,b,c,r,intervals,f)
                if abs(f-.95)<.1:
                    error=abs(f-float(conditional_cdf(lr,s,k,256)))
                    max_refinement=max(max_refinement,error)
                    assert error<2e-9
    examples={}
    for name,a,b,c,k in (
        ('observed_zero',0.,0.,0.,3),('strong_aligned',64.,16.,32.,3),
        ('exact_anti_aligned',7.,28.,-14.,3),('strong_opposed',100.,100.,-90.,3),
        ('orthogonal',100.,100.,0.,3),('many_trials_mean_gram',3584.,2048.,1024.,1536)):
        examples[name]={kind:confidence_set(a,b,c,k,kind) for kind in ('omnibus','score','split')}
        examples[name]['conditional_lr']=conditional_set(a,b,c,k)
    return {'examples':examples,'max_checked_cdf_refinement':max_refinement}


def run(replicates=4096,seed=20260917):
    result=checks()
    rng=np.random.default_rng(seed)
    rows=[]
    # N latent groups, J repeated readings per group, standardized true norm t.
    scenarios=[(1,1,16.,.25),(8,1,2.,1.),(64,1,2.,1.),(512,1,2.,1.),
               (64,8,2.,1.),(8,64,2.,1.),(1,512,2.,1.),
               (1,1,2.,0.),(64,1,2.,0.),(512,1,2.,0.),(4096,1,2.,0.)]
    for groups,repeats,r,t in scenarios:
        k=3*groups
        flags={key:[] for key in ('omnibus','score','split','conditional_lr','fixed_lr_chi1','strength_lr_chi1')}
        fixed_lrs,strength_lrs=[],[]
        max_refinement=0.
        for start in range(0,replicates,64):
            count=min(64,replicates-start)
            x,y=rng.normal(size=(2,count,groups,3))
            # Sufficient within-group averages, standardized by sigma/sqrt(J).
            for i in range(groups):
                x[:,i,i%3]+=r*t*math.sqrt(repeats)
                y[:,i,i%3]+=t*math.sqrt(repeats)
            a,b,c=(x*x).sum((1,2)),(y*y).sum((1,2)),(x*y).sum((1,2))
            current=coverage_flags(a,b,c,k,r)
            residual,z,w,fitted=statistics(a,b,c,r)
            lr=lr_statistic(residual,fitted,z)
            conditional=conditional_cdf(lr,fitted,k)
            _,maximum=eigenvalues(a,b,c)
            strength_lr=k*np.maximum(0.,strength_gain(maximum/k)-strength_gain(fitted/k))
            current['conditional_lr']=conditional<=.95
            current['fixed_lr_chi1']=lr<=chi_quantile(1,.95)
            current['strength_lr_chi1']=strength_lr<=chi_quantile(1,.95)
            for key in flags: flags[key].extend(current[key].tolist())
            fixed_lrs.extend(lr.tolist())
            strength_lrs.extend(strength_lr.tolist())
            if start==0:
                high=conditional_cdf(lr,fitted,k,256)
                max_refinement=float(np.max(abs(high-conditional)))
                assert max_refinement<1e-7,(groups,repeats,max_refinement)
                threshold=fitted_threshold(maximum,k)
                assert np.all((fitted>=threshold)==current['conditional_lr'])
        row={'latent_groups':groups,'readings_per_group':repeats,'total_readings':groups*repeats,
             'mass_over_s':r,'acceleration_snr_per_reading':t,
             'coverage':{key:float(np.mean(value)) for key,value in flags.items()},
             'fixed_lr_quantiles':np.quantile(fixed_lrs,[.5,.95]).tolist(),
             'strength_lr_quantiles':np.quantile(strength_lrs,[.5,.95]).tolist(),
             'cdf_refinement_first64':max_refinement}
        print(json.dumps(row),flush=True)
        rows.append(row)
    result.update({'replicates_per_scenario':replicates,'seed':seed,'calibration':rows})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--replicates',type=int,default=4096)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.replicates<64: parser.error('Use at least 64 replicates')
    result=checks() if args.checks_only else run(args.replicates)
    result.update({'source_commit':SOURCE_COMMIT,
                   'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'joint_checks_sha256':hashlib.sha256(Path(__file__).with_name('joint_checks.py').read_bytes()).hexdigest()})
    print('All assertions passed.',flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(json_safe(result),indent=2,allow_nan=False)+'\n',encoding='utf-8')

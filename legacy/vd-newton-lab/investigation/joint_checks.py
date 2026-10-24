"""Joint common-mass posterior and exact Gaussian confidence-set benchmarks.

Research implementation: independent isotropic channels with common uncertainty
ratio, after standardization. Three confidence sets retain unbounded/disconnected
shapes. No set filters the mass point. NumPy and Python standard library only.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np

SOURCE_COMMIT = '8782a293297330f7a354d87ebe62792018fa8866'
KINDS = ('omnibus', 'score', 'split')


@lru_cache(None)
def nodes(n):
    z,w = np.polynomial.legendre.leggauss(n)
    return (z+1)/2,w/2


def gamma_p(a,x):
    """Regularized lower incomplete gamma, series/continued fraction."""
    if x <= 0:
        return 0.
    factor = math.exp(a*math.log(x)-x-math.lgamma(a))
    if x < a+1:
        term = total = 1/a
        for j in range(1,10000):
            term *= x/(a+j)
            total += term
            if abs(term) < abs(total)*2e-15:
                return min(1.,max(0.,factor*total))
    else:
        tiny = 1e-300
        b = x+1-a
        c,d = 1/tiny,1/b
        h = d
        for j in range(1,10000):
            an = -j*(j-a)
            b += 2
            d = b+an*d
            c = b+an/c
            if abs(d) < tiny: d = tiny
            if abs(c) < tiny: c = tiny
            d = 1/d
            change = d*c
            h *= change
            if abs(change-1) < 2e-15:
                return min(1.,max(0.,1-factor*h))
    raise ArithmeticError('Gamma evaluation did not converge')


@lru_cache(None)
def chi_quantile(df,p):
    if df < 1 or not 0 < p < 1:
        raise ValueError('Positive degrees of freedom and interior probability required')
    if df == 1:
        return NormalDist().inv_cdf((1+p)/2)**2
    lo,hi = 0.,float(df+10*math.sqrt(2*df)+20)
    while gamma_p(df/2,hi/2) < p: hi *= 2
    for _ in range(70):
        mid = (lo+hi)/2
        if gamma_p(df/2,mid/2) < p: lo = mid
        else: hi = mid
    return (lo+hi)/2


def h_chi(k,noncentral_squared,order=96):
    """Mean norm of N(mu,I_k), with |mu|^2 supplied; one-dimensional integral."""
    t,w = nodes(order)
    beta = t*math.pi/2
    sine2 = np.sin(beta)**2
    logarithm = .5*k*np.log1p(-sine2)-.5*np.asarray(noncentral_squared)[...,None]*sine2
    return math.sqrt(2/math.pi)*math.pi/2*np.sum(
        w*(-np.expm1(logarithm))/sine2,axis=-1)


def posterior(a,b,c,k,true_ratio,order=64,h_order=96):
    """Jeffreys joint posterior: CDF at true_ratio and a declared stacked readout.

    a,b,c are Gram summaries of standardized stacked force/acceleration. The
    readout uses expected norms of the full stacks; this is an explicit choice
    for multiple trials. Physical mass = returned ratio * sigma_F/sigma_a.
    """
    a,b,c = (np.asarray(v) for v in (a,b,c))
    t,w = nodes(order)
    split = math.atan(true_ratio)
    theta = np.r_[split*t,split+(math.pi/2-split)*t]
    weights = np.r_[split*w,(math.pi/2-split)*w]
    sn,cs = np.sin(theta),np.cos(theta)
    residual2 = a[...,None]*cs**2+b[...,None]*sn**2-2*c[...,None]*sn*cs
    fitted2 = a[...,None]*sn**2+b[...,None]*cs**2+2*c[...,None]*sn*cs
    residual2,fitted2 = np.maximum(0.,residual2),np.maximum(0.,fitted2)
    # Centering exponent is a numerical factor shared by all theta at each datum.
    likelihood = np.exp(-.5*(residual2-residual2.min(axis=-1,keepdims=True)))
    density = likelihood*h_chi(k,fitted2,h_order)*weights
    cdf = density[...,:order].sum(-1)/density.sum(-1)
    first_moments = likelihood*(k+fitted2)*weights
    readout = (first_moments*sn).sum(-1)/(first_moments*cs).sum(-1)
    return readout,cdf


def statistics(a,b,c,r):
    a,b,c=(np.asarray(v,dtype=float) for v in (a,b,c))
    if r<=1:
        cs=1/math.sqrt(1+r*r)
        sn=r*cs
    else:
        sn=1/math.sqrt(1+1/r/r)
        cs=sn/r
    # Reconstruct a two-coordinate representative of the Gram matrix. This
    # avoids cancellation of |S|^2 near exactly antiparallel data.
    xs=np.sqrt(a)
    ys=np.zeros_like(xs)
    np.divide(c,xs,out=ys,where=xs>0)
    rank_one=a*b==c*c
    yp=np.sqrt(np.where(rank_one,0.,np.maximum(0.,b-ys*ys)))
    rp,sp=xs*cs-ys*sn,xs*sn+ys*cs
    rt,st=-yp*sn,yp*cs
    # Preserve an exact S=0 candidate from the supplied Gram relation. A
    # square-root round trip can otherwise introduce a spurious perpendicular
    # component, e.g. a=7,b=28,c=-14,r=2.
    null_ratio=np.full_like(xs,np.nan)
    np.divide(-c,a,out=null_ratio,where=a>0)
    singular=rank_one&(r==null_ratio)
    sp,st=np.where(singular,0.,sp),np.where(singular,0.,st)
    t,s=rp*rp+rt*rt,sp*sp+st*st
    cross=rp*sp+rt*st
    score = np.zeros_like(np.asarray(t,dtype=float))
    np.divide(cross*cross,s,out=score,where=s>0)
    score = np.minimum(t,np.maximum(0.,score))
    orthogonal = np.maximum(0.,t-score)
    # Explicit null-event convention at S=0: retain it for score/split sets.
    return t,score,orthogonal,s


def coverage_flags(a,b,c,k,r,level=.95):
    t,z,w,s = statistics(a,b,c,r)
    p = math.sqrt(level) if k > 1 else level
    split_ok = z <= chi_quantile(1,p)
    if k > 1: split_ok &= w <= chi_quantile(k-1,p)
    return {'omnibus':t<=chi_quantile(k,level),
            'score':(z<=chi_quantile(1,level))|(s==0),
            'split':split_ok|(s==0)}


def positive_roots(coefficients):
    coefficients = np.asarray(coefficients,dtype=float)
    while len(coefficients)>1 and coefficients[-1]==0: coefficients=coefficients[:-1]
    if np.max(abs(coefficients))==0: return []
    roots = np.polynomial.polynomial.polyroots(coefficients/np.max(abs(coefficients)))
    selected = sorted(float(v.real) for v in roots
                      if v.real>0 and abs(v.imag)<1e-8*(1+abs(v.real)))
    unique=[]
    for value in selected:
        if not unique or abs(value-unique[-1])>1e-7*(1+value): unique.append(value)
    return unique


def confidence_set(a,b,c,k,kind='score',level=.95):
    """Prototype polynomial inversion; finite endpoints included, 0/infinity open.

    Returns dimensionless mass intervals. Root accuracy is checked on fixtures
    and generic sampled Gram matrices below; extreme or unresolved multiple
    roots are outside that numerical validation. The coverage proofs concern
    the exact inequalities, not a floating-point root routine.
    """
    if kind not in KINDS: raise ValueError('Unknown confidence construction')
    if min(a,b)<0 or c*c>a*b+1e-10*max(1.,a*b): raise ValueError('Invalid Gram matrix')
    p = math.sqrt(level) if kind=='split' and k>1 else level
    q = chi_quantile(1,p)
    if kind=='omnibus' or a*b==c*c:
        threshold=chi_quantile(k,level) if kind=='omnibus' else q
        roots=positive_roots([a-threshold,-2*c,b-threshold])
        if kind!='omnibus' and a>0 and c<0: roots.append(-c/a)
    else:
        roots=positive_roots([c*c-q*b,2*c*(a-b-q),(a-b)**2-2*c*c-q*(a+b),
                              -2*c*(a-b+q),c*c-q*a])
        if kind=='split' and k>1:
            qw=chi_quantile(k-1,p)
            determinant=max(0.,a*b-c*c)
            roots+=positive_roots([determinant-qw*b,-2*qw*c,determinant-qw*a])
        # An exactly antiparallel stack has an S=0 candidate with a declared convention.
        if a>0 and c<0 and a*b==c*c: roots.append(-c/a)
    root_map={math.atan(r):r for r in roots}
    angles=sorted(set([0.,math.pi/2]+list(root_map)))
    accepted=[]
    for left,right in zip(angles[:-1],angles[1:]):
        ratio=math.tan((left+right)/2)
        if bool(coverage_flags(a,b,c,k,ratio,level)[kind]): accepted.append([left,right])
    # Roots of a nonpositive polynomial may be isolated, notably S=0 conventions.
    for theta in angles[1:-1]:
        r=root_map[theta]
        t,z,w,s=statistics(a,b,c,r)
        tolerance=1e-8*max(1.,a+b)
        if kind=='omnibus': okay=t<=chi_quantile(k,level)+tolerance
        elif kind=='score': okay=z<=q+tolerance or s<=tolerance**2
        else: okay=(z<=q+tolerance and (k==1 or w<=chi_quantile(k-1,p)+tolerance)) or s<=tolerance**2
        if okay and not any(lo<=theta<=hi for lo,hi in accepted): accepted.append([theta,theta])
    merged=[]
    for interval in sorted(accepted):
        if merged and interval[0]<=merged[-1][1]+1e-12: merged[-1][1]=max(merged[-1][1],interval[1])
        else: merged.append(interval)
    return [[0. if lo==0 else root_map.get(lo,math.tan(lo)),
             math.inf if hi==math.pi/2 else root_map.get(hi,math.tan(hi))]
            for lo,hi in merged]


def contains(intervals,r):
    return any(lo<=r<=hi for lo,hi in intervals)


def checks():
    assert abs(chi_quantile(3,.95)-7.814727903251179)<2e-12
    for df in (1,2,3,24,192,1536):
        for probability in (.5,.95,math.sqrt(.95)):
            q=chi_quantile(df,probability)
            assert abs(gamma_p(df/2,q/2)-probability)<2e-12
    for k in (1,3,24,192,1536):
        exact=math.sqrt(2/math.pi) if k%2 else math.sqrt(math.pi/2)
        for j in range(1 if k%2 else 2,k,2): exact *= (j+1)/j
        assert abs(float(h_chi(k,0.))/exact-1)<2e-12
        low=h_chi(k,np.array([0.,1.,100.,4096.]),96)
        high=h_chi(k,np.array([0.,1.,100.,4096.]),192)
        assert np.max(abs(low/high-1))<2e-10
        point,cdf=posterior(0.,0.,0.,k,4.)
        assert abs(point-1)<1e-12 and abs(cdf-2*math.atan(4)/math.pi)<1e-12
    from gaussian_conditioning import summary
    for force,acceleration in (([4,0,0],[1,1,0]),([8,0,0],[4,0,0])):
        x,y=np.array(force,float),np.array(acceleration,float)
        point,_=posterior(x@x,y@y,x@y,3,2.)
        assert abs(point/summary(x,y,'tube',256)['mass_over_s']-1)<1e-11
    rng=np.random.default_rng(20260916)
    for a,b,c,r in ((7.,28.,-14.,2.),(28.,7.,-14.,.5),(2.,8.,-4.,2.)):
        assert statistics(a,b,c,r)[3]==0
        for kind in KINDS:
            assert contains(confidence_set(a,b,c,3,kind),r)==bool(coverage_flags(a,b,c,3,r)[kind])
    for k in (1,3,24):
        for _ in range(30):
            x,y=rng.normal(size=(2,k))
            a,b,c=float(x@x),float(y@y),float(x@y)
            for kind in KINDS:
                intervals=confidence_set(a,b,c,k,kind)
                for r in np.geomspace(1e-5,1e5,101):
                    assert contains(intervals,r)==bool(coverage_flags(a,b,c,k,r)[kind]),(k,kind,intervals,r)
    examples={}
    for name,a,b,c,k in (('observed_zero',0.,0.,0.,3),('strong_aligned',64.,16.,32.,3),
                         ('exact_anti_aligned_singular',100.,100.,-100.,3),
                         ('strong_opposed',100.,100.,-90.,3),
                         ('orthogonal',100.,100.,0.,3),('many_trials_mean_gram',3584.,2048.,1024.,1536)):
        examples[name]={kind:confidence_set(a,b,c,k,kind) for kind in KINDS}
    return {'confidence_set_examples':examples}


def run(replicates=2048,seed=20260916,null_study=False):
    result=checks()
    rng=np.random.default_rng(seed)
    rows=[]
    scenarios=([(n,2.,0.) for n in (1,8,64,512,4096)] if null_study else
               [(1,16.,.25),(1,2.,1.),(8,2.,1.),(64,2.,1.),(512,2.,1.)])
    for trials,r,t in scenarios:
        k=3*trials
        flags={kind:[] for kind in KINDS}
        credible,points=[],[]
        max_cdf_difference=0.
        for start in range(0,replicates,32):
            count=min(32,replicates-start)
            x,y=rng.normal(size=(2,count,trials,3))
            # Distinct fixed excitation directions, not supplied to the estimator.
            for i in range(trials):
                x[:,i,i%3]+=r*t
                y[:,i,i%3]+=t
            a,b,c=(x*x).sum((1,2)),(y*y).sum((1,2)),(x*y).sum((1,2))
            current=coverage_flags(a,b,c,k,r)
            for kind in KINDS: flags[kind].extend(current[kind].tolist())
            point,cdf=posterior(a,b,c,k,r)
            points.extend(point.tolist())
            credible.extend(((cdf>=.025)&(cdf<=.975)).tolist())
            if start==0:
                point_hi,cdf_hi=posterior(a,b,c,k,r,128,192)
                max_cdf_difference=float(np.max(abs(cdf-cdf_hi)))
                assert max_cdf_difference<1e-7
                assert np.max(abs(point/point_hi-1))<1e-7
        row={'trials':trials,'stack_dimension':k,'mass_over_s':r,'acceleration_snr_per_trial':t,
             'coverage':{kind:float(np.mean(flags[kind])) for kind in KINDS},
             'joint_jeffreys_credible_coverage':float(np.mean(credible)),
             'median_stacked_readout_over_truth':float(np.median(points)/r),
             'cdf_refinement_difference_first32':max_cdf_difference}
        rows.append(row)
        print(json.dumps(row),flush=True)
    result.update({'seed':seed,'replicates_per_scenario':replicates,
                   'study':'true_zero_excitation' if null_study else 'positive_excitation',
                   'calibration':rows})
    return result


def json_safe(value):
    if isinstance(value,dict): return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value,list): return [json_safe(v) for v in value]
    if isinstance(value,float) and math.isinf(value): return 'infinity' if value>0 else '-infinity'
    return value


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--null-study',action='store_true',help='Repeat trials with exactly zero true excitation')
    parser.add_argument('--replicates',type=int,default=2048)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.replicates<64: parser.error('Use at least 64 replicates')
    result=checks() if args.checks_only else run(args.replicates,null_study=args.null_study)
    result['source_commit']=SOURCE_COMMIT
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print('All assertions passed.',flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(json_safe(result),indent=2,allow_nan=False)+'\n',encoding='utf-8')

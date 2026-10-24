"""Independent scalar/vector crossover identity checks; scratch only."""
from math import exp, erfc, erf, pi, sqrt
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator'))
from estimator import positive_normal_mean, estimate, isotropic_covariance

def phi(t): return exp(-t*t/2)/sqrt(2*pi)
def Phi(t): return .5*erfc(-t/sqrt(2))
def correction(t): return phi(t)/(t*Phi(t))
def H(t): return positive_normal_mean(t,1.)
def bisect(fn,lo,hi):
    sign=fn(lo)
    assert sign*fn(hi)<0
    for _ in range(70):
        mid=(lo+hi)/2
        if fn(mid)*sign>0: lo=mid
        else: hi=mid
    return (lo+hi)/2

def R3(t,n=128):
    if t==0: return sqrt(2/pi)
    x,w=np.polynomial.legendre.leggauss(n)
    s=t*(x+1)/2
    D=float((w @ np.exp((s*s-t*t)/2))*t/2)
    return erf(t/sqrt(2))/D

def main():
    c=sqrt(2/pi)
    k=(1-c*c)/c
    tol=[dict(tolerance=tol,threshold=bisect(lambda t:correction(t)-tol,.0001,10))
         for tol in [.1,.05,.01,.001,.0001]]
    sampled=[]
    for t in [0,.25,.5,1,2,3,4,5,10]:
        sampled.append(dict(snr=t,H=H(t),R3=R3(t),
            relative_channel_correction=correction(t) if t else None))
    # Near-zero Taylor expansion and bounded log sensitivity.
    for t in [-1e-3,1e-3]:
        approx=c+(1-c*c)*t+.5*c*(2*c*c-1)*t*t
        assert abs(H(t)-approx)<1e-10
    for t in np.linspace(0.01,12,121):
        lam=phi(t)/Phi(t)
        hp=1-lam*H(t)
        assert 0<hp<=1
        assert 0<t*hp/H(t)<=1
    # Exact scalar shrinkage bracket and common-SNR cancellation.
    for tf in [.1,.3,1,2,5,10]:
        for ta in [.1,.3,1,2,5,10]:
            m=H(tf)/H(ta)
            assert min(1,tf/ta)-1e-13<=m<=max(1,tf/ta)+1e-13
            assert abs(m/(tf/ta)-1)<=correction(min(tf,ta))+1e-12
    vector=[]
    for t in [0,1,3,5]:
        res=estimate([0,0,0],[t,0,0],isotropic_covariance(3,1,1),
                     direction_order=48,ratio_order=768)
        target=c/R3(t)
        vector.append(dict(snr=t,full_mass=res.mass,closed_form_mass=target,
                           abs_error=abs(res.mass-target)))
        assert abs(res.mass-target)<2e-7
    # An incompatible antiparallel pair selects a boundary as noise shrinks.
    contrad=[]
    for sig in [.5,.25,.125]:
        res=estimate([2,0],[-1,0],isotropic_covariance(2,sig,sig),
                     direction_order=320,ratio_order=1024)
        contrad.append(dict(sigma=sig,mass=res.mass,leading=2/sig**2,
                             ratio_to_leading=res.mass*sig**2/2))
    result=dict(c=c,zero_log_sensitivity=k,thresholds=tol,curves=sampled,
                full_vector_checks=vector,contradiction=contrad)
    print(json.dumps(result,indent=2))
    Path(__file__).with_name('math_results.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__': main()

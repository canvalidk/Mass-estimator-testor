"""Deterministic one-dimensional risks in the force-SNR -> infinity limit."""
from functools import lru_cache
import json
import math
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator'))
from estimator import positive_normal_mean

def phi(x):
    return math.exp(-x*x/2)/math.sqrt(2*math.pi)
def Phi(x):
    return .5*math.erfc(-x/math.sqrt(2))
def H(x):
    return positive_normal_mean(x,1.)
def inverse_H(v):
    lo,hi=-max(10.,2/v),max(10.,2*v)
    for _ in range(65):
        mid=(lo+hi)/2
        if H(mid)<v:lo=mid
        else:hi=mid
    return (lo+hi)/2
def norm_cdf(r,a):
    return Phi(r-a)+Phi(r+a)-1+(phi(r+a)-phi(r-a))/a
def success_gap(a):
    ours=Phi(inverse_H(2*a)-a)-Phi(inverse_H(a/2)-a)
    norm=norm_cdf(2*a,a)-norm_cdf(a/2,a)
    return norm-ours
@lru_cache(None)
def nodes(n):
    return np.polynomial.legendre.leggauss(n)
def risk_gap(a,n=256):
    x,w=nodes(n)
    z=a+12*x
    weights=12*w*np.exp(-.5*(z-a)**2)/math.sqrt(2*math.pi)
    ours=float(weights @ np.array([math.log(H(v)/a)**2 for v in z]))
    r=(a+12)*(x+1)/2
    pdf=r/(a*math.sqrt(2*math.pi))*np.exp(-.5*(r-a)**2)*(-np.expm1(-2*r*a))
    norm=float((a+12)/2*w @ (pdf*np.log(r/a)**2))
    return ours-norm
def root(fn,lo,hi):
    flo=fn(lo)
    assert flo*fn(hi)<0
    for _ in range(50):
        mid=(lo+hi)/2
        if fn(mid)*flo>0:lo=mid
        else:hi=mid
    return (lo+hi)/2
if __name__=='__main__':
    out=dict(factor2_lower=root(success_gap,.1,.3),factor2_upper=root(success_gap,1.5,3.),
             squared_log_upper_256=root(risk_gap,2.,3.),
             squared_log_upper_512=root(lambda a:risk_gap(a,512),2.,3.),
             squared_log_upper_1024=root(lambda a:risk_gap(a,1024),2.,3.))
    print(json.dumps(out,indent=2))
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')

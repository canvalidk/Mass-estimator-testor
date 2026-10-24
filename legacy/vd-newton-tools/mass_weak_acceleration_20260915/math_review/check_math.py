"""Independent exact-known-force, weak-acceleration derivation checks.

Internal scratch for parent integration; no claim of novelty.
"""
import json
from math import exp, log, sqrt, pi, erfc
from pathlib import Path
from statistics import NormalDist
import numpy as np

def ndtr(x):
    return .5*erfc(-x/sqrt(2))

def log_ndtr(x):
    if x > -10:
        return log(ndtr(x))
    term = total = 1.0
    for j in range(1,100):
        nxt = term * -(2*j-1)/(x*x)
        if abs(nxt) >= abs(term):
            break
        total += nxt
        term = nxt
        if abs(nxt)<1e-16*abs(total):
            break
    return -.5*x*x-.5*log(2*pi)-log(-x)+log(total)

class norm:
    @staticmethod
    def pdf(x): return exp(-x*x/2)/sqrt(2*pi)
    @staticmethod
    def logpdf(x): return -x*x/2-.5*log(2*pi)
    @staticmethod
    def sf(x): return ndtr(-x)
    @staticmethod
    def ppf(p): return NormalDist().inv_cdf(p)

def brentq(fn,lo,hi,xtol=1e-12):
    # Bisection is enough for these monotone equations.
    flo=fn(lo)
    fhi=fn(hi)
    assert flo*fhi<=0,(lo,hi,flo,fhi)
    for _ in range(150):
        mid=(lo+hi)/2
        fm=fn(mid)
        if fm==0 or hi-lo<=xtol:
            return mid
        if fm*flo>0:
            lo,flo=mid,fm
        else:
            hi=mid
    return (lo+hi)/2

def quad(fn,lo,hi):
    x,w=np.polynomial.legendre.leggauss(200)
    return (sum(float(ww)*fn(lo+(hi-lo)*(float(xx)+1)/2) for xx,ww in zip(x,w))*(hi-lo)/2,)

P = 0.95
DELTA = 1 - P

def H(t):
    return t + exp(norm.logpdf(t) - log_ndtr(t))

def alpha_q(t, p=P):
    ld = log(1-p)
    return brentq(lambda a: log_ndtr(t-a) - log_ndtr(t) - ld,
                  0, max(10, t+10), xtol=1e-13)

def sampling_cutoff(theta, p=P):
    ld = log(1-p)
    lo = -max(10, 10/theta)
    return brentq(lambda t: log_ndtr(t-theta)-log_ndtr(t)-ld,
                  lo, theta+10, xtol=1e-12)

def A(x):
    return x*ndtr(x)+norm.pdf(x)

def L_bin(theta, k):
    return ndtr(k-theta)-ndtr(-k-theta)

def G_bin(theta, k):
    return (A(k+theta)-A(k-theta)-theta)/k

def mean_bin(k):
    return ((k*k+1)*(ndtr(k)-.5)+k*norm.pdf(k))/k

def q_bin(k, p=P):
    return brentq(lambda a: G_bin(a,k)-p, 0, k+10, xtol=1e-13)

def hybrid_coverage(theta, k, p=P):
    # Reports only E within [-k,k], precise signed t outside.
    # Endpoint theta=q_bin is included in the confidence set.
    tb=sampling_cutoff(theta,p)
    qb=q_bin(k,p)
    low_success=max(0, ndtr(-k-theta)-ndtr(tb-theta))
    high_success=norm.sf(max(k,tb)-theta)
    bin_success=L_bin(theta,k) if theta<=qb else 0
    return low_success+high_success+bin_success

def main():
    exact=[]
    for theta in [.01,.1,.25,.5,1,2,3,5,10]:
        tc=sampling_cutoff(theta)
        cov=ndtr(theta-tc)
        ident=1-DELTA*ndtr(tc)
        assert abs(cov-ident)<2e-13
        exact.append(dict(theta=theta,cutoff=tc,coverage=cov))
    observations=[]
    for t in [-5,-3,-2,-1,0,.5,1,2,3,5]:
        q=alpha_q(t)
        assert abs(exp(log_ndtr(t-q)-log_ndtr(t))-DELTA)<1e-13
        observations.append(dict(t=t,mean_alpha_over_sigma=H(t),
                                 mass_point_scale=1/H(t),
                                 q95_alpha_over_sigma=q,
                                 mass_lower_scale=1/q))
    bins=[]
    for k in [.1,.5,1,2,3,5]:
        mb=mean_bin(k)
        qm=quad(lambda a:a*L_bin(a,k)/k,0,k+12)[0]
        assert abs(mb-qm)<1e-10
        qb=q_bin(k)
        qconf=brentq(lambda a:L_bin(a,k)-DELTA,0,k+10)
        right=ndtr(qb-k)
        eps=1e-7
        assert abs(hybrid_coverage(qb+eps,k)-ndtr(qb+eps-k))<1e-12
        bins.append(dict(k=k,mean_alpha_over_sigma=mb,mass_point_scale=1/mb,
                         q95_alpha_over_sigma=qb,mass_lower_scale=1/qb,
                         binary_confidence_upper=qconf,
                         binary_only_right_coverage=1-L_bin(qb,k),
                         hybrid_right_coverage=right,
                         hybrid_left_coverage=ndtr(qb+k)))
    results=dict(exact_y_zero=dict(mean=H(0),point=1/H(0),q95=alpha_q(0),
                    lower=1/alpha_q(0),gaussian_q95=norm.ppf(P)),
                 observed=observations,exact_coverage=exact,censored=bins)
    print(json.dumps(results,indent=2))
    Path(__file__).with_name('results.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':
    main()

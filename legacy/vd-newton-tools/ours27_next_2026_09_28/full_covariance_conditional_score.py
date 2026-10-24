"""Exact full-covariance conditional score: scratch uncertainty benchmark.

Not a mass-law replacement. VD-docs snapshot:
b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

At tested m, y=(F,a)~N(B_m a*, Sigma), any fixed aligned candidate a*.
r=F-m*a, V=Var(r). The GLS candidate mean is
b=a-(Sigma_aF-m*Sigma_aa)V^-1 r.
Then b~N(a*,K), r~N(0,V), and Cov(b,r)=0; they are independent.
The flat log-likelihood score is b.T V^-1 r. Conditional on b this
is exactly N(0,b.T V^-1 b). Thus Z below is exactly standard normal
whenever b is nonzero, without any candidate-size prior. The zero set
has probability zero under positive-definite measurement covariance.
For the literal all-zero observed vector the benchmark is undefined;
that point is assigned no directional evidence rather than a mass law.

Independent readings can be stacked with full block covariance. This
also allows arbitrary correlations in the stacked data if supplied.
This is an exact confidence-test pivot, NOT a normalized law of mass:
Z(m) need not be monotone and can vanish at the wrong stationary branch.
Ratios of conditional p-values are not posterior odds.
"""
from pathlib import Path
import json
import numpy as np


def terms(y, sigma, m):
    d=len(sigma)//2
    r=y[..., :d]-m*y[..., d:]
    c=sigma[d:,:d]-m*sigma[d:,d:]
    v=sigma[:d,:d]-m*(sigma[:d,d:]+sigma[d:,:d])+m*m*sigma[d:,d:]
    vi=np.linalg.inv(v)
    b=y[...,d:]-r@vi@c.T
    bv=b@vi
    score=np.sum(bv*r,axis=-1)
    variance=np.sum(bv*b,axis=-1)
    return score/np.sqrt(variance),b,r,v


def main():
    rng=np.random.default_rng(202609282)
    d=4
    c=rng.normal(size=(2*d,2*d))
    sigma=c@c.T/4+.25*np.eye(2*d)
    m=2.3
    candidates=[np.zeros(d),np.array([1.,-.3,2.,.2]),np.array([10.,-3.,20.,2.])]
    rows=[]
    errors={}
    chol=np.linalg.cholesky(sigma)
    for a in candidates:
        y=np.r_[m*a,a]+rng.normal(size=(250000,2*d))@chol.T
        z,b,r,v=terms(y,sigma,m)
        rows.append({'candidate':a.tolist(),'mean':float(z.mean()),'variance':float(z.var()),'coverage_95':float(np.mean(abs(z)<=1.959963984540054)),'quantiles':np.quantile(z,[.025,.5,.975]).tolist()})
        assert abs(z.mean())<.01
        assert abs(z.var()-1)<.02
        assert abs(np.mean(abs(z)<=1.959963984540054)-.95)<.003
    y=np.r_[2.,-1.,.5,4., .3,.2,-.6,.1]
    z,b,r,v=terms(y,sigma,m)
    rot,_=np.linalg.qr(rng.normal(size=(d,d)))
    full=np.zeros((2*d,2*d));full[:d,:d]=rot;full[d:,d:]=rot
    zr=terms(full@y,full@sigma@full.T,m)[0]
    errors['rotation']=float(abs(z-zr))
    scaling=np.diag(np.r_[np.full(d,3.),np.full(d,.7)])
    zu=terms(scaling@y,scaling@sigma@scaling.T,m*3/.7)[0]
    errors['units']=float(abs(z-zu))
    swap=np.block([[np.zeros((d,d)),np.eye(d)],[np.eye(d),np.zeros((d,d))]])
    zs=terms(swap@y,swap@sigma@swap.T,1/m)[0]
    errors['swap_sign']=float(abs(z+zs))
    eps=1e-5
    def loglik(x):
        _,_,rr,vv=terms(y,sigma,x)
        return -.5*rr@np.linalg.solve(vv,rr)
    derivative=(loglik(m+eps)-loglik(m-eps))/(2*eps)
    score=b@np.linalg.solve(v,r)
    errors['score_derivative']=float(abs(score-derivative))
    assert max(errors.values())<1e-8,errors
    # A literally zero independent axis adds exactly zero score and variance.
    base=np.array([2.,.3]);zb=terms(base,np.eye(2),1.4)[0]
    null_errors=[]
    for sd in [1.,.1,1e-3]:
        zz=terms(np.array([2.,0.,.3,0.]),np.diag([1.,sd*sd,1.,sd*sd]),1.4)[0]
        null_errors.append(float(abs(zz-zb)))
    assert max(null_errors)<1e-12
    out={'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8','trials_per_cell':250000,'conditional_score_cells':rows,'invariance_max_errors':errors,'null_axis_max_error':max(null_errors),'scope':'Exact test pivot only; not a joint candidate law or final mass estimator.'}
    Path(__file__).with_name('full_covariance_conditional_score_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()

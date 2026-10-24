"""A continuous instrument-only scale coupling for diagonal noise pairs.

For centered Gaussian noise densities n_g on the same F,a channel space:
  rho_gh = integral sqrt(n_g n_h)
         = |S_g|^(1/4)|S_h|^(1/4) / |(S_g+S_h)/2|^(1/2).
It is a positive semidefinite Gram kernel, invariant under a common invertible
channel recoding (including units and F/a swap). Draw Z~N(0,rho), then
  beta_g = sin^2[(pi/2) Phi(Z_g)].
Each beta has the arcsine prior. Equal noise laws share the same beta; widely
separated absolute precisions decorrelate continuously.

Learn this joint state across readings:
 p(m,Z|D) proportional pi(m) phi_rho(Z)
   prod_g (1-beta_g)^(M_g/2) exp(beta_g H_g(m)/2).

This is a constructive candidate, not a uniqueness claim. Choosing the
Gaussian copula built from square-root noise overlap is an explicit extra
principle. It removes an arbitrary grouping threshold and a correlation
length parameter. General full-covariance block selection is not solved here.
"""
from pathlib import Path
from functools import lru_cache
import json
import math
import numpy as np


def overlap(a,b):
    sa,la=np.linalg.slogdet(a)
    sb,lb=np.linalg.slogdet(b)
    sm,lm=np.linalg.slogdet((a+b)/2)
    assert sa==sb==sm==1
    return math.exp(.25*(la+lb)-.5*lm)


def cdf(z):
    z=np.asarray(z)
    return .5*(1+np.fromiter(
        (math.erf(float(t)/math.sqrt(2)) for t in z.flat),float).reshape(z.shape))


@lru_cache(maxsize=None)
def gaussian_rule(order):
    x,w=np.polynomial.hermite.hermgauss(order)
    return np.sqrt(2)*x,w/np.sqrt(np.pi)


def density(epsilon, order=160, grid=1001):
    theta=np.linspace(0,np.pi/2,grid)
    H=(3*np.sin(theta)+np.cos(theta))**2
    rho=2*epsilon/(1+epsilon**2) if epsilon>0 else 0.
    x,w=gaussian_rule(order)
    z1=x[:,None]+np.zeros((order,order))
    z2=rho*x[:,None]+np.sqrt(max(0,1-rho*rho))*x[None,:]
    b1=np.sin(np.pi/2*cdf(z1))**2
    b2=np.sin(np.pi/2*cdf(z2))**2
    weight=w[:,None]*w[None,:]*np.sqrt(1-b1)*(1-b2)
    # Active axis has M1=1; the two zero axes share their own noise pair, M2=2.
    law=np.array([np.sum(weight*np.exp(.5*h*b1)) for h in H])
    law/=np.trapezoid(law,theta)
    c=np.r_[0,np.cumsum((law[:-1]+law[1:])*np.diff(theta)/2)]
    qs=np.tan(np.interp([.025,.5,.975],c,theta))
    return theta,law,{"epsilon":epsilon,"rho":rho,
                     "mass_quantiles":qs.tolist(),
                     "mass_cdf_at_one":float(np.interp(np.pi/4,theta,c))}


def exact_one_scale(d,order=400,grid=1001):
    z,w=np.polynomial.legendre.leggauss(order)
    phi=(z+1)*np.pi/4
    beta=np.sin(phi)**2
    weight=w*np.pi/4*np.cos(phi)**d
    theta=np.linspace(0,np.pi/2,grid)
    H=(3*np.sin(theta)+np.cos(theta))**2
    law=np.exp(.5*H[:,None]*beta)@weight
    return law/np.trapezoid(law,theta)


def main():
    rows=[]
    _,active,base=density(0)
    _,isotropic,iso=density(1)
    assert np.max(abs(active-exact_one_scale(1)))<2e-8
    assert np.max(abs(isotropic-exact_one_scale(3)))<2e-8
    for e in (1, .999, .9,.5,.1,.01,.001,0):
        _,law,row=density(e)
        row["max_density_difference_from_1D"]=float(np.max(abs(law-active)))
        rows.append(row)
    _,law160,_=density(.5,160)
    _,law200,_=density(.5,200)
    refinement=float(np.max(abs(law160-law200)))
    assert refinement<2e-8
    rng=np.random.default_rng(20260929)
    s=[]
    for i in range(9):
        a=rng.normal(size=(2,2))
        s.append(a@a.T+.1*np.eye(2))
    gram=np.array([[overlap(a,b) for b in s] for a in s])
    assert np.linalg.eigvalsh(gram).min()>-1e-12
    transform=np.array([[0.,3.],[.2,0.]])
    changed=[transform@a@transform.T for a in s]
    transformed=np.array([[overlap(a,b) for b in changed] for a in changed])
    err=float(np.max(abs(gram-transformed)))
    assert err<1e-12
    out={"source_commit":"592d144531fa4a30d1c2229a4c3b528e84043d77",
         "null_axis_continuation":rows,
         "refinement_160_to_200":refinement,
         "one_dimensional_endpoint_error":float(np.max(abs(active-exact_one_scale(1)))),
         "three_dimensional_endpoint_error":float(np.max(abs(isotropic-exact_one_scale(3)))),
         "gram_min_eigenvalue":float(np.linalg.eigvalsh(gram).min()),
         "units_swap_error":err,
         "status":"Proper smooth coupling candidate for supplied diagonal/independent channel-pair groups; not uniquely forced or full-covariance selected."}
    Path(__file__).with_name("noise_overlap_candidate_results.json").write_text(
        json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()

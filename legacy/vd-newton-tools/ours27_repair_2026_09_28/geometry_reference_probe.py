"""Exploratory algebra checks, not a replacement estimator or a VD record.

Managed context: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
Extends the local September 27 projector/reference-model probes.
The mass-first scale measure below requires a NEW reference-model premise.
No simulation scores are used to choose a prior or fit a hyperparameter.
"""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss


def gl(a, b, n):
    x, w = leggauss(n)
    return (a+b)/2+(b-a)*x/2, (b-a)*w/2


def geometry(sigma, m):
    d = len(sigma)//2
    eig, u = np.linalg.eigh(sigma)
    whitening = (u / np.sqrt(eig)) @ u.T
    b = np.vstack([m*np.eye(d), np.eye(d)])
    h = whitening @ b
    derivative = whitening @ np.vstack([np.eye(d), np.zeros((d, d))])
    k = np.linalg.inv(h.T @ h)
    p = h @ k @ h.T
    t = (np.eye(2*d)-p) @ derivative @ k @ h.T
    dp = t+t.T
    v = sigma[:d,:d]-m*(sigma[:d,d:]+sigma[d:,:d])+m*m*sigma[d:,d:]
    g = float(np.trace(np.linalg.solve(v,k)))
    return g, float(.5*np.trace(dp @ dp)), p, dp, k, v


def reference_law(f, a, theta, scale_order=180):
    """One pair, s=1; prior uniform theta, size reference d beta/(1-beta).

    Law marginal = integral_0^1 (1-beta)^(M/2-1) exp(beta H2/2) d beta.
    Substitute t=(1-beta)^(M/2), making base quadrature uniform in t.
    """
    M = len(f)
    h2 = np.sum((f[:,None]*np.sin(theta)+a[:,None]*np.cos(theta))**2,axis=0)
    t,w = gl(0,1,scale_order)
    beta = 1-t**(2/M)
    exponents = .5*beta[:,None]*h2[None,:]
    shift = float(np.max(exponents))
    integrated = w @ np.exp(exponents-shift)
    density = integrated/np.trapezoid(integrated, theta)
    return density, h2


def quantiles(density, theta):
    cdf=np.r_[0,np.cumsum((density[1:]+density[:-1])*np.diff(theta)/2)]
    return np.tan(np.interp([.025,.5,.975],cdf,theta)).tolist()


def main():
    rng=np.random.default_rng(280927)
    c=rng.normal(size=(6,6))
    sigma=c@c.T+np.eye(6)
    rot,_=np.linalg.qr(rng.normal(size=(3,3)))
    transform=np.zeros((6,6));transform[:3,:3]=rot;transform[3:,3:]=rot
    swap=np.block([[np.zeros((3,3)),np.eye(3)],[np.eye(3),np.zeros((3,3))]])
    scaling=np.diag([3.,3.,3.,.7,.7,.7])
    maxima={k:0. for k in ['trace_identity','rotation','reciprocal','units','derivative','diagonal_reduction']}
    for m in np.geomspace(.05,30,25):
        g,gp,p,dp,k,v=geometry(sigma,m)
        maxima['trace_identity']=max(maxima['trace_identity'],abs(g/gp-1))
        maxima['rotation']=max(maxima['rotation'],abs(geometry(transform@sigma@transform.T,m)[0]/g-1))
        maxima['reciprocal']=max(maxima['reciprocal'],abs(geometry(swap@sigma@swap.T,1/m)[0]/(m**4*g)-1))
        maxima['units']=max(maxima['units'],abs(geometry(scaling@sigma@scaling.T,3/.7*m)[0]/g*(3/.7)**2-1))
        eps=1e-5*m
        fd=(geometry(sigma,m+eps)[2]-geometry(sigma,m-eps)[2])/(2*eps)
        maxima['derivative']=max(maxima['derivative'],float(np.max(np.abs(fd-dp))))
        ratios=np.array([5.,2.,.8])
        sd=np.diag(np.r_[np.ones(3),1/ratios**2])
        expected=np.sum((ratios/(ratios**2+m*m))**2)
        maxima['diagonal_reduction']=max(maxima['diagonal_reduction'],abs(geometry(sd,m)[0]/expected-1))
    assert max(maxima.values()) < 2e-8,maxima

    # Correlated one-axis normalization has a closed form, independent of data.
    one_axis=[]
    theta,w=gl(0,np.pi/2,1000)
    for rho in [-.9,-.4,0.,.6,.95]:
        sf,sa=2.,.7
        s=sf/sa
        m=s*np.tan(theta)
        sigma1=np.array([[sf*sf,rho*sf*sa],[rho*sf*sa,sa*sa]])
        density=np.sqrt(np.linalg.det(sigma1))/(sf*sf-2*rho*sf*sa*m+sa*sa*m*m)
        integral=float(w@(density*s/(np.cos(theta)**2)))
        exact=float(np.arccos(-rho))
        assert abs(integral-exact)<1e-10
        one_axis.append({'rho':rho,'arc_length':integral,'exact':exact})

    # Appending a null learns shared size even though all-null datasets leave
    # the mass law flat. This rejects a candidate, not a defect in flat ours27.
    theta=np.linspace(0,np.pi/2,4001)
    cases=[]
    for zeros in [0,1,5,50]:
        f=np.r_[3.,np.zeros(zeros)];a=np.r_[1.,np.zeros(zeros)]
        d,h2=reference_law(f,a,theta)
        dfine,_=reference_law(f,a,theta,scale_order=360)
        qerr=float(np.max(np.abs(d-dfine)))
        assert qerr<2e-5,qerr
        flat=np.exp(.5*(h2-np.max(h2)));flat/=np.trapezoid(flat,theta)
        cases.append({'appended_nulls':zeros,'reference_95_median':quantiles(dfine,theta),'flat_95_median':quantiles(flat,theta),'refinement_max_density_difference':qerr})
    null_errors=[]
    for M in [1,2,6,63]:
        d,_=reference_law(np.zeros(M),np.zeros(M),theta)
        null_errors.append(float(np.max(np.abs(d-2/np.pi))))
    assert max(null_errors)<1e-12
    out={'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8','geometry_max_errors':maxima,'correlated_1d_arc':one_axis,'mass_first_reference_candidate_null_counterexample':cases,'all_null_density_max_error':max(null_errors),'status':'Geometry verified; reference-size candidate is not adopted and fails marginal appended-null invariance.'}
    Path(__file__).with_name('geometry_reference_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()

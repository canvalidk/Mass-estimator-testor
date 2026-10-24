"""A diagnostic candidate outside ordinary candidate-prior integration.

NOT an adopted mass equation. No claim of a data-independent joint prior.
One common isotropic instrument pair; theta=atan(m/s), standardized x,y.

At the unrestricted projective TLS direction let h>=q be the Gram eigenvalues,
Delta=h-q. The exact conditional score C/sqrt(H) has derivative -Delta/sqrt(h)
at that direction, giving observed local information Delta^2/h. Flat profile
curvature is Delta. Matching only these local curvatures selects beta=Delta/h.

Candidate mass law: uniform-angle prior times exp(-beta*Q(theta)/2), beta=0
if all data are zero. To make a full candidate law one could keep ours27's
conditional N(A_theta,I), then use its ratio-of-mean-sizes readout; that is
an additional assignment, not derived here. This is NOT another size prior:
the data-dependent exponent cannot be obtained by ordinary Gaussian updating
from a fixed candidate measure while retaining null-law neutrality.

Passes same-coordinate rotations, regrouping, units, swap, literal measured
null-law neutrality and fixed-candidate regular width calibration. It is NOT
an exact confidence distribution, and local curvature does not select its
global shape. The global law fails multiplication of separately reduced
batch laws; raw Gram-state combination still works. Correct instrument noise
is retained; this diagnostic reweights the mass law rather than declaring a
new measurement covariance.

Source context: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def summary(x, y):
    p=np.sum(x*x,axis=-1); q=np.sum(y*y,axis=-1); d=np.sum(x*y,axis=-1)
    gap=np.hypot(p-q,2*d); h=(p+q+gap)/2
    beta=np.divide(gap,h,out=np.zeros_like(h),where=h>0)
    return p,q,d,gap,h,beta


def log_kernel(s,theta,adjusted=True):
    p,q,d,gap,h,beta=s
    if not adjusted: beta=np.ones_like(beta)
    theta=np.asarray(theta)
    H=p[...,None]*np.sin(theta)**2+q[...,None]*np.cos(theta)**2+d[...,None]*np.sin(2*theta)
    return .5*beta[...,None]*(H-h[...,None])


def cdf_at(s,theta0,adjusted,nodes=96):
    t,w=leggauss(nodes)
    lo=(t+1)*theta0/2
    hi=theta0+(t+1)*(np.pi/2-theta0)/2
    a=np.exp(log_kernel(s,lo,adjusted))@w*theta0/2
    b=np.exp(log_kernel(s,hi,adjusted))@w*(np.pi/2-theta0)/2
    return a/(a+b)


def main():
    rng=np.random.default_rng(28092831)
    n=4000
    cells=[]
    for M,kappa in [(1,1.),(6,1/3),(63,1/3),(63,1.),(400,1.),(400,3.)]:
        for theta0 in [np.pi/4,np.arctan(.25)]:
            v=np.full(M,np.sqrt(kappa))
            x=np.sin(theta0)*v+rng.normal(size=(n,M))
            y=np.cos(theta0)*v+rng.normal(size=(n,M))
            s=summary(x,y)
            base=cdf_at(s,theta0,False)
            changed=cdf_at(s,theta0,True)
            refined=cdf_at(tuple(a[:50] for a in s),theta0,True,160)
            error=float(np.max(abs(refined-changed[:50])))
            assert error<2e-9,error
            cells.append({'M':M,'R_squared_per_component':kappa,'mass_over_s':float(np.tan(theta0)),
                'trials':n,'flat_95_coverage':float(np.mean((base>=.025)&(base<=.975))),
                'local_information_95_coverage':float(np.mean((changed>=.025)&(changed<=.975))),
                'mean_beta':float(np.mean(s[-1])),'regular_beta_target':kappa/(1+kappa),
                'cdf_refinement_error':error})

    theta=np.linspace(0,np.pi/2,4001)
    x=np.array([3.,-1.]); y=np.array([1.,2.])
    batch=log_kernel(summary(x,y),theta)
    separate=log_kernel(summary(x[:1],y[:1]),theta)+log_kernel(summary(x[1:],y[1:]),theta)
    # A theta-independent difference could be removed by normalization; this cannot.
    discrepancy=float(np.ptp(batch-separate))
    assert discrepancy>.1
    null=log_kernel(summary(np.r_[x,0.],np.r_[y,0.]),theta)
    assert np.max(abs(null-batch))<1e-14
    swap=log_kernel(summary(y,x),np.pi/2-theta)
    assert np.max(abs(swap-batch))<1e-14
    rot,_=np.linalg.qr(rng.normal(size=(2,2)))
    rotated=log_kernel(summary(rot@x,rot@y),theta)
    assert np.max(abs(rotated-batch))<1e-13
    zeros=log_kernel(summary(np.zeros(3),np.zeros(3)),theta)
    assert np.max(abs(zeros))==0
    # Curvature matching and regular calibration do not determine global shape.
    # L=Q-q>=0 vanishes quadratically at the unrestricted TLS direction.
    # Adding -c L^2/(1+T) changes neither curvature nor the regular local
    # limiting normal law, for any fixed c>=0. It preserves the same 1-pair
    # rotations, units, channel swap, and literal null-law neutrality.
    # This demonstrates freedom; c is NOT tuned or recommended.
    tt,ww=leggauss(500); tt=(tt+1)*np.pi/4;ww*=np.pi/4
    ambiguity=[]
    for f,a in [(1.,0.),(2.,0.),(2.,1.),(-2.,1.)]:
        st=summary(np.array([f]),np.array([a]))
        base_log=log_kernel(st,tt)
        Avec=f*np.sin(tt)+a*np.cos(tt)
        mu=np.sqrt(2/np.pi)*np.exp(-Avec*Avec/2)+Avec*np.array([math.erf(v/np.sqrt(2)) for v in Avec])
        L=-2*base_log/st[-1] # M=1 has beta=1 at nonzero data
        T=f*f+a*a
        for c in [0.,1.,4.]:
            lk=base_log-c*L*L/(1+T)
            measure=ww*np.exp(lk-np.max(lk))
            point=float((measure@(np.sin(tt)*mu))/(measure@(np.cos(tt)*mu)))
            ambiguity.append({'force':f,'acceleration':a,'c':c,'ratio_of_mean_sizes':point})
    out={'status':'Local-information diagnostic, not derived full candidate law or adopted equation.',
        'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8','cells':cells,
        'batch_vs_separate_log_kernel_difference_range':discrepancy,
        'global_shape_counterfamily':ambiguity,
        'null_rotation_swap_zero_checks':'passed'}
    Path(__file__).with_name('local_information_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()

"""Audit of curvature-matched mass laws, not an adopted definition.

ONE PAIR: LET h>=q BE DATA GRAM EIGENVALUES, Delta=h-q, I=Delta^2/h.
If u=(h-H(theta))/Delta=sin^2(theta-theta_TLS), then both
 ell_1=-I*u/2 and ell_2=-I*(u+u^2)/2
have the same unrestricted mode and local curvature -I. Both preserve measured
null neutrality, unit changes, channel swap, rotations and coordinate regrouping.
Both are proper on the finite angle interval, and their regular concentrating
limits agree since the extra term is I*O(delta_theta^4)=O(1/I).
At Delta=0 define both kernels as constant, as H itself is then constant.
Their finite-data laws differ. More generally, every increasing g with g(0)=0,
g'(0)=1 gives this same local information. Thus scalar profile power is uniquely
fixed only AFTER selecting the power family, not as the entire probability law.

FULL COVARIANCE LOCAL EXTENSION
Let ell=.5 sum z_i'P_i(m)z_i, s=ell', t_i=P_i z_i.
Conditional on candidate-m tangent projections, the exact Gaussian score has
variance Vc=sum ||P_i'(m)t_i||^2=sum b_i' V_i^-1 b_i. At a stationary mhat,
the normalized score Z=s/sqrt(Vc) has derivative -j/sqrt(Vc), j=-ell''.
Hence matching local curvature of a scalar profile power suggests beta=j/Vc.
At stationary points this ratio is invariant under reparameterization of mass.
It is still only a local diagnostic; multiple/degenerate extrema need care.

EXACT DEGENERATE COUNTEREXAMPLE TO THE FULL-COV SCALAR EXTENSION
Two independent 1D readings with ratios s1=sqrt(2)+1, s2=sqrt(2)-1, and
standardized readings x_i=y_i=k. For u=log m and a=log(sqrt(2)+1),

 ell(u)=k^2[1+.5 sech(u-a)+.5 sech(u+a)]
       =k^2[1+sqrt(2) cosh(u)/(1+cosh(u)^2)].

Its unique unrestricted maximum is u=0 (negative m gives lower values).
ell''(0)=0, but ell''''(0)=-3k^2/sqrt(2)<0, while
Vc(0)=k^2*(1+1/sqrt(2))/4>0. Thus beta=j/Vc=0 makes the scalar-power proposal
throw away every reading for every k, even though the nonconstant profile
concentrates quartically. This does NOT criticize the original one-pair formula,
where Delta=0 implies a genuinely constant profile. It identifies an obstacle
to extending that formula solely by local second-order information.

A POSSIBLE EXPLICIT SHAPE AXIOM (NOT YET PHYSICALLY DERIVED)
Require log(p/pi) to contain only constant, sin(2theta) and cos(2theta)
harmonics, with preferred direction equal to TLS. This forces the Bingham
angular family; conditional local information then fixes its concentration.
Gaussian measurement errors alone do not imply this shape after candidate
integration. Maximum entropy or minimum change subject merely to a derivative
constraint at one point does not determine the global density.
Attaching the original Gaussian candidate conditional to the modified mass
law defines a positive mathematical joint law, but does not establish that
it arose by ordinary Bayesian updating from a data-independent candidate prior.

LARGER EXACT DOMAIN FOR THE ONE-PAIR CONSTRUCTION
If the full measurement error covariance is separable C_channel tensor
S_component, whiten channel and component spaces separately. Every candidate
still has one common whitened channel direction u(m) times an M-vector, so
the profile remains a 2x2 Gram quadratic. Thus beta=Delta/h remains in [0,1]
and Delta=0 means the entire profile is constant. This allows spatial or
reading covariance and common channel correlation. Positive mass maps to
an instrument-defined angular interval, not necessarily length pi/2.
If the starting law is taken uniform in geometric angle, it is uniform on
that interval. Different channel-noise ratios generally destroy separability,
which is the situation exploited by the quartic counterexample above.
"""

import json
import numpy as np

def profile(theta,s,x,y):
    S=np.sin(theta);C=np.cos(theta)
    den=np.sqrt((S/s)**2+C*C)
    su=S/s/den;cu=C/den
    t=x*su+y*cu;n=x*cu-y*su
    d=s/(s*s*C*C+S*S)
    dd=-2*s*(1-s*s)*S*C/(s*s*C*C+S*S)**2
    ell=.5*np.sum(t*t,axis=-1)
    score=np.sum(d*t*n,axis=-1)
    j=np.sum(d*d*(t*t-n*n)-dd*t*n,axis=-1)
    vc=np.sum(d*d*t*t,axis=-1)
    return ell,score,j,vc

def maximum(s,x,y):
    grid=np.arange(1024)*np.pi/1024
    ell=profile(grid[:,None],s,x,y)[0]
    i=int(np.argmax(ell));step=np.pi/1024
    lo=grid[i]-step;hi=grid[i]+step
    ratio=(np.sqrt(5)-1)/2
    a=hi-ratio*(hi-lo);b=lo+ratio*(hi-lo)
    fa=profile(a,s,x,y)[0];fb=profile(b,s,x,y)[0]
    for _ in range(65):
        if fa<fb:
            lo=a;a=b;fa=fb;b=lo+ratio*(hi-lo);fb=profile(b,s,x,y)[0]
        else:
            hi=b;b=a;fb=fa;a=hi-ratio*(hi-lo);fa=profile(a,s,x,y)[0]
    theta=(lo+hi)/2
    for _ in range(3):
        _,score,j,_=profile(theta,s,x,y)
        theta+=score/j
    vals=profile(theta,s,x,y)
    return theta,vals

def random_counterexample():
    rng=np.random.default_rng(28092845)
    for trial in range(4000):
        s=np.exp(rng.uniform(-4,4,size=3))
        x=rng.normal(size=3)*3;y=rng.normal(size=3)*3
        theta,(ell,score,j,vc)=maximum(s,x,y)
        if j/vc>1.05 and 0<theta<np.pi/2:
            return {'trial':trial,'s':s.tolist(),'x':x.tolist(),'y':y.tolist(),
                'theta_at_projective_max':float(theta),'score':float(score),
                'curvature':float(j),'conditional_score_variance':float(vc),
                'beta':float(j/vc)}
    return None

def shape_comparison():
    th=np.linspace(0,np.pi/2,40001)
    h,q,mode=8.,2.,.7
    I=(h-q)**2/h
    u=np.sin(th-mode)**2
    result=[]
    for power in [0,1]:
        w=np.exp(-.5*I*(u+power*u*u))
        c=np.r_[0,np.cumsum((w[:-1]+w[1:])/2*np.diff(th))]
        c/=c[-1]
        qs=np.interp([.025,.5,.975],c,th)
        result.append({'quartic_coefficient':power,'local_information':I,
          'mass_quantiles':np.tan(qs).tolist()})
    return result

if __name__=='__main__':
    print(json.dumps({'same_curvature_different_global_laws':shape_comparison(),
      'full_covariance_beta_greater_than1':random_counterexample(),
      'quartic_maximum':{'s':[float(np.sqrt(2)+1),float(np.sqrt(2)-1)],
        'ell_second_derivative':0.,'ell_fourth_derivative_at_k1':float(-3/np.sqrt(2)),
        'conditional_score_variance_at_k1':float((1+1/np.sqrt(2))/4)}},indent=2))

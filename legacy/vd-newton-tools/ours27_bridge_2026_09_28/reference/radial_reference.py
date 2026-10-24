"""Ordered reference measure in the invariant rank-one Gaussian experiment.

SCOPE: one noise pair, M standardized component pairs, all sharing theta.
z_parallel=Z' u_theta, z_perp=Z' u_theta_perp, H=||z_parallel||, T=||Z||_F^2.
For a fixed candidate radius R and Haar direction on S^(M-1),

 L(theta,R) proportional exp[-(T+R^2)/2] A_M(R H),
 A_M(z)=0F1(M/2;z^2/4), b_M(z)=A_M'(z)/A_M(z).

The same is the likelihood of the maximal invariant scatter Z Z' under the
right O(M) action. Let H~noncentral-chi_M(R). Then

 s_theta = R b_M(RH) (z_parallel . z_perp)/H,
 s_R     = H b_M(RH)-R,
 I_theta_theta = R^2 E[b_M(RH)^2],
 I_RR = E[(H b_M(RH)-R)^2], I_theta_R=0.

An equivalent curvature expression is
 I_RR=1-E[H^2 b_M'(RH)],
 b_M'=1-b_M^2-(M-1)b_M/z.

The usual compact-set ordered reference construction, with theta as the
interest parameter and R as the nuisance parameter IN THIS REDUCED EXPERIMENT,
therefore gives pi(theta,R)dtheta dR proportional sqrt(I_RR(R)) dtheta dR.
The theta factor is constant because all information entries are independent
of theta. Near R=0, I_RR~2R^2/M and I_theta_theta~R^4/M. At large R with fixed
M, I_RR->1 and I_theta_theta=R^2-(M-1)+o(1). The radial weight is linear at
the origin and asymptotically flat. Every finite-data posterior is proper,
including M=1, and has finite radial moments.

REGULAR GROWING-M LIMIT
For R^2=kappa M with fixed kappa>0, H^2/M->1+kappa and
b_M(RH)->sqrt(kappa/(1+kappa)). Thus
 I_theta_theta ~ R^4/(R^2+M),
 I_RR -> 2*kappa/(1+2*kappa).
The first expression is exactly the flat estimator's Godambe information.
The second stays positive, so this reference prior has no singular effect in
that regular limit. Integrating R has the correct first-order width in this
invariant experiment; there is no Gaussian random-size assumption.

REFERENCE-PRINCIPLE CAVEAT
Taking the reference construction before reducing the fixed-direction normal
model is different. In full parameters (theta,R,xi), Fisher is orthogonal:
I_theta_theta=R^2, I_RR=1, I_xi_xi=R^2*g_sphere. Ordered scalar R then Haar
direction nuisance groups give uniform dR; grouping all v nuisance coordinates
together gives flat d^M v. Reduction first then reference gives the weight
derived here. Theta priority ALONE does not choose among these constructions.
It needs a commitment that the invariant experiment, rather than the full
direction-labelled experiment, is the reference-design experiment.

The new selecting principle can be stated precisely: count radii uniformly in
their invariant-data distinguishability coordinate eta_M(R)=int_0^R sqrt(I_RR)
dR, and prioritize mass angle over this nuisance. This is not a fitted prior.
It is nevertheless an extra information principle, not yet a consequence of
Newtonian composition. Joint Jeffreys would additionally multiply by
sqrt(I_theta_theta), another distinct priority choice.

KNOWN-EMPTY CHECK / REJECTION AS A COMPLETE SOLUTION
For measured F=2,a=0, unit noises, M=1, the derived law gives point2.323798594
and P(m<1)=.187754. Append two zero axes, each with force and acceleration SD
epsilon. In the invariant experiment their standardized zeros and dimension
are epsilon-independent. For every epsilon>0 the new law has P(m<1)=.314946;
the physical ratio-of-mean point tends to1.912131618 as epsilon->0 rather than
the M=1 answer. Thus it fails known-empty continuity despite passing the
regular-width and all-null-symmetry checks. This is a structural counterexample,
not a finite quadrature artefact; changing epsilon never enters its mass law.

Primary background: Bernardo's Reference Analysis explains ordered nuisance
reference construction; Bernardo & Rueda's Bayesian Hypothesis Testing gives
the uniform-radius reference after retaining full normal-model nuisance
structure (Section on multivariate normal means):
https://www.uv.es/~bernardo/IntStatRev.pdf . The formulas above are derived here.
"""

import json
from functools import lru_cache
import numpy as np

def bessel_ratio(z,M):
    z=np.asarray(z)
    if M==1:
        return np.tanh(z)
    if M==3:
        small=np.abs(z)<1e-3
        safe=np.where(small,1.,z)
        out=1/np.tanh(safe)-1/safe
        return np.where(small,z/3-z**3/45+2*z**5/945,out)
    # Continued fraction I_(nu+1)(z)/I_nu(z), nu=M/2-1.
    depth=int(max(60,np.max(np.abs(z))+20))
    r=np.zeros_like(z,dtype=float)
    for j in range(depth,0,-1):
        r=z*z/(M+2*j+r)
    return z/(M+r)

@lru_cache(None)
def quadrature(M,n):
    x,w=np.polynomial.hermite.hermgauss(n)
    z=np.sqrt(2)*x; w=w/np.sqrt(np.pi)
    if M==1:
        return z,0.,w
    else:
        shape=(M-1)/2
        k=np.arange(n,dtype=float)
        diag=2*k+shape
        off=np.sqrt(np.arange(1,n)*(np.arange(1,n)+shape-1))
        Jac=np.diag(diag)+np.diag(off,1)+np.diag(off,-1)
        lag,Q=np.linalg.eigh(Jac)
        wl=Q[0,:]**2
        return z[:,None],2*lag[None,:],w[:,None]*wl[None,:]

def fisher(R,M,n=128):
    z,off,weights=quadrature(M,n)
    H=np.sqrt((R+z)**2+off)
    b=bessel_ratio(R*H,M)
    sr=H*b-R
    return {'I_RR':float(np.sum(weights*sr*sr)),
       'I_theta_theta':float(R*R*np.sum(weights*b*b)),
       'R_score_mean':float(np.sum(weights*sr))}

def angular_reference(M,nr=180,nt=180):
    x,w=np.polynomial.legendre.leggauss(nr)
    R=(x+1)*6;wr=w*6
    xt,wt=np.polynomial.legendre.leggauss(nt)
    theta=(xt+1)*np.pi/4;wt=wt*np.pi/4
    ir=np.array([fisher(r,M)['I_RR'] for r in R])
    radial=wr*np.sqrt(ir)*np.exp(-R*R/2)
    z=R[:,None]*(2*np.sin(theta))[None,:]
    if M==1:
        A=np.cosh(z)
        first_magnitude_factor=np.ones_like(z)
    elif M==3:
        A=np.sinh(z)/z
        # If the two added axes' physical noise goes to zero, the retained
        # physical magnitude is R E[|xi_1| | R,H], for xi uniform on S^2.
        first_magnitude_factor=1-np.tanh(z/2)/z
    else:
        raise ValueError('This check uses M=1 and M=3 closed forms only.')
    L=radial@A
    p_below1=float(np.sum(wt[theta<np.pi/4]*L[theta<np.pi/4])/np.sum(wt*L))
    stacked= (radial*R)@A
    retained=(radial*R)@(A*first_magnitude_factor)
    def ratio(B):
        return float((wt*np.sin(theta))@B/((wt*np.cos(theta))@B))
    return {'M':M,'P_mass_below1':p_below1,
        'stacked_point':ratio(stacked),
        'point_when_hidden_physical_noise_tends_zero':ratio(retained)}

if __name__=='__main__':
    table=[]
    for M in [1,3,9]:
        for R in [.02,.1,.5,1.,3.,10.]:
            table.append({'M':M,'R':R,**fisher(R,M)})
    asymptotic=[]
    for M in [100,400]:
        for kappa in [.25,1.,4.]:
            R=np.sqrt(kappa*M)
            f=fisher(R,M,n=32)
            asymptotic.append({'M':M,'kappa':kappa,**f,
              'theta_information_over_Godambe_limit':f['I_theta_theta']/(R**4/(R*R+M)),
              'radial_information_limit':2*kappa/(1+2*kappa)})
    print(json.dumps({'fisher':table,'growing_dimension_checks':asymptotic,
        'known_empty_reference_checks':[angular_reference(1),angular_reference(3)]},indent=2))

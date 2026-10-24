"""Exact size-direction coupling for isotropic candidate measures.

Scratch mathematical audit; no final mass definition is asserted. Managed
context is VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

THE GENERAL IDENTITY
--------------------
Normalize the common noise pair and put
  A_theta=x*sin(theta)+y*cos(theta), B_theta=x*cos(theta)-y*sin(theta),
  H=||A||^2, Q=||B||^2, C=A.B, T=H+Q.
For a fixed isotropic positive candidate measure mu_M, assume its Gaussian-
damped transform and required moments are finite. Then

  F_M(H)=integral exp(-||v||^2/2 + A.v) mu_M(dv),
  evidence(theta)=exp(-T/2) F_M(H),
  kappa(H)=2*d log(F_M)/dH.

Differentiating the transform in A gives the EXACT posterior identities

  E[v|A] = kappa A,
  Cov(v|A) = kappa I + 2*kappa'(H) A A',
  E||v||^2 = H*kappa^2 + M*kappa + 2H*kappa'.

Thus for M>=2 every unit vector perpendicular to A has posterior candidate
variance kappa. The parallel variance is kappa+2H*kappa'. Positivity gives
kappa>=0 and kappa+2H*kappa'>=0, but does not determine either function.

Writing ell=log evidence,
  ell' = kappa C,
  ell'' = kappa(Q-H)+2*kappa' C^2.
At an unrestricted TLS maximum with eigenvalues h>=q and C=0,
  -ell'' = kappa(h)*(h-q).
For aligned data q=0, this is kappa(h)*h=sqrt(h)*E[v_parallel].

The same result is the missing-information identity: conditional complete
angular information E[A.v]=kappa H, minus score variance Var(B.v)
=kappa Q+2*kappa' C^2. In particular, the flat candidate posterior N(A,I)
does retain directional uncertainty; the precise issue is that its radial
measure fixes kappa=1 and E||v||^2=H+M. It would be inaccurate to say it ignores
candidate-direction uncertainty altogether.

RADIUS AND DIRECTION
--------------------
Write mu_M(dv)=rho_M(dr)*sigma_{M-1}(dn), v=r*n, using normalized spherical
area. At fixed r the posterior direction is von Mises-Fisher with concentration
a*r, a=sqrt(H). If L_M(t) is its mean cosine,

  kappa = E[r*L_M(a*r)]/a
        = E[r^2*(1-n_parallel^2)]/(M-1), M>=2.

The second equality follows by integrating the spherical coordinate density
exp(a*r*t)*(1-t^2)^((M-3)/2) by parts. Thus direction uncertainty and size
weighting cannot be assigned separately without changing the same law.

Integrating n first gives exp(-r^2/2) Z_M(a*r) rho_M(dr). Integrating r first
gives integral exp(-r^2/2+a*r*n_parallel) rho_M(dr), followed by spherical
integration. Tonelli gives identical answers for this SAME nonnegative
measure. Flat dv instead of flat dr changes rho_M from dr to r^(M-1)dr;
that discrepancy is a change of measure, not an order-of-integration effect.

ASYMPTOTIC RADIAL VOLUME MECHANISM
--------------------------------
Let H=h*M with h>1 and r=sqrt(M)*s. If the radial measure has no exponential
in M weighting near its interior saddle, its leading log integrand is
  M[-s^2/2+s*sqrt(h)*t + (1/2)*log(1-t^2)].
The stationary equations give t^2=1-1/h, s^2=h-1, kappa=1-1/h. With h=lambda+1
this is lambda/(lambda+1), the regular angular information reduction.

More generally a radial power rho_M(dr)=r^(gamma*M+o(M))dr contributes
gamma*log(s) to the action. The positive saddle obeys
  h*kappa^2+(1-h)*kappa-gamma=0.
Flat Cartesian measure has gamma=1 and kappa=1. Fixed radial powers have
gamma=0 and the corrected leading kappa. This is a conditional Laplace
argument, requiring an interior dominant saddle and regular radial weights;
it does not select one finite-data radial law. Many measures share this limit.

OBSERVABLE-VARIANCE HAAR CANDIDATE
--------------------------------
Take t>=1, v|t ~ N(0,(t-1)I), and starting measure dt/t. Let beta=1-1/t.
Then, with a=M/2 and z=H/2,
  p(beta|A) proportional (1-beta)^(a-1)*exp(z*beta), 0<beta<1,
  v|beta,A ~ N(beta*A,beta*I).
Consequently kappa=E[beta], kappa'=Var(beta)/2 exactly. Equivalently u=1-beta
has a Gamma(a,rate=z) law truncated to (0,1). For H>0,
  kappa=1-gamma_lower(a+1,z)/(z*gamma_lower(a,z));
at H=0, kappa=2/(M+2). This gives a proper posterior for every finite datum
and M>0, even though dt/t is improper. Its large-H limit is 1-M/H.

FULL COVARIANCE
---------------
Whiten the whole observation z. Let E_m have orthonormal columns spanning
the whitened codirectional subspace, P_m=E_m E_m', rank M, A_m=E_m'z.
The hierarchical observed covariance is I+(t-1)P_m. Its evidence is
  t^(-M/2)*exp(-Q_m/2-H_m/(2t)), H_m=z'P_m z, Q_m=||z||^2-H_m.
Integrating dt/t again gives exp(-||z||^2/2) F_M(H_m). It is continuous
through profile ties and preserves higher-order peaks: at the exact quartic
example, the fourth derivative is multiplied by positive kappa, not erased.
In physical candidate coordinates a, K_m=(B_m' Sigma^-1 B_m)^-1, b_m=K_m B_m'
Sigma^-1 y, the conditional law is N(beta*b_m,beta*K_m).

This construction is coherent but NOT forced merely by isotropy or by Gaussian
instruments. dt/t is the conditional Jeffreys scale measure for fixed m, since
I_tt=M/(2t^2). Joint Jeffreys for varying m is different:
  I_mm=(t-1)^2/(2t)*tr[(P_m')^2], I_mt=0,
  sqrt(det I) dm dt proportional sqrt(tr[(P_m')^2])*(t-1)/t^(3/2) dm dt.
Choosing fixed pi(m) times conditional dt/t is a substantive choice of which
parameter measure to hold fixed, although it needs no object-specific size.

It also inherits a known-empty-axis discontinuity in its angular law: adding
a measured-zero axis of equal force/acceleration SD epsilon leaves H_m and
the original whitened data unchanged, but changes M in F_M. The resulting
law depends on the added dimension independently of epsilon>0. Exact zero
variance with the empty axis removed is a different limiting experiment.
"""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss


def beta_moments(M,H,order=256):
    """Remove only the M=1 endpoint singularity; retain the large-M interior."""
    nodes,w=leggauss(order);v=(nodes+1)/2;w=w/2
    if M==1:
        u=v*v;beta=1-u
        logs=H*beta/2
    else:
        u=v;beta=1-u
        logs=H*beta/2+(M/2-1)*np.log(u)+np.log(M/2)
    p=w*np.exp(logs-max(logs));p=p/sum(p)
    mean=float(p@beta);second=float(p@(beta*beta))
    return {'M':M,'H':H,'kappa':mean,'beta_variance':second-mean*mean,
            'candidate_radius_squared_mean':H*second+M*mean,
            'log_F_up_to_M_only_constant':float(max(logs)+np.log(sum(w*np.exp(logs-max(logs)))))}


def radial_direction_quadrature(M,H,radial_power,order=260):
    """Explicit radius/direction quadrature for polynomial radial measures."""
    nodes,w=leggauss(order)
    rmax=np.sqrt(H+M)+12
    r=(nodes+1)*rmax/2;wr=w*rmax/2
    t=nodes;wt=w
    logp=(-r[:,None]**2/2 + np.sqrt(H)*r[:,None]*t[None,:]
          +radial_power*np.log(r[:,None])
          +(M-3)/2*np.log1p(-t[None,:]**2))
    shift=float(np.max(logp))
    mass=wr[:,None]*wt[None,:]*np.exp(logp-shift)
    norm=float(np.sum(mass));p=mass/norm
    para=r[:,None]*t[None,:]
    radius2=r[:,None]**2
    kappa=float(np.sum(p*para)/np.sqrt(H))
    transverse=float(np.sum(p*radius2*(1-t[None,:]**2))/(M-1))
    radius_mean=float(np.sum(p*radius2))
    parallelvar=float(np.sum(p*para*para)-(np.sum(p*para))**2)
    # Two elimination orders are explicitly performed on the same weighted matrix.
    radius_first=float(np.dot(wt,np.sum(wr[:,None]*np.exp(logp-shift),axis=0)))
    direction_first=float(np.dot(wr,np.sum(wt[None,:]*np.exp(logp-shift),axis=1)))
    assert abs(radius_first/direction_first-1)<2e-14
    assert abs(transverse-kappa)<5e-10
    return {'M':M,'H':H,'radial_power':radial_power,'kappa':kappa,
            'transverse_variance':transverse,'parallel_variance':parallelvar,
            'candidate_radius_squared_mean':radius_mean,
            'elimination_order_relative_difference':radius_first/direction_first-1}


def main():
    radial=[]
    for M in (3,20,80,320):
        for power in (0,M-1):
            row=radial_direction_quadrature(M,2*M,power)
            if power==M-1:
                assert abs(row['kappa']-1)<1e-11
                assert abs(row['candidate_radius_squared_mean']-3*M)<1e-8
            radial.append(row)
    beta=[]
    for M in (1,3,20,80,320):
        for H in (0.,4.,2.*M):
            row=beta_moments(M,H,512)
            if H==0:
                assert abs(row['kappa']-2/(M+2))<1e-4
            # Independently check derivative of kappa against Var(beta)/2.
            step=1e-3
            derivative=(beta_moments(M,H+step,512)['kappa']-beta_moments(M,H-step,512)['kappa'])/(2*step)
            assert abs(derivative-row['beta_variance']/2)<1e-8
            row['kappa_derivative_finite_difference']=derivative
            beta.append(row)
    # Exact quartic data: H(0)=k^2*(2+sqrt(2)), ell_profile''''=-3k^2/sqrt(2).
    quartic=[]
    for k in (1.,10.,100.):
        H=k*k*(2+np.sqrt(2))
        row=beta_moments(2,H,512)
        quartic.append({'k':k,'H_at_peak':H,'kappa':row['kappa'],
                        'marginal_log_evidence_fourth_derivative':-3*k*k*row['kappa']/np.sqrt(2)})
    # Angular odds at theta=pi/2 versus theta=0 for base F=2,a=0.
    # Endpoint values stand for their interior limits; any fixed positive pi cancels
    # from the ratio of the two dimension-dependent likelihood-odds ratios.
    odds=[]
    for M in (1,2):
        value=beta_moments(M,4.,1024)['log_F_up_to_M_only_constant']
        odds.append({'M':M,'likelihood_odds_H4_vs_H0':float(np.exp(value))})
    out={'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
         'radial_direction_identity_checks':radial,'observable_variance_mixture':beta,
         'quartic_survival':quartic,'known_empty_axis_likelihood_odds':odds,
         'status':'Exact coupling and specified joint-law construction; uniqueness and known-empty continuity remain unresolved/failing respectively.'}
    Path(__file__).with_name('size_direction_identity_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()

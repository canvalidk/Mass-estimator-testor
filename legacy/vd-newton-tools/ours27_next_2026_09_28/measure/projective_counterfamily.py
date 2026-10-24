"""Projective, rotation-invariant candidate-size priors still leave shape free.

Use conventional inverse-Wishart notation:
    Omega_d ~ IW_d(nu_d=delta+d-1, Psi=(delta-2) I_d), delta>2.

Then E[Omega_d]=I_d, and each k x k principal block is
IW_k(delta+k-1,(delta-2)I_k). This is exactly dimension-projective and puts
positive density on every SPD covariance in every finite dimension. All
members are proper, centered in the same noise covariance, rotation-invariant,
unit-covariant after whitening, and reciprocal under channel swap because
the candidate-subspace frame changes orthogonally. They differ only in shape.

Source for the underlying marginalization identity, using Dawid's alternative
dimension-free parameter delta: Shah, Wilson, Ghahramani (2014), section 2,
https://www.cs.cmu.edu/~andrewgw/tprocess.pdf . Our use here is an independent
counterexample, not a proposal to adopt that prior.

The inverse-Wishart construction is over candidate covariance, NOT noise.
Measurement errors remain Gaussian with their declared known covariance.
For a proper starting law pi(m), posterior totality follows because the
Gaussian integrated likelihood is positive and bounded, and the conditional
candidate covariance Beta=Omega(I+Omega)^-1 is between 0 and I. Conditional
candidate means and first norms are uniformly bounded in whitened units;
physical acceleration candidates decay O(1/m) as m grows, as in ours27.

The numerical check uses one 1D reading x=2,y=0, s=1, uniform angle prior,
Gaussian measurement noise, and the unchanged ratio-of-mean magnitudes.
Every delta uses E[Omega]=1, so changed answers prove unresolved *shape*,
even after the default candidate covariance is fixed equal to noise.

WARNING: exact dimension marginalization is stronger and distinct from merely
regrouping known axes. This family passes the stronger prior restriction but
does not claim that dropping an actually measured axis preserves its evidence.
Shared covariance learning still violates exact appended-null invariance.
"""

import json
import math
import numpy as np

def solve(delta, n=256):
    x,w=np.polynomial.legendre.leggauss(n)
    theta=(x+1)*np.pi/4
    wt=w*np.pi/4
    beta=(x+1)/2
    wb=w/2
    alpha=delta/2
    scale=(delta-2)/2
    # Prior under omega variance -> beta, with one marginal data determinant.
    logbase=-(alpha+1)*np.log(beta)+(alpha-.5)*np.log1p(-beta)-scale*(1-beta)/beta
    logbase-=np.max(logbase)
    bh=beta[:,None]
    h=2*np.sin(theta)[None,:]
    joint=np.exp(logbase[:,None]+.5*bh*h*h)*wb[:,None]*wt[None,:]
    erf=np.vectorize(math.erf)(np.sqrt(bh/2)*h)
    absmean=np.sqrt(2*bh/np.pi)*np.exp(-.5*bh*h*h)+bh*h*erf
    mass=float(np.sum(joint*absmean*np.sin(theta)[None,:])/np.sum(joint*absmean*np.cos(theta)[None,:]))
    endpoint_odds=float(np.sum(np.exp(logbase+2*beta)*wb)/np.sum(np.exp(logbase)*wb))
    return {'delta':delta,'prior_E_Omega':1,'point':mass,'mass_law_endpoint_odds':endpoint_odds}

if __name__=='__main__':
    print(json.dumps({'reading':[2,0],'s':1,'results':[solve(d) for d in [3,4,8,32]],
      'quadrature_repeat_delta4':solve(4,384)},indent=2))

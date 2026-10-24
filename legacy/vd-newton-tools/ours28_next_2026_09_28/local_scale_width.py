"""Analytic large-N width of a per-reading radial mass evidence.

At the tested common-noise mass angle:
  A ~ N(0,t I_d), B ~ N(0,I_d), independent; H=|A|^2.
For ell=log F(H), kappa=2 F'/F:
  score=kappa(H) A.B
  score variance=E[H kappa^2]
  expected negative curvature=(1-1/t) E[H kappa].

The Gaussian integration-by-parts identity supplies the second equality.
For exact equality for every t>1, uniqueness of Laplace transforms implies
  2 kappa' + (d/H-1) kappa + kappa^2=0
or 2 H F''+(d-H)F'=0. A finite positive Gaussian scale-mixture evidence is
analytic at H=0 with kappa(0)>0, and cannot obey this ODE.

This concerns a fixed, nondegenerate per-reading scale law as N grows.
It is NOT an impossibility result for a scale distribution learned jointly
from the readings, or for intervals calibrated under their own prior world.
"""
from pathlib import Path
import json
import math
import numpy as np
from functools import lru_cache


@lru_cache(maxsize=None)
def beta_rule(a, b, order=240):
    """Normalized Gauss-Jacobi nodes for Beta(a,b), via Golub-Welsch."""
    # Jacobi weight on [-1,1]: (1-x)^alpha (1+x)^beta.
    alpha, beta = b - 1, a - 1
    n = np.arange(order, dtype=float)
    ab = alpha + beta
    diag = np.empty(order)
    diag[0] = (beta-alpha)/(ab+2)
    k = n[1:]
    diag[1:] = (beta*beta-alpha*alpha)/((2*k+ab)*(2*k+ab+2))
    k = np.arange(1, order, dtype=float)
    off = np.empty(order-1)
    off[0] = math.sqrt(4*(1+alpha)*(1+beta) /
                      ((2+ab)**2*(3+ab)))
    k = k[1:]
    off[1:] = np.sqrt(4*k*(k+alpha)*(k+beta)*(k+ab) /
                     ((2*k+ab)**2*(2*k+ab+1)*(2*k+ab-1)))
    matrix = np.diag(diag)+np.diag(off, 1)+np.diag(off, -1)
    x, v = np.linalg.eigh(matrix)
    return (x+1)/2, v[0]**2


@lru_cache(maxsize=None)
def gamma_rule(shape, order=240):
    n = np.arange(order, dtype=float)
    diag = 2*n+shape
    k = np.arange(1, order, dtype=float)
    off = np.sqrt(k*(k+shape-1))
    x, v = np.linalg.eigh(np.diag(diag)+np.diag(off,1)+np.diag(off,-1))
    w = v[0]**2
    return x, w/w.sum()


def kernel(h, d, a=.5, c=.5, order=240):
    beta, weight = beta_rule(a, c+d/2, order)
    h = np.asarray(h)
    exponent = .5*h[...,None]*beta
    maximum = np.max(exponent, axis=-1)
    weights = np.exp(exponent-maximum[...,None])*weight
    norm = weights.sum(axis=-1)
    kappa = (weights*beta).sum(axis=-1)/norm
    second = (weights*beta*beta).sum(axis=-1)/norm
    return kappa, .5*(second-kappa*kappa)


def cell(d, omega, a=.5, c=.5, order=240):
    t = 1+omega*omega
    x, w = gamma_rule(d/2+1, order)
    h = 2*t*x
    k, kp = kernel(h, d, a, c, order)
    v = d*t*np.dot(w,k*k)
    j = d*(t-1)*np.dot(w,k)
    ratio = v/j
    coverage = math.erf(1.959963984540054/math.sqrt(2*ratio))
    # Check direct curvature on Gamma(d/2), not its size-biased law.
    xx, ww = gamma_rule(d/2, order)
    hh = 2*t*xx
    kk, kkp = kernel(hh,d,a,c,order)
    direct_j = np.dot(ww,kk*(hh-d)-2*kkp*hh)
    return {
        "d":d,"omega":omega,"beta_prior":[a,c],
        "score_variance":v,"expected_curvature":j,
        "curvature_identity_error":float(abs(j-direct_j)),
        "sampling_variance_over_law_variance":ratio,
        "asymptotic_95_coverage":coverage
    }


def main():
    rows=[]
    for d in (1,3):
        for omega in (.25,.5,1.,2.,4.):
            for a in (.25,.5,1.,2.):
                rows.append(cell(d,omega,a,a))
    check=cell(3,.5,order=400)
    coarse=cell(3,.5)
    error=abs(check["sampling_variance_over_law_variance"]-
              coarse["sampling_variance_over_law_variance"])
    assert error<1e-8,error
    assert max(x["curvature_identity_error"] for x in rows)<2e-6
    out={"source_commit":"592d144531fa4a30d1c2229a4c3b528e84043d77",
         "cells":rows,"refinement_240_to_400":error,
         "scope":"Regular fixed per-reading scale-law widths. Test candidates, not a mass-definition axiom."}
    Path(__file__).with_name("local_scale_width_results.json").write_text(
        json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps({"arcsine_cells":[r for r in rows if r["beta_prior"][0]==.5],
                      "refinement":error},indent=2))


if __name__=="__main__":
    main()

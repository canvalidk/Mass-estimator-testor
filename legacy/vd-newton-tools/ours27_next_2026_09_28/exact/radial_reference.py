"""Exact radial invariant reduction and its reference-prior limitations.

Scratch derivation, not a proposed final mass equation.
Managed snapshot: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

One noise pair, M scalar components, standardized x,y. At angle theta let
A=x sin(theta)+y cos(theta), H=||A||^2, and candidate v=R eta with eta
uniform on S^(M-1). The compact orthogonal group supplies the unique uniform
orientation measure. Integrating that measure gives

 L(theta,R) proportional to exp(-R^2/2) 0F1(M/2; R^2 H/4).

This is the rank-one noncentral-Wishart invariant likelihood. It also follows
directly from E_eta exp(R eta'A)=0F1(M/2;R^2 H/4). The omitted factor
exp(-(||x||^2+||y||^2)/2) depends only on observations. For M=1, the sphere
is the two signs and 0F1(1/2;R^2 H/4)=cosh(R sqrt(H)).

If one additionally chooses UNIFORM RADIUS dR (not forced by rotations),
the exact radial integrals are

 I_k(H)=integral_0^infinity R^k exp(-R^2/2) 0F1(M/2;R^2 H/4)dR
       =2^((k-1)/2) Gamma((k+1)/2)
          1F1((k+1)/2;M/2;H/2), k>-1.

Thus the mass-angle law is proportional to I_0(H(theta)), and the ratio of
mean sizes is
 mhat=s* integral sin(theta) I_1(H(theta))dtheta
           / integral cos(theta) I_1(H(theta))dtheta.
This posterior and its size means exist for every finite data vector,
including M=1 and all-zero data. It preserves one-pair swap, force/accel
units, O(M) component rotations, and regrouping the SAME M coordinates.
The full-vector prior is ||v||^(1-M)dv dtheta, so it is different from
independent flat candidate components. The old flat-v law corresponds to
k=M-1, for which 1F1(M/2;M/2;H/2)=exp(H/2).

Why this is not a finished unique replacement:
1. Rotational symmetry fixes eta but permits arbitrary radial weights.
2. The improper family ||v||^(1-M) has compatible one-coordinate marginals
   for M>1: integral (r^2+z^2)^(-M/2)dz
        =sqrt(pi)Gamma((M-1)/2)/Gamma(M/2) * r^(1-M).
   But its M=2 -> M=1 marginal diverges logarithmically. Directly defining
   the proper M=1 posterior does not repair this prior projectivity failure.
3. An isotropic global radius does not directly settle full heterogeneous
   covariance geometry or conditional candidate composition across designs.

A specific reference-rule alternative (not a new uniqueness theorem):
Under tested theta, H~noncentral-chi-square_M(R^2). The radial Jeffreys
weight is J_M(R)=sqrt(E[(d_R log f(H|R))^2]). The score is
 -R + sqrt(H) I_(M/2)(R sqrt(H))/I_(M/2-1)(R sqrt(H)).
Expansion gives J_M(R)~sqrt(2/M)*R near zero; J_M(R)->1 for large R.
Both statements concern the invariant reduced experiment. Choosing the
ordered reference rule (theta of interest, R nuisance) is extra input;
joint Jeffreys, or using the unreduced experiment, need not agree. This
radial weight has a proper posterior, but its large-radius behavior retains
the M2 -> M1 marginal divergence. It therefore does not close the problem.

The code independently checks radial integration identities and the
Jeffreys limiting behavior for M=1 and M=3 using Gaussian/Laguerre
quadrature. No simulation coverage or uniqueness claim is made here.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.laguerre import laggauss


def hypergeom(a, b, z):
    """0F1 if a=None, otherwise 1F1; positive arguments in this probe."""
    z = np.asarray(z, dtype=float)
    term = np.ones_like(z)
    total = term.copy()
    for k in range(1, 1000):
        term *= z / (k * (b + k - 1))
        if a is not None:
            term *= a + k - 1
        total += term
        if np.max(abs(term / total)) < 2e-15:
            return total
    raise AssertionError("Hypergeometric series did not converge")


def numeric_radial(dimension, h, power, nodes=240):
    x, w = leggauss(nodes)
    upper = np.sqrt(h) + np.sqrt(dimension + abs(power)) + 12
    radius = (x + 1) * upper / 2
    weights = w * upper / 2
    integrand = radius**power * np.exp(-radius*radius/2)
    integrand *= hypergeom(None, dimension/2, radius*radius*h/4)
    return float(weights @ integrand)


def radial_jeffreys(dimension, radius, nodes=100):
    x, wx = hermgauss(nodes)
    parallel = radius + np.sqrt(2) * x
    wx = wx / np.sqrt(np.pi)
    if dimension == 1:
        score = -radius + parallel * np.tanh(radius * parallel)
        return float(np.sqrt(wx @ score**2))
    if dimension == 3:
        z, wz = laggauss(nodes)
        norm = np.sqrt(parallel[:, None]**2 + 2*z[None, :])
        kappa = radius * norm
        small = kappa < .01
        ratio = np.empty_like(kappa)
        kk = kappa[small]
        ratio[small] = kk/3 - kk**3/45 + 2*kk**5/945
        kk = kappa[~small]
        ratio[~small] = 1/np.tanh(kk) - 1/kk
        score = -radius + norm * ratio
        return float(np.sqrt(wx @ (score**2) @ wz))
    raise ValueError(dimension)


def main():
    integration = []
    for dimension in (1, 2, 6, 21):
        for h in (0., 1., float(dimension), float(dimension+25)):
            for power in (0, 1, dimension-1):
                got = numeric_radial(dimension, h, power)
                expected = float(2**((power-1)/2)*math.gamma((power+1)/2)
                                 * hypergeom((power+1)/2, dimension/2, h/2))
                error = abs(got/expected-1)
                assert error < 2e-11, (dimension,h,power,error)
                integration.append({"M":dimension,"H":h,"power":power,"relative_error":error})
    jeffreys = []
    for dimension in (1, 3):
        for radius in (.0001, .01, .1, .5, 1., 2., 5., 10.):
            value = radial_jeffreys(dimension, radius)
            refined = radial_jeffreys(dimension, radius, 140)
            jeffreys.append({"M":dimension,"radius":radius,"J_M":value,
                             "relative_to_small_R_limit":value/(np.sqrt(2/dimension)*radius),
                             "quadrature_refinement_absolute_error":abs(refined-value)})
            assert abs(refined-value)<1e-5
            if radius==.0001:
                assert abs(value/(np.sqrt(2/dimension)*radius)-1)<1e-6
            if radius==10.:
                assert abs(value-1)<.02
    divergent_marginal = [{"cutoff":cutoff,"integral_M2_prior_over_extra_coordinate_at_v1_1":2*math.asinh(cutoff)}
                          for cutoff in (10.,100.,1000.,1e6)]
    out={"source_commit":"b4bfc2040031ef010416334d0c5f4aa2ad05cbd8",
         "max_radial_identity_relative_error":max(r["relative_error"] for r in integration),
         "radial_identity_checks":integration,"radial_Jeffreys_checks":jeffreys,
         "M2_to_M1_marginal_divergence":divergent_marginal,
         "scope":"Exact one-pair radial candidates; unresolved radial premise, M1 projectivity, and full covariance extension."}
    Path(__file__).with_name("radial_reference_results.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k!="radial_identity_checks"},indent=2))


if __name__=="__main__":
    main()

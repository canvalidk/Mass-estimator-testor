"""Tensor-Fisher boundary cutoff audit. Scratch only, 2026-09-28.

Starting measure q(R)dR=|I+R|^{-(d+1)/2}dR, R>0.  With N readings,
beta=R(I+R)^(-1) has posterior proportional

 |I-beta|^{(N-d-1)/2} exp(tr(beta S)/2) d beta, 0<beta<I.

S=sum A_i A_i'. Conditional candidates v_i|beta,A_i are
N(beta A_i,beta). Consider N=1, d=2 or3, and a proper fixed mass measure.
Every finite cutoff gives a proper joint law. The following are EXPLICIT
different normalized limiting laws, so propriety of their individual
pushforwards does not define a cutoff-independent tensor measure.

d=3,N=1: write X=I-beta with ordered eigenvalues 0<x1<x2<x3<1.
At S=0 its eigenvalue density is
 (x2-x1)(x3-x1)(x3-x2)(x1*x2*x3)^(-3/2).

1. Eigenvalue ceiling R<L I is x1>delta=1/(1+L). Integrating x1 near0
   has leading factor 2 delta^(-1/2). The remaining eigenvalue density is
   (x3-x2)/sqrt(x2*x3), 0<x2<x3<1. There is one beta eigenvalue1;
   its perpendicular 2x2 block B has matrix-Beta parameters (3/2,1),
   with density |I-B|^(-1/2). Thus E tr(beta)=1+2*(3/5)=11/5.

2. Determinant ceiling |I+R|<L is x1*x2*x3>1/L. The lower x1 bound is
   delta/(x2*x3), delta=1/L. Its integral contributes
   2 delta^(-1/2)*sqrt(x2*x3). The remaining density is x3-x2.
   Now B is uniform in matrix Lebesgue measure 0<B<I, matrix-Beta
   parameters (3/2,3/2). E tr(beta)=1+2*(1/2)=2.

Both cutoffs are rotation invariant, dimensionless, nested, and exhaust the
same SPD cone. They give different candidate covariances even at zero data:
E(vv')=(11/15)I versus (2/3)I. For nonzero data their mass evidences also
differ: F(H)=E exp(H s/2), s=u'beta u, has F'(0)=11/30 versus1/3.
The finite-data limits are these baseline boundary measures tilted by
exp(tr(beta S)/2), which is bounded between1 and exp(tr(S)/2).
They induce proper weak limits for the actual candidates and mass: the
Gaussian conditional kernel has mean beta*A and covariance beta bounded
on 0<=beta<=I, and is weakly continuous there. Thus the different limits
cannot be dismissed as divergent nuisance coordinates with one common
observable pushforward. Finite-data angular differences below verify this.

The singular-Wishart extension with A~Wishart_d(d+1,I), C~Wishart_d(N,I),
beta=(A+C)^(-1/2)A(A+C)^(-1/2) is a THIRD rule. For d3,N1 it has TWO
beta eigenvalues1 (rank(I-beta)=1), unlike both cutoff limits above.
Its mean is ((d+1)/(d+1+N))I=(4/5)I, trace12/5. It is not the normalized
eigenvalue- or determinant-ceiling limit of this improper interior measure.

d=2,N=1: the divergence is logarithmic. Both ceilings above give one beta
eigenvalue1, an independent uniform remaining eigenvalue b in(0,1), and
uniform orientation. E tr(beta)=3/2. This agrees with the singular-Wishart
extension in this particular dimension. Homothetic cutoff perturbations
change log L only by bounded terms, so this borderline limit has robustness.
It is NOT universal over invariant exhaustions: the rotation-invariant
cutoff log(1+lambda_max(R)) < [1+b] log L, b=lambda_min/(1+lambda_min),
weights the leading logarithmic divergence by 1+b. The remaining b density
becomes (2/3)(1+b), and E tr(beta)=14/9 instead of3/2. Both exhaustions are
nested and eventually contain every fixed R. No direction labels are used.

The executable checks the d3 finite-cutoff means by exact integration of x1
and 2D Gaussian quadrature, and evaluates both limiting angular evidences.
No Monte Carlo or calibration target is used. This is a rejection of a
cutoff-independent extension, not a claim that no additional axiom could
select a particular boundary rule.
"""

from __future__ import annotations

import json
import math
import numpy as np


def gl(n, lo=0.0, hi=1.0):
    z, w = np.polynomial.legendre.leggauss(n)
    return lo + (hi - lo) * (z + 1) / 2, w * (hi - lo) / 2


def finite_cutoff_trace(delta, cutoff, n=160):
    # Integrate x1 exactly; use yj=sqrt(xj) for the remaining smooth limit.
    if cutoff == "eigen":
        ymin = math.sqrt(delta)
    elif cutoff == "determinant":
        ymin = delta ** (1 / 6)
    else:
        raise ValueError(cutoff)
    y3, w3 = gl(n, ymin, 1)
    r2, wr2 = gl(n)
    y2min = np.full_like(y3, math.sqrt(delta)) if cutoff == "eigen" else delta**0.25 / np.sqrt(y3)
    y2 = y2min[:, None] + (y3 - y2min)[:, None] * r2[None, :]
    x2, x3 = y2**2, y3[:, None] ** 2
    lower = np.full_like(x2, delta) if cutoff == "eigen" else delta / (x2 * x3)
    # Integral of x1^(-3/2) (x2-x1)(x3-x1), lower<x1<x2.
    i0 = (
        2 * x2 * x3 * (lower**-0.5 - x2**-0.5)
        - 2 * (x2 + x3) * (x2**0.5 - lower**0.5)
        + (2 / 3) * (x2**1.5 - lower**1.5)
    )
    i1 = (
        2 * x2 * x3 * (x2**0.5 - lower**0.5)
        - (2 / 3) * (x2 + x3) * (x2**1.5 - lower**1.5)
        + (2 / 5) * (x2**2.5 - lower**2.5)
    )
    weight = (
        w3[:, None] * wr2[None, :] * (y3 - y2min)[:, None]
        * 4 * y3[:, None] * y2 * (x3 - x2) * (x2 * x3) ** -1.5
    )
    z = np.sum(weight * i0)
    return float(np.sum(weight * ((3 - x2 - x3) * i0 - i1)) / z)


def limiting_evidence(h_values, kind, n=64):
    # For a uniform direction on S2, |u dot e| is uniform on[0,1].
    t, tw = gl(n)
    phi, pw = gl(n, 0, math.pi / 2)
    b = np.sin(phi) ** 2
    if kind == "eigen":
        bw = pw * 3 * np.sin(phi) ** 2 * np.cos(phi)  # Beta(3/2,1)
        s = t[:, None] ** 2 + (1 - t[:, None] ** 2) * b[None, :]
    elif kind == "determinant":
        bw = pw * (16 / math.pi) * np.sin(phi) ** 2 * np.cos(phi) ** 2
        s = t[:, None] ** 2 + (1 - t[:, None] ** 2) * b[None, :]
    elif kind == "singular_wishart":
        bw = pw * 3 * np.sin(phi) ** 2 * np.cos(phi)
        s = 1 - t[:, None] ** 2 * b[None, :]
    else:
        raise ValueError(kind)
    weights = tw[:, None] * bw[None, :]
    return np.array([np.sum(weights * np.exp(h * s / 2)) for h in h_values])


def angular_case(kind):
    # One reading d3, x=(2,0,0),y=0, unit noise, uniform mass angle.
    t1, w1 = gl(64, 0, math.pi / 4)
    t2, w2 = gl(64, math.pi / 4, math.pi / 2)
    theta, tw = np.concatenate([t1, t2]), np.concatenate([w1, w2])
    evidence = limiting_evidence(4 * np.sin(theta) ** 2, kind)
    return {
        "evidence_H4_relative_H0": float(limiting_evidence([4], kind)[0]),
        "P_mass_below_one": float(np.dot(w1, evidence[:len(t1)]) / np.dot(tw, evidence)),
    }


def main():
    result = {
        "d3_N1_expected_trace_limits": {"eigen": 11/5, "determinant": 2, "singular_wishart": 12/5},
        "d3_N1_finite_cutoffs": {
            kind: {str(delta): finite_cutoff_trace(delta, kind) for delta in (1e-2, 1e-4, 1e-6, 1e-8)}
            for kind in ("eigen", "determinant")
        },
        "quadrature_check_n240": {kind: finite_cutoff_trace(1e-8, kind, 240) for kind in ("eigen", "determinant")},
        "d3_N1_limiting_angular_laws": {kind: angular_case(kind) for kind in ("eigen", "determinant", "singular_wishart")},
        "d2_N1_trace_limits": {"ordinary_ceiling": 3/2, "curved_ceiling": 14/9},
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Proper tensor half-Cauchy constructions: scratch derivation, 2026-09-28.

Source snapshot VD-docs@592d144531fa4a30d1c2229a4c3b528e84043d77.
Read the four local staged sources in ../sources: ours28_candidate_push_scale_
as_nuisance, ours28_size_measure_jeffreys_check, ours28_known_empty_vs_measured_
empty, and ours27_c3_from_betweenness_and_the_scale_fork (all 2026-09-28).
No construction below is adopted. Candidate covariances are measures over
hypothetical aligned vectors, not assumed ontic distributions of pushes.

1. PROPER GRASSMANN TENSOR
-------------------------
In standardized candidate coordinates, take R>0 and
 q_d(R) proportional |R|^(-1/2)|I+R|^(-d).
Equivalently W=R(I+R)^(-1) has real matrix-Beta(d/2,d/2) law:
 q(W) proportional |W|^(-1/2)|I-W|^(-1/2), 0<W<I.
It is proper for every d>=1. Its eigenvalue density has arcsine weights
times the real Vandermonde product. It is orthogonally invariant, and
R and R^(-1) have the same distribution. In d1, sqrt(R) is half-Cauchy.
It comes from two independent square Gaussian/Wishart blocks, or the
principal block of a Haar-uniform d-plane projector in R^(2d).

Given one shared R and N candidate-coordinate readings A_i, let
 S=sum A_i A_i'. The joint posterior has W density
 |W|^(-1/2)|I-W|^((N-1)/2) exp(tr(W S)/2).
Thus it is a tilted matrix-Beta(d/2,(d+N)/2), proper for every N>=0.
Conditional candidates are v_i|W,A_i ~ N(W A_i,W).
For full instrument covariance use a_i*=K_i^(1/2)v_i in a specified frame,
F_i*=m a_i*, and multiply by the declared proper mass measure pi(m).
Evidence is finite since the determinant likelihood factor is <=1 and
the exponential factor is bounded by exp(sum||A_i||^2/2). Conditional
candidate moments are bounded by the observed projections and instruments.
Mass/size readout is therefore total. Units and channel-swap equivariance
follow when K and the candidate frames are transformed consistently; full
matrix inversion symmetry is a DIFFERENT, additional push/noise premise.

For N=1, the marginal mass evidence is the ordinary Kummer function
  g_d(H) proportional 1F1(d/2; d+1/2; H/2), H=||A||^2.
The same actual one-vector prior is obtained from the SCALAR mixture
  v|lambda ~ N(0,lambda I_d), lambda~BetaPrime(d/2,1/2).
Reason: R11 has that beta-prime law, and isotropy plus the Gaussian
characteristic function determines the entire one-vector law.
Thus a per-reading tensor changes the radial shape; it supplies no additional
directional flexibility to one isotropic candidate cloud. Its effective
scalar shrinkage prior has E[lambda/(1+lambda)]=d/(d+1), not1/2 for d>1.

This family is not projective: a principal k block of R_d has matrix-F
parameters (d/2,k/2), not (k/2,k/2). Already R11 is BetaPrime(d/2,1/2),
whereas the intended d1 marginal is BetaPrime(1/2,1/2).
Known-empty reduction also fails with a local noise-unit metric. Adding
zero axes of noise SD epsilon in both channels leaves the enlarged H
unchanged at EVERY epsilon>0, but replaces g_1 by g_d. It cannot approach
g_1 as epsilon->0. This remains distinct from a declaration restricting
the candidate space before assigning the measure.

2. WHY FULL INVERSION SYMMETRY CANNOT REPAIR PROJECTIVITY
-------------------------------------------------------
Suppose an SPD random R is orthogonally invariant, R has the same law as
R^(-1), and R11 has the same law as 1/R11. Schur complements give
  (R^(-1))11 >= 1/R11,
with equality iff the off-diagonal first row vanishes. Inversion and scalar
reciprocity give the two sides equal distributions. Apply a bounded strictly
increasing function to conclude equality almost surely (no moments needed).
Repeat for every row: R is diagonal almost surely. Orthogonal invariance,
including 45-degree rotations of coordinate planes, then forces R=lambda I.

This also applies when only ACTUAL CANDIDATE marginals are required to be
projective: the normal variance mixture v1|R11 identifies the distribution
of R11 by uniqueness of its Laplace transform. A scalar half-Cauchy candidate
marginal therefore forces the self-reciprocal BetaPrime(1/2,1/2) variance.
Consequently a nontrivial projective tensor extension must relax full matrix
inversion symmetry, even while preserving scalar half-Cauchy marginals.

3. CONSTRUCTIVE PROJECTIVE EXTENSION
-----------------------------------
Let T~chi-square_1, independent of C_d~InverseWishart_d(nu=d, scale=I_d),
and set R_d=T C_d. Every distribution is proper. The principal k block of
C_d is InverseWishart_k(k,I_k), so R_d is projective in d. Its density is
  q_d(R) proportional
  |R|^(-(2d+1)/2) [1+tr(R^(-1))]^(-(d^2+1)/2).
It is orthogonally invariant and independent of mass in the declared
standardized frame. R11 is chi-square_1 / chi-square_1, exactly the intended
half-Cauchy-squared law. It cannot have full tensor inversion symmetry by2.

For V=[v_1,...,v_N], integrate C to obtain the proper matrix-t law
 p(V|T) proportional T^(-dN/2)|I+VV'/T|^(-(d+N)/2).
The construction is projective both in spatial coordinates and reading
columns. For N=1, v|T is a multivariate Cauchy, so v=sqrt(T) z/|g| with
z~N_d(0,I),g~N(0,1). Since sqrt(T)=|g0|, this is EXACTLY the original
shared scalar half-Cauchy mixture for one reading. Its one-reading posterior,
mass law, and ratio-of-means point coincide with ours28, including the
known-empty difficulty. Only N>=2 can use the extra covariance structure.

FRAME REQUIREMENT: the matrix-t density depends on cross-reading inner
products. It is invariant under one common row rotation, but not independent
rotations of the readings. To share a tensor across readings while retaining
their own coordinate conventions, one must carry declared transport maps
identifying their tensor frames. In an isotropic instrument experiment no
such alignment is supplied intrinsically. Treating all noise-coordinate
arrays as one numerical frame silently adds that identification. A scalar
shared scale requires none. This is why the positive tensor construction
is useful only when a covariance-sharing group and its transport are part
of the experimental structure, not merely inferred from a common mass.

4. REGULAR WIDTH OF A LEARNED TENSOR
----------------------------------
For a group with one common instrument projector family, let U have
orthonormal columns and choose a local horizontal frame U^T dot(U)=0. Set
D=(I-UU^T) dot(U), Tgeom=D^T D. If fixed candidate vectors have empirical covariance
C_N=N^(-1)sum v_i v_i' -> C>0, the covariance likelihood learns R=C.
Any proper prior with positive continuous density at C, including3, leaves
this regular leading result unchanged. Write B=C(I+C)^(-1). With independent
standard tangent/normal noise e_i,f_i, the mass score is
  sum f_i' D B(v_i+e_i).
Its variance and expected negative curvature both equal
  N tr[D C(I+C)^(-1) C D'].
The expected cross Hessian and score covariance with the tensor nuisance
vanish. This is a fixed-candidate second-moment statement, not an assumption
of ontic Gaussian pushes; regularity still requires a CLT, interior C,
local identification, and a sufficiently informative growing group.

Heterogeneous projectors do NOT automatically inherit that equality. With
horizontal frames U_i and T_i=D_i'D_i, one shared tensor still learns the
unweighted C_N, while
  H = sum tr[T_i B v_i v_i'],
  V = sum tr[T_i B(v_i v_i'+I)B],
  V-H = -sum tr[T_i B(v_i v_i'-C_N)(I-B)].
Cross nuisance terms still vanish in these horizontal coordinates. Equality
holds for common T_i or asymptotic factorization of candidate covariance and
instrument geometry, not arbitrary fixed candidate/design correlations.
Different replicated instrument groups, each learning its own R_g, restore
the groupwise equality. For a scalar in each group, no orientation transport
is needed. Continuous sharing across instrument designs needs an additional
rule. Arbitrary full-covariance frames also have within-subspace angular
velocity; the horizontal transport used here must be specified rather than
silently treating independent frame changes as one shared R coordinate.

The numerical section checks closed-form angular differences, the simple
heterogeneous scalar width identity, and the independent-frame dependence.
It does not fit a calibration correction or adopt a new equation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import numpy as np


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return lo + (hi-lo)*(x+1)/2, w*(hi-lo)/2


def beta_evidence(h, a, b, n=96):
    phi, pw = gl(n, 0, math.pi/2)
    v = np.sin(phi)**2
    weights = pw * np.sin(phi)**(2*a-1) * np.cos(phi)**(2*b-1)
    weights /= weights.sum()
    return np.exp(np.outer(np.asarray(h)/2, v)) @ weights


def angular_probability(d, construction):
    theta1, tw1 = gl(96, 0, math.pi/4)
    theta2, tw2 = gl(96, math.pi/4, math.pi/2)
    theta = np.concatenate([theta1,theta2])
    a = d/2 if construction == 'grassmann' else 0.5
    b = (d+1)/2
    evidence = beta_evidence(4*np.sin(theta)**2, a, b)
    probability = np.dot(tw1,evidence[:len(theta1)]) / np.dot(np.concatenate([tw1,tw2]),evidence)
    return {'P_mass_below_one':float(probability), 'null_mean_scalar_beta':a/(a+b)}


def width_identity():
    v2 = np.array([9.0,1.0])
    sensitivity = np.array([1.0,9.0])
    c = v2.mean()
    b = c/(1+c)
    h = np.sum(sensitivity*b*v2)
    variance = np.sum(sensitivity*b*b*(v2+1))
    bg = v2/(1+v2)
    grouped_h = np.sum(sensitivity*bg*v2)
    grouped_variance = np.sum(sensitivity*bg**2*(v2+1))
    return {'shared_beta':float(b),'H':float(h),'V':float(variance),
            'V_over_H':float(variance/h),'grouped_H':float(grouped_h),
            'grouped_V':float(grouped_variance)}


def main():
    result = {
        'source_commit':'592d144531fa4a30d1c2229a4c3b528e84043d77',
        'one_reading_x2_y0': {
            construction:{str(d):angular_probability(d,construction) for d in (1,2,3)}
            for construction in ('grassmann','projective_half_cauchy')
        },
        'grassmann_prior_component_shrinkage_mean':{str(d):d/(d+1) for d in (1,2,3)},
        'projective_prior_component_shrinkage_mean':{str(d):0.5 for d in (1,2,3)},
        'heterogeneous_scalar_example':width_identity(),
        'independent_frame_change': {
            'd2_N2_T1_parallel_columns_det':3,
            'd2_N2_T1_orthogonal_columns_det':4,
            'conditional_prior_density_orthogonal_over_parallel':(4/3)**(-2),
        },
    }
    output = Path(__file__).with_name('proper_tensor_results.json')
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()

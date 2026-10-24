"""Scope of instrument-only covariance constructions and a nonunique escape.

This is a scratch proof, not an adopted equation or a project record.

A. WHAT UNITS AND DIRECT-SUM LOCALITY ACTUALLY FORCE
---------------------------------------------------
Let W(m,Sigma) be a covariance for physical acceleration candidates. Require
physical rotations, force/acceleration unit changes, channel swap, and
direct-sum LOCALITY across independent axes. Taking both channel unit factors
equal to epsilon leaves m fixed and forces

    W(m,epsilon^2 Sigma)=epsilon^2 W(m,Sigma).

Locality then applies this identity to any independent added noise block:
its W block scales epsilon^2 and the retained W block stays unchanged. Thus
the candidate covariance in that block's own noise units does not change.

For a maximum-entropy Gaussian family this also fixes the entire standardized
candidate law. Any mass-law or shared-size-learning effect of a zero coordinate
is therefore independent of epsilon and cannot fade as epsilon -> 0. For a
general non-Gaussian measure, covariance scaling alone does NOT fix the full
standardized law: higher-order shape could depend on global instrument ratios.
To extend the result, locality and scaling must constrain the full measure.

This DOES NOT force W proportional to K=(B' Sigma^-1 B)^-1. For example,

    W_H=K,
    W_A=(Sigma_aa + Sigma_FF/m^2)/4

both satisfy the stated rotations, units, swap and direct-sum requirements.
At m=1, Sigma_FF=diag(1,4), Sigma_aa=I, cross covariance zero,
W_H=diag(.5,.8), while W_A=diag(.5,1.25). They are not proportional.

B. A STRONGER EXPERIMENT-RECODING AXIOM DOES FORCE c K
----------------------------------------------------
Suppose W is constructed only from the Gaussian linear experiment (B,Sigma),
and is invariant under every invertible observation recoding
 (B,Sigma)->(A B,A Sigma A'),
as well as covariant under arbitrary parameter-coordinate changes B->B T.
Whiten Sigma and map the normalized columns of B to a canonical embedding.
At that embedding covariance must commute with every orthogonal parameter
rotation, so it is c I. Undoing the transformations gives W=c K. Direct sums
make c dimension-independent. This is stronger than physical rotations,
separate channel units and swapping; it must be named as an extra premise.

Maximum entropy conditional on W=c K selects a Gaussian, but it does not
choose c. A fixed c gives beta=c/(1+c) and reported local variance 1/(beta J).
Actual repeated-error variance remains (J+G)/J^2. Matching them requires
c=J/G, which depends on candidate strength. No fixed instrument-only c can
match all admissible candidate sequences. Learning c introduces the shared
scale's already-proved null-feedback problem.

C. AN ELEMENTARY GAUSSIAN CHARACTERIZATION, WITH EXTRA PREMISES
--------------------------------------------------------------
Suppose proper radial candidate densities are simultaneously projective under
dimension deletion and zero-section compatible: conditioning an exactly zero
coordinate returns the lower-dimensional candidate density. The latter makes
their shape the same function h(||v||^2), up to normalization, in all dimensions.
Projectivity says T h=C h, where T h(t)=int_R h(t+x^2)dx. Applying T twice and
using polar coordinates gives

    C^2 h(t)=pi int_t^infinity h(s)ds.

Differentiating yields h'=-(pi/C^2)h, so the density is Gaussian. This is a
statement about WHOLE candidate laws, not merely the mass output's known-empty
continuity. The latter alone is weaker, so this is not a no-go for the original
task unless these extra candidate-law premises are independently justified.

D. KNOWN-EMPTY CONTINUITY ALONE HAS A CONSTRUCTIVE ESCAPE
--------------------------------------------------------
Direct-sum locality is stronger than regrouping the SAME full experiment.
A global design functional can respect regrouping while comparing noise
levels across coordinates. An explicit example is

    D(Sigma)=.5*(Sigma_aa/lambda_max(Sigma_aa)
                 +Sigma_FF/lambda_max(Sigma_FF)).

It is positive-definite, rotation-covariant, invariant under channel unit
changes, and unchanged by channel swap. It fails direct-sum locality.
For an independent zero axis whose two SDs are both epsilon times fixed
values, D_hidden=O(epsilon^2), while D_retained is unchanged for sufficiently
small epsilon, because both largest eigenvalues remain in the retained block.

For ANY fixed p>0 define standardized candidates v=D^(p/2)u. With a proper,
projective prior on u, the hidden standardized candidate tends to zero, its
zero-reading likelihood tends a constant, and shared-scale feedback disappears.
Physical acceleration candidates a=K^(1/2)v shrink as epsilon^(1+p) on that
axis. A compatible starting-law example is

    pi_p(m) proportional sqrt(tr[T(m) D^p]),
    T=U_m'^T(I-U_mU_m')U_m'.

The hidden term in this law also disappears as epsilon -> 0. It is proper,
unit-covariant and reciprocal. Dominated convergence applies to the marginal
likelihood for a proper latent prior because the Gaussian integrated likelihood
is bounded. Conditional Gaussian candidate posterior first moments have the
usual K^(1/2) bound, giving totality and disappearance of the hidden physical
component. This exhibits a real escape from a UNIVERSAL known-empty no-go.

But D versus D^2 already produces two different admissible design rules.
Maximum entropy does not select p, the spectral comparison functional D, or
the mixing measure. They are new mathematical choices requiring a physical
premise; this construction is a counterexample, not the sought equation.

E. REGROUPING AND SEQUENTIAL-DESIGN LIMITS
-----------------------------------------
For regrouping the SAME M standardized coordinates, take D over the full
declared stacked design and use a scalar shared variance u|tau~N(0,tau I_M).
Every regrouping carries the same D and same tau, so the law is unchanged.
This can evade the continuity obstruction, but scalar-size learning still does
not solve arbitrary anisotropic width calibration.

A full tensor Omega on an original d-dimensional physical group can retain
the regular width benefit, but it introduces correlation/group labels which
must survive axis splitting. It is NOT automatically invariant under replacing
that model by independent scalar readings or by a new arbitrary tensor on a
different chunking. Claiming all properties at once would be incorrect.

The global D and pi_p also depend on the declared instrument design. Adding a
new noise level can change earlier candidate weighting, even before its data
are read. The T20 requirement that the design be specified in advance therefore
still matters. No assertion of a freely extendable sequential prior is made.
"""

import json
import numpy as np

if __name__=='__main__':
    ff=np.diag([1.,4.]); aa=np.eye(2); m=1.
    K=np.linalg.inv(m*m*np.linalg.inv(ff)+np.linalg.inv(aa))
    A=(aa+ff/(m*m))/4
    print(json.dumps({'harmonic_conditional_K':K.tolist(),
      'arithmetic_covariance_counterexample':A.tolist(),
      'axiswise_ratio':(np.diag(A)/np.diag(K)).tolist()},indent=2))

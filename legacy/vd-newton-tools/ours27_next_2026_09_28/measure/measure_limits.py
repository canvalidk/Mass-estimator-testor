"""Scratch derivations: totality, projection, and compact candidate geometry.

1. THERE IS NO MINIMAL DETERMINANT-POWER REPAIR
--------------------------------------------
For d-dimensional standardized candidate covariance Omega, take
 q_eta(Omega) = |I+Omega|^(-(d+1)/2-eta), Omega>0.
With N Gaussian readings, its posterior integral is a tilted matrix-beta
integral. It is finite iff eta+N/2>(d-1)/2. For the required N=1,d=3 domain,
eta>1/2, while eta=1/2 diverges logarithmically. Consequently no *smallest*
proper correction exists in this family. Choosing one exponent is an extra
premise, not a consequence of propriety.

2. A TOTAL PROJECTIVE IMPROPER MEASURE DOES EXIST, BUT IS NOT FORCED
----------------------------------------------------------------
The delta -> 0 boundary of the Dawid inverse-Wishart process with fixed
positive scale tau I has density
 q_0(Omega) = |Omega|^(-d) exp[-tau tr(Omega^-1)/2].
This is improper as a prior but proper after every N>=1, any d. To see it,
put A=Omega^-1. The symmetric-matrix Jacobian gives dOmega=|A|^(-(d+1))dA.
The posterior's determinant factor becomes
 |A|^(N/2-1) |I+A|^(-N/2) exp[-tau tr(A)/2],
times a positive data-dependent factor bounded above and below. Near the
zero-eigenvalue boundary, matrix-gamma integrability requires N>0. At infinity
the exponential controls the tail. tau>0 is essential; tau=0 reintroduces a
nonintegrable boundary as Omega -> 0, since the observation likelihood there
is a positive constant.

The family is projective as a sigma-finite limit of normalized, projective
inverse-Wishart measures; the dimension-dependent normalization factors have
finite positive ratios as delta -> 0. It yields an explicit candidate-space
measure on V in R^(d x N):
 dmu(V) proportional det(tau I+VV')^(-(N+d-1)/2) dV.
This is finite at the origin and improper only at infinity. Gaussian
observations then make the candidate posterior proper.

But neither tau=1 nor this inverse-Wishart functional form follows from
rotation, units, projection consistency, or totality. The counterfamily in
projective_counterfamily.py demonstrates the remaining freedom even after
fixing E[Omega]=I and requiring proper priors.

3. MATRIX CAUCHY / GRASSMANN GEOMETRY
------------------------------------
The proper delta=1,tau=1 inverse-Wishart process gives
 dmu(V) proportional det(I+VV')^(-(N+d)/2)dV.
This is the canonical volume in the graph chart [I_d;V'] of the real
Grassmannian G(d,d+N): its metric determinant is
 det(I+VV')^(-(N+d)). This selects a parameter-free matrix-Cauchy law *if*
uniform Grassmann geometry is independently stipulated. It is not a theorem
that the Gaussian observation metric calls for this compactification.

4. ONE JOINT MASS/SIZE GRAPH VOLUME
----------------------------------
For the actual aligned standardized candidate Z=U_m V, take its graph
[I_N;Z] in G(N,N+2d). Pull back the canonical Grassmann metric. Let
 T(m)=U_m'^T(I-U_mU_m')U_m'. The joint volume is

 dmu(m,V) proportional
 sqrt(tr[T VV'(I+VV')^-1]) det(I+VV')^(-(N+d)/2) dm dV.

The cross terms from a nonhorizontal frame cancel by the metric's Schur
complement; graph_volume_checks.py verifies this with correlated Sigma.
This measure simultaneously specifies the mass and size weights. It is
proper: the square-root factor is at most sqrt(tr T), whose integral over
positive m is finite for positive-definite instrument covariance, while
the matrix-Cauchy candidate factor is a proper density. Rotation, units and
channel swap respect it, with the correct coordinate Jacobians.

HOWEVER TWO STRUCTURAL FAILURES PREVENT ITS ADOPTION:

(a) It is not projective across readings. Already d=1, its candidate density
is proportional ||v||(1+||v||^2)^(-(N+2)/2). With N=1 it vanishes at v=0.
The first-coordinate marginal at zero for N=2 is proportional to
 int_R |w|/(1+w^2)^2 dw = 1, strictly positive. Thus the N=1 law is not the
marginal of the N=2 law. A fixed declared design could define it, but it is
not one coherent prior over an indefinitely extendable sequence.

(b) It fails the known-empty continuity requirement. One 1D reading F=2,a=0,
noise=(1,1), has ratio-of-mean readout 1.9898912729. Add a second zero axis,
both noise SDs epsilon, same ratio. The marginal mass law is independent of
epsilon>0 because standardization removes epsilon. Its point tends to
1.838380705 as epsilon -> 0, rather than the 1D answer. P(m<1) changes from
approximately .261275 to .312012, independently of epsilon. This is not a
metric/Jacobian omission: the full metric determinant passes independent
direct evaluation to 7e-16 relative error.

The all-null one-noise-pair test still gives mass s by channel symmetry. That
test is therefore insufficient to detect either failure.

5. ORDINARY HAUSDORFF VOLUME DOES NOT FIX THE WIDTH
------------------------------------------------
For one noise pair, aligned candidates in R^(2M) are (v sin(theta),v cos(theta)).
Their Euclidean induced volume is ||v|| d^Mv dtheta. Its marginal law is
flat evidence times E[||v|| | m,data]. In a regular growing-M experiment,
log of this factor changes curvature only O(1), whereas the missing width
term is O(M). It cannot supply the required first-order correction.
"""

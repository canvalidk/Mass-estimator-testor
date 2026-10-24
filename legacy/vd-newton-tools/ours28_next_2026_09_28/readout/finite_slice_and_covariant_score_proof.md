# Finite-law flat slices and covariant scores

Scratch proof, 2026-09-28. No production or managed-document edits.

Source scope: the four newly managed September 28 records supplied in
`../sources/`, at VD-docs@592d144531fa4a30d1c2229a4c3b528e84043d77, especially
`ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md`.
The September 6 composition uniqueness source was also fetched through the
GitHub connector at this same commit. Its blob is unchanged,
1bb77de5c16b1d1102b2e74866e5a4219ed10c53.

## Result

The finite-support flat-slice step can be proved directly for arbitrary
candidate spaces. No compact prize interval or monotonicity in a prize order
is needed. Continuous strict betweenness gives a fixed score at each mass
level. Readout covariance under reflections then makes that score exactly
invariant, without any assumption that the candidate law is isotropic.

With full stacked orthogonal recoding this derives C3 on finite laws. With
Wasserstein-1 continuity it derives C3 on all finite-first-moment laws.
The existing two composition axioms then force the ratio of mean norms.

Two qualifications are essential:

1. Per-reading frames/signs alone supply only block orthogonal covariance,
   leaving a per-reading size-profile invariant. They do not yet force C3.
2. The older continuity-at-constants condition C5 is insufficient for the
   extension from finite laws. A counterexample below proves this explicitly.

## 1. Abstract finite-law theorem

Let C be any set of candidates and let m:C->(0,infinity) be the mass of a
single candidate. Let P_fin(C) contain every finitely supported probability
law on C. Assume a readout R:P_fin(C)->(0,infinity) with:

* Point calibration: R(delta_c)=m(c).
* Strict betweenness: for 0<t<1, R((1-t)P+tQ) lies strictly between R(P)
  and R(Q) when those differ, and equals their common value otherwise.
* Finite-simplex continuity: for each fixed finite set of candidates, R is
  continuous in their probability weights, including boundary weights.

Fix r with candidates a,b satisfying m(a)<r<m(b).

**Theorem.** There exists a function phi_r:C->R, unique up to multiplication
by a positive constant, such that for every finitely supported P,

    R(P)<r iff integral phi_r dP<0,
    R(P)=r iff integral phi_r dP=0,
    R(P)>r iff integral phi_r dP>0.

In particular sign(phi_r(c))=sign(m(c)-r). The score is defined by the
readout's behavior on hypothetical laws, rather than by a particular observed
data set or a particular fitted candidate distribution.

### Proof on a finite simplex

Take a finite set S containing a,b and write Delta_S for its probability
simplex. Define U={P:R(P)<r}, W={P:R(P)>r}, and L={P:R(P)=r}.
Betweenness makes all three convex; continuity makes U and W relatively open.
Both have points in the relative interior of Delta_S: perturb the point
masses at a and b by sufficiently small fully supported mixtures.

Finite-dimensional separation of the disjoint convex relative interiors
provides a nonconstant affine function ell and a separating level, absorbed
into ell, with ell<=0 on U and ell>=0 on W. The inequalities are actually
strict everywhere on U and W. For example, if P is in U with ell(P)=0,
mix it slightly with an interior Q in W having ell(Q)>0. Continuity keeps
the mixture in U for sufficiently small mixing weight, whereas ell is
positive, a contradiction. The reverse argument handles W.

For P in L, every nontrivial mixture toward a is in U by strict betweenness,
and every mixture toward b is in W. Taking limits gives ell(P)<=0 and
ell(P)>=0. Thus ell(P)=0. Conversely a point with ell=0 cannot be in U or W.
Therefore L=Delta_S intersect {ell=0}, with the stated signs on either side.

The hyperplane cuts the relative interior of the simplex, since it separates
interior points of U and W. Its zero set consequently contains an open subset
of a codimension-one affine hyperplane. Any other separating affine function
with the same zeros is a positive multiple. An affine function on a
probability simplex is the expectation of its values at the vertices.

### Patching the scores across finite supports

Normalize every score by phi_r(a)=-1. Its value at b is fixed by the unique
t in (0,1) with R((1-t)delta_a+t delta_b)=r:

    phi_r(b)=(1-t)/t.

Existence follows from continuity. Uniqueness follows from strict
betweenness, which makes the readout strictly increasing along that segment.

For any candidate c, construct the score on {a,b,c}. On a larger finite set
containing these candidates, restriction has the same zero set and the same
normalization, so its value at c must agree. This defines phi_r(c) globally.
Every finitely supported P is covered by its support together with a,b,
proving existence and consistency of the representation.

For an elementary uniqueness check, compare any two proposed global scores.
The r-mixture of a and b fixes their relative scales. Mixing each additional
candidate with an anchor on the opposite side fixes its score uniquely;
candidates with mass exactly r have score zero. This also proves uniqueness
without an infinite-dimensional separation theorem.

### Why the strict form matters

The ordinary median is not a counterexample or an example of this theorem:
median(0.9 delta_1+0.1 delta_3)=1 violates strict betweenness. Quantile
conditions are generally sign inequalities, not a single expectation-zero
condition on every finite law. Weak betweenness alone does not give the flat
slice argument just proved.

## 2. Covariance of the readout forces invariance of its score

Let a group G act on C, preserving each candidate's mass. Suppose

    R(g#P)=R(P) for every g in G and every P in P_fin(C).

The function c->phi_r(g c) represents exactly the same r-level readout as
phi_r(c). Score uniqueness therefore gives a positive number a_r(g) with

    phi_r(g c)=a_r(g) phi_r(c) for every c.

These multipliers form a multiplicative character:
a_r(gh)=a_r(g)a_r(h). In particular, if g is an involution, then
a_r(g)^2=1. Positivity forces a_r(g)=1.

Consequently if G is generated by involutions, phi_r is exactly invariant
under G. Orthogonal groups are generated by reflections. Thus orthogonal
readout covariance forces an orthogonally invariant score. This proof needs
neither group averaging nor any continuity assumption on g->a_r(g).

This is the needed bridge from one common frame transformation of a whole
law to a per-candidate score that cannot exploit direction. The original
law-dependent direction weight is not a fixed score; strict betweenness is
what prevents that escape.

## 3. Application to aligned candidate pairs

For fixed stacked dimension M>=1, let

    C_M={(X u,Y u): X,Y>0, u in S^(M-1)}.

Here X and Y are physical force and acceleration norms in declared reference
units, and m(c)=X/Y. Suppose the group in section 2 is the whole O(M), acting
on both channels by the same orthogonal transformation.

Its orbits are exactly the pairs (X,Y). Therefore

    phi_r(X u,Y u)=psi_r(X,Y).

For a finite law P, each atom can now be moved to the same chosen direction
e1 WITHOUT changing its score. Although this is not one common rotation of
the entire law, score invariance makes the operation valid atom by atom.
It follows that

    R(P)=R(T#P),  T(X u,Y u)=(X e1,Y e1).

Any two laws with the same magnitude-pair law have identical images under T.
This proves C3 for every finite law.

Only a covariance requirement on the READOUT was used. At no point was
P=g#P, nor was a starting candidate measure required to be invariant.

## 4. Covariant recoding is not isotropy of the starting measure

Let the Gaussian experiment have observation y, candidate c, covariance
Sigma, and starting measure mu. For g in O(M), let G=diag(g,g). Recoding the
same experiment sends

    y -> Gy,   Sigma -> G Sigma G',   mu -> G#mu.

The posterior then sends P->G#P. This is true for an anisotropic or grouped
measure too: its full dependence structure is transported. Requiring the
scalar readout to agree before and after this recoding does not say that mu
equals G#mu. In particular it does not force one shared Gaussian scale over
all coordinates, nor invalidate a tensor or per-reading measure merely
because its canonical-coordinate density is not O(M)-invariant.

There is still a physical/domain premise to identify. The user has accepted
each reading's own frame and signs, not explicitly every rotation mixing
coordinates from different readings. If only block rotations are allowed,
section 2 yields block-invariant scores, not scores depending only on the
total norms.

Also, the readout must really be a function of the candidate-pair law alone.
If it receives partition projectors or other provenance as EXTRA inputs and
those too transform under G, a profile-dependent score can remain covariant.
The accepted instrument-independence statement excludes external instrument
matrices at the readout stage. Whether it also excludes a separate reading
partition is the precise remaining question. It need not be answered by
requiring the GENERATING MEASURE to forget that partition.

Thus the proposed operational route is:

    construct the law with all genuine experiment dependence preserved;
    then read the mass from its intrinsic aligned-pair geometry, covariantly
    under orthogonal recoding, without a second partition/provenance input.

This is weaker than assuming C3 outright: direction-law functionals were
initially allowed, and betweenness plus covariance have now excluded them.
It is stronger than only per-reading own frames; that gap is not hidden.

A concrete sufficient added operation is an invertible orthogonal digital
recombination of two scalar channels from different readings:

    (F_j,F_k) -> (cos(t) F_j+sin(t) F_k, -sin(t) F_j+cos(t) F_k),

with the identical recombination on the acceleration channels and the full
candidate measure/covariance transported. Force channels are only combined
with force channels, so independent force/acceleration obtainability is
preserved; F*=m a* is preserved candidate by candidate. Such plane operations,
together with the already accepted within-reading changes and signs,
generate full O(M). Calling this a representation change of the same
experiment, rather than a new set of independently fitted readings, supplies
the exact operational premise used by the theorem. If its original reading
partition must remain an additional readout input, this operation does not
erase that input and C3 still does not follow.

## 5. Exact survivor under only per-reading rotations

Partition u into reading blocks u_i and put

    h(u)=1+sum_i ||u_i||^4,
    R_h(P)=E[h(u) X]/E[h(u) Y].

Since sum_i ||u_i||^2=1, h lies between 1 and 2. It is unchanged by every
reading's own rotation or sign convention and by reading permutations.
It is homogeneous of degree zero, instrument-independent, and continuous
on the candidate space. Its numerator and denominator are linear in P,
so it obeys strict betweenness. On a common candidate direction, the same
h appears in C1 and C2, so both composition identities hold exactly. It fixes
all constant-mass candidate laws and has the required constant-limit
continuity. Nevertheless it distinguishes associations with per-reading
size profiles and therefore violates C3.

This proves that the user's currently listed symmetries plus betweenness
cannot alone close C3. Full stacked recoding at the readout stage, or another
principle with equivalent orbit-erasing power, is genuinely needed.

## 6. Extension to continuous candidate laws

Use the metric

    d(c,c')=||F-F'||+||a-a'||

in fixed reference units, and let P_1(C_M) be laws with finite first moments.
Assume R is continuous for Wasserstein-1 convergence on this domain. This
means continuity when candidate clouds converge weakly and their first
size moments converge, not continuity under arbitrary weak convergence with
uncontrolled outliers.

Every P in this domain has finite-support approximations P_n converging in
Wasserstein-1. The canonicalization map T of section 3 is Lipschitz: each norm
map is 1-Lipschitz. Thus T#P_n -> T#P in Wasserstein-1 too. Finite-law C3 and
continuity give

    R(P)=lim R(P_n)=lim R(T#P_n)=R(T#P).

Consequently C3 extends to every positive, finite-first-moment candidate law.
Zero candidate vectors, where their mass is not encoded by the pair, remain
outside this theorem and need the estimator's separate exact/null typing.

This continuity assumption also supplies finite-simplex continuity and the
older C5: L1 convergence of a positive magnitude to a constant can be coupled
on a fixed canonical direction to give Wasserstein-1 convergence.

## 7. Why older C5 cannot replace the extension assumption

Here is an exact counterexample, for M>=2. Given P, let nu be its direction
marginal on S^(M-1), and split nu into its atomic and nonatomic parts. Give
candidates whose direction belongs to an atom of nu weight 1, and give the
remaining direction component weight 2. Define

    R_atom(P)=E[w_P(u) X]/E[w_P(u) Y].

Although written with a P-dependent weight, both weighted integrals are
LINEAR in P. The atomic part of a mixture's direction marginal is the sum
of its components' atomic parts: a nonatomic measure assigns zero mass to
the countable union of their atoms. Therefore R_atom satisfies strict
betweenness and mixture continuity.

Every bijective direction transformation preserves the atomic/nonatomic
split. Hence it is invariant under the entire orthogonal group and under
the accepted frames/signs. It is instrument-independent. Common-direction
composition retains the same direction marginal, so C1 and C2 hold. It is
calibrated on constants and indeed on every constant-mass law. Finally,

    |E[w X]/E[w]-c| <= 2 E|X-c|

proves the old constant-limit C5, uniformly over changing direction laws.

Yet C3 fails. Give probability 1/2 to X=Y=1,u=e1 and probability 1/2 to
X=2,Y=1,u uniform on the circle. Then R_atom=(1/2+2)/(1/2+1)=5/3. If both
magnitude cases instead have direction e1, their magnitude law is identical
but the readout is 3/2.

Approximating the uniform circle by finitely many directions makes every
direction atomic, so each approximating law has readout 3/2, while its limit
has readout 5/3. The first magnitude moments are unchanged. Thus exactly the
Wasserstein/weak-cloud continuity required in section 6 fails.

This counterexample shows that the continuity gap is substantive. Finite-law
and even arbitrary-mixture continuity, C5, composition, and full readout
rotation covariance do not alone establish C3 for all continuous laws.

## 8. Completion: the unique norm score and readout

Once sections 3 and 6 establish C3, apply C1 and C2 to positive magnitude
variables. For one arbitrary joint law define

    r=R(X,Y), x=R(X,1), y=1/R(1,Y).

Compute R(X+1,Y+1) by splitting force first and acceleration first. The two
expressions differ by

    -(r y-x)^2 / [(x+r)(y+1)(y r+1+x y+y)].

The denominator is positive, so r y=x. Consequently R(X,Y)=T(X)/S(Y),
where T(X)=R(X,1) and S(Y)=1/R(1,Y). C1/C2 make T and S additive under
arbitrary couplings; C3 makes them depend only on each magnitude's law.

For iid copies X_i, additivity and rational homogeneity give
T((X_1+...+X_n)/n)=T(X). Their average converges in L1 to E[X], so C5 and
calibration give T(X)=E[X]. The same argument gives S(Y)=E[Y]. Therefore

    R(P)=E_P ||F*|| / E_P ||a*||.

The fixed level-r score is, up to a positive factor depending only on r,

    phi_r(c)=||F*||-r ||a*||.

Indeed this score and the score constructed in section 1 have identical
finite-law zero sets, so their proportionality follows from the proved score
uniqueness. Linear norm-score dependence is a consequence of composition
after the covariant-score step, not an initial assumed weighted-mean family.

## Exact advance and exact remaining premises

Proved here: the finite-support flat-slice theorem; score invariance from
covariance and involutions; finite-law C3 for full stacked recoding; extension
under first-moment cloud continuity; and the resulting unique norm score.

Still not supplied by the cited provisional rulings alone: full stacked
recoding as a symmetry of the LAW-ONLY READOUT, and the stated continuity
for continuous clouds. Neither premise requires isotropy of the constructed
candidate law. Neither chooses a size prior or solves its sharing problem.

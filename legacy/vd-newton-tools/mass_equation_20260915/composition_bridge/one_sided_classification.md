# One-sided additive mismatch: classification and the norm-exchange theorem

Independent proof audit, 2026-09-15. Scratch only; no production or managed
edits. All vectors below use a fixed common standardized Euclidean coordinate
convention. Exchanging force and acceleration is an additional symmetry of
this mathematical residual; standardizing their units does not itself
justify that symmetry physically.

## 1. Assumptions and complete one-sided classification

Let B:R^3 x R^3 -> R^3 be jointly continuous and suppose:

    B(F1+F2,a)=B(F1,a)+B(F2,a),                                  (A)
    B(QF,Qa)=Q B(F,a),  Q in SO(3),                              (R)
    B(t a,a)=0 for a!=0 and t>0.                                 (C)

Thus additivity is assumed only in the force argument, and vanishing is
initially required only for positive collinear pairs.

For a!=0, put r=||a||, u=a/r and P_u=I-u u^T. The complete class is

    B(F,a)=A(r) P_u F+D(r)(F cross u).                            (1)

Here A and D are arbitrary continuous real functions on (0,infinity)
subject to

    A(r)->0 and D(r)->0 as r->0.                                 (2)

At a=0 one must have B(F,0)=0. Conditions (1)-(2), together with that
definition, are also sufficient for all three assumptions and joint
continuity. In particular, one-sided additivity does not by itself force
the residual to have only a cross-product term.

### Proof

Continuity plus the additive Cauchy equation makes B linear over the reals
in its first argument. Therefore B(F,a)=M(a)F, with M(a) a continuous
matrix-valued function. Equivariance becomes

    M(Qa)=Q M(a)Q^T.

Fix a=r u. Every rotation about u fixes a, so M(a) commutes with all such
rotations. Decompose R^3 as the u axis and its perpendicular plane. There
can be no off-diagonal blocks: the plane has no nonzero vector fixed by
every rotation. The axis block is a scalar. The plane block commuting
with all planar rotations has the form A I+D J, where J sends F to
F cross u. These coefficients depend only on r because proper rotations
act transitively on the unit sphere.

Condition (C) gives M(a)a=0, so the axis block vanishes. This proves (1).
Force linearity also extends collinearity vanishing from positive to
negative multiples and to F=0.

For any fixed nonzero F, the pairs (F,tF), t>0, are positive collinear.
Continuity as t->0 therefore gives B(F,0)=0. Conversely, evaluating
M(r e3) on e1 shows that continuity at zero requires A(r),D(r)->0.
These conditions are sufficient because

    ||B(F,a)||<=sqrt(A(r)^2+D(r)^2) ||F||,

which tends uniformly to zero for bounded F as r->0.

### Residual norm and detection of collinearity

The two terms P_u F and F cross u are perpendicular and have the same
length. Define rho(r)=sqrt(A(r)^2+D(r)^2). Then

    ||B(F,a)||=rho(r)||P_u F||
             =[rho(r)/r] ||F cross a||.                         (3)

The residual is not identically zero exactly when rho is positive
somewhere. To detect every noncollinear pair with a!=0 requires the
stronger condition rho(r)>0 for every r>0. A merely nonzero residual
may vanish identically for some acceleration magnitudes.

## 2. Main theorem: exchanging residual norms selects area exactly

Add the norm-exchange condition

    ||B(F,a)||=||B(a,F)|| for all F,a.                            (N)

Then there is one constant k>=0 such that

    rho(r)=k r for every r>0,
    ||B(F,a)||=k ||F cross a|| for every F,a.                     (4)

If B is not identically zero, k>0. Thus the fixed-ball conditioning event
is already identical at every tolerance, not only asymptotically:

    {F dot a>0, ||B(F,a)||<epsilon}
      ={F dot a>0, ||F cross a||<epsilon/k}.                     (5)

For a supplied nonsingular Gaussian ambient probability law, the proven
area-conditioning limit consequently selects the flat-magnitude measure
up to this constant scale. Neither acceleration additivity nor a unique
oriented residual vector is needed for that conclusion.

### Proof

Choose orthogonal F and a with arbitrary positive lengths f and r.
Equation (3) and norm exchange give

    rho(r) f=rho(f) r.

Thus rho(r)/r is the same constant k at every positive radius. Substitute
back into (3) to prove (4). If k=0, (1) is identically zero; therefore
nontriviality is exactly what excludes the vacuous case.

This proof uses only orthogonal pairs from the full norm-exchange axiom.
Requiring (N) for all pairs is a natural sufficient statement, rather than
a claim that its domain is logically minimal.

### The remaining vector freedom

For k>0, the whole continuous norm-exchange family can be written

    A(r)=k r cos(phi(r)),
    D(r)=k r sin(phi(r)),                                        (6)

for a continuous angle lift phi on (0,infinity). Equivalently it is any
continuous pair A,D satisfying A^2+D^2=k^2 r^2; the angle description is
not needed for the theorem. Every such choice is continuous at a=0
because its norm is bounded by k||a||||F||, even when phi oscillates
without a limit as r->0.

The residual can rotate within the plane perpendicular to a while its
norm, zero set, and all fixed Euclidean-ball events remain identical.

## 3. Explicit nonbilinear, norm-symmetric examples

### Transverse projection example

For a!=0 define

    B_P(F,a)=||a||F-(F dot a)a/||a||,

and set B_P(F,0)=0. This is (1) with A(r)=r,D(r)=0. It is jointly
continuous, force-additive, proper-rotation equivariant, zero on collinear
pairs, and

    ||B_P(F,a)||^2=||F||^2||a||^2-(F dot a)^2
                 =||F cross a||^2.

It therefore satisfies norm exchange. It even has separate positive
degree-one scaling:

    B_P(sF,t a)=s t B_P(F,a),  s,t>0.

Nevertheless it is not acceleration-additive. For orthonormal e1,e2,

    B_P(e1,e1)+B_P(e1,e2)=e1,
    B_P(e1,e1+e2)=(e1-e2)/sqrt(2).

It is even in its second argument, B_P(F,-a)=B_P(F,a), and is not
channel-exchange antisymmetric. Thus neither norm exchange nor independent
positive degree-one scaling can be mistaken for full bilinearity.

### Radius-dependent rotation example

For r>0 take A(r)=r cos(log r), D(r)=r sin(log r). This gives a
continuous-at-zero residual with the same norm as F cross a, but changing
the acceleration magnitude changes its transverse rotation. It need not
obey positive scaling covariance in the second argument, despite satisfying
all assumptions of the norm-exchange theorem.

## 4. Stronger exchange antisymmetry recovers the cross product itself

Instead of (N), assume the vector identity

    B(F,a)=-B(a,F).                                              (S)

Force additivity then automatically yields acceleration additivity:

    B(F,a1+a2)=-B(a1+a2,F)
      =-B(a1,F)-B(a2,F)=B(F,a1)+B(F,a2).

Continuity makes B bilinear. The elementary SO(3) classification of
bilinear vector maps then gives

    B(F,a)=c F cross a                                           (7)

with one constant c. In this route, collinearity vanishing also follows
from antisymmetry and force linearity, rather than needing to be assumed.

One can alternatively see (7) inside the one-sided classification.
B(a,F) is perpendicular to F. Antisymmetry therefore forces B(F,a) to
be perpendicular to F as well, eliminating A(r)P_u F in (1). The
remaining expression is h(r)F cross a. Exchanging inputs at arbitrary
nonparallel lengths forces h(r)=h(f), so h is constant.

Thus exchange antisymmetry can replace a separate assumption of acceleration
additivity. It is still an additional modeling requirement, and is stronger
than the norm exchange sufficient to select the conditioning event.

## 5. Independent positive scaling: what it does and does not force

Give this phrase the explicit degree-one meaning

    B(sF,t a)=s t B(F,a),  s,t>0.                                (H)

The first factor is already supplied by force linearity. The second means
A(tr)=t A(r), D(tr)=t D(r), so

    A(r)=alpha r, D(r)=delta r                                  (8)

for constants alpha,delta. Conversely these forms satisfy (H). Therefore

    B(F,a)=alpha[||a||F-(F dot a)a/||a||]+delta F cross a

for a!=0, with value zero at a=0. Its norm is

    sqrt(alpha^2+delta^2) ||F cross a||.

Thus independent positive degree-one scaling is another sufficient added
condition for exactly the area tolerance and the flat limit. It does not
force the oriented residual to be a cross product: the projection example
is a counterexample.

Adding oddness in the second input, B(F,-a)=-B(F,a), removes alpha in
(8). Full signed degree-one scaling therefore gives (7). A scaling rule
with an unspecified degree or unspecified output factor is a different
assumption and does not have this conclusion.

## 6. Reflections and axial parity

For polar inputs F and a, the projection term in (1) transforms as a polar
vector, while the cross term transforms as an axial vector. Under all
orthogonal transformations the respective transformation rules are

    polar: B(QF,Qa)=Q B(F,a),
    axial: B(QF,Qa)=det(Q) Q B(F,a).

Proper rotations alone do not distinguish these parities.

- Axial O(3) equivariance forces A=0, leaving
  B(F,a)=h(||a||)F cross a. The function h can still vary; continuity
  at a=0 requires r h(r)->0, not that h itself have a finite limit.
- Polar O(3) equivariance forces D=0, leaving the transverse projection
  family. In particular, the norm-exchange projection example is a valid
  polar residual.

The inversion Q=-I proves these claims directly. In (1), sending both
inputs to their negatives changes the sign of the projection term and
leaves the cross term unchanged. Polar covariance demands the former
behavior; axial covariance demands the latter.

Axial parity together with norm exchange forces |h(r)|=k. Continuity
on the connected positive-radius domain and k>0 then force its sign to
be constant as well, giving (7). Axial parity together with (H) also
gives (7). Neither axial parity nor force additivity alone makes h constant.

## 7. Counterexamples locating the added assumptions

These demonstrate specific remaining freedoms; they are not a claim that
every premise above has been proved logically indispensable.

1. Without norm exchange or degree-one acceleration scaling,
   B(F,a)=(1+||a||^2)F cross a is continuous, force-additive, axial
   equivariant, and detects collinearity, but its residual norm is not
   a constant multiple of area. Thus even axial parity does not close
   the radial freedom.
2. Without force additivity,
   B(F,a)=(1+||F||^2+||a||^2)F cross a remains continuous, equivariant,
   collinearity-zero, and norm-symmetric, but again has a variable
   area multiplier.
3. Without rotational equivariance, L(F cross a) for a fixed nonconformal
   linear map L is bilinear and exchange-antisymmetric, with symmetric
   norm, but generally its norm is not a constant times Euclidean area.
4. Without collinearity vanishing,
   B(F,a)=(F dot a)a/||a|| for a!=0, extended by zero at a=0,
   is continuous, force-additive, SO(3)-equivariant and norm-symmetric.
   Its norm is |F dot a|, which is nonzero on positive collinear pairs.

For the general class (1), a fixed residual ball corresponds locally to
an area tolerance epsilon r/rho(r). The resulting interior reference
factor relative to flat positive magnitudes is [r/rho(r)]^2 in three
dimensions. Establishing a normalized global limit for arbitrary rho
requires positivity and integrability; no such unrestricted claim is
made here. Under the norm-exchange conclusion (4), this factor is simply
the constant k^(-2), so
the previously proved flat-limit theorem applies directly and globally.

## 8. Boundary of the claimed bridge

The strongest useful conclusion here is about the scalar discrepancy:

    continuous force-additive residual
    + proper-rotation equivariance
    + positive-collinearity vanishing
    + norm exchange
    + nontriviality
    => a fixed positive multiple of parallelogram area as its norm.

With a fixed Euclidean-ball tolerance and positive alignment, this selects
the original flat conditioning construction. The vector residual itself
need not be unique or bilinear.

Residual additivity is not automatically the existing readout composition
axiom. Norm exchange is not automatically implied by reciprocity of a
final mass point. A precise operational reason to impose either must still
be supplied. Gaussianity does not supply those premises, a fixed residual
ball, or an ambient probability law when the starting object is only a
likelihood. Additional directional arguments such as anisotropic covariance
tensors also change the equivariant classification and must be included
explicitly if they belong to the residual's inputs.

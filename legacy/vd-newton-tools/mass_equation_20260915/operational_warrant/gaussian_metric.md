# Gaussian metric, residual re-expression, and what it can warrant

Independent bounded proof audit, 2026-09-15. Scratch only; no production or
managed edits. Local context read: investigation/why_flat_conditioning.md.
The calculations concern one vector pair and introduce no population of
excitations across trials.

## 1. Distinguish the supplied cloud from a family of sampling experiments

Two probability objects must not be conflated:

1. A supplied ambient uncertainty law g over possible vector values,
   for example (X,Y)~N((U,V),I_6), centered at observed U,V.
2. A sampling family in which noisy observations at candidate true means
   (x,y) are (x+xi,y+eta), with xi,eta independent N(0,I_3).

The covariance of X cross Y under the supplied cloud in item 1 is one fixed
matrix determined by U,V. Evaluating a sampling covariance at every candidate
pair (x,y) instead defines a covariance field. Using that field to normalize
a discrepancy inside the cloud is an additional policy. It is not merely
recomputing the covariance of the same supplied cloud at each point.

The same distinction applies to a local Fisher metric: a sampling family can
define a metric at candidate states, but metric-based compatibility tolerance
is a specified construction, not an automatic conditional probability rule.

## 2. Linearized cross-product covariance at compatibility

At compatible true standardized means x u,y u, with ||u||=1, put

    r^2=x^2+y^2,   P_u=I-u u^T.

The first-order perturbation of B=X cross Y is

    delta B=x u cross eta+xi cross y u
           =u cross(x eta-y xi).

It is a Gaussian vector in the two-dimensional plane u-perpendicular, with

    Cov(delta B)=r^2 P_u.                                       (1)

Consequently its covariance-whitened norm is ||delta B||/r, and its squared
norm has a chi-square_2 distribution under this linearized sampling model.
This last statement is exact for the linearized random variable, not for
the nonlinear measured cross product at finite noise.

In physical coordinates with independent isotropic channel errors,

    Cov(delta(F cross a))
      =(f^2 sigma_a^2+alpha^2 sigma_F^2)P_u.                     (2)

Dividing the channels by sigma_F and sigma_a returns (1). The covariance has
rank two; whitening concerns its transverse image, not a nonsingular
three-dimensional Gaussian.

A symmetric extension of the resulting standardized scalar defect is

    D_lin(X,Y)=||X cross Y||/sqrt(||X||^2+||Y||^2),               (3)

with value zero at the common origin. This choice has the correct local
normal scale at every nonzero compatible pair. It is not an exact finite-
noise chi-square pivot. It is also not exactly Euclidean distance away from
compatibility, although it has the same first-order normal sensitivity.

### Conditional on one observed channel is a different covariance

For independent channels, condition instead on a realized Y. Then X cross Y
is exactly a rank-two Gaussian image of X, with covariance ||Y||^2 P_Y
in standardized coordinates. This conditions away Y's error rather than
propagating both errors near a compatible candidate. It is not (1), and
its conditional mean need not vanish under a compatible true pair when
the realized Y points away from the true direction. These two meanings
of 'conditional covariance' must be kept distinct.

## 3. What the local Gaussian metric selects

At a compatible pair (xu,yu), the Euclidean surface element of the
positive-collinearity cone is

    dVol_C=r^2 dx dy dOmega(u).                                  (4)

The derivative of B on the normal space is r times an isometry into
u-perpendicular. Its two singular values are r,r, giving normal Jacobian
r^2. Therefore:

- fixed raw cross-product balls have local normal area pi epsilon^2/r^2;
  multiplying (4) leaves the flat reference;
- fixed balls after the local covariance whitening (1) have local normal
  area pi epsilon^2; multiplying (4) retains r^2.

The scalar defect (3) gives the global limiting law

    g(xu,yu) r^2 dx dy dOmega(u),                                (5)

after normalization. Thus the Gaussian normal metric provides a route to
the previously studied tube measure. It does not provide a second proof
of flat weighting.

This limit is global, including the axes and origin. For a general radial
width h(x,y), an area event

    X dot Y>0, ||X cross Y||<epsilon h(||X||,||Y||)

has an allowed angular cap with argument epsilon h/(xy). The uniform
cap bound A_3(q)<=C min(q^2,1) gives, after multiplying by x^2y^2 and
dividing by epsilon^2, an envelope proportional to

    h(x,y)^2 exp[-c(x^2+y^2)].

For h=r this is integrable. Dominated convergence then supplies (5), not
merely a calculation away from the cone's singular strata.

## 4. Coarea and changes of independent defining functions

For a smooth codimension-p constraint C, let b take values in R^p and be a
regular local defining function: b=0 on C, and its derivative J=Db has
rank p there. Under an ambient density g and ordinary Euclidean balls,
the coarea/local-normal calculation gives the surface density

    g/J_b dH_C,
    J_b=sqrt(det(J J^T)).                                       (6)

Suppose another defining function is b_tilde=A(z)b, with A(z) invertible
on C. On C, Db_tilde=A J, because the derivative of A is multiplied by
b=0. Resetting the tolerance to a unit Euclidean ball in the new residual
coordinates changes (6) by

    1/|det A(z)|.                                               (7)

Thus state-dependent re-expression can change the conditioned law even
though it does not change the exact constraint.

There are several distinct invariances:

- A true coordinate change carries the tolerance metric or set along with
  it; then the underlying events, and hence the conditional laws, are
  identical.
- A fixed nonsingular transformation of p independent residual coordinates
  has constant determinant, so even if a fresh Euclidean ball is used, the
  limiting factor in (7) is constant and cancels in normalization.
- A fixed smooth residual-space diffeomorphism b_tilde=phi(b), with phi(0)=0
  and nonsingular derivative at zero, has the same constant-factor property.
- A common strictly increasing scalar transformation of a nonnegative norm,
  such as replacing area by squared area, merely relabels the small-radius
  events and has the same limit. This does not require a nonzero derivative
  at zero.
- A state-dependent scalar multiplier generally changes (7). Re-expressing
  a constraint and silently resetting its tolerance are different actions.

### Whitening restores local re-expression invariance

Let Sigma be a fixed positive-definite ambient error covariance. The
linearized residual covariance at a candidate point of C is

    C_b=J Sigma J^T.

The quadratic discrepancy b^T C_b^(-1)b is invariant on the linearized
normal space under b->A b and C_b->A C_b A^T. A positive-definite covariance
therefore fixes a whitened norm up to an orthogonal output rotation.

Small balls in this norm give the surface density

    g sqrt(det(J Sigma J^T))/sqrt(det(J J^T)) dH_C.              (8)

Both determinants acquire the same factor |det A|, so (8) is invariant
under the local re-expression. For Sigma=I it is simply g dH_C: the
Euclidean tube measure. For general constant Sigma it is, up to a constant,
g times surface volume in the ambient Mahalanobis metric Sigma^(-1).
Indeed

    sqrt(det(J Sigma J^T))/sqrt(det(J J^T)) dH_C
      =sqrt(det Sigma) dH_(Sigma^-1).

If the re-expression varies off C, recomputing its full derivative includes
terms involving b times derivatives of A. Exact off-constraint invariance
requires explicit metric transport; covariance propagation gives the
stated invariance on C and of the limiting law.

For the physical compatible parameters (f,alpha,u), the Mahalanobis surface
element under independent isotropic channel covariances is proportional to

    (f^2/sigma_F^2+alpha^2/sigma_a^2) df d alpha dOmega

in three dimensions. The omitted factor 1/(sigma_F sigma_a) is constant.
This agrees with the standardized tube weight (5).

### The cross product is a redundant, critical-valued residual

B maps R^6 into R^3, but at nonzero compatible pairs its derivative has
rank two. One must not insert det(DB DB^T)=0 into a regular three-dimensional
coarea formula. Use two independent residual coordinates in a local frame
of u-perpendicular, or the two-dimensional normal Jacobian as above.

The zero cross product is a critical value of this three-component map.
For example, its small-ball probability is of order epsilon^2, not
epsilon^3. A regular conditional law at the single null value B=0 is not
uniquely specified by an almost-everywhere conditional-density convention;
the shrinking sets still need to be stated.

This redundancy also matters when a fixed nonorthogonal 3-by-3 output map
is applied to B and a fresh ordinary ball is used. Its determinant on the
varying transverse plane can depend on u, despite the map being constant
in R^3. That need not preserve the limit. Transporting the output metric
does preserve the events; transformed covariance whitening preserves the
normal norm as described above.

## 5. Exact centered sampling covariance of the cross product

For arbitrary fixed standardized means X,Y, let

    B_obs=(X+xi) cross (Y+eta),
    xi,eta independent N(0,I_3),
    R^2=||X||^2+||Y||^2.

Then E[B_obs]=X cross Y. Expanding the centered error gives

    B_obs-E[B_obs]=X cross eta+xi cross Y+xi cross eta.

The first two terms have covariances ||X||^2 I-X X^T and
||Y||^2 I-Y Y^T. The last has covariance 2I: each component is a
difference of two independent products of standard normals. Cross
covariances vanish by independence and centered odd Gaussian moments.
Therefore the exact centered sampling covariance field is

    C(X,Y)=(R^2+2)I-X X^T-Y Y^T.                                (9)

The first-order part is positive semidefinite, so C>=2I everywhere.
At a compatible pair this is

    C(xu,yu)=(r^2+2)P_u+2u u^T.                                 (10)

It differs from the rank-two linearized covariance (1), because the
quadratic error term produces a nonzero fluctuation in every direction.
In particular, the exact nonlinear cross product is not a Gaussian
vector merely because its input errors are Gaussian.

For independent physical errors with isotropic variances, the counterpart is

    C_phys(F,a)=sigma_a^2(||F||^2 I-F F^T)
               +sigma_F^2(||a||^2 I-a a^T)
               +2 sigma_F^2 sigma_a^2 I.

## 6. Using the exact covariance as a candidate field is another law

Suppose one explicitly chooses (9), evaluated at each candidate pair,
as a covariance field for measuring its cross product B=X cross Y.
Since B is perpendicular to both X and Y,

    C(X,Y)B=(R^2+2)B,
    B^T C(X,Y)^(-1) B=||B||^2/(R^2+2).                         (11)

These identities are exact off the compatible cone as well as on it.
Thus a fixed ball in that field norm is exactly the event

    X dot Y>0,
    ||X cross Y||<epsilon sqrt(||X||^2+||Y||^2+2).               (12)

For any supplied nonsingular Gaussian density g, the global limit of (12)
has density relative to flat magnitude coordinates

    g(xu,yu)(x^2+y^2+2) dx dy dOmega(u),                         (13)

after normalization. In particular, its normalizing probability is

    pi epsilon^2 integral (x^2+y^2+2)g(xu,yu)dx dy dOmega
      +o(epsilon^2).

The proof uses the cap domination in Section 3 with h=sqrt(x^2+y^2+2).
The envelope is a polynomial times Gaussian decay, so it is integrable
over all positive magnitudes. The axes and origin supply no additional
limiting component. This is a bounded global argument, not just a formal
surface weight.

Equation (13) is distinct from both flat and tube conditioning. The
added 2 remains at zero excitation. With standardized error variance
tau^2 rather than one, the covariance is

    tau^2(R^2 I-X X^T-Y Y^T)+2 tau^4 I,

and the analogous weight, after removing constant factors, is
R^2+2 tau^2 on the cone. The local metric weight R^2 is recovered as
the leading small-noise construction at nonzero compatible points.

Using exact covariance does not make (11) an exact chi-square statistic.
It supplies second moments of a non-Gaussian residual. Moreover, using
candidate means to define the field is an added policy. Neither exact
normality, exact coverage, nor statistical optimality has been proved.

## 7. Using the covariance of the supplied cloud instead

If the input cloud is N((U,V),I_6), its cross-product covariance is the
single matrix

    C0=(||U||^2+||V||^2+2)I-U U^T-V V^T.

It does not vary as candidate (X,Y) values are integrated. Conditioning
on B^T C0^(-1)B<epsilon^2 generally produces a direction-dependent
flat-coordinate weight, not (13).

Let E_u be any orthonormal basis matrix for u-perpendicular. The small
residual ellipse in that plane has area multiplier

    rho(u)=1/sqrt(det(E_u^T C0^(-1) E_u))
          =sqrt(det C0)/sqrt(u^T C0 u).                          (14)

The second equality is the principal-minor/Schur-complement identity in
the orthonormal basis [E_u,u]. Thus the limiting density is

    g(xu,yu) rho(u) dx dy dOmega(u).

All eigenvalues of C0 are positive, so the event is bounded between two
constant area tolerances and the global Gaussian domination applies.
At U=V=0, C0=2I and rho is constant: this fixed-covariance construction
reduces to the flat law. It does not reduce to the candidate-field law
r^2+2. This is a concrete distinction between the two probability uses.

## 8. Consequences for the additive-balance bridge

Choosing a linear elimination map L_Y with L_Y Y=0 makes L_Y(X-mY)=L_Y X
mass-free and additive in X at fixed Y. That is a coherent way to preserve
the balance error's additive structure. The earlier norm-exchange theorem
can then select raw area as the residual norm under its other premises.

It does not uniquely select a physical tolerance. In particular:

- covariance whitening can leave the raw additive residual unchanged while
  assigning a state-dependent tolerance to its norm;
- the normalized residual need not itself remain additive in X, because
  its local covariance or reference radius depends on both inputs;
- insisting on additivity of the normalized discrepancy is an additional
  requirement, different from preserving additivity of the raw balance
  error before statistical scaling;
- a fixed tolerance in algebraic residual units, a fixed tolerance in
  local Gaussian normal units, and a fixed tolerance under an exact
  covariance field are different policies.

The local covariance metric has a definite mathematical advantage if the
chosen requirement is invariance under equivalent regular defining
functions: covariance and derivative Jacobians cancel in (8). That is a
conditional structural warrant, not a proof that the resulting mass
uncertainty is universally best. The Gaussian model can support several
explicit constructions; choosing which error quantity should have a fixed
tolerance remains part of the operational specification.

# Does physical mass composition force an additive off-constraint residual?

Independent scratch audit, 2026-09-15. No production or managed edits.

## Verdict

Even granting the project's composition requirements in their full stated
strength, they constrain the readout of an already supplied probability law.
They do not entail bilinearity of a diagnostic evaluated on incompatible
vectors, and do not select the probability reference used to construct that
law. Two additional bridges would be needed:

1. A reason that the off-constraint diagnostic must preserve both additive
   vector operations as an oriented residual.
2. A reason that its small-tolerance neighborhoods have a state-independent
   scale in a fixed norm.

The cross-product theorem is sound conditional on its assumptions. The issue
is whether physical/readout composition independently requires those new
assumptions. The explicit countermodels below show that it does not do so
as presently stated.

## 1. What the managed composition axioms actually say

For every appropriate propagated joint probability law, the readout satisfies

    R(f1+f2,alpha)=R(f1,alpha)+R(f2,alpha),
    1/R(f,alpha1+alpha2)=1/R(f,alpha1)+1/R(f,alpha2).

The rigorous uniqueness note expressly restricts these statements to one
propagated law. It does not identify that law with a fresh fit using another
measurement or another reference assignment. Its conclusion is

    R_P(f,alpha)=E_P[f]/E_P[alpha].

This formula obeys both requirements for every allowed P, regardless of
whether P arose from flat magnitudes, metric-tube weighting, angular
conditioning, or another normalizable construction. Expectation additivity
does the work; the density of P cancels out of the proof of the identities.

The source's force-lumping interpretation assumes the parts share an
acceleration and hence a compatible direction. The inverse-mass operation
likewise uses its specified common force and direction/sign convention.
Those operations do not automatically extend to arbitrary incompatible
observed vectors.

## 2. Physical composition does not constrain a residual away from its zero set

If F1=m1 a and F2=m2 a, then F1+F2=(m1+m2)a. Any diagnostic that vanishes
on the compatible set obeys

    D(F1+F2,a)=0=D(F1,a)+D(F2,a)

on these exactly compatible inputs. This is true whether D is linear,
quadratic, normalized by magnitude, or otherwise nonlinear off the set.
The common-force statement has the same limitation.

Thus a theorem about exact physical composition only sees zeros of the
diagnostic. Inferring additivity on all incompatible vectors is an additional
extension, not an algebraic consequence of those zeros.

## 3. Explicit alternative measure and residual

In fixed dimensionless Euclidean coordinates, define

    B_tube(X,Y)=(X cross Y)/sqrt(||X||^2+||Y||^2),

with value zero at the origin. This map is continuous everywhere because
||X cross Y||<= (||X||^2+||Y||^2)/2. It is rotation equivariant, has the
appropriate axial transformation under reflections, changes sign on channel
exchange, and has the same collinearity zero set as the cross product.
The separate positive-dot condition distinguishes the positive branch.

It is not additive. With X1=X2=e1 and Y=e2,

    B_tube(X1+X2,Y)=2 e3/sqrt(5),
    B_tube(X1,Y)+B_tube(X2,Y)=sqrt(2) e3.

At a regular compatible pair with r=sqrt(x^2+y^2), an orthonormal normal
displacement n satisfies area=r||n|| and total ambient squared norm
r^2+||n||^2. Consequently

    ||B_tube||=r||n||/sqrt(r^2+||n||^2) approximately ||n||.

A constant small tolerance therefore has the same local weighting as a
metric tube: rho proportional to r^(d-1). The separately proved literal
distance-tube construction already supplies a global alternative law

    dP_tube proportional to g(xu,yu) r^(d-1) dx dy dOmega.

For the normalized-cross-product band itself, away from the small-origin
region its allowed normal radius is

    epsilon*r/sqrt(r^2-epsilon^2),

which tends to epsilon and is at most sqrt(2)epsilon when r>=sqrt(2)epsilon.
The remaining origin region has negligible ambient probability at the tube's
epsilon^(d-1) scale. This gives the same limit under the Gaussian envelope.

Now retain R_P=E_P[f]/E_P[alpha]. Both readout composition axioms still hold
for all propagated joint laws, and exact fixed-mass inputs still return that
mass. The alternative measure and nonadditive residual therefore do not
violate physical mass composition or the readout theorem.

This is a counterexample to the proposed implication, not a recommendation
to choose the tube policy. A comparison of independently reconstructed laws
for parts and lump is a different claim and cannot be substituted for C1.

## 4. Scalar squared discrepancy exposes a representation issue

Let B0=X cross Y and A=||B0||. The scalar discrepancy A^2 is nonnegative,
continuous, and nonadditive. Nevertheless,

    {A^2<epsilon}={A<sqrt(epsilon)}.

It produces exactly the same limiting flat law. The vector recoding
B1=||B0|| B0 is also nonadditive but preserves the same tolerance events
after that scalar reparameterization. Thus residual additivity is not even
necessary for the resulting conditioning law; it is a property of a chosen
oriented representation of the mismatch.

More strongly, a nonnegative scalar discrepancy additive in each argument
on the whole signed vector space must be identically zero:

    D(X,Y)+D(-X,Y)=D(0,Y)=0.

The bilinear axiom must therefore refer to a signed/oriented vector residual,
not generically to error magnitude, squared error, or statistical severity.
This distinction prevents a false transfer from Gaussian quadratic error
scoring to additive residual values.

## 5. Bilinearity alone leaves a tolerance-scale freedom

Keep the bilinear B0 exactly, but consider

    ||B0(X,Y)|| < epsilon h(X,Y),   X dot Y>0.

At a regular compatible pair, for h continuous and positive there, the local
normal radius is epsilon h(xu,yu)/r. Multiplication by the surface Jacobian
r^(d-1) gives local magnitude weight

    rho_h(x,y,u) proportional to h(xu,yu)^(d-1).

Global convergence and normalizability still require control of h near axes
and at infinity; they must not be assumed for arbitrary h. The already
established examples have the required control:

    h=1       -> flat magnitude weight,
    h=r       -> metric-tube weight,
    h=xy      -> angular weight.

Hence the cross-product uniqueness theorem selects the flat reference only
when paired with the independent choice of a constant tolerance in its
fixed Euclidean residual norm. Allowing a state-dependent error scale
reintroduces measure freedom even while the residual remains bilinear.

Equal absolute cross-product tolerances are not automatically equal Gaussian
significance tolerances: the sampling area variance changes with excitation,
as the previous single-pair contribution derives. That fact does not select
the tube either; it shows why the neighborhood rule needs its own warrant.

## 6. Explicit lumping reference calculation

Suppose the parts have positive force magnitudes f1,...,fN, one common
acceleration magnitude alpha, and one common direction u. Consider the
locally finite reference measure

    dPi_parts=(product_i df_i) d(alpha) dOmega(u).

The shared acceleration and direction appear once. This is not a product
of independently normalized single-pair priors or a probability law.

For N=2, put s=f1+f2 and t=f1/s, with s>0 and 0<t<1. The Jacobian is

    df1 df2 = s ds dt.

Integrating out the unobserved split t yields s ds, not flat ds. Generally,
put f_i=s t_i on the positive unit simplex. Its coordinate volume is
1/(N-1)!, and the transformation has Jacobian s^(N-1), so

    pushforward(dPi_parts)
      = [s^(N-1)/(N-1)!] ds d(alpha) dOmega(u).                      (1)

This is a valid reference-measure identity even though the total measure is
improper: preimages of compact ranges of s have finite simplex volume.

If the likelihood depends only on (s,alpha,u), marginalizing the part split
adds exactly the factor in (1). With part-specific evidence, the likelihood
must additionally be integrated over the simplex; it cannot be dropped.

Thus flat references at the parts do not push forward to a flat reference
for a freshly fitted lump. Larger totals have more reference volume of
possible decompositions. This is a counting property of the selected
reference, not additional measured information about the mass.

The result reinforces the source's distinction: C1 holds for the probability
law actually propagated from the parts, while a new flat-reference fit for
the total starts with another law. Demanding equality between them would be
an additional projective-consistency requirement. It is not already C1,
and (1) shows that a naive flat-at-every-level prescription would fail it.

## 7. A precise statement the investigation can support

The strongest justified current statement is:

> If the reconciliation diagnostic is additionally required to be a
> continuous bilinear rotation-equivariant oriented vector, it is a scaled
> cross product. If compatibility is additionally approached through fixed
> Euclidean balls in that residual, the flat magnitude law follows. The
> existing physical and readout composition commitments do not, by themselves,
> establish either additional requirement.

The next useful investigation is to specify an operation on uncertain,
possibly incompatible inputs that would require those extra properties.
It must distinguish propagating an already supplied law from refitting new
evidence, and residual-vector algebra from the tolerance used to assign
probability volume. No new across-trial excitation model is needed.

## Provenance and coverage

Managed sources reread through the connector at
canvalidk/VD-docs@8782a293297330f7a354d87ebe62792018fa8866:

* Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md,
  especially section 4's physical interpretation and propagated-law caution,
  and section 5's explicit remaining measure choice.
* Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md,
  especially its domain, C1/C2, and conditional-on-P scope.

Local context was the existing investigation README and
why_flat_conditioning.md, plus the already-audited tube and area constructions.
This is a bounded logical-bridge audit, not a new corpus or literature review.

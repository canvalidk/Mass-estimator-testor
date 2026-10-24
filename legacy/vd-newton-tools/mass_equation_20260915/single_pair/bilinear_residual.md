# A bilinear oriented mismatch selects area conditioning in three dimensions

Independent bounded proof audit, 2026-09-15. Scratch only. No other files
edited. This is an additional residual-and-tolerance model; it is not a
consequence of the existing mass-readout composition theorem or Gaussianity.

## 1. Precise assumptions and conclusion

Let B:R^3 x R^3 -> R^3 satisfy:

1. Continuity.
2. Separate additivity:

       B(X1+X2,Y)=B(X1,Y)+B(X2,Y),
       B(X,Y1+Y2)=B(X,Y1)+B(X,Y2).

3. Proper-rotation equivariance, with no additional directional arguments:

       B(QX,QY)=Q B(X,Y),  Q in SO(3).

Then for one fixed scalar c,

       B(X,Y)=c (X cross Y).                                     (1)

If B is not identically zero, c!=0. A nontriviality requirement is essential:
the zero map satisfies all three assumptions and also vanishes on every
collinear pair, but cannot impose collinearity through a small tolerance.

The condition that B vanish on collinear pairs is automatic under assumptions
1-3. It may be stated to explain the intended interpretation, but is not an
independent axiom of this characterization.

## 2. Elementary proof

### Separate additivity plus continuity gives real bilinearity

Fix Y. Additivity first gives B(nX,Y)=nB(X,Y) for integers n, then
B(qX,Y)=qB(X,Y) for rationals q. Approximate any real scalar by rationals
and use continuity. Do the same in the second argument. Thus B is a real
bilinear map.

### Proper rotations already force the diagonal to vanish

For X!=0, every rotation fixing X must also fix B(X,X). The vectors fixed
by all rotations about the X axis form exactly that axis, so B(X,X) is
parallel to X.

Choose a half-turn Q about an axis perpendicular to X. Then QX=-X.
Bilinearity and equivariance give

       B(X,X)=B(-X,-X)=Q B(X,X)=-B(X,X).

Consequently B(X,X)=0. This also holds at X=0. Polarizing the identity,

       0=B(X+Y,X+Y)=B(X,Y)+B(Y,X),

shows that B is alternating. Therefore B(X,tX)=0 for every real t, including
the antiparallel case.

### One basis value determines the map

Take a positively oriented orthonormal basis e1,e2,e3. The half-turn
Q=diag(-1,-1,1) sends (e1,e2) to (-e1,-e2), so

       B(e1,e2)=B(-e1,-e2)=Q B(e1,e2).

Its fixed subspace is the e3 axis. Write B(e1,e2)=c e3.
The proper cyclic rotation e1->e2->e3->e1 then gives

       B(e2,e3)=c e1,
       B(e3,e1)=c e2.

Alternation fixes the reversed pairs, and diagonal values are zero.
Bilinearity now gives

       B(X,Y)=c[(X2 Y3-X3 Y2)e1
                +(X3 Y1-X1 Y3)e2
                +(X1 Y2-X2 Y1)e3],

which is (1). Conversely, every constant multiple of the cross product
satisfies the assumptions. No representation-theory result is needed.

## 3. Reflection and axial-vector qualification

Under proper rotations, polar vectors and axial vectors have the same
transformation rule. Under arbitrary orthogonal transformations Q,

       (QX) cross (QY)=det(Q) Q(X cross Y).                       (2)

For polar input vectors such as force and acceleration, the nonzero residual
in (1) is therefore an axial vector: it represents an oriented plane via
the three-dimensional cross-product convention.

If one additionally demands ordinary polar-vector covariance under all of
O(3), meaning B(QX,QY)=Q B(X,Y) even for reflections, only the zero map
survives. Set Q=-I: bilinearity gives B(-X,-Y)=B(X,Y), while that proposed
transformation rule gives -B(X,Y).

The consistent full-orthogonal rule for a nonzero residual is instead

       B(QX,QY)=det(Q) Q B(X,Y).

Its Euclidean norm is reflection-invariant. Thus the scalar tolerance
||B||<epsilon has no handedness preference even though an oriented vector
representation must keep the axial qualification.

## 4. Consequence for the conditioning measure

Add two explicit conditions to the residual characterization:

- B is nonzero;
- use a fixed Euclidean ball ||B(X,Y)||<epsilon as the mismatch tolerance,
  with its radius independent of X and Y.

Equation (1) then makes the event exactly

       ||X cross Y||<epsilon/|c|.

To select the positive-mass branch, also impose X dot Y>0. Bilinearity
cannot distinguish parallel from antiparallel collinearity: both give
zero residual. This sign restriction is a separate necessary choice for
the positive compatible family.

For an explicitly supplied nonsingular Gaussian ambient density g, the
area-conditioning theorem in `area_conditioning_proof.md` gives

       P(X dot Y>0, ||B||<epsilon)
          ~ pi (epsilon/|c|)^2 Z_flat,

and the normalized limiting law

       g(fu,alpha u) df d alpha dOmega(u)/Z_flat,
       f,alpha>0.                                                (3)

The nonzero scalar c changes the tolerance's units or scale, and possibly
the orientation sign; it cancels from the normalized limit. Therefore these
additional residual and acceptance assumptions select the flat-magnitude
reference within the specified three-dimensional class.

## 5. What the characterization warrants

This can provide a concrete reason to consider the original reference:
if the mismatch is required to add separately when either input is added,
to transform properly under spatial rotations, and to be compared against
a fixed ball in its own output space, then the mismatch is an oriented
area and the resulting limit is flat in the magnitudes.

Each requirement in that argument is an extra modeling assumption.
In particular:

- Separate additivity of a residual B(X,Y) is not the mass-readout
  composition theorem's additivity of latent summaries or probability
  measures. A connection between them would need a separate argument.
- Gaussian measurement errors do not require an additive mismatch
  residual or a constant tolerance in that residual.
- Dividing the cross product by a magnitude-dependent scale, as in angular
  discrepancy or normalized normal distance, generally breaks separate
  additivity. Such choices are excluded by this added residual axiom,
  rather than disproved as conditioning constructions.
- The map is assumed to have only X and Y as directional inputs. If an
  anisotropic covariance, external tensor, or preferred direction is an
  additional transforming argument, other equivariant constructions can
  exist. For example, a covariance-dependent linear map applied to
  X cross Y need not be a scalar multiple at a fixed anisotropic design.
- A variable radius or a different acceptance geometry is another
  assumption and can change the induced measure. The fixed-ball premise
  is part of the result.

The conclusion is conditional: these specified additional axioms select
area conditioning up to scale. It establishes neither universal optimality
nor that a physical measurement protocol must use those axioms. If g began
only as a measurement likelihood, the earlier distinction between a
likelihood and a supplied ambient probability law also remains in force.

# Flat-magnitude law from shrinking parallelogram area: rigorous proof

Independent proof audit, 2026-09-15, by the proof-audit subagent. Scratch only.
This uniquely named file preserves the proof that previously occupied
`area_conditioning.md`; that shared filename now contains a separate audit
and is not edited here. The shorter exact-coordinate proof below also audits
the root agent's proposed simplification.

## 1. Statement and scope

Let (X,Y) in R^d x R^d have an explicitly supplied nonsingular joint Gaussian
probability density g, with arbitrary finite mean and positive-definite joint
covariance. Independence or isotropic covariance is not needed. The Euclidean
coordinates used in the constraint are fixed. For d>=2 define

    E_epsilon = {X dot Y>0, ||X wedge Y||<epsilon},
    ||X wedge Y|| = sqrt(||X||^2 ||Y||^2-(X dot Y)^2).

Write p=d-1 and let kappa_p be the volume of the unit ball in R^p. Put

    Z_flat = integral_S^(d-1) integral_0^infinity integral_0^infinity
                 g(xu,yu) dx dy dOmega(u).

Then 0<Z_flat<infinity and

    P(E_epsilon) / epsilon^p -> kappa_p Z_flat.                    (1)

The conditional ambient laws converge weakly to the compatible probability
measure

    g(xu,yu) dx dy dOmega(u) / Z_flat,  x,y>0.                     (2)

Thus the limiting reference is flat in the two positive magnitudes, with
ordinary sphere-area direction measure. In d=2 the leading factor in (1)
is 2 epsilon; in d=3 it is pi epsilon^2.

The result is a precise conditioning construction. It does not show that
the Gaussian law, or exact compatibility alone, uniquely chooses an area
tolerance over another approach to the same probability-zero event.

## 2. Short proof by an exact Cartesian slice

The event X=0 has probability zero. Write

    X=xu,   x>0, u in S^(d-1),
    Y=a u+w,   a real, w in u-perpendicular.

Here a is the component of Y parallel to X, not its magnitude before taking
the limit. Orthogonal decomposition gives the exact identities

    X dot Y=xa,
    ||X wedge Y||=x||w||,
    dX dY=x^p dx dOmega(u) da dw.                                (3)

The measure dw is Euclidean Lebesgue measure in the normal fiber
u-perpendicular. It needs no preferred orientation of that fiber.

Therefore E_epsilon is exactly the domain

    x>0, a>0, ||w||<epsilon/x.

Substitute w=(epsilon/x)z with ||z||<1. The factor x^p from dX and the
factor (epsilon/x)^p from dw cancel exactly. For any bounded continuous
ambient test function h,

    epsilon^(-p) E[h(X,Y) 1_E_epsilon]
      = integral_(x,a>0) integral_S^(d-1) integral_(z in B_1(u-perp))
          h(xu,a u+(epsilon/x)z)
          g(xu,a u+(epsilon/x)z) dz dOmega(u) dx da.               (4)

Every finite nonsingular Gaussian admits constants C,c>0 with

    g(X,Y)<=C exp[-c(||X||^2+||Y||^2)].                           (5)

For example, bound the inverse covariance below by a positive multiple
of the identity and absorb the finite mean into C. Since a u is
orthogonal to z, the absolute integrand in (4) is bounded by

    C ||h||_infinity
      exp[-c(x^2+a^2+epsilon^2 ||z||^2/x^2)]
      <= C ||h||_infinity exp[-c(x^2+a^2)].                       (6)

This envelope is integrable in x,a, sphere direction, and the unit normal
ball. For every x,a>0, the integrand in (4) tends pointwise to

    h(xu,a u)g(xu,a u).

Dominated convergence applies globally, including x arbitrarily close to
zero; it is not necessary to make a local tube approximation near the axes.
Integrating the limiting integrand over z supplies kappa_p. Taking h=1
proves (1); dividing the general result by (1) proves (2).

The reduced-coordinate densities of (x,a,u) also converge in L1, hence in
total variation, relative to dx da dOmega. They are the z integrals in
(4), divided by their convergent positive normalizing constants; the same
domination proves the L1 claim. Section 3 separately establishes this for
the original magnitude coordinates (||X||,||Y||,X/||X||).

Ambient total-variation convergence is impossible: each finite-epsilon
conditional law is absolutely continuous in R^(2d), while (2) is
concentrated on a lower-dimensional cone. The ambient assertion is weak
convergence, not convergence in total variation.

The proof also permits continuous tests with polynomial growth in the
ambient norms. A polynomial times the Gaussian envelope can be absorbed
into a Gaussian envelope with a smaller exponent. Thus both magnitude means
converge, and so does their ratio-of-means readout. Expectations singular
at an axis, such as X's magnitude divided by Y's magnitude, require a
separate integrability argument and do not follow from this observation.

## 3. Original polar proof, retaining the magnitude coordinates

Write X=xu, Y=yv, with x,y>0 and u,v on S^(d-1). The ambient volume is

    dX dY=x^p y^p dx dy dOmega(u) dOmega(v).

For the angle theta between u and v, positive dot product means
theta<pi/2, and the area condition is xy sin(theta)<epsilon. For
q=epsilon/(xy), the allowed v form a cap of area

    A_d(q)=|S^(d-2)| integral_0^arcsin(min(1,q)) sin^(d-2)(theta)dtheta.

When q>=1 this is an entire positive hemisphere. The two needed bounds are

    A_d(q)/q^p -> kappa_p as q->0,
    A_d(q)<=C_d min(q^p,1) for every q>0.                         (7)

The elementary inequalities sin(theta)<=theta and
arcsin(q)<=(pi/2)q give a finite dimension-dependent C_d. In particular,

    (xy/epsilon)^p A_d(epsilon/(xy))<=C_d                         (8)

uniformly, even when xy<epsilon. No small-angle approximation is applied
in that region.

For bounded continuous h, the scaled expectation is

    integral dx dy dOmega(u) (xy/epsilon)^p
      integral_allowed_v h(xu,yv)g(xu,yv)dOmega(v).

After the v integral, (5) and (8) bound its absolute integrand by

    C_d C ||h||_infinity exp[-c(x^2+y^2)].

At each fixed x,y>0 and u, the cap shrinks to u. Its pointwise limit is
kappa_p h(xu,yu)g(xu,yu). Dominated convergence gives the same theorem.
Applying it to the nonnegative reduced-coordinate densities gives L1
convergence, before and after normalization. This establishes the stated
total-variation convergence of (||X||,||Y||,X/||X||).

## 4. Axes, origin, and the saturated-cap region

Both proofs show that the axes introduce no hidden component in the limit.
For example, uniformly in small epsilon,

    epsilon^(-p) P(E_epsilon, ||X||<delta or ||Y||<delta)=O(delta).

In the slice proof, ||X||<delta means x<delta, and ||Y||<delta implies
a<delta; the envelope (6) therefore proves this bound directly. A ball
of radius delta about the common origin gives an O(delta^2) bound by
restricting both x and a. After division by (1), taking delta to zero
rules out a limiting atom on either axis or at the origin.

For the region where the original magnitudes satisfy xy<epsilon, the
polar bound gives

    epsilon^(-p) P(E_epsilon, xy<epsilon)
      <=C' integral_(xy<epsilon) exp[-c(x^2+y^2)]dx dy
      =O(epsilon[1+|log epsilon|]).                              (9)

The last estimate follows by splitting x into (0,epsilon),
(epsilon,1), and (1,infinity). Their contributions are respectively
O(epsilon), O(epsilon|log epsilon|), and O(epsilon). The slice proof
also gives (9), since ||Y||>=a and x||Y||<epsilon implies xa<epsilon.

Thus the saturated-cap region has absolute probability
O(epsilon^d[1+|log epsilon|]), smaller than the leading epsilon^(d-1)
term. This makes explicit why the small-magnitude region does not undo
the magnitude-factor cancellation.

## 5. Exact normal-coordinate interpretation

At a positive compatible point (fu,alpha u), put r=sqrt(f^2+alpha^2).
For w perpendicular to u define the ambient pair

    X=fu-(alpha/r)w,
    Y=alpha u+(f/r)w.                                            (10)

This is a normal displacement of Euclidean norm ||w||. Direct calculation
gives exact identities

    ||X wedge Y||=r||w||,
    X dot Y=f alpha[1-||w||^2/r^2].                               (11)

When ||w||<r, the singular values of [X Y] are r and ||w||. The compatible
point is the nearest rank-one projection, with positive aligned magnitudes.
Thus r in (11) is the projection's norm, or the larger singular value,
not the ambient radius sqrt(||X||^2+||Y||^2)=sqrt(r^2+||w||^2).

The induced surface volume on the compatible cone is

    dVol_cone=r^p df d alpha dOmega(u).

Away from the singular strata, area<epsilon is a normal ball of radius
epsilon/r. Its leading volume kappa_p epsilon^p r^(-p) cancels the r^p
surface factor and leaves flat magnitudes. The global normal chart needs
||w||<r and degenerates at its boundaries; Sections 2-4 justify the full
limit without relying on a uniform normal-chart approximation there.

## 6. Dimension one

At d=1, the wedge norm is identically zero. For every positive epsilon,
E_epsilon is just XY>0. With counting measure on S^0={-1,+1},

    P(XY>0)=sum_(u=+-1) integral_(x,y>0) g(xu,yu)dx dy.

The flat-magnitude compatible law therefore holds exactly after this sign
conditioning; epsilon does nothing. There is no lower-dimensional limit
or small angular cap to analyze. Formally p=0 and kappa_0=1, but the
polar proof involving S^(d-2) must not be applied at d=1. The Cartesian
slice proof can be interpreted with a zero-dimensional normal fiber.

## 7. Independent Gaussian normalization check

For independent X,Y~N(0,I_d), the area is the product of independent
chi_d and chi_(d-1) lengths: condition on the direction of X and take
Y's perpendicular length. The positive-dot-product sign independently
retains one half of the parallel Gaussian signs.

Gamma duplication gives the area's nth moment as

    2^n Gamma((d+n)/2) Gamma((d-1+n)/2)
      /[Gamma(d/2) Gamma((d-1)/2)]
      =Gamma(d-1+n)/Gamma(d-1).

The product has a moment-generating function near zero because
RQ<=(R^2+Q^2)/2, so these moments identify a Gamma(d-1,1) distribution.
Consequently P(E_epsilon) is one half of that gamma CDF. In particular,

    d=2: P(E_epsilon)=(1-exp(-epsilon))/2 ~ epsilon/2,
    d=3: P(E_epsilon)=[1-(1+epsilon)exp(-epsilon)]/2 ~ epsilon^2/4.

Direct evaluation gives

    Z_flat=2^(-d) pi^(1-d/2)/Gamma(d/2),
    kappa_(d-1) Z_flat=1/[2 Gamma(d)],

agreeing with the exact small-epsilon probabilities in all d>=2.

## 8. Interpretation and remaining choices

This supplies the original flat reference with an explicit Gaussian
conditioning interpretation: small absolute parallelogram area, together
with positive alignment. Dropping positive alignment would retain the
antiparallel branch as well.

It does not make this the only possible conditioning law. In d>=2 exact
collinearity has probability zero. Constant angular tolerance, fixed
Euclidean normal distance, and fixed area specify different neighborhoods
of that event. Their respective magnitude weights are (xy)^p, r^p, and
1. Coherence of the area construction does not prove unique optimality
or derive it from the composition/readout axioms.

Fixed positive channel rescalings X=F/sigma_F and Y=a/sigma_a transform
dx dy into df d alpha by a constant factor, so flat standardized magnitudes
are flat physical magnitudes after normalization. Common orthogonal spatial
rotations preserve area. Arbitrary joint whitening or a changed spatial
metric must be declared, because it need not preserve the physical
collinearity constraint or the chosen tolerance.

Finally, g is assumed to be an ambient probability density. A Gaussian
measurement likelihood, viewed only as a function of latent true vectors,
does not by itself uniquely specify that probability law. A normalized
Gaussian likelihood with an ambient Lebesgue reference is one possible
choice; it has not become a uniquely determined posterior merely by
applying this conditioning theorem.

# Area, angle, and normal distance

Bounded independent mathematical audit, 2026-09-15. Scratch only.

## The same zero set can have different normal sensitivities

Work in standardized Cartesian coordinates X,Y in R^d, d>=2. Put

\[
x=\|X\|,\quad y=\|Y\|,\quad
r=\sqrt{x^2+y^2},\quad
\mathcal A=\|X\wedge Y\|=xy\sin\gamma.
\]

Positive compatibility is `X=x u,Y=y u,x,y>0`. A tolerance construction also
needs positive alignment, such as X dot Y>0, to exclude negative mass.
Let c=d-1 be the codimension of this manifold in R^(2d).

At a compatible point define normal displacement

\[
n={y\,\delta X_\perp-x\,\delta Y_\perp\over r}\in u^\perp.
\]

The wedge vector has the first-order expansion

\[
X\wedge Y=-r\,u\wedge n+O(\|\delta\|^2),
\qquad \mathcal A=r\|n\|+O(\|\delta\|^2).
\]

The manifold's Euclidean surface volume is `r^c dx dy dOmega`. Thus:

| Tolerance residual | First-order normal sensitivity | Induced density relative to dx dy dOmega |
|---|---|---|
| Area A | r ||n|| | 1: flat magnitudes |
| Angle, or A/(xy) | r ||n||/(xy) | (xy)^c |
| Normal distance, or A/r to first order | ||n|| | r^c |

For example, area<epsilon allows a normal ball of radius epsilon/r. Its
volume contributes r^(-c), cancelling the surface factor r^c. Angular
tolerance allows a normal ball of radius epsilon xy/r and retains the
Cartesian magnitude factors (xy)^c. A fixed Euclidean tube retains surface
volume directly.

The scalar norm A is not differentiable at its zero. The rigorous first-order
statement concerns the wedge vector's derivative, or equivalently its normal
sensitivity, rather than an ordinary scalar gradient at compatibility.

Multiplying the ambient Gaussian density by these limiting measures gives
the three different conditional laws already investigated. The argument
applies to finite nonzero observed centers too. The zero set alone does not
select how a vanishing neighborhood is weighted.

### A compact general family

For fixed real a,b, use the defect

\[
\delta_{a,b}={\mathcal A\over(xy)^b r^a}.
\]

Its normal sensitivity is `r^(1-a)(xy)^(-b)||n||`, so conditioning on a small
positive-alignment tolerance induces

\[
\boxed{\rho_{a,b}(x,y)=(xy)^{b(d-1)}r^{a(d-1)}.}
\]

Write beta=b(d-1), eta=a(d-1). In magnitude polar coordinates its null law
has weight

\[
r^{2\beta+\eta+1}e^{-r^2/2}
(\sin\theta\cos\theta)^\beta\,dr\,d\theta.
\]

The normalized Gaussian conditional law on the strict positive interior is
proper precisely when

\[
\beta>-1,\qquad 2\beta+\eta>-2.
\]

The first condition controls the axes, and the second controls the common
origin. Gaussian decay handles infinity. For d=3 these conditions are
`b>-.5` and `a+2b>-1`. If they fail, the proposed interior density cannot be
normalized; one cannot simply use the formal local expression as a proper
conditional distribution. A boundary-dominated limit would need separate
analysis.

These conditions also justify the global thin-tolerance limit. At fixed x,y
the accepted positive angular cap has angular radius controlled by
`epsilon r^a (xy)^(b-1)`. Its probability is bounded by a constant times that
quantity to power c, even when the cap saturates. After multiplication by
the ambient radial factors, an integrable Gaussian envelope times rho_(a,b)
dominates the rescaled density. Dominated convergence then upgrades the local
normal calculation to a normalized limit.

All of these are homogeneous, exchange-symmetric and rotation-invariant.
They do not exhaust such measures: generally a homogeneous symmetric density
can be written `r^p h(theta)` with h(theta)=h(pi/2-theta). At null it is proper
when p>-2 and h is nonnegative, integrable, and not zero almost everywhere.
Homogeneity and symmetry leave an entire angular function open.

### The axis exponent controls mass tails

For the power family, the null radial law is chi_(2 beta+eta+2), independently
of the ratio angle. The ratio R=M/(sigma_F/sigma_a) has density

\[
p_R(z)={2\over B((\beta+1)/2,(\beta+1)/2)}
{z^\beta\over(1+z^2)^{\beta+1}}.
\]

The null ratio-of-means point equals sigma_F/sigma_a for every admissible
(a,b). The mass law depends on b through beta; a changes the radial law but
not this null ratio distribution. In particular, area and distance-tube
conditioning share the half-Cauchy null mass law while having different full
latent magnitude laws.

For every finite Gaussian dataset, not only null data, the power family's
ratio density is proportional to m^beta near zero and m^(-beta-2) in the far
tail, with positive finite asymptotic coefficients. Thus precisely the power
moments with `-(beta+1)<p<beta+1` exist. The radial exponent changes the tail
coefficients but not these exponents. All fixed-order log moments exist in
this power family. More general integrable angular functions h need not have
all log moments: propriety alone does not determine tail regularity.

## Exact Gaussian sampling law of squared area

Now take **random observations** with compatible true standardized means
x u,y u and independent unit Gaussian errors. No observed dot-positive
filter is applied in the following sampling identity.

For r=sqrt(x^2+y^2)>0 rotate the two channel columns:

\[
P={xX+yY\over r},\qquad Q={-yX+xY\over r}.
\]

Then P~N(r u,I_d), Q~N(0,I_d), independently, and the rotation preserves the
wedge: P wedge Q=X wedge Y. Conditionally on P, the squared perpendicular
length of Q is chi-square with d-1 degrees of freedom, with a law independent
of P. Consequently

\[
\boxed{\mathcal A^2=UV,\quad
U\sim\chi^2_d(r^2),\quad V\sim\chi^2_{d-1},\quad U\perp V.}
\]

The same statement holds at r=0 with any channel rotation. In particular,

\[
\boxed{E[\mathcal A^2]=(d-1)(r^2+d).}
\]

The extra d(d-1) term is a genuine noise floor at weak excitation. For d=3,
the expectation is 2r^2+6. A constant area threshold therefore is not a
constant Gaussian-significance threshold across different true magnitudes.

With W the squared component of Q parallel to P, U,V,W are independent,
W~chi-square_1, and the observed squared total radius is U+V+W. Therefore

\[
{\mathcal A^2\over\|X\|^2+\|Y\|^2}={UV\over U+V+W}.
\]

This is not an exact chi-square pivot at finite signal. At true null its mean
is `d(d-1)/(2(d+1))`, which is 0.75 in d=3. At strong signal it tends in
distribution to chi-square_(d-1), whose mean is 2 in d=3. The null expectation
follows by writing the three independent central gamma variables as a common
sum times Dirichlet proportions.

## What warrant the Gaussian geometry supplies

At strong compatible signals, the normal displacement n has d-1 independent
unit Gaussian coordinates to first order. Thus area/r is locally measured in
normal-noise units; angular discrepancy has local noise scale r/(xy). This
explains why a fixed noise-metric tube has a local Gaussian/Fisher geometric
warrant.

It does not prove that the tube is the uniquely correct reference measure.
Constant area, constant angle and constant normal displacement express
different choices about what is held equal while approaching the same exact
compatibility relation. A physical protocol might warrant one of those
defects; Gaussian measurement errors alone do not make them interchangeable.

Nor should the strong-signal normalization be advertised as an exact weak-
signal significance test. Observed radius is noisy, and A/r_observed is only
asymptotically the Euclidean normal distance. The earlier exact mass-confidence
pivots have a different, explicitly candidate-wise calibration argument.

No new weighting policy, rejection gate, or ranking is installed by this audit.

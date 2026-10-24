# Exact conditional likelihood-ratio confidence construction

Working derivation, 2026-09-15. Scratch only; no production or managed edits.
The Gaussian model, positive reported mass domain, and larger comparison family
are specified separately below.

## 1. Candidate likelihood ratio

Let X,Y in R^K be independent standardized Gaussian vectors with means
(r v,v), where r>0 is the mass ratio and v is an unrestricted fixed nuisance
vector. For each candidate r define

    R=(X-rY)/sqrt(1+r^2), S=(rX+Y)/sqrt(1+r^2),
    t=||R||^2, s=||S||^2, H=R dot S.

At the true r, R~N(0,I_K) independently of S. Profiling v out at a specified
candidate leaves residual t. Profiling the direction of the line over the
entire projective channel plane, allowing signed slopes and the infinite-slope
limit, gives the smaller eigenvalue lambda_min of Gram(X,Y). The twice-log
likelihood-ratio statistic against that larger family is

    Lambda=t-lambda_min
      = [t-s+sqrt((t-s)^2+4H^2)]/2.                                  (1)

It is nonnegative. The same statistic can be computed in the (R,S) basis
because the channel rotation preserves the Gram eigenvalues.

For any c>0, Lambda is the nonnegative root of
lambda^2-(t-s)lambda-H^2=0. Hence

    Lambda<=c iff c^2-c(t-s)-H^2>=0.

## 2. Exact conditional calibration

For s>0, condition on S and put

    Z=H/sqrt(s), U=Z^2, W=t-U.

Then U~chi-square_1 and W~chi-square_(K-1), independently; their joint law
does not depend on S or the nuisance vector. Substituting H^2=sU gives

    Lambda<=c iff U/c+W/(c+s)<=1.                                   (2)

Thus the conditional CDF is

    F_K(c;s)=P[U/c+W/(c+s)<=1].                                      (3)

For K>1 it can be evaluated by a one-dimensional integral, for example

    F_K(c;s)=2 integral_0^sqrt(c) phi(z)
       F_chi-square_(K-1)((c+s)(1-z^2/c)) dz.                         (4)

For K=1, W=0 and F_1(c;s)=F_chi-square_1(c) for every s. At s=0, use the
continuous calibration extension: Lambda=t and F_K(c;0)=F_chi-square_K(c).
The event S=0 has probability zero at each fixed true candidate, but this
extension also supplies an explicit rule for exactly degenerate supplied data.

For any interior probability p, F_K(c;s) is continuous and strictly increasing
from zero to one as c goes from zero to infinity. Define c_K(s) by

    F_K(c_K(s);s)=p, with p=.95 below.

Accepting Lambda(r)<=c_K(s(r)) therefore has exactly 95% probability at every
fixed true r and every fixed nuisance v, conditionally on S and unconditionally.
Restricting the returned set to r>0 preserves that coverage for positive truths.

### Endpoints, monotonicity and a useful three-dimensional check

For fixed c, increasing s reduces W/(c+s), so F_K(c;s) increases. Thus c_K(s)
decreases continuously, strictly for K>1, with

    c_K(0)=chi-square_(K,.95),
    lim_(s->infinity) c_K(s)=chi-square_(1,.95).

For K=3 the W CDF is exponential, allowing the independent formula

    F_3(c;s)=erf(sqrt(c/2))
      -sqrt(2c/pi) exp(-c/2) integral_0^1 exp[-s(1-u^2)/2]du.          (5)

It reduces to the chi-square_3 CDF at s=0 and tends to the chi-square_1 CDF
as s grows. The bounded integral in (5) avoids separately evaluating growing
erfi factors.

## 3. Large-dimension calibration

If K tends to infinity and the conditioning values satisfy s/K->kappa>1,
then W/K->1. For fixed positive c,

    F_K(c;s) -> F_chi-square_1(c(1-1/kappa)).

Consequently

    c_K(s) -> chi-square_(1,.95)/(1-1/kappa),

or c_K(s) approximately chi-square_(1,.95)/(1-K/s) in that regime. The
ratio s/K must stay above one; this approximation cannot be used at the
true-null transition s=K+O(sqrt(K)), and is meaningless when it turns negative.

For context, if s/K->kappa<1, the leading critical scale is K-s, rather than
a constant. At the transition s=K+gamma sqrt(K), write c=a sqrt(K). Then

    Lambda/sqrt(K) =>
      [sqrt(2)N-gamma+sqrt((sqrt(2)N-gamma)^2+4U)]/2,

with N standard normal independently of U~chi-square_1. Equivalently its CDF
at a>0 is E Phi((gamma+a-U/a)/sqrt(2)). This identifies the sqrt(K) critical
scale in weak/null conditions without applying the strong-signal formula there.

## 4. Exact inversion reduces to one quadratic

For fixed observed data put

    A=||X||^2, B=||Y||^2, C=X dot Y,
    L=lambda_max(Gram(X,Y)).

The trace and eigenvalues are invariant under the candidate rotation. Since
Lambda=t-lambda_min and t+s=A+B, one has the exact identity

    Lambda+s=L.

Because Lambda<=c_K(s) iff F_K(Lambda;s)<=.95, candidate acceptance is

    P[U/(L-s)+W/L<=1] <= .95,                                      (6)

with its limiting interpretation when L-s=0. For a fixed L>0 the expression
in (6) is strictly decreasing in s over (0,L). At s=0 it equals the
chi-square_K CDF at L; at s=L its limit is zero.

If L<=chi-square_(K,.95), every projective candidate passes. Otherwise there
is a unique b in (0,L) such that

    F_K(L-b;b)=.95.

The complete accepted set is simply s(r)>=b, or

    (A-b)r^2+2Cr+(B-b)>=0,   r>0.                                 (7)

Thus inversion needs one scalar calibration root per dataset and then one
quadratic inequality, not a critical-curve calculation at every candidate.
L=0 means X=Y=0 and all positive candidates pass.

Equation (7) can give a bounded interval, one-sided set, two rays, all positive
values, or an empty set. It does not create an isolated singleton: that would
require b to equal an interior maximum of s, but the interior maximum is L
and calibration gives b<L. Zero and infinity are parameter-boundary
symbols, not finite mass values. Finite endpoints satisfying equality belong
to the set. Even if L exceeds its whole-projective cutoff, every positive
candidate can still pass when the quadratic holds on the positive domain.

## 5. What it changes relative to the directional score

The score vanishes whenever H=0. If t<s this is the unrestricted profile
minimum and Lambda=0 too. If t>s, however, it is the profile maximum and
Lambda=t-s. Such a candidate is no longer accepted automatically merely
because its score vanishes; it must pass a calibrated likelihood improvement
comparison.

For strongly opposed Gram summaries A=B=100,C=-90, the positive stationary
score root is r=1 with t=190,s=10 and Lambda=180. The positive CLR set is
empty: the best signed line is substantially better than every positive line.

For A=B=100,C=0, all candidate lines have the same residual, Lambda=0 everywhere,
and the CLR set is the whole positive domain. This is an informative limit:
CLR measures relative evidence about the mass direction, not the absolute size
of the full residual. The full-residual and split tests can reject these data
at every mass, whereas CLR recognizes that no particular mass fits better.

In a growing-trial, positive-average-excitation regime, the CLR local critical
inflation is matched by the corresponding inflation of Lambda relative to Z^2.
It therefore has the same leading root-N local width as the exactly calibrated
score procedure, including the nuisance-direction noise contribution. That
does not establish a universal finite-sample power ordering.

## 6. The larger alternative must be stated accurately

The physical reported hypotheses remain positive mass. The statistic compares
their fit to a larger geometric family allowing signed slopes and the vertical
line. A calibrated test may use this larger benchmark without interpreting
negative mass as a physical possibility. Restriction of accepted candidates
to the positive domain does not change coverage at any positive true mass.

It is a different test from profiling only over positive masses. In the
original X,Y coordinates the positive-family best residual is

    lambda_min when C>=0,
    min(A,B) when C<0,

using boundary suprema when needed. The resulting positive-alternative LR is
no greater than (1). Reusing the critical curve derived above for that smaller
statistic would therefore be conservative, rather than generally exact.

For exact recalibration of the positive-alternative LR, the conditional
distribution also depends on the candidate's angular position: allowable
relative rotations are (-theta,pi/2-theta). The sign of Z matters as well as
Z^2. It cannot be replaced by the current c_K(s) without another argument.

This conditional LR construction is a mass-confidence procedure under the
assumed Gaussian law. An empty set does not itself suppress a separately
defined point or constitute a policy for rejecting Newton II. Those uses,
and the comparison family's operational meaning, remain separate decisions.

## Status

All displayed identities and calibration claims are direct derivations for
the stated model. Numerical inversion still requires bounded root and CDF
validation. No claim of uniquely optimal uncertainty, point-estimator
superiority, or historical novelty is made.

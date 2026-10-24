# Exact mass-confidence sets without a law-rejection gate

Working mathematical audit, 2026-09-15. Scratch only.

## 1. Residual test inversion is exact

For independent Gaussian observations `Y_i=(F_i,a_i)` with fixed true means
`(m v_i,v_i)` and known positive definite covariances Sigma_i, put

\[
r_i(m)=F_i-ma_i,\qquad
V_i(m)=\Sigma_{FF,i}+m^2\Sigma_{aa,i}
-m(\Sigma_{Fa,i}+\Sigma_{aF,i}).
\]

At the true m, the residuals are independent mean-zero Gaussians with these
covariances. Hence, with K=Nd,

\[
T(m)=\sum_i r_i(m)^TV_i(m)^{-1}r_i(m)\sim\chi^2_K.
\]

Therefore the random set

\[
\mathcal C_T(Y)=\{m>0:T(m)\le\chi^2_{K,.95}\}
\]

has exactly 95% coverage at **every** fixed positive m and every collection
of fixed nuisance vectors v_i. No prior or excitation distribution is needed.
The nuisance dimension has already been eliminated by the contrasts; no
further subtraction of one degree of freedom is appropriate while testing a
specified candidate m. The test is not based on the minimized statistic over m.

A candidate-dependent statistic is legitimate for test inversion. Under
each fixed null candidate its distribution is known. This does not turn the
statistic's Gaussian residual density into a likelihood for the original
fixed observed data, and does not contradict the earlier likelihood warning.

The set may be disconnected, unbounded, or empty. An empty set has sampling
probability at most .05 at any true parameter, since it necessarily excludes
that true mass. Reporting this confidence-set outcome need not reject Newton
II or suppress a separately defined point estimate. It records the result of
a particular 95% mass-confidence construction under the assumed law.

Known Gaussian covariances and independence are substantive assumptions.
With a known full covariance across trials, stack the residuals and use their
full covariance to obtain the same chi-square result. Estimated covariance
or non-Gaussian noise requires a different calibration argument.

## 2. Homogeneous covariance gives a Fieller quadratic

For independent isotropic channels with common sigma_F,sigma_a, let

\[
A=\sum_i\|a_i\|^2-q\sigma_a^2,\quad
B=-2\sum_iF_i\cdot a_i,\quad
C=\sum_i\|F_i\|^2-q\sigma_F^2,
\qquad q=\chi^2_{K,.95}.
\]

Then `C_T={m>0:A m^2+B m+C<=0}`. This includes empty, singleton, bounded,
one-sided, all-positive and two-ray possibilities. Finite roots satisfying
equality are included; zero and infinity are parameter boundaries, not masses.
For K=1 this is the known-variance Fieller construction.

Isotropic within-trial channel correlation rho changes B to
`-2 sum F_i dot a_i + 2 q rho sigma_F sigma_a`; the common denominator is
`sigma_F^2-2m rho sigma_F sigma_a+m^2 sigma_a^2`. General anisotropic
covariance does not automatically reduce to a scalar quadratic.

## 3. Exact null shape probabilities for K=3

At true zero excitation, stack the standardized observations as independent
X,Y~N(0,I_K). Put r=m/(sigma_F/sigma_a), a=||X||^2, b=||Y||^2, c=X dot Y.
Then

\[
T(r)={a-2cr+br^2\over1+r^2}.
\]

For c>=0, its supremum over positive r is max(a,b). For c<0, its supremum
is the largest eigenvalue of the 2x2 Gram matrix of X,Y. Sign symmetry and
independence of a,b therefore give

\[
P_0(\mathcal C_T=(0,\infty))
=\tfrac12(.95)^2+\tfrac12P(\lambda_{\max}(W_2(K,I))\le q).
\]

For K=3, the ordered Wishart eigenvalue density is

\[
\tfrac14(\lambda_1-\lambda_2)e^{-(\lambda_1+\lambda_2)/2},
\qquad\lambda_1>\lambda_2>0.
\]

Integrating it gives

\[
P(\lambda_{\max}\le q)=1-e^{-q}-q e^{-q/2},\qquad
P(\lambda_{\min}>q)=e^{-q}.
\]

With q=7.814727903251173, the exact null shape probabilities are:

| Shape | Probability |
|---|---:|
| Whole positive domain | 0.8725359068 |
| Two positive rays | 0.0299640932 |
| One-sided, both orientations combined | 0.0950000000 |
| Bounded interval | 0.0010481277 |
| Empty | 0.0014518723 |

This does not contradict the previous 90% whole-domain necessity for
**connected intervals**. The new object is a possibly disconnected set.
It accepts both endpoint limits with probability `.95^2=.9025`, but sometimes
excludes an intermediate band. Taking its interval hull would fill that
band and give a conservative confidence interval; it changes the output.

## 4. An exact conditional directional pivot

For homogeneous independent isotropic channels, stack standardized vectors
as X,Y in R^K and let r=m/(sigma_F/sigma_a). For each candidate r define

\[
R(r)={X-rY\over\sqrt{1+r^2}},\qquad
S(r)={rX+Y\over\sqrt{1+r^2}}.
\]

At the true r, with stacked w=(v_1/sigma_a,...,v_N/sigma_a),

\[
R\sim N(0,I_K),\qquad
S\sim N(\sqrt{1+r^2}\,w,I_K),\qquad R\perp S.
\]

Consequently, conditionally on every nonzero S,

\[
Z(r)={R\cdot S\over\|S\|}\sim N(0,1).
\]

The conditional distribution does not depend on S or w, so inversion of
`Z(r)^2<=chi^2_(1,.95)` gives another exact 95% confidence set for all fixed
true masses and all nuisance vectors. It is not the same confidence set as
the omnibus residual construction.

### Quartic shape and the stationary maximum

Write a=||X||^2, b=||Y||^2, c=X dot Y and D=a-b. Away from S=0,

\[
Z(r)^2={\{c(1-r^2)+rD\}^2
\over(1+r^2)(ar^2+2cr+b)}.
\]

Its inversion is a quartic inequality. With threshold q its coefficients,
from r^4 down to r^0, are

    c^2-q a,
    -2c(D+q),
    D^2-2c^2-q(a+b),
    2c(D-q),
    c^2-q b.

The numerator also vanishes at stationary points of T(r). If c>0 the
positive stationary root minimizes T; if c<0 it maximizes T. Thus for K>=2
and nondegenerate data, Z inversion is nonempty and always includes that
stationary root, even when it is a very poor full-residual reconciliation.

For example X=(10,0), Y=(-10,1) has a positive stationary maximum near
r=1.0050125. Z is zero there, while T is about 200.5. The directional test has
ignored the large residual orthogonal to S. This is a power/shape limitation,
not a failure of its exact coverage guarantee. In particular, it gives no
reason to select that accepted root as a point mass estimate.

### S=0 requires an explicit convention

For each fixed true candidate the event S=0 has probability zero, so changing
the test there does not change pointwise coverage. But when K=1 and observed
X,Y have opposite signs, S(r)=0 occurs at a random positive candidate. Blindly
multiplying by ||S||^2 admits an undefined zero of the quartic.

For K=1, cancel the factor first: Z^2=R^2 wherever S!=0, with a continuous
extension at S=0. This recovers the ordinary Fieller statistic. Exact negative
collinearity in K>1 similarly has a continuous squared-statistic extension
equal to ||R||^2 at its S=0 point. If X=Y=0, defining Z=0 for every candidate
is natural. These conventions should be stated before numerically inverting.

## 5. Efficiency: fixed nuisance dimension versus new nuisances per trial

For S!=0 the exact derivative is

\[
{dZ\over dr}={\|R\|^2-\|S\|^2-Z^2\over(1+r^2)\|S\|}.
\]

If the nuisance direction is learned accurately, meaning
`(1+r^2)||w||^2 >> K`, Z has the efficient leading local slope
`-||w||/sqrt(1+r^2)`. This includes fixed dimension with strong signal, or
repeated measurements of one unchanged pair **aggregated first** into its
Gaussian sample mean.

For N different nuisance vectors, let K=Nd and
`L_N=sum ||v_i/sigma_a||^2`, with `L_N/N -> ell>0`. Then at the true r,

\[
{1\over\sqrt N}{dZ\over dr}
\ \xrightarrow{P}\ -{\ell\over\sqrt{d+(1+r^2)\ell}}.
\]

Local inversion therefore has root-N width, with leading half-width

\[
z_{.975}{\sqrt{d+(1+r^2)\ell}\over\ell\sqrt N}.
\]

Compared with the fixed-parameter Fisher benchmark this is inflated by
`sqrt(1+d/((1+r^2)ell))`. Estimating a growing number of nuisance vectors
from one noisy pair each is not the same asymptotic problem as fixed-
dimensional nuisance estimation. An unqualified efficiency claim would be
wrong in this regime.

By contrast, the omnibus chi-square threshold is of order K+O(sqrt K).
At fixed positive average excitation its typical local slack is O(sqrt N),
while mass displacement changes the residual mean-square by O(N delta^2).
Its usual width is therefore of order N^(-1/4). This describes this inversion
and regime, not an impossibility theorem for exact confidence procedures.

## 6. Optional orthogonal-residual refinement

For K>=2, conditionally on S, decompose

\[
W(r)=\|R\|^2-Z(r)^2.
\]

Then Z^2~chi^2_1 and W~chi^2_(K-1) are independent, with their joint law also
independent of S. Choose fixed acceptance probabilities p_Z,p_W with
`p_Z p_W=.95`. Invert the intersection

\[
Z(r)^2\le\chi^2_{1,p_Z},\qquad
W(r)\le\chi^2_{K-1,p_W}.
\]

This is exactly 95% pointwise calibrated. A symmetric choice is
`p_Z=p_W=sqrt(.95)`. It allocates some error probability to the orthogonal
residual instead of ignoring it. The strong anti-aligned example above is
excluded because its W is huge despite Z=0.

The tradeoff is explicit. At a strong compatible signal W is almost fixed
throughout the local mass interval, while Z changes rapidly with mass. When
W passes, the local root-N interval remains, with a wider Z cutoff. When W
fails, the local set disappears. In the strong-signal limit the probability
of an empty set tends to `1-p_W`; for the symmetric split this is about
2.532%. This remains a mass-confidence output, not a law-rejection gate.

Useful identities explaining the behavior are

\[
W(r)={ab-c^2\over\|S(r)\|^2},\qquad
{dW\over dr}=-{2ZW\over(1+r^2)\|S\|}.
\]

For strong compatible data c>0 and W has its global minimum at the good
stationary root. Its minimum approaches the true-candidate orthogonal
chi-square component. Dropping only the unwanted stationary maximum from
the original Z confidence set, without a calibrated additional test, would
not retain the exact coverage argument.

For K=1, W is identically zero, so its acceptance probability is 1. Use the
Z test at .95; using sqrt(.95) would instead give 97.468% coverage.

These exact constructions illustrate that a calibration requirement still
does not uniquely specify an uncertainty set. Their shape, power, nuisance
treatment and intended use remain part of the investigation.

# Independent mathematical review: known force, weak acceleration

Internal contribution for integration by the lead investigator, 2026-09-15.
Scope: exact known positive force and direction, a signed scalar Gaussian
acceleration reading, known standard uncertainty, and the original flat
positive-acceleration reference measure. No novelty claim. All results below
are derived directly; `check_math.py` checks numerical identities, quadrature,
root inversion, and the coarsened-report coverage formula.

## 1. Model and exact posterior

Write `F>0`, true acceleration `alpha>0`, mass `m=F/alpha`, and

\[
Y\sim N(\alpha,\sigma^2),\qquad d\Pi=d\alpha\quad(\alpha>0).
\]

Set `t=y/sigma`, `theta=alpha/sigma`, and let `phi,Phi` be the standard-normal
density and CDF. The reference measure is improper but the posterior is proper:

\[
p(\alpha\mid y)=\frac{\phi((\alpha-y)/\sigma)}{\sigma\Phi(t)},
\quad\alpha>0.
\]

For an isotropic vector acceleration measurement with exact known force
direction, the signed projection along that direction supplies this `Y`; the
orthogonal components contribute a parameter-independent likelihood factor.
For a known general positive-definite acceleration covariance `S`, use
`y=(u^T S^{-1} a)/(u^T S^{-1}u)` and
`sigma²=1/(u^T S^{-1}u)`. Force uncertainty is not covered by these formulas.

The original ratio-of-means readout reduces to

\[
H(t)=t+\frac{\phi(t)}{\Phi(t)},\qquad
\widehat m=\frac{F}{\sigma H(t)}.
\]

`H` is strictly increasing, positive, and maps the real line onto `(0,infinity)`.
Its derivative is the variance of the dimensionless truncated normal, so lies
strictly between zero and one. The point is the reciprocal of the posterior
mean inverse mass. It is the posterior squared-error-optimal prediction of
acceleration (equivalently inverse mass), transformed back to mass; it is not
the posterior mean mass.

Indeed

\[
p(m\mid y)=\frac{F}{\sigma\Phi(t)m^2}
\phi\!\left(\frac{F}{\sigma m}-t\right),\qquad m>0.
\]

The `m^{-2}` upper tail makes the posterior mean mass infinite at every finite
reading. The flat `d alpha` reference corresponds to a mass reference density
proportional to `m^{-2}`. This assumption must accompany any interpretation.

## 2. A posterior lower mass bound with an exact coverage theorem

For posterior level `p`, let `delta=1-p`. The acceleration upper quantile is

\[
Q_p(t)=t+\Phi^{-1}(1-\delta\Phi(t)),\qquad q_p=\sigma Q_p(t),
\]

and the posterior lower mass bound is `L_p(y)=F/q_p`. This is exactly a
posterior statement `P(m >= L_p | y)=p`.

For numerics, avoid subtracting nearly equal values by solving

\[
\log\Phi(t-Q_p)-\log\Phi(t)=\log\delta.
\]

The family of posteriors has a monotone likelihood ratio in `t`, so `Q_p(t)`
is strictly increasing. Given a fixed true `theta>0`, let `t_*` be its unique
inverse: `Q_p(t_*)=theta`. Since the sampling observation is
`T=theta+Z`, `Z~N(0,1)`,

\[
\begin{aligned}
C_p(\theta)
&=P_\theta\{L_p(Y)\le m\}
=P_\theta\{T\ge t_*\}=\Phi(\theta-t_*)\\
&=1-\delta\Phi(t_*).
\end{aligned}
\]

Consequently `p<C_p(theta)<1` for every finite positive `theta`, with limits
`1` as `theta -> 0` and `p` as `theta -> infinity`. It decreases with `theta`.
Thus the continuous-reading posterior lower bound is conservative under
repeated sampling throughout this restricted experiment.

For `p=.95`, the coverage is approximately:

| True alpha/sigma | Coverage |
|---:|---:|
| 0.5 | 0.9999999994 |
| 1 | 0.99928205 |
| 2 | 0.97380084 |
| 3 | 0.95478111 |
| 5 | 0.95001984 |

Values displayed as `1` at smaller signal are rounded, not exact.

At the precise reading `y=0`,

\[
E[\alpha\mid0]=\sigma\sqrt{2/\pi},\qquad
\widehat m=1.25331414\,F/\sigma,\qquad
L_{.95}=0.51021346\,F/\sigma.
\]

## 3. Same-information confidence benchmark

The classical one-sided Gaussian upper confidence endpoint for acceleration is
`U(y)=y+z_p sigma`, with exact sampling coverage `p`. Its positive-parameter
confidence set is `(0,U]` if `U>0`, and empty otherwise. Transforming gives the
mass set `[F/U,infinity)` if `U>0`, empty otherwise. Setting `U=0` after clipping
does not turn that empty set into a finite-mass determination.

Since `Q_p(t)>t+z_p`, the posterior lower mass bound is always weaker wherever
the Gaussian mass endpoint is finite. At `y=0` the Gaussian coefficient is
`1/1.64485363 = .60795683`, versus posterior `.51021346`. The posterior procedure
always supplies a nonempty one-sided set; the exact Gaussian procedure can
produce an empty set with probability `Phi(-z_p-theta)<delta`. This is a useful
benchmark for this specific problem, not an alternative to the full estimator
or its composition requirements.

## 4. A finite point does not become an accurate high-mass estimate

At fixed `F,sigma`, take `theta -> 0`, hence true mass
`m=F/(sigma theta) -> infinity`. Under a common standard-normal `Z`,

\[
\frac{\widehat m}{m}
=\frac{\theta}{H(\theta+Z)}\longrightarrow0
\quad\text{almost surely}.
\]

The estimator itself remains of order `F/sigma`; its relative underestimation
becomes arbitrarily severe in this limit. This does not invalidate the useful,
conservative lower bound. It separates a successful one-sided inference from
an accurate finite central determination. Other finite posterior quantiles
also remain of order `F/sigma` in this limit.

For increasingly negative exact readings, `H(t) ~ 1/(-t)` and
`mhat ~ F(-y)/sigma²`. The likelihood favors the boundary `alpha=0` whenever
`y<=0` (equivalently unbounded mass), whereas the chosen reference measure and
point readout give a finite output.

For repeated readings of one unchanged positive acceleration, replacing
`y,sigma` by `mean(Y),sigma/sqrt(n)` yields consistency at every fixed finite
positive mass. The weak-signal limiting result above is not inconsistency at
a fixed interior parameter.

## 5. A nondetection report is a different observation

Suppose the instrument reports only

\[
E_k=\{|Y|\le k\sigma\},\qquad k>0.
\]

Its likelihood is

\[
L_k(\theta)=\Phi(k-\theta)-\Phi(-k-\theta).
\]

Since this function is even on the extended real line and its full integral
is `2k`, the positive-half integral is `k`. Thus

\[
p(\theta\mid E_k)=L_k(\theta)/k,\quad\theta>0.
\]

The posterior mean and CDF are explicit:

\[
E[\theta\mid E_k]=
\frac{(k^2+1)(\Phi(k)-1/2)+k\phi(k)}k,
\]

\[
G_k(x)=\frac{A(k+x)-A(k-x)-x}{k},\qquad
A(v)=v\Phi(v)+\phi(v).
\]

The point is `(F/sigma)/E[theta|E_k]`; the mass lower bound is
`F/(sigma B_p(k))`, where `G_k(B_p)=p`. The precise zero-reading formulas are
recovered only as `k -> 0`. For broad bins, `E[theta|E_k] ~ k/2+1/(2k)`.

| Threshold k | Mean alpha/sigma | Point / (F/sigma) | B_.95(k) | Lower bound / (F/sigma) |
|---:|---:|---:|---:|---:|
| 0 (limit) | .79788456 | 1.25331414 | 1.95996398 | .51021346 |
| .5 | .83072148 | 1.20377289 | 2.03940908 | .49033811 |
| 1 | .92466022 | 1.08147834 | 2.25413708 | .44362874 |
| 2 | 1.24711564 | .80185026 | 2.90234586 | .34454887 |
| 3 | 1.66659885 | .60002441 | 3.67111354 | .27239691 |

### Coarsening can remove conservative sampling coverage

Consider a fully specified reporting procedure: within `[-k,k]` report only
nondetection and use its posterior bound; outside, retain the signed exact
reading and use its posterior bound. Let `B=B_p(k)` and `t_*(theta)` be as above.
Because the bin posterior is a mixture of the exact-reading posteriors,
`Q_p(-k)<B<Q_p(k)`. The exact sampling coverage is

\[
C^{\rm bin}_p(\theta)=
\begin{cases}
C_p(\theta),&\theta\le Q_p(-k),\\
\Phi(\theta+k),&Q_p(-k)<\theta\le B,\\
\Phi(\theta-k),&B<\theta<Q_p(k),\\
C_p(\theta),&\theta\ge Q_p(k).
\end{cases}
\]

At `theta=B` the lower mass endpoint is included, hence the left-hand value.
Immediately above it, coverage approaches `Phi(B-k)`. For nominal `.95`:

| k | Coverage just above B |
|---:|---:|
| .5 | .93814777 |
| 1 | .89510391 |
| 2 | .81656341 |
| 3 | .74892590 |

These are exact analytic counterexamples to carrying the continuous-reading
coverage guarantee over to a censored instrument report. The posterior bound
still has its declared posterior probability under the stated reference measure.
The failed guarantee is repeated-sampling coverage.

If one observes *only* a binary detection/nondetection indicator and reports
all positive masses when detected, the corresponding right-hand coverage at
`B` is `1-L_k(B)` (for k=1: `.89567260`). A sampling-valid binary benchmark
can use the acceleration endpoint solving `L_k(theta_c)=delta`, when
`L_k(0)>delta`; when `L_k(0)<=delta`, even every positive acceleration can be
rejected on the rare nondetection event without exceeding the allowed error
rate. This illustrates the information cost and discreteness of binary data;
it is not a unique required reporting policy.

## Scope and integration recommendations

1. The exact-Y coverage theorem is a positive practical result under explicit
   assumptions. It does not require Monte Carlo calibration.
2. The point limitation and lower-bound guarantee should be presented together.
3. Keep precise signed readings separate from rounded values, thresholded
   nondetections, and magnitude-only observations. Their likelihoods differ.
4. Do not extrapolate the coverage theorem to uncertain force, unknown variance,
   force-direction uncertainty, vector magnitudes without direction, or arbitrary
   correlations. Those require fresh analysis.
5. The normal-tail bound is an available benchmark for this subproblem. Nothing
   here establishes a competitor meeting the full original composition brief.

## Numerical verification

`check_math.py` uses only the standard library and NumPy, with stable normal
log-tail evaluation and bisection. It checks the coverage identity, exact
quantile inversion, the bin mean against Gaussian quadrature, and the hybrid
coverage jump against a direct sum of the three observable regions. Output:
`results.json`. No stochastic simulation was needed for these results.

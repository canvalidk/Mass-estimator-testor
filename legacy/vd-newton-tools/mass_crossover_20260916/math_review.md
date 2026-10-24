# Independent mathematical contribution: the mass crossover

16 September 2026. Scratch contribution for integration by the lead investigator.
This examines the original Gaussian likelihood with flat positive-magnitude
measure `df d alpha dOmega`; it does not select a new measure or assert novelty.
Managed-context snapshot supplied for this task:
`canvalidk/VD-docs@8782a293297330f7a354d87ebe62792018fa8866`.
Local implementation inspected: `estimator.py` and `practical_formulas.py` in
`03_trace_and_evaluator/mass_estimator/`.

## Principal contributions

1. An exact scalar interpolation and shrinkage theorem, with explicit tolerance
   thresholds for approaching the familiar ratio.
2. A global logarithmic sensitivity bound, including signed readings at and
   below zero, which formalizes the removal of the raw ratio's singularity.
3. A proof that the full fixed-covariance estimator is real analytic in the
   readings, and an exact three-dimensional crossover when one observed vector
   is zero and the common direction is unknown.
4. A distinction between smooth interpolation, added information, and physical
   compatibility. Each requires a different claim.

## 1. Exact known-direction scalar reduction

Assume an externally supplied physical orientation, independent Gaussian
force and acceleration errors with known positive standard deviations, and
positive true magnitudes. Use signed projections `x,y` on that orientation.
Then the posterior factorizes into two separate positive-truncated Gaussians.
Writing

\[
\lambda(t)=\phi(t)/\Phi(t),\quad H(t)=t+\lambda(t),\quad
h_\sigma(x)=\sigma H(x/\sigma),\quad c=\sqrt{2/\pi},
\]

the original ratio-of-means readout is exactly

\[
\boxed{\widehat m=h_{\sigma_F}(x)/h_{\sigma_a}(y).}
\]

It is **identical** to taking the positive-truncated Gaussian mean of each
channel separately and dividing. There is no additional scalar information
extraction over that exact baseline. Possible advantages over raw ratios,
hard thresholds, or other methods need those methods and a loss criterion
specified. This equality does not remove the direction coupling of the full
vector problem or extend to correlated scalar errors.

### Limits and the important constants

\[
h_\sigma(0)=c\sigma,\qquad h_\sigma(x)\sim x
\quad(x/\sigma\to+\infty).
\]

Thus

| Regime | Leading mass readout |
|---|---|
| Zero force reading, positive resolved acceleration | `c sigma_F / y` |
| Positive resolved force, zero acceleration reading | `x / (c sigma_a)` |
| Both readings zero | `sigma_F / sigma_a` |
| Both positive and resolved | `x/y` |

The coefficients are `c=.79788456` and `1/c=1.25331414`. The transition is not
exactly between `sigma_F/y` or `x/sigma_a` and `x/y`; those expressions specify
scales, with coefficients determined by the construction.

For increasingly negative signed input, `H(t) ~ -1/t`; the model continues
smoothly rather than interpreting a negative reading as a positive magnitude.
This continuation is conditional on positive true quantities and the error
model. It does not establish that an extreme contradictory reading is likely.

## 2. Monotonicity and exact shrinkage toward the instrument ratio

For the truncated variable `X_t ~ N(t,1)` restricted to `(0,infinity)`,

\[
E[X_t]=H(t),\qquad H'(t)=\operatorname{Var}(X_t)
=1-\lambda(t)H(t),\qquad 0<H'(t)<1.
\]

Therefore the point increases strictly with signed force and decreases
strictly with signed acceleration. Define the positive-reading correction

\[
r(t)=\frac{H(t)}t-1=\frac{\phi(t)}{t\Phi(t)},\qquad t>0.
\]

It is strictly decreasing, since

\[
\frac{d}{dt}\log r(t)=-t-1/t-\lambda(t)<0.
\]

Let `t_F=x/sigma_F`, `t_a=y/sigma_a`, both positive, and
`m_0=sigma_F/sigma_a`. Then

\[
\frac{\widehat m}{x/y}=
\frac{1+r(t_F)}{1+r(t_a)},
\]

and the exact bracket is

\[
\boxed{\min(m_0,x/y)\le\widehat m\le\max(m_0,x/y).}
\]

Proof: `H` increases, but `H(t)/t` decreases, so if `t_F>t_a`,
`1<H(t_F)/H(t_a)<t_F/t_a`; reverse the roles otherwise.

In particular, equal signed signal-to-noise ratios give `mhat=x/y=m_0`
exactly at every signal strength. Agreement with `x/y` therefore cannot by
itself establish that the observations resolve a mass well.

### A quantitative resolved-reading criterion

If **both** positive signed signal-to-noise ratios exceed `T`, then

\[
\left|\frac{\widehat m}{x/y}-1\right|\le r(T).
\]

This is a deterministic approximation bound on the reported point, not a
confidence statement or bound on error relative to the unknown true mass.

| Guaranteed maximum point difference from raw ratio | Required minimum signed SNR in both channels |
|---:|---:|
| 10% | 1.46622942 |
| 5% | 1.76110605 |
| 1% | 2.37845044 |
| 0.1% | 3.11555104 |
| 0.01% | 3.73465910 |

For example `r(3)=.00147928`, or approximately `.148%`. These thresholds apply
to the known-direction independent scalar formula. Direction uncertainty can
create much larger, algebraically decaying corrections in the full estimator.

## 3. Smooth local response and a global sensitivity bound

At the crossover origin,

\[
H(t)=c+(1-c^2)t+\tfrac12c(2c^2-1)t^2+O(t^3).
\]

Consequently, near `t_F=t_a=0`,

\[
\frac{\widehat m}{m_0}
=1+\kappa(t_F-t_a)+O(t_F^2+t_a^2),\quad
\kappa=\frac{1-c^2}c=.4554295765.
\]

The log sensitivity is

\[
\ell'(t)=\frac{d\log H(t)}{dt}=\frac{H'(t)}{H(t)}.
\]

It satisfies the simple global bound

\[
\boxed{0<\ell'(t)\le1\quad\text{for every finite signed }t.}
\]

Here is a self-contained proof of the non-obvious part. The positive-truncated
Gaussian has log-concave survival function `S`; the Gaussian upper-tail hazard
is increasing. Since `S(0)=1`, log concavity gives
`S(x+y)<=S(x)S(y)`. Hence

\[
E[X_t^2]/2=\int_0^\infty\int_0^\infty S(x+y)\,dx\,dy
\le (E[X_t])^2.
\]

Thus `Var(X_t)<=H(t)^2`. The explicit variance formula also gives
`Var(X_t)<1`, and therefore
`ell'<=min(H,1/H)<=1`.

For fixed standard deviations and any two finite pairs of signed readings,

\[
\boxed{
|\log\widehat m(x',y')-\log\widehat m(x,y)|
\le\frac{|x'-x|}{\sigma_F}+\frac{|y'-y|}{\sigma_a}.}
\]

This bounds numerical response to changes in standardized inputs. In contrast,
the raw positive ratio has dimensionless log derivative `1/t`, unbounded near
zero. The bound does not assert statistical accuracy, resistance to an incorrect
noise model, or an optimal point estimate under mass error.

## 4. The full flat-measure estimator is analytic at fixed positive covariance

For `z=(f u,alpha u)`, `f,alpha>0`, `u` on the compact unit sphere, fixed SPD
covariance `Sigma`, precision `A=Sigma^{-1}`, and observed vector `y`, define

\[
N_g(y)=\int g(z)e^{-(z-y)^TA(z-y)/2}\,df\,d\alpha\,d\Omega,
\qquad\widehat m(y)=N_f(y)/N_\alpha(y).
\]

The quadratic form bounds the tails by an integrable Gaussian in `(f,alpha)`.
On compact sets of readings, every derivative in `y` is a polynomial in `z`
times an integrable Gaussian bound. More strongly, after factoring out
`exp(-y^TAy/2)`, the remaining integral is a Gaussian-dominated Laplace
transform. Thus both numerator and denominator are real analytic. Both are
strictly positive, so their ratio is real analytic at every finite reading.
The same reasoning gives smoothness as `Sigma` varies in the open SPD cone.

An exact local derivative under the normalized compatible-pair law `P` is

\[
\nabla_y\log\widehat m
=A\left(
\frac{\operatorname{Cov}_P(z,f)}{E_Pf}
-\frac{\operatorname{Cov}_P(z,\alpha)}{E_P\alpha}
\right).
\]

No threshold or mathematical discontinuity is required. This result concerns
the integral, not errors caused by finite numerical quadrature. It does not
imply a smooth or unique limit on the singular boundary `Sigma=0`.

### Local vector response at the noise-scale origin

For independent isotropic channels in dimension `d`, set
`X=F_observed/sigma_F`, `Y=a_observed/sigma_a`. Expanding the angular integral
and using `E[uu^T]=I/d` gives

\[
\boxed{
\frac{\widehat m}{m_0}
=1+\frac{\|X\|^2-\|Y\|^2}{2d}
+O((\|X\|+\|Y\|)^4).}
\]

Thus the full unknown-direction response is quadratic near the origin, unlike
the linear signed-projection response with known direction. The angle enters
the mass point only at fourth order. The second-order numerator and denominator
cross terms cancel in their ratio; this says nothing about all uncertainty
features of the compatible-pair distribution.

If `||X||=||Y||`, channel exchange combined with an orthogonal transformation
swapping `X,Y` leaves the model invariant. Consequently `mhat=m_0` **exactly**
for every angle and signal strength, including opposed observations.

## 5. An exact three-dimensional vector crossover on either zero-reading axis

Continue with independent isotropic channels and unknown common direction.
When the observed force vector is zero, the force magnitude factorizes as a
half-normal, so `E[f]=c sigma_F` exactly. Define

\[
t=\|\mathbf a\|/\sigma_a,\qquad
D(t)=e^{-t^2/2}\int_0^t e^{s^2/2}\,ds,\qquad
R_3(t)=\frac{2\Phi(t)-1}{D(t)},\quad R_3(0)=c.
\]

Then

\[
\boxed{E[\alpha]=\sigma_aR_3(t),\qquad
\widehat m=m_0\frac{c}{R_3(t)}.}
\]

`D(t)` is `sqrt(2)` times Dawson's integral at `t/sqrt(2)`. No special-function
dependency is needed for the verification: bounded Gaussian quadrature of the
displayed integral is adequate at the tested values.

**Derivation.** The three-dimensional angular average is
`sinh(tr)/(tr)`, where `r=alpha/sigma_a`. Its radial normalizer is

\[
Z(t)=\int_0^\infty e^{-r^2/2}\frac{\sinh(tr)}{tr}\,dr
=\frac{\sqrt{\pi/2}}t\int_0^t e^{s^2/2}\,ds.
\]

The first radial numerator is
`sqrt(pi/2) exp(t²/2)(2Phi(t)-1)/t`; taking their ratio yields `R_3`.
The observed-acceleration-zero case follows by exchanging channels:
`mhat=m_0 R_3(||F||/sigma_F)/c`.

The response is increasing because the radial family has a monotone likelihood
ratio in `t`: for angular moment-generating function `M(q)`,
`d log M(tr)/dt=r E_q[u_1]` increases with positive `r`, since
`E_q[u_1]>=0` and its derivative is a variance. Thus the two mass-axis curves
are monotone in their respective resolved companion magnitude.

The useful expansions are

\[
R_3(t)=c(1+t^2/6+t^4/72+O(t^6)),\quad t\to0,
\]

\[
R_3(t)=t-1/t-2/t^3+O(t^{-5}),\quad t\to\infty.
\]

Therefore the expected noise-scale formulas are recovered smoothly as the
companion resolves, but unknown direction creates an `O(t^{-2})` relative
correction. At `t=5`, `R_3=4.77906668`, appreciably below `5`, while the
known-direction `H(5)=5.00000149` is nearly `5`.

| Companion SNR | Known-direction H(t) | Unknown-direction R_3(t) |
|---:|---:|---:|
| 0 | .79788456 | .79788456 |
| 1 | 1.28759997 | .94192851 |
| 2 | 2.05524786 | 1.49143363 |
| 3 | 3.00443784 | 2.53658245 |
| 5 | 5.00000149 | 4.77906668 |
| 10 | 10.00000000 | 9.89789180 |

When one companion vector pins a nonzero direction while the other remains
at bounded signal-to-noise ratio, the conditional direction concentrates at
the companion orientation; the other channel reduces to the signed-projection
`h` formula. The exact axis calculation quantifies one important finite-SNR
correction to that statement.

## 6. Compatibility and the high-signal limit

For positive, nonzero, codirectional readings and covariance
`epsilon² Sigma_0`, the unique zero-residual compatible pair is the observed
pair. Ordinary interior Laplace concentration gives `mhat -> ||F||/||a||`,
with regular local corrections of order `epsilon²`. This statement assumes
the magnitudes stay away from zero.

For fixed incompatible readings, shrinking covariance concentrates the law
at nearest compatible pairs in the declared covariance metric, including
boundary points or multiple minimizers. It does not turn their norm ratio
into a verified Newtonian mass.

Two useful examples with equal isotropic channel variances:

- Equal-length opposed vectors give `mhat=1` in equal-unit normalization at
  every uncertainty by symmetry, while their minimum discrepancy diverges as
  uncertainty shrinks. The point's apparent agreement with the norm ratio is
  uninformative about compatibility.
- For normalized observations `F=(2,0)`, `a=(-1,0)` and channel standard
  deviation `epsilon`, the stronger force observation selects direction
  `(1,0)`. The weaker acceleration is explained near its positivity boundary:
  `E alpha ~ epsilon²`, `E f -> 2`, hence `mhat ~ 2/epsilon² -> infinity`.
  The raw norm ratio stays `2`. Direct laboratory calculations at
  `epsilon=.5,.25,.125` give points `9.8836,34.9296,131.3283`, approaching the
  leading values `8,32,128`.

## 7. Verification and useful scope

`check_math.py` and `math_results.json` are adjacent scratch files. The script
checks the Taylor expansion, scalar monotonic derivatives, shrinkage brackets,
tolerance inversion, and the exact three-dimensional zero-force formula
against the full laboratory quadrature. At companion SNRs `0,1,3,5`, the latter
agrees to less than `3e-11` in the mass point. The contradiction example is an
asymptotic diagnostic, not asserted to equal its leading formula at finite noise.

The strongest practical statement established here is that the declared model
provides a continuous, bounded-sensitivity scalar readout that joins exact
noise-scale limits to the resolved ratio, with calculable transition widths.
The vector law is smooth as well, but direction uncertainty materially changes
the transition. Showing lower error, more informative uncertainty, or practical
superiority still requires a specified experiment, comparator, and objective.

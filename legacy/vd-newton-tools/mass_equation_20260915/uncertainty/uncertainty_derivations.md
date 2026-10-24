# Mass uncertainty: a bounded first audit

Working contribution, 2026-09-15. This is scratch for the current investigation,
not a new estimator policy or a managed record.

## Scope and sources

Managed sources were read through the GitHub connector at
`canvalidk/VD-docs@8782a293297330f7a354d87ebe62792018fa8866`:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md`

The latter is an exploratory mathematical contribution; its optimizer and
uncertainty deductions are not assumed to be ratified policy. The present user
steering assumes Newton II holds; the older specification's rejection behavior
does not govern this investigation.

Local implementation inspected:

- `03_trace_and_evaluator/mass_estimator/README.md`
- `03_trace_and_evaluator/mass_estimator/estimator.py`
- `03_trace_and_evaluator/mass_estimator/practical_formulas.py`
- `03_trace_and_evaluator/mass_estimator/uncertainty_walkthrough.py`

These implement the flat positive-magnitude measure. The practical formulas
already distinguish induced mass-law intervals, local approximate standard
uncertainty, and exact logarithmic spread. No production or managed file was
changed, and no literature search was undertaken.

## 1. There are different laws here

Write `P` for the normalized law on compatible latent states `(f, alpha, u)`,
`M=f/alpha`, `A=E_P[alpha]`, and `B=E_P[f]`. The reported point is `B/A`.

The mass uncertainty implemented by `Estimate.quantile` is the pushforward of
**P** through `M`. It is not the sampling law of the estimator across repeated
measurement datasets. It also is not the law under which the point is a mean.

Define a second law by

\[
dQ_\alpha={\alpha\over A}\,dP.
\]

Then exactly

\[
E_{Q_\alpha}[M]={E_P[f]\over E_P[\alpha]}=\widehat m.
\]

For inverse mass, `1/m_hat=E_Qf[1/M]`, where `dQ_f=(f/B)dP`.
These are two differently weighted laws on the same latent states. The two mean
identities do not establish that mass uncertainty should be reported using
either tilted law. That is an additional decision.

### Explicit default-null example

Set `s=sigma_F/sigma_a` and `R=M/s`. Under independent isotropic Gaussian errors,
flat positive magnitudes, and both measured central vectors zero:

\[
p_R(r)={2\over\pi(1+r^2)},\quad
q_{\alpha,R}(r)={1\over(1+r^2)^{3/2}},\qquad r>0.
\]

Both correspond to the same reported point `m_hat=s`, but:

| Law | Mass mean | Mass median / s | Central 95% interval / s |
|---|---:|---:|---:|
| P | infinite | 1 | [0.0392901, 25.4517] |
| Q_alpha | s | 1/sqrt(3) = 0.577350 | [0.0250078, 4.38784] |

The Q_alpha CDF is `r/sqrt(1+r^2)`, so its quantile is
`p/sqrt(1-p^2)`. Its second moment is infinite. Indeed,

\[
\int_0^T {r^2\over(1+r^2)^{3/2}}dr
=\operatorname{asinh}T-{T\over\sqrt{1+T^2}}
\sim\log(2T)-1.
\]

Thus identifying a probability law under which the point is a mean does not
automatically supply a finite standard deviation for it.

## 2. The default infinite moments are not restricted to zero readings

This follows from an elementary boundary calculation. After integrating the
common direction, suppose the magnitude density has the form

\[
g(f,a)={1\over Z}f^p a^q H(f,a),\qquad f,a>0,\quad p,q>-1.
\]

For the finite-data positive definite Gaussian model, `H` is strictly positive
and smooth at each boundary away from the origin, with Gaussian decay. This
includes full covariance and anisotropy, and fixed or integrated direction.
The ratio density is

\[
p_M(m)={m^p\over Z}\int_0^\infty a^{p+q+1}H(ma,a)da.
\]

Changing variables to `f=ma` gives

\[
p_M(m)={m^{-(q+2)}\over Z}
\int_0^\infty f^{p+q+1}H(f,f/m)df.
\]

Dominated convergence therefore yields

\[
p_M(m)\sim C_0m^p\quad(m\downarrow0),\qquad
p_M(m)\sim C_\infty m^{-(q+2)}\quad(m\to\infty),
\]

where both constants are finite and strictly positive:

\[
C_0={1\over Z}\int_0^\infty a^{p+q+1}H(0,a)da,
\qquad
C_\infty={1\over Z}\int_0^\infty f^{p+q+1}H(f,0)df.
\]

Consequently,

\[
E_P[M^t]<\infty\quad\Longleftrightarrow\quad -(p+1)<t<q+1.
\]

At either equality the divergence is logarithmic. In particular, for the
default flat measure `p=q=0`, **every finite dataset with positive definite
Gaussian uncertainty has infinite `E_P[M]`, `E_P[1/M]`, and ordinary mass second
moment**. The acceleration-weighted law has exponents `p=0,q=1`, so its first
mass moment is finite and its second is infinite.

This is about exact tail behavior. A strongly informative dataset can have an
extremely small tail coefficient and a very narrow central interval. Infinite
variance by itself does not tell us that practical central uncertainty is large.
Nor does a finite-sample numerical average certify that an exact moment exists.

All fixed-order absolute log-mass moments remain finite in these power models,
because the log density has exponential tails. Hard lower bounds, singular
covariance, exact branches, or measures with different boundary behavior need
separate analysis; this theorem does not cover them automatically.

## 3. An exact family with the same null point and different uncertainty

Consider the declared magnitude measures

\[
d\Pi_k\propto(f\alpha)^k\,df\,d\alpha\,d\Omega,
\qquad k>-1.
\]

They are homogeneous under separate unit rescalings and symmetric under channel
exchange. They are improper before multiplication by the Gaussian likelihood;
the null normalized law is proper for precisely the displayed range of k.
These symmetry properties do not prove that the measures are physically
appropriate or uniquely required.

At measured zero-zero with independent isotropic noise, define
`nu=k+1`, `U=f/sigma_F`, and `V=alpha/sigma_a`. Then U and V are independent with

\[
U^2/2,V^2/2\sim\operatorname{Gamma}(\nu/2,1).
\]

Their means are identical, so for every k:

\[
\widehat m=s={\sigma_F\over\sigma_a}.
\]

But the density of `R=M/s=U/V` is

\[
p_R(r)={2\over B(\nu/2,\nu/2)}
{r^{\nu-1}\over(1+r^2)^\nu},\qquad r>0.
\]

Equivalently, `R^2` has the equal-degrees-of-freedom F law. With `L=log R`,

\[
p_L(\ell)={2^{1-\nu}\over B(\nu/2,\nu/2)}
\operatorname{sech}^{\nu}\ell,
\]

\[
E[e^{tL}]={\Gamma((\nu+t)/2)\Gamma((\nu-t)/2)\over\Gamma(\nu/2)^2},
\qquad |t|<\nu,
\]

\[
E[L]=0,\qquad
\operatorname{Var}(L)={1\over2}\psi_1(\nu/2).
\]

Here psi_1 is trigamma, the derivative of the logarithmic derivative of Gamma.
Thus the whole log spread is exactly known:

| k | Point / s | Var(log(M/s)) | Central 95% interval / s |
|---:|---:|---:|---:|
| -0.5 | 1 | 8.59866 | [0.00214850, 465.442] |
| 0 | 1 | 2.46740 | [0.0392901, 25.4517] |
| 1 | 1 | 0.822467 | [0.160128, 6.24500] |
| 2 | 1 | 0.467401 | [0.254500, 3.92928] |
| 9 | 1 | 0.110661 | [0.518698, 1.92791] |
| 99 | 1 | 0.0101007 | [0.821080, 1.21791] |

As k approaches -1 from above, log variance diverges. As k tends to infinity,
log variance is asymptotic to `1/(k+1)` and central intervals collapse around s.
Nothing has been added to the measured data. This is concentration imposed by
the declared weighting. Large k also favors increasingly large latent
magnitudes, of order `sigma*sqrt(k)`, despite central readings of zero.

For every fixed k, shrinking both noise scales by the same factor leaves the
entire null law of M unchanged. This passes the existing null-scaling test.
The test therefore does not determine the null spread.

The symmetric family has finite mass mean iff k>0 and finite mass variance iff
k>1. Choosing a k that makes moments exist is a substantive weighting choice,
not a numerical repair. Under P, even when the mass mean exists, it is generally
different from s; s remains the ratio-of-magnitude-means point.

## 4. Independent check of the bounded angular family

The main investigation proposes another family. In standardized polar
coordinates `f=sigma_F*r*sin(theta)`, `alpha=sigma_a*r*cos(theta)`, set

\[
\rho_\lambda=\exp(\lambda\sin2\theta),\qquad\lambda\in\mathbb R.
\]

At null, radius and angle remain independent, radius keeps its default law, and

\[
p_\theta(\theta)={e^{\lambda\sin2\theta}\over Z_\theta},
\quad Z_\theta=\int_0^{\pi/2}e^{\lambda\sin2\theta}d\theta.
\]

Symmetry about pi/4 keeps `m_hat=s` exactly. Direct variable changes give

\[
p_R(r)={e^{2\lambda r/(1+r^2)}\over Z_\theta(1+r^2)},
\qquad
p_L(\ell)={e^{\lambda\operatorname{sech}\ell}\over2Z_\theta\cosh\ell}.
\]

For each finite lambda, the measure is bounded above and below by positive
constants. Its ratio density retains positive constant density near zero and
the `1/r^2` high tail, so both the P mass mean and inverse-mass mean still
diverge. Nevertheless, Laplace expansion at `theta=pi/4` gives

\[
\operatorname{Var}(\log R)\sim1/\lambda\quad(\lambda\to+\infty).
\]

This demonstrates arbitrary central concentration while preserving both the
point and the default tail exponents. Boundary moment existence alone cannot
settle the uncertainty freedom.

Technical caveat: rho_lambda is smooth in the positive interior and at each
axis away from the origin, but has no direction-independent limit at the
origin when lambda is nonzero. The origin has zero integration measure, and
the current construction does not require continuity there. A globally smooth
alternative `exp(lambda*2xy/(1+x^2+y^2))`, with standardized x and y, preserves
exchange symmetry and boundary exponents while also changing the radial law.
No preference between these policies is asserted here.

## 5. What this first contribution establishes

1. A probability law and the readout functional need separate justification.
2. The mass-law mean need not be the reported point, even when it exists.
3. Exact moment existence and practical central spread are different questions.
4. The observed null point and the null scaling requirement leave large freedom
   in the uncertainty law. Channel exchange and unit covariance do too.
5. A narrowing caused by a weighting change is not evidence that the same
   measurements have acquired more information about mass.

The unresolved question is what experimental or operational commitment, if
any, selects the flat law, an alternative measure, or a constrained family.
That question should be resolved before treating one uncertainty width as a
data-determined consequence of the readout theorem.

## Verification

`verify_uncertainty.py` uses only NumPy and the standard library. It checks
normalization and exact log variance for six powers by independent numerical
quadrature, two closed CDFs, the non-null high-mass boundary coefficient by
direct magnitude integration, and the Q_alpha mean and second-moment
divergence. `verification.json` contains the numerical results. These verify
the displayed deductions; they do not calibrate coverage across repeated
experiments or establish a preferred measure.

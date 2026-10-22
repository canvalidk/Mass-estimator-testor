# Mass estimator: practical limits and uncertainty

**Date:** 2026-09-07. **Status:** local derivation and tested implementation,
awaiting transfer to VD-docs; not submitted. Newton II is assumed correct,
following the user's September 7 steering. No discrepancy rejection is applied.

## The result a user should receive

Supply measured force and acceleration vectors and their measurement covariance.
Return the agreed point and an interval:

$$
\boxed{\widehat m=\frac{E[f]}{E[\alpha]},\qquad
        [m_{.025},m_{.975}]\text{ for a 95\% interval}.}
$$

Here $f,\alpha>0$ are compatible latent magnitudes and $M=f/\alpha$.
Expectations and quantiles use the same law: Gaussian measurement likelihood,
flat $df\,d\alpha$, and uniform unknown common direction unless supplied.
The point is a ratio of means; the interval describes the possible latent
ratio. It need not be symmetric about the point. Its probability is conditional
on this model; repeated-experiment coverage is a separate property.

For a single finite measure of spread, use

$$
\boxed{s_{\log m}=\operatorname{SD}[\log(M/m_{\rm ref})].}
$$

Any fixed positive reference mass gives the same spread. This number is
dimensionless and unchanged on switching from mass to inverse mass. A spread
of 0.03 means approximately 3% uncertainty when the distribution is narrow.
The factor $e^{s_{\log m}}$ describes a logarithmic scale of variation;
it is not automatically a 68% interval or a claim of a lognormal distribution.

The calculator should perform the integration. Users should not need to
select special formulas for zeros or anti-alignment. The formulas below explain
what that calculation must do and provide useful independent checks.

## 1. Recovery of F/a: exactly what converges

Write the aligned measured vectors as
$\widetilde{\mathbf F}=F\mathbf u_0$,
$\widetilde{\mathbf a}=A\mathbf u_0$, with $F,A>0$.
For any fixed positive definite covariance shape $\Sigma_0$, set
$\Sigma=\varepsilon^2\Sigma_0$. Then

$$
\boxed{\widehat m\longrightarrow F/A\quad\text{as }\varepsilon\to0.}
$$

Every fixed interior mass quantile also tends to $F/A$, and
$s_{\log m}\to0$. This holds with unknown common direction and full
covariance. It is a limit about increasingly precise, nonzero aligned inputs.

**Alignment alone is insufficient.** Reducing the angle to zero while keeping
measurement uncertainty fixed gives the noisy aligned answer. It need not
give the raw magnitude quotient.

Example: 2D unknown direction, $F=2,A=1$, independent isotropic noise with
the same numerical standard uncertainty $\sigma$ in both channels:

| $\sigma$ | $\widehat m$ | Conditional 95% mass interval | $s_{\log m}$ |
|---:|---:|---:|---:|
| 1 | 1.520356 | [0.15595, 25.38838] | 1.22913 |
| 0.5 | 1.937850 | [0.73631, 13.07936] | 0.73587 |
| 0.25 | 1.999925 | [1.22141, 4.06331] | 0.30784 |
| 0.125 | 2.000000 | [1.54782, 2.70847] | 0.14245 |

An accurate central quotient does not by itself establish small uncertainty.

### Proof of the full-vector limit

Let $z(f,\alpha,\mathbf u)=(f\mathbf u,\alpha\mathbf u)$,
$y=(F\mathbf u_0,A\mathbf u_0)$, and
$R=(y-z)^T\Sigma_0^{-1}(y-z)$. The normalized law has density proportional
to $e^{-R/(2\varepsilon^2)}$.

The unique zero of $R$ is $(F,A,\mathbf u_0)$. Outside any neighborhood
of this point, $R$ has a strictly positive lower bound: direction space is
compact, and $R$ grows quadratically as either magnitude goes to infinity.
In local coordinates near the zero, the map $z$ has full column rank
$d+1$; the local normalizing integral is of order $\varepsilon^{d+1}$.
The probability outside that neighborhood is exponentially small compared
with this polynomial normalization.

The same bound with factors $f$ or $\alpha$, together with the Gaussian
tail bound at infinity, establishes convergence of their means to $F,A$.
Their ratio therefore tends to $F/A$. Continuity of $f/\alpha$ near
$(F,A)$ establishes convergence in probability and of all interior quantiles.
The logarithmic singularities at the axes are integrable; their neighborhoods
are separated from the unique zero and exponentially suppressed. Applying the
same bounds to $|\log(f/\alpha)|$ and its square gives convergence of the
logarithmic variance to zero. This argument does not apply at $F=A=0$.

## 2. A compact exact formula when direction is known

With independently known common orientation $\mathbf u$ and independent
isotropic Gaussian channels, use signed projections
$x=\mathbf u\cdot\widetilde{\mathbf F}$,
$y=\mathbf u\cdot\widetilde{\mathbf a}$. Define

$$
h_\sigma(x)=x+\sigma\frac{\phi(x/\sigma)}{\Phi(x/\sigma)}.
$$

Here $\phi,\Phi$ are the standard normal density and CDF. Then exactly

$$
\boxed{\widehat m=\frac{h_{\sigma_F}(x)}{h_{\sigma_a}(y)}.}
$$

Proof: perpendicular residuals are constant when the direction is fixed and
cancel on normalization. Each positive magnitude is independently normal
with its measured projection as location, conditioned to be positive.
Integrating that one-dimensional density gives $h$.

This formula also handles zero and negative projections. Do not substitute
measured norms and call the direction known merely because the observed arrows
align. The full unknown-direction estimator integrates over direction.

For positive projections, put $z_F=x/\sigma_F$, $z_a=y/\sigma_a$, and
$\delta(z)=\phi(z)/(z\Phi(z))$. Then

$$
\frac{\widehat m}{x/y}=\frac{1+\delta(z_F)}{1+\delta(z_a)}.
$$

If both standardized signals are at least $k>0$,

$$
\boxed{\left|\frac{\widehat m}{x/y}-1\right|
       \le\frac{2\phi(k)}{k}.}
$$

Indeed $\Phi(z)\ge1/2$, $\phi(z)/z$ decreases for $z>0$, and
$|\delta_F-\delta_a|/(1+\delta_a)\le\max(\delta_F,\delta_a)$.
The bound is 0.296% at $k=3$, and less than 0.000060% at $k=5$.
This bounds the point correction in the known-direction case, not interval
width or the error of an uncertainty approximation.

There is also an exact finite-noise aligned case with unknown direction:
if $F/\sigma_F=A/\sigma_a$, normalized magnitudes are exchangeable.
Thus $E[f]/\sigma_F=E[\alpha]/\sigma_a$ and
$\widehat m=\sigma_F/\sigma_a=F/A$. This includes arbitrarily weak
positive signals, so exact equality is not evidence of precise mass information.

## 3. The familiar uncertainty formula emerges

For nonzero aligned signals that are well resolved, under independent isotropic
Gaussian errors the leading local uncertainty is

$$
\boxed{s_{\log m}\simeq\frac{u_{m,\mathrm{local}}}{\widehat m}
\simeq\sqrt{\left(\frac{\sigma_F}{F}\right)^2+
             \left(\frac{\sigma_a}{A}\right)^2}.}
$$

The $\sigma$'s are component standard uncertainties along the common
direction. The derivation is simply
$d\log m=df/F-d\alpha/A$. At alignment, under isotropy, direction
perturbations are transverse and decouple from these first-order magnitude
perturbations. With isotropic cross-covariance $cI$, add
$-2c/(FA)$ inside the square root. The covariance propagation step is the
usual first-order calculation described in
[JCGM 100:2008, sections 5.1–5.2](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf).
Its applicability to this constrained estimator is derived here.

If the log distribution is sufficiently close to normal, a convenient local
95% approximation is

$$
[\widehat m\,e^{-1.96s_{\log m}},\quad
 \widehat m\,e^{+1.96s_{\log m}}].
$$

This uses the central estimator as an asymptotic center; it is not an exact
finite-noise quantile identity. For routine output the full quantiles are
already available and avoid choosing a universal signal-strength cutoff.

Example: $F=20\,\mathrm N,A=10\,\mathrm{m/s^2}$,
$\sigma_F=0.6\,\mathrm N,\sigma_a=0.2\,\mathrm{m/s^2}$.
The point is 2 kg. The simple relative spread is 0.036056; full numerical
$s_{\log m}$ is 0.036095. The approximate interval is [1.86354, 2.14645] kg;
the full interval is [1.86158, 2.14458] kg.

### Full covariance, without assuming away directional correlations

Let $B$ contain an orthonormal basis perpendicular to the aligned direction
$\mathbf u$. In coordinates $(df,d\alpha,d\eta)$, the measurement
Jacobian and mass gradient are

$$
J=\begin{pmatrix}\mathbf u&0&FB\\0&\mathbf u&AB\end{pmatrix},
\qquad g=(1/A,-F/A^2,0,\ldots,0)^T.
$$

The general leading local expression is

$$
\boxed{u_{m,\mathrm{local}}^2=g^T(J^T\Sigma^{-1}J)^{-1}g.}
$$

It follows by expanding the measurement map to first order and normalizing
the resulting Gaussian in the $d+1$ tangent coordinates. A smooth flat
magnitude/direction measure changes higher-order terms. Correlations between
longitudinal and transverse errors can affect this result; simply inserting
longitudinal variances in the scalar formula is not sufficient for arbitrary
covariance. The helper implements the matrix expression.

## 4. Why logarithmic spread and intervals remain usable at zero

For any finite data and positive definite Gaussian covariance under this
flat-magnitude construction, let $q(f,\alpha)$ be the normalized compatible
magnitude density after integrating direction. It is positive on the axes.
Changing variables gives

$$
p_M(m)=\int_0^\infty\alpha q(m\alpha,\alpha)\,d\alpha
      =\frac1{m^2}\int_0^\infty f q(f,f/m)\,df.
$$

Hence $p_M(m)\sim K/m^2$ as $m\to\infty$, with
$K=\int_0^\infty f q(f,0)df>0$. The first and second raw moments of
$M$ diverge: an exact ordinary standard deviation does not exist. This is
true even when the tail's weight is tiny and almost all probability is near
$F/A$. The local $u_m$ above describes the narrow central approximation;
it must not be labelled the exact SD of the full ratio law.

At the other end $p_M(0^+)$ is finite and positive. Consequently
$\log(M/m_{\rm ref})$ has exponentially decreasing tails on both sides
and finite moments. Its spread and every interior mass quantile remain usable.
This uncertainty concerns the latent mass ratio given these readings. The
sampling spread of the reported point over repeated experiments is another
quantity and is not calculated by calling it $\operatorname{SD}(M)$.

At the measured double zero, with independent isotropic channels, put
$s=\sigma_F/\sigma_a$. The two magnitudes are independent half-normals,
so

$$
\boxed{\widehat m=s,\quad
 p_M(m)=\frac{2s}{\pi(s^2+m^2)},\quad
 m_p=s\tan\frac{\pi p}{2},\quad s_{\log m}=\frac\pi2.}
$$

Thus the exact conditional 95% interval is
$[0.0392901s,25.4517s]$. The log density for $L=\log(M/s)$ is
$1/(\pi\cosh L)$, not a normal density. One derivation of its variance
uses $f/\sigma_F=|Z_1|$, $\alpha/\sigma_a=|Z_2|$: the Mellin transform
of $|Z|$ is $2^{t/2}\Gamma((1+t)/2)/\sqrt\pi$, whose log second
derivative at zero is $\pi^2/8$. Subtracting the two independent logs
doubles the variance to $\pi^2/4$.

Shrinking both instrument uncertainties by a common factor at measured zero
leaves this ratio distribution unchanged. The central number tracks their
ratio; it is not a uniquely identified physical mass from a null trial.
Exactly supplied zero vectors also leave every positive mass compatible.

For a single measured zero and a well-resolved other channel, independent
isotropic errors give the useful leading limits

$$
F=0:\quad \widehat m\sim\sqrt{2/\pi}\,\sigma_F/A;
\qquad A=0:\quad \widehat m\sim F/(\sqrt{2/\pi}\,\sigma_a).
$$

For known direction these follow immediately from $h_\sigma(0)$.
For unknown direction the resolved channel concentrates the common direction,
giving the same leading terms as both noise scales shrink together. These
are boundary limits to zero or infinity, rather than an arbitrary replacement
of a division by zero. The existing behavioral suite checks them.

## 5. Implementation and verification

The user-facing call is:

```python
from estimator import estimate, isotropic_covariance
from practical_formulas import summarize

result = summarize(estimate(
    [20, 0], [10, 0], isotropic_covariance(2, .6, .2),
    direction_order=768, ratio_order=8192))
# mass, lower, upper, probability, log_standard_deviation
```

Use the measured full covariance in place of the example isotropic covariance
when appropriate. Neither this call nor its summary applies a rejection gate.
Arbitrary misalignment is still reconciled conditional on the law and model.

For the existing ratio-angle representation $M=s\tan\theta$, the exact
mathematical expression used for logarithmic spread is

$$
\bar\ell=\int_0^{\pi/2}\log(\tan\theta)p_\theta(\theta)d\theta,
\qquad
s_{\log m}^2=\int_0^{\pi/2}(\log\tan\theta-\bar\ell)^2p_\theta(\theta)d\theta.
$$

The new helper integrates the same piecewise linear density used by the
quantile routine. Fourth-power substitutions at both endpoints handle the
integrable logarithmic singularities without introducing a mass cutoff.
Original direction/angle quadrature remains numerical: refinement is required
for new extreme inputs. The test at the resolved example checks refinement
from 512/4096 to 768/8192; null log spread is checked against $\pi/2$.

Ten new tests cover: alignment versus zero noise; full-vector convergence;
the known-direction formula and explicit error bound; exact equality at equal
standardized signals; null intervals/log spread; reciprocal and unit symmetry;
correlated scalar propagation; general-covariance conditioning and rotation;
the resolved-signal approximation and refinement; invalid approximation inputs.
Proofs establish the stated limits; finite test cases check their numerical
implementation. No claim of a universally simplest closed full-vector formula
or of certified numerical error for all inputs is made.

Verification on September 7: **67 tests passed, no failures, errors, or skips**,
including all ten new tests. Existing gate tests remain historical contract
checks, as recorded in the steering note; they do not gate the new readout.

Files: `03_trace_and_evaluator/mass_estimator/practical_formulas.py`,
`test_practical_formulas.py`, `practical_examples.py`, and saved
`results/practical_examples.json`. Full suite status is in `results/tests.json`.

## Provenance and scope

Managed context: `canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`,
resolved in this task on September 7. Relevant discovery used `catalog.yaml`
and the mass-estimator collection. Basis: the already read
`Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
and `mass_estimator_additional_derivations_2026-09-05.md` in that directory.
The working spec's earlier rejection policy is superseded for this investigation
by the user's local `mass_estimator_law_assumed_correct_steering_2026-09-07.md`.

Local active scope: `estimator.py`, existing behavioral tests and uncertainty
walkthrough, plus the new files listed above. This is a scoped derivation from
the agreed law, not an exhaustive search of VD materials or a reanalysis of
the separate published-data trials. Mathematical limit/bound and uncertainty
derivations here are new local analysis. JCGM is cited only for the conventional
first-order propagation step. Managed context has not been edited.

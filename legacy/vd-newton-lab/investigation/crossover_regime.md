# Mass estimation between unresolved and resolved measurements

**Ninth investigation contribution — 16 September 2026.**

## Result and practical meaning

The proposed equation supplies a mathematically smooth transition from
noise-scale mass estimates to the force/acceleration ratio for nonzero,
well-resolved compatible (codirectional) measurements. This
property holds for the full compatible-vector construction, including at
measured zeros. In the known-direction case we can prove both how the
transition regularizes the ratio and a bound on its sensitivity to readings.

The first sampling comparisons support a **bounded candidate application**:
one well-resolved vector and one partially resolved vector. For example, with
true force signal-to-noise 8 and acceleration signal-to-noise 1, the full 3D
equation recovers mass within a factor of two in 79.22% of simulated trials;
the measured-magnitude ratio succeeds in 59.07%. At acceleration signal-to-noise
3 the ranking reverses: 93.64% versus 95.89%. Some other smooth rules also
outperform the equation in parts of the scalar transition region.

Thus smoothness is established, and an advantage over simple division is
demonstrated in specific settings. General superiority, a unique continuous
solution, and experimental readiness are not established. These conclusions
concern the declared model and readout; they do not select its reference measure.

## 1. Scope and distinctions

The latent compatible pair is

$$\mathbf F^*=f\mathbf u,\qquad \mathbf a^*=\alpha\mathbf u,
\qquad f,\alpha>0,\quad \|\mathbf u\|=1.$$

The investigated law and point are

$$dP\propto e^{-Q/2}\,df\,d\alpha\,d\Omega,
\qquad \widehat m=E_P[f]/E_P[\alpha].$$

Newton II is assumed correct. No observation is rejected. Measurement errors
are Gaussian with supplied, known positive-definite covariance; the two
channels address the same object and time in an inertial frame. The numerical
comparisons specialize to independent errors and, for vectors, isotropic
noise. Correlated or anisotropic performance is outside this comparison.

Keep three questions separate:

1. **Continuity in observed data:** does a small reading change cause a small
   output change at fixed instrument uncertainties?
2. **Sampling accuracy:** across repeated noisy experiments at fixed true
   force and acceleration, how often is the mass estimate useful?
3. **Conditional uncertainty:** what distribution of latent mass accompanies
   this point under the chosen model, and how does its interval cover truth?

An observed signal-to-noise ratio and a true signal-to-noise ratio are different.
The analytical transition thresholds below use observed readings; the simulation
tables are indexed by true signals, which practitioners do not know exactly.

## 2. Exact scalar transition

Assume the common physical orientation is supplied independently of these
readings. Let $F,a$ denote their **signed projections** on that orientation.
With independent errors, put

$$H(t)=t+\frac{\phi(t)}{\Phi(t)},\quad
c=\sqrt{2/\pi},\quad s=\frac{\sigma_F}{\sigma_a}.$$

Then the exact equation is

$$\boxed{\widehat m=s\frac{H(F/\sigma_F)}{H(a/\sigma_a)}.}$$

This is exactly the ratio of the two separately computed positive-truncated
normal means. In this restricted model a joint implementation cannot extract
additional information over that exact factorized calculation. Substituting
observed vector lengths for the projections is not the same calculation.

### Limits and response near zero

$$H(0)=c,\qquad H(t)\sim t\ (t\to+\infty),\qquad
H(t)\sim1/|t|\ (t\to-\infty).$$

For positive resolved companion readings, the limiting mass formulas are
$c\sigma_F/a$ at zero force and $F/(c\sigma_a)$ at zero acceleration.
At two zero readings the point is $s$. These constants follow from the
Gaussian/flat-positive law; the more general phrases $\sigma_F/a$ and
$F/\sigma_a$ describe scales, not universally mandated numerical values.

The expansion

$$H(t)=c+(1-c^2)t+\tfrac12c(2c^2-1)t^2+O(t^3)$$

shows that partly resolved signal is retained immediately; it is not replaced
by a fixed zero-channel value up to a threshold.

### What the smoothing does

For **positive observed** $t_F=F/\sigma_F,t_a=a/\sigma_a$,

$$\boxed{\min(s,F/a)\le\widehat m\le\max(s,F/a).}$$

Proof: $H$ is increasing, while $H(t)/t=1+\phi(t)/(t\Phi(t))$ is decreasing
for $t>0$. Apply both inequalities to the larger and smaller standardized
readings. Thus the correction pulls the positive raw ratio toward the
instrument-uncertainty ratio. This controls noise sensitivity but can also
understate a genuinely very large mass. The bracket does not apply to negative
raw readings, although the estimator remains positive there.

Equal positive standardized readings give $\widehat m=F/a=s$ even when both
are tiny. Agreement with division therefore does not demonstrate identification.

### A global scalar stability bound

At fixed uncertainties, for any two finite signed reading pairs,

$$\boxed{|\Delta\log\widehat m|
\le |\Delta F|/\sigma_F+|\Delta a|/\sigma_a.}$$

Let $X_t$ have density proportional to $e^{-(x-t)^2/2}$ on $x>0$.
Differentiating its exponential-family mean gives
$H'(t)=\operatorname{Var}(X_t)=1-\lambda(t)H(t)$, where
$\lambda=\phi/\Phi$; hence $0<H'<1$. Its survival function has increasing
hazard and is log-concave, so $S(x+y)\le S(x)S(y)$. Integrating gives
$E[X_t^2]/2\le E[X_t]^2$, or $\operatorname{Var}(X_t)\le H(t)^2$.
Therefore

$$0<\frac{d\log H}{dt}\le\min(H,1/H)\le1.$$

Integrate along the two coordinates to obtain the bound. At zero the log
sensitivity is $(1-c^2)/c=0.455430$ per standard uncertainty. The raw positive
ratio has log sensitivity $1/t$, which diverges as its reading approaches zero.
This theorem concerns numerical/data sensitivity, not a bound on mass error.

### When does the scalar correction become small?

If both positive observed standardized readings are at least $T$, then

$$\left|\frac{\widehat m}{F/a}-1\right|
\le r(T),\qquad r(T)=\frac{\phi(T)}{T\Phi(T)}.$$

This follows by writing the quotient of the two corrections as
$(1+r(t_F))/(1+r(t_a))$. These are sufficient bounds, not detection cutoffs:

| Maximum relative difference from raw division | Sufficient observed SNR in both channels |
|---|---:|
| 10% | 1.4663 |
| 5% | 1.7612 |
| 1% | 2.3785 |
| 0.1% | 3.1156 |

The displayed thresholds are rounded upward. Closeness to a noisy raw ratio
is not a guarantee of closeness to true mass.

## 3. The full vector transition

### Smoothness for general fixed covariance

Write the measured pair as $y$ and a compatible latent pair as
$v=(f\mathbf u,\alpha\mathbf u)$. Apart from a factor depending only on $y$,
each defining integral has integrand

$$w(f,\alpha)e^{-v^T\Sigma^{-1}v/2+y^T\Sigma^{-1}v},$$

with $w=1,f,$ or $\alpha$. Positive definiteness supplies a Gaussian tail bound;
on every compact set of readings it dominates the integrand and all its
derivatives. Expanding the reading-dependent exponential about any finite
reading gives a locally absolutely integrable power series, with termwise
integration justified by the same Gaussian bound. Thus the integrals are real
analytic, finite and positive. Consequently
their ratio is real analytic at every finite reading, including measured
zeros. This statement fixes a positive-definite covariance; it does not include
singular zero-uncertainty limits. No global vector Lipschitz constant is asserted.
At fixed incompatible observed vectors, reducing the uncertainties does not
generally make the answer approach the ratio of their measured lengths.

### A useful exact 3D reduction

For independent isotropic errors define standardized observed vectors
$X=\widetilde{\mathbf F}/\sigma_F$, $Y=\widetilde{\mathbf a}/\sigma_a$,
and $p=\|X\|^2,q=\|Y\|^2,b=X\cdot Y$. Integrating direction and radial
amplitude yields

$$\boxed{\frac{\widehat m}{s}=
\frac{\int_0^{\pi/2}\sin\theta\,e^{h(\theta)^2/2}\,d\theta}
{\int_0^{\pi/2}\cos\theta\,e^{h(\theta)^2/2}\,d\theta},\quad
h^2=p\sin^2\theta+q\cos^2\theta+2b\sin\theta\cos\theta.}$$

This uses the radial reduction in the seventh contribution. It makes the
direction dependence explicit and permits fast repeated-sampling checks.
For example, at observed standardized lengths 4 and 1, changing their angle
from 0 to 90 to 180 degrees changes $\widehat m/s$ from 2.993 to 4.456 to
6.765, while the measured-length ratio remains 4. The exact values are in the
reproduction output. Retaining their relative direction changes the inference.

### Exact transition along a zero-vector axis

At $\widetilde{\mathbf F}=0$, let $t=\|Y\|$ and define

$$D(t)=e^{-t^2/2}\int_0^t e^{z^2/2}\,dz,\qquad
R_3(t)=\frac{2\Phi(t)-1}{D(t)},\quad R_3(0)=c.$$

Then

$$E[f]=c\sigma_F,\quad E[\alpha]=\sigma_a R_3(t),\qquad
\boxed{\widehat m=s\,c/R_3(t).}$$

Swapping the channels gives the reciprocal axis formula. To derive it,
integrate the direction average $\sinh(rt)/(rt)$ against $e^{-r^2/2}$.
Its normalization $Z$ obeys
$tZ(t)=\sqrt{\pi/2}\int_0^t e^{z^2/2}dz$, while its first-moment numerator is
$\sqrt{\pi/2}e^{t^2/2}\operatorname{erf}(t/\sqrt2)/t$.

| Observed nonzero-channel SNR $t$ | Known-direction mean $H(t)$ | Unknown-direction axis mean $R_3(t)$ |
|---|---:|---:|
| 0 | .797885 | .797885 |
| .5 | 1.009160 | .831823 |
| 1 | 1.287600 | .941929 |
| 2 | 2.055248 | 1.491434 |
| 3 | 3.004438 | 2.536582 |
| 5 | 5.000001 | 4.779067 |

In particular,
$R_3(t)=c[1+t^2/6+t^4/72+O(t^6)]$ near zero and
$R_3(t)=t-1/t+O(t^{-3})$ at large positive $t$.
The unknown-direction correction persists much longer than the scalar
positive-boundary correction. The scalar 1% threshold cannot be exported
to this vector axis. Known direction is additional information, so this
table is explanatory, not a same-information performance comparison.

Both weak standardized vectors give
$\widehat m/s=1+(\|X\|^2-\|Y\|^2)/6+O((\|X\|+\|Y\|)^4)$ in 3D.
Equal standardized lengths give the point $s$ exactly at every angle;
the accompanying uncertainty can still change greatly.

## 4. Fair comparisons in repeated measurements

### Design and metrics

Scalar study: true standardized signals independently range over
$\{.25,.5,1,2,3,5,8\}$, 49 settings, 500,000 IID Gaussian reading pairs
per setting. Vector study: 16 settings, 32,768 IID 3D reading pairs per
setting, true vectors codirectional. The same error draws are reused across
methods and settings. There is no pooling of settings into a claimed overall
winner and no selection of trials by their observed signal.

Primary comparison: the fraction of **all** estimates in $[m/2,2m]$.
Nonpositive, infinite and undefined outputs count as failures. We also report
mean squared log error for everywhere-positive rules, median estimate/truth,
and uncertainty coverage where applicable. Factor-of-two success is an
illustrative operational criterion, not a universal definition of usefulness.

The scalar JSON adds capped squared-log loss with cap $(\log10)^2$, assigning
the maximum to invalid outputs; that cap is an explicit comparison convention.
It labels positive-subset raw-ratio summaries as conditional. No ordinary
mean squared error is claimed for Gaussian raw division: its denominator
has density at zero and its squared reciprocal is nonintegrable.

### Scalar comparisons

Besides the existing ratio of means, compare the ratio of channel transforms:

- **Floor:** $g(t)=\max(t,c)$, continuous with a derivative kink.
- **Algebraic:** $g(t)=(t+\sqrt{t^2+4c^2})/2$, smooth and positive.
- **Posterior geometric:** $\exp(E[\log f]-E[\log\alpha])$ under the same
  flat-positive posterior. This minimizes posterior expected squared log
  loss; it is a different readout and need not satisfy the composition theorem.
- **Raw:** signed division, with failures retained.

The first two are declared illustrative controls with matching zero and
resolved limits, not claims of published complete mass estimators. They were
specified before examining results and their constants were not fitted.

With true force SNR 8, percentages within a factor of two are:

| True acceleration SNR | Equation | Floor | Algebraic | Geometric | Raw | Raw invalid |
|---|---:|---:|---:|---:|---:|---:|
| .25 | 8.99 | 0 observed | 15.85 | 35.10 | 14.94 | 40.12 |
| .5 | 48.51 | 65.91 | 44.54 | 68.08 | 28.95 | 30.88 |
| 1 | 79.34 | 83.39 | 70.86 | 74.74 | 52.51 | 15.89 |
| 2 | 89.45 | 80.82 | 89.72 | 79.98 | 80.26 | 2.26 |
| 3 | 94.18 | 92.17 | 95.94 | 89.25 | 92.17 | .13 |

Individual proportion Monte Carlo standard errors are at most .071 percentage
point. The output also gives paired difference standard errors. Zero observed
successes for the floor is not a proof of zero success probability.

The criterion matters: at (8,2), algebraic smoothing has slightly higher
factor-of-two success, while the equation has smaller squared-log error
(.18704 versus .19033). At (8,.25), the raw ratio has higher factor-of-two
success but much worse capped log loss (3.23990 versus 1.90256), including its
40.12% invalid outputs. No single ranking describes both properties.

### Full 3D comparisons

The competitors receive the same two observed vectors and instrument sigmas:

- **Norm ratio:** $\|\widetilde{\mathbf F}\|/\|\widetilde{\mathbf a}\|$.
- **Tube:** the same likelihood with radial reference weight
  $(f^2/\sigma_F^2+\alpha^2/\sigma_a^2)df\,d\alpha\,d\Omega$, and the same
  ratio-of-means readout. This is the previously investigated alternative
  measure, a sensitivity comparator rather than the user's selected law.
- **Positive profile:** minimize Gaussian discrepancy over positive mass and
  the common latent vector. This is a positivity-constrained TLS-type fit.
  Its unattained zero/infinity optima count as failures; no arbitrary finite
  replacement is supplied.

The JSON additionally reports a fixed floor applied to the two measured norms.
No tuned family of floors was optimized. The profile objective in standardized
mass $r=m/s$ is $(p-2br+qr^2)/(1+r^2)$; its positive interior stationary point
obeys $br^2+(q-p)r-b=0$. When $b\le0$, the optimum is at an open boundary,
apart from the nonidentified tie $p=q,b=0$. Full-support Gaussian noise gives
boundary outcomes positive probability, so a boundary-valued profile rule has
infinite unconditional squared-log loss even when a finite high-signal sample
misses those outcomes. Its finite-positive frequency and factor-of-two success
remain useful comparisons. With true force SNR 8:

| True acceleration SNR | Equation | Norm ratio | Tube | Positive profile | Equation minus norm ratio, percentage points ± MC standard error |
|---|---:|---:|---:|---:|---:|
| .25 | 8.94 | 3.31 | 9.73 | 14.82 | +5.63 ± .19 |
| .5 | 47.99 | 18.89 | 49.08 | 28.53 | +29.09 ± .30 |
| 1 | 79.22 | 59.07 | 79.18 | 51.81 | +20.15 ± .24 |
| 2 | 89.00 | 89.13 | 88.68 | 79.35 | −.13 ± .16 |
| 3 | 93.64 | 95.89 | 93.54 | 91.41 | −2.25 ± .10 |

At (8,1), the equation also reduces mean squared log error from .49607 to
.28816. At (8,3), its squared-log error is worse (.13767 versus .11952).
The clear advantage over the norm ratio at (8,.5) and (8,1) does not distinguish
the flat reference from the nearby tube alternative.

These results also prevent an overly optimistic weak-signal claim: at (8,.25),
only 8.94% of equation outputs fall within a factor of two, and its median
estimate is .2740 of truth. Improving on the norm ratio does not make this a
well-measured mass. The eighth contribution addresses useful lower-bound
information in such cases.

The equal-signal diagonal needs care: when both true SNRs are equal, true
mass equals the instrument ratio, favouring rules that shrink toward it.
For instance the scalar floor beats the equation at (.25,.25), but that
cannot establish recovery of arbitrary masses from two weak channels.

### Uncertainty is still a separate obligation

The full-vector central 95% conditional interval covered truth in 93.98%
of trials at (8,.25), 97.35% at (8,1), and 95.39% at (8,3). The first value's
Monte Carlo standard error is .131 percentage point, so this study does not
support universal 95% sampling coverage. These intervals remain conditional
posterior statements. Finite positive points do not certify accurate masses
or calibrated intervals.

### Follow-up: convergence and the physical degrees of freedom

The table holding force SNR at 8 while changing acceleration SNR is **not a
vanishing-noise experiment at fixed physical mass**. It changes the true mass
as well as acceleration resolution. Its finite-noise rankings do not establish
an accuracy gap that persists in the exact-measurement limit.

A direct follow-up holds the physical vectors at `(8,0,0)` and `(3,0,0)`, so
the true mass remains `8/3`, and decreases both component noise standard
deviations from 1 to `1/k`. Each setting uses 65,536 independent Gaussian
vector pairs; methods and precision settings reuse the same standardized
noise draws. Neither method is supplied the true direction or mass.

| Precision multiplier k | Force / acceleration SNR | Our squared-log error | Norm-ratio squared-log error | Our excess error relative to norm ratio |
|---:|---:|---:|---:|---:|
| 1 | 8 / 3 | .137940 | .119557 | 15.38% |
| 2 | 16 / 6 | .034738 | .031688 | 9.63% |
| 4 | 32 / 12 | .008138 | .007959 | 2.25% |
| 8 | 64 / 24 | .002003 | .001991 | .586% |
| 16 | 128 / 48 | .000498707 | .000497888 | .164% |

These finite-sample results support a shrinking gap, not a persistent
high-precision disadvantage. A strict finite-noise ranking can remain while
the numerical difference tends to zero.

There is also a local analytical explanation. Write the true pair as
`F0 u, A0 u`, with positive fixed magnitudes, and the observed errors as
`epsilon sigma_F xi, epsilon sigma_a eta`. For independent isotropic noise,
both estimators have the same first-order log-error expansion:

$$
\log\widehat m-\log m_0
=\epsilon\left(\frac{\sigma_F}{F_0}\xi\cdot u
 -\frac{\sigma_a}{A_0}\eta\cdot u\right)+O_P(\epsilon^2).
$$

For the norm ratio this follows by differentiating log vector length.
For the compatible-pair point, the exact angular integral concentrates around
its interior maximum. Its exponent is even about that maximum; the sine and
cosine moment integrals have the same local multiplicative correction. Their
ratio therefore approaches the positive profile/TLS mass, whose first-order
expansion is the expression above. Boundary terms are exponentially small
along this fixed, nonzero compatible-pair limit.

The angular dependence is a second-order difference in that limit. With
observed standardized squared lengths `p,q` and cosine `c` close to 1, the
positive profile solution satisfies

$$
\log\frac{m_{\rm profile}}{m_{\rm norm}}
=\frac{p-q}{p+q}(1-c)+O((1-c)^2).
$$

Since the angle is of first order in noise, `1-c` is of second order.
This explains how directional corrections can affect finite-noise accuracy
without preventing convergence. Fixed incompatible observed vectors with
shrinking *declared* errors are a different limit: they need not recover the
observed length ratio.

**No extra physical direction is introduced.** The latent model uses one
common unit vector `u`, with two directional degrees of freedom, and two
positive magnitudes. Its four parameters are equivalent to an unknown mass
plus a three-component true acceleration. Measured vectors can disagree in
direction because their measurement errors have transverse components; the
latent physical vectors are constrained to be parallel. Fixing the generating
direction to an axis loses no generality for this isotropic, rotation-invariant
comparison and does not supply that direction to either estimator.

This checks a particular suspected structural error. It does not establish
that flat magnitude weighting is uniquely warranted or that the composition
readout minimizes repeated-experiment error. Those inferential choices remain
separate from the physical parameter count.

The follow-up script and JSON are in
`.tools/mass_crossover_20260916/convergence_audit.py` and
`convergence_audit.json`. Seed: 2026091601. The script uses the exact angular
moment integrals, split at the **observed** maximum, with 96 nodes per segment.
Checks at 192 nodes and against the previous integration implementation
agreed within `6.0e-12` relative error on the checked subsets. No estimator
or earlier simulation code was changed.

## 5. What could responsibly be offered to a user?

A supported description is:

> A positive, smooth Newton-II-constrained mass readout, using declared
> measurement uncertainties, which bridges noise-dominated and resolved
> observations without an arbitrary detection switch. Under independent
> isotropic Gaussian noise, initial comparisons show substantial improvement
> over measured-magnitude division in some partially resolved regimes.

Its distinctive practical case would need to include the whole required
package: shared direction, covariance handling, positive finite output,
the composition-motivated readout, and an appropriate uncertainty output.
Smoothness alone neither selects the flat measure nor proves that its point
is optimal for a specified use. Its scalar boundary function is established
positive-normal-mean mathematics: [Chang, Shinozaki and Strawderman,
equation 1.4](https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf)
explicitly reproduce it and attribute it to Katz (1961). Our quotient is a
mass construction using that ingredient, not a mass formula attributed to Katz.

The next practical validation must name an experiment and acceptable error,
then test its actual force/acceleration calibration, correlations, offsets,
direction uncertainty and reporting protocol. This contribution supplies a
testable domain claim and counterexamples to broader claims; it does not yet
justify deployment or identify a best rule for every experiment.

## 6. Reproduction, verification and provenance

Reproduction scripts:
[scalar comparisons](crossover_scalar_checks.py) and
[vector comparisons and analytical checks](crossover_vector_checks.py).
They require Python and NumPy only. From the workspace root, for example:

```powershell
python -B 03_trace_and_evaluator/mass_estimator/investigation/crossover_scalar_checks.py --output .tools/mass_crossover_20260916/scalar_results.json
python -B 03_trace_and_evaluator/mass_estimator/investigation/crossover_vector_checks.py --output .tools/mass_crossover_20260916/vector_results.json
```

The actual run used bundled Python 3.12.14 and NumPy 2.3.5. Seeds are 20260916
(scalar) and 2026091601 (vector). Output includes script and baseline hashes,
all settings, metrics and numerical checks. There are 24,500,000 scalar and
524,288 vector setting-trials; error draws are shared across settings, so these
are not counts of independent experiments across the entire grid.

Checks include scalar identities and transition thresholds; the exact 3D axis
formula; reciprocity; six fixtures against the baseline direction/radial
implementation (worst relative point difference $3.74\times10^{-10}$); and
angular quadrature refinement on 128 observations per vector setting (worst
point/CDF difference $2.98\times10^{-14}$). The posterior-geometric scalar
comparison was checked against higher-order integration (maximum log-mean
difference $8.61\times10^{-8}$) and the exact half-normal log moment.
Independent analytical review checked the radial formula, profile boundaries,
axis formula and scalar sensitivity proof. These scoped checks do not certify
all extreme floating-point inputs.

Managed context was resolved from `canvalidk/VD-docs/main` to
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
The catalog inventory and the following discussion sources at that commit
ground this contribution (paths under `Newton-analysis/_dscn_mass_estimator/`):

- `README.md`: current discussion scope and retained/historical distinctions.
- `inertial_mass_estimator_full_working_spec_d1.md`: integral construction;
  its rejection wording is superseded for this investigation by the steering below.
- `mass_estimator_law_assumed_correct_steering_2026-09-07.md`: retain conditional
  inference under the law, including discrepant observations.
- `mass_readout_composition_uniqueness_2026-09-06.md`: readout characterization.
- `mass_estimation_zero_zero_two_cases.md` and
  `mass_estimator_originality_research_2026-09-05.md`: zero distinctions and
  scalar mathematical ancestry, with the limitations examined in later local work.

Local active sources: `estimator.py`, `practical_formulas.py`, and investigation
contributions 1–8, particularly [operational warrant](operational_warrant.md)
and [weak acceleration information](weak_acceleration_information.md).
These active local contributions are not assumed to have been ingested into
the managed collection. No managed files were changed. This was a focused
transition study, not an exhaustive resurvey of the managed or inbox material.

External primary-source checks were limited to the Chang et al. manuscript
and the abstract of [Golub and Van Loan, *An Analysis of the Total Least Squares
Problem*](https://epubs.siam.org/doi/10.1137/0717073). The latter supports the
general errors-in-both-channels comparison family; the positivity-constrained
profile specialization here is derived directly from the stated likelihood.
No complete published predecessor claim follows from this search.

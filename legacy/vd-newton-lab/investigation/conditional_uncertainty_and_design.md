# Calibrated mass uncertainty and learning excitation strength

**Working contribution 5 — 15 September 2026.** Continues
[joint inference and uncertainty](joint_inference_and_uncertainty.md).
Newton II remains a premise. The mass point is still a separate readout;
this contribution develops a confidence benchmark and explains two ways the
experiment's treatment of unknown excitation affects uncertainty.

## What has advanced

1. **An exact conditional likelihood-ratio confidence construction is now
   derived and implemented.** It retains 95% coverage at every fixed positive
   mass and excitation under the stated Gaussian model, including true zero
   excitation. Computing its mass set requires one scalar calibration root
   followed by a quadratic inequality.
2. **Learning excitation strength explains the missing uncertainty factor.**
   An integrated Gaussian excitation model with unknown strength has the same
   informative interior mass optimum as the fixed-vector profile, but changes
   its curvature. In our earlier example it supplies precisely the missing
   factor that changes local mass variance from \(5/N\) to \(8/N\).
3. **The way measurements are grouped matters.** Repeated readings of a stable
   vector reduce the cost of estimating its unknown direction and magnitude.
   Treating those readings as unrelated latent vectors discards a real equality
   supplied by the experimental protocol.

The new confidence construction is an adaptation of established conditional
likelihood-ratio inference, not a claim of a new statistical principle. The
results narrow the investigation without selecting a final probability measure
or replacing the ratio-of-means equation.

## 1. The sampling model and the guarantee

There is one positive mass and \(N\) fixed unknown acceleration vectors:

$$
\mathbf F_i\sim N(m\mathbf v_i,\sigma_{F,i}^2I_d),\qquad
\mathbf a_i\sim N(\mathbf v_i,\sigma_{a,i}^2I_d).
$$

Measurement errors are independent across channels and trials; the positive
variances are known. For the formulas here, the instrument ratio
\(s_0=\sigma_{F,i}/\sigma_{a,i}\) is common. Stack the standardized observations
into \(X,Y\in\mathbb R^K\), where \(K=Nd\), and set \(r=m/s_0\). For a candidate,

$$
R_r=\frac{X-rY}{\sqrt{1+r^2}},\qquad
S_r=\frac{rX+Y}{\sqrt{1+r^2}}.
$$

At the true ratio, \(R\sim N(0,I_K)\) independently of \(S\), whose mean contains
all the unknown excitation. Thus conditioning on \(S\) leaves a known residual
distribution even when the excitation is arbitrarily weak or zero.

The target guarantee is

$$
P_{m,V}\{m\in\mathcal C(X,Y)\}=.95
$$

for every fixed positive \(m\) and fixed \(V\). It is a sampling guarantee for
a random set. It does not assign 95% posterior probability to the observed set
and does not specify a probability law for the composition theorem's readout.

## 2. A likelihood-ratio statistic with exact conditional calibration

Let \(t=\|R\|^2\), \(s=\|S\|^2\), and \(H=R\cdot S\). The symbol \(s\) in this
section is squared fitted length, not the instrument ratio \(s_0\).

Compare the candidate's profile likelihood with the best fit over the larger
family of all signed slopes, including the limiting vertical line. Twice the
log-likelihood difference is

$$
\Lambda=t-\lambda_{\min}
=\frac12\left[t-s+\sqrt{(t-s)^2+4H^2}\right], \tag{1}
$$

where \(\lambda_{\min}\) is the smaller eigenvalue of the two-column Gram
matrix. Its value is unchanged by rotating \((X,Y)\) into \((R,S)\).

For fixed nonzero \(S\), define

$$
U=H^2/s\sim\chi^2_1,\qquad W=t-U\sim\chi^2_{K-1}.
$$

These variables are independent. Direct algebra gives, for \(c>0\),

$$
\boxed{\Lambda\le c\quad\Longleftrightarrow\quad
\frac Uc+\frac W{c+s}\le1.} \tag{2}
$$

Consequently the conditional CDF is

$$
F_K(c;s)=P\!\left\{\frac Uc+\frac W{c+s}\le1\right\}.
$$

Choose \(c_K(s)\) to satisfy \(F_K(c_K(s);s)=.95\). Accepting a mass when
\(\Lambda\le c_K(s)\) has exactly 95% conditional coverage and therefore exactly
95% unconditional coverage. The critical value depends on the observed fitted
length, with no unknown excitation parameter left in its calibration.

The limiting cases are useful checks:

$$
c_K(0)=\chi^2_{K,.95},\qquad
c_K(s)\longrightarrow\chi^2_{1,.95}\quad(s\to\infty).
$$

For \(K>1\) the curve is decreasing. For \(K=1\), \(W=0\) and the critical value
is \(\chi^2_{1,.95}\) everywhere. At exactly \(S=0\), use \(\Lambda=t\) and the
\(K\)-degree threshold directly; no arbitrary accepted singleton is needed.

This is the known-covariance conditional likelihood-ratio construction in
[Moreira, *A Conditional Likelihood Ratio Test for Structural Models*
(2003)](https://drphilipshaw.com/Moreira-2003.pdf), specialized to our standardized
Gaussian vector model. Theorem 1 gives the conditioning principle; Section 3,
equation (3), gives the corresponding eigenvalue statistic. The mass
interpretation, positive-domain restriction, and experimental comparisons here
are stated separately from that antecedent.

### The comparison family is an explicit choice

The reported candidates are positive masses. The statistic uses signed slopes
as a larger geometric reference family; it does not assert that negative mass
is physically possible. Restricting the accepted set to \(r>0\) preserves
coverage at every positive true mass.

Profiling the reference fit over **positive masses only** gives a different,
smaller statistic. Reusing the threshold above would generally be conservative.
Exact recalibration of that alternative would also depend on the candidate's
position relative to the positive-mass boundaries. That variation has not been
implemented here.

## 3. The complete mass set is a quadratic inequality

Write

$$
A=\|X\|^2,\quad B=\|Y\|^2,\quad C=X\cdot Y,\quad
L=\lambda_{\max}=\frac{A+B+\sqrt{(A-B)^2+4C^2}}2.
$$

Since \(t+s=A+B\), equation (1) implies \(\Lambda+s=L\). Along this curve,

$$
F_K(L-s;s)=P\!\left\{\frac U{L-s}+\frac WL\le1\right\}
$$

decreases with \(s\). If \(L\le\chi^2_{K,.95}\), every mass passes. Otherwise
there is a unique \(b\in(0,L)\) satisfying \(F_K(L-b;b)=.95\), and the set is

$$
\boxed{\mathcal C_{\rm CLR}/s_0
=\{r>0:(A-b)r^2+2Cr+(B-b)\ge0\}.} \tag{3}
$$

Thus calibration needs one scalar root per dataset, followed by ordinary
quadratic inversion. A large \(L\) alone does not imply a narrow mass set;
the full inequality still determines what directions are distinguishable.

The possible sets are bounded intervals, one-sided intervals, two rays, the
whole positive domain, or empty. Finite endpoints satisfying equality are
included; zero and infinity are parameter boundaries. An isolated singleton
cannot occur: it would require \(b\) to equal an interior maximum \(L\), whereas
calibration gives \(b<L\).

### How it compares with the earlier constructions

The earlier score test accepts a candidate automatically when \(H=0\). At
a stationary **maximum** of the full residual, that can be a poor mass fit.
Here the likelihood ratio is \(t-s\) when \(H=0,t>s\), so the candidate must
still pass a calibrated comparison.

The following are deterministic examples in dimensionless mass, with \(K=3\):

| Data summaries | Score set | Split-residual set | Conditional LR set |
|---|---|---|---|
| Aligned: \(A=64,B=16,C=32\) | \([1.2251,4.0387]\) | \([1.1485,4.6706]\) | \([1.2177,4.0905]\) |
| Opposed: \(A=B=100,C=-90\) | \([0.93265,1.07222]\) | Empty | Empty |
| Orthogonal: \(A=B=100,C=0\) | All positive | Empty | All positive |

The orthogonal example explains what CLR measures: every mass direction fits
equally poorly, so the likelihood comparison supplies no preference among them.
The split-residual construction additionally responds to absolute residual size.
These are different uncertainty choices despite their shared coverage level.

None of these sets suppresses the ratio-of-means point or acts as a rule for
rejecting Newton II. An empty confidence set is an outcome of the specified
mass procedure; its operational use remains a separate decision.

## 4. Learning excitation strength changes the likelihood width

One possible physical model for independently generated excitations is
\(\mathbf v_i\stackrel{iid}{\sim}N(0,\tau^2I_d)\), with unknown \(\tau\).
For homogeneous channel scales, integrating the vectors gives each standardized
coordinate pair covariance

$$
I_2+\gamma u_ru_r^T,\qquad
u_r=\frac{(r,1)^T}{\sqrt{1+r^2}},\qquad
\gamma=\frac{\tau^2}{\sigma_a^2}(1+r^2)\ge0.
$$

Here \(\gamma\) measures total standardized excitation strength. Gaussian
measurement errors do not imply this distribution of true excitations.

Let \(M=K^{-1}\begin{pmatrix}A&C\\C&B\end{pmatrix}\), the empirical second-moment
matrix **about the specified zero mean**, and \(q_r=u_r^TMu_r=\|S_r\|^2/K\).
It is not an automatically centered sample covariance. Relative to the
zero-excitation likelihood,

$$
\ell(r,\gamma)-\ell_0
=\frac K2\left[-\log(1+\gamma)+\frac{\gamma}{1+\gamma}q_r\right].
$$

Maximizing over the unknown strength gives

$$
\widehat\gamma(r)=(q_r-1)_+,\qquad
\boxed{\ell_{\rm strength}(r)-\ell_0=
\begin{cases}
\frac K2(q_r-1-\log q_r),&q_r>1,\\
0,&q_r\le1.
\end{cases}} \tag{4}
$$

The fixed-vector profile instead has \(\ell_{\rm fixed}(r)=Kq_r/2+\text{constant}\).
At an informative interior optimum both select the same mass direction. Their
curvatures differ exactly by

$$
\frac{-\ell_{\rm strength}''}{-\ell_{\rm fixed}''}
=\frac{q_{\widehat r}-1}{q_{\widehat r}}. \tag{5}
$$

For the earlier \(d=3,m=2,\sigma_F=\sigma_a=1,\|\mathbf v_i\|=1\) example,
\(q_{\widehat r}\to8/3\); this factor tends to \(5/8\). Inverse curvature therefore
changes from \(5/N\) to \(8/N\), matching the actual leading sampling variance.
Learning excitation strength produces the correction from a declared
likelihood, without inserting that example's factor by hand.

### Why the true-null behavior improves

At true zero excitation, \(M=I_2+O_P(K^{-1/2})\). Thus
\(q_r-1=O_P(K^{-1/2})\) uniformly over directions. Fixed-vector profile
differences grow as \(\sqrt K\); equation (4) instead gives

$$
\ell_{\rm strength}(r)-\ell_0
=\frac K4(q_r-1)_+^2+o_P(1)=O_P(1).
$$

The learned-strength model does not accumulate an increasing preference from
those null fluctuations. It still cannot identify mass at \(\tau=0\).

This is an explanation of likelihood behavior, not a newly calibrated
posterior. A posterior would additionally require a reference or prior for
mass and strength. At the zero-strength boundary, ordinary regular-model
likelihood-ratio and Bayesian calibration arguments do not apply.

Under the fixed-vector experiment used in our simulations, (4) is evaluated as
a **working second-moment likelihood**. Its informative interior point and
leading curvature have the useful agreement above, without claiming that the
fixed vectors were actually sampled from a zero-mean Gaussian population.
If the Gaussian excitation distribution is to have literal physical meaning,
its mean, variability, and independence need support from the protocol.

## 5. Numerical comparison at positive and zero excitation

The new experiment uses 4,096 independent datasets per scenario, seed 20260917,
with \(d=3\). Known measurement variances are used throughout. Each fixed vector
has the declared magnitude and its direction cycles through the coordinate axes;
those vectors are not revealed to the estimator.

The table contrasts three likelihood-ratio procedures. **Only the conditional
column has the exact finite-sample guarantee derived above.** The other columns
use the ordinary \(\chi^2_{1,.95}\) cutoff as explicit comparisons. These are not
the Bayesian quantile intervals reported in contribution 4.

| Fixed latent vectors | Acceleration signal/noise | Fixed-vector LR, ordinary cutoff | Learned-strength LR, ordinary cutoff | Conditional LR |
|---:|---:|---:|---:|---:|
| 8 | 1 | 87.99% | 95.53% | 95.24% |
| 64 | 1 | 88.04% | 94.95% | 94.87% |
| 512 | 1 | 87.60% | 94.87% | 94.87% |
| 1 | 0 | 82.89% | 97.39% | 94.95% |
| 64 | 0 | 26.25% | 95.53% | 94.97% |
| 512 | 0 | 13.84% | 94.78% | 94.87% |
| 4,096 | 0 | 8.96% | 94.78% | 95.00% |

The reference true mass ratio is 2. At zero excitation its value does not affect
the data distribution, but it remains the fixed candidate whose coverage is
measured. Near 95% coverage the Monte Carlo standard error is 0.34 percentage
points per scenario. The full output also includes the full-residual, score,
and split methods, a weak-acceleration ratio-16 case, and repeated-group designs.

The close-to-95% learned-strength entries are not proof of calibration. Its
null limiting statistic is not \(\chi^2_1\); it has an atom at zero. An independent
large simulation of that limiting distribution gives about **94.65%** coverage
for the unrestricted reference alternative. This illustrates why numerical
proximity alone is insufficient. The conditional construction has a proof that
does not require a strong-signal approximation.

In the null rows, the fixed-vector LR's 95th percentile grows from 6.61 at one
vector to 411.11 at 4,096 vectors. The learned-strength LR's corresponding
percentile stays between 2.75 and 3.95 in this run. The ordinary cutoff is 3.84.
This directly exhibits the difference between increasing spurious evidence
and bounded null likelihood variation.

## 6. Repeating a stable vector changes the uncertainty cost

Suppose each of \(N\) latent vectors has \(J\) independent repeated readings.
Combine those readings into their sufficient group means first:

$$
\overline{\mathbf F}_i\sim N(m\mathbf v_i,\sigma_F^2I_d/J),\qquad
\overline{\mathbf a}_i\sim N(\mathbf v_i,\sigma_a^2I_d/J).
$$

With known covariance, the remaining within-group scatter supplies no further
likelihood information about the means. This reduction uses the equality of
the true vectors within each group; it is not available for genuinely changing
excitations.

Put \(e=N^{-1}\sum_i\|\mathbf v_i\|^2>0\) and
\(V_0=\sigma_F^2+m^2\sigma_a^2\). In the regular local regime, the sampling
variance of the mass center is

$$
\frac{V_0}{JNe}+\frac{d\sigma_F^2\sigma_a^2}{J^2Ne^2}.
$$

Its ratio to the joint Jeffreys posterior's smaller local variance is

$$
\boxed{\kappa_J=1+
\frac{d\sigma_F^2\sigma_a^2}{J V_0e}.} \tag{6}
$$

These are leading local normal variances, not assertions about finite exact
second moments of every possible mass statistic.

For the running example, \(\kappa_J=1+0.6/J\). Keeping \(J\) fixed and adding
new vectors gives the following limiting coverage for the joint posterior's
nominal 95% interval:

| Repeated readings per vector \(J\) | Local variance ratio | Limiting coverage |
|---:|---:|---:|
| 1 | 1.6000 | 87.87% |
| 2 | 1.3000 | 91.44% |
| 4 | 1.1500 | 93.24% |
| 8 | 1.0750 | 94.13% |
| 16 | 1.0375 | 94.57% |
| 64 | 1.0094 | 94.89% |

This is an asymptotic calculation, not a universal recommendation for a repeat
count. If \(N\) stays fixed and \(J\) grows, the gap vanishes when at least one
vector is nonzero. At actual zero excitation, the standardized group means
have the same distribution for every \(J\); repetition creates no mass
information and does not repair the growing-\(N\) null pathology.

### What this says about information use

With a fixed reading budget \(B_0=NJ\) and comparable excitation,

$$
I_{m\cdot V}=\frac{B_0e}{V_0},\qquad
V_{\rm local}=I_{m\cdot V}^{-1}
+\frac{d\sigma_F^2\sigma_a^2}{J B_0e^2}.
$$

The Fisher benchmark is unchanged by this allocation, but the extra uncertainty
from estimating unknown vectors decreases with repeats. For the example, the
ratio of inverse local variance to that Fisher benchmark is 62.5% at \(J=1\)
and about 93.0% at \(J=8\).

This is a **benchmark ratio for the nuisance-estimation cost**, not a percentage
of raw information discarded by the ratio-of-means equation. The growing-
nuisance experiment need not attain the fixed-dimensional Fisher bound. It
becomes an actionable information comparison when the protocol really allows
the same latent vector to be measured repeatedly, or supplies other structure
that helps estimate it.

For a fair estimator comparison, supply the same observations, their correct
covariance, and the same warranted sharing structure to every method. Comparing
algorithms that have been given different equality constraints would confound
information supplied with information used.

## 7. Where this returns us to the original equation

The composition theorem still selects a readout only after a joint probability
law is supplied. CLR supplies a calibrated comparison target for mass
uncertainty, not that joint law. Equation (4) identifies a promising way to
model and learn excitation strength, but does not choose its prior or prove
that a particular physical experiment follows it.

The next direct question for the ratio-of-means equation is therefore precise:
construct its joint law with an explicit, learnable excitation model; retain
the declared aggregation operation; and compare both its point and uncertainty
against the calibrated benchmark. The comparison should include reciprocity,
dependence on the strength prior, genuine repeated readings, new excitations,
and the zero-signal boundary. This is where tuning can acquire a defensible
purpose rather than merely improving a few familiar examples.

## Reproduction and scope

[conditional_checks.py](conditional_checks.py) contains the conditional CDF,
the one-root quadratic inversion, the learned-strength likelihood, and the
simulation. It uses NumPy and the Python standard library, importing the
preceding contribution's checks. From the workspace root:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/conditional_checks.py' --replicates 4096 --output '.tools/mass_equation_20260915/conditional_lr/results_4096.json'
```

The run covers **45,056 datasets across 11 scenarios**. Checks include gamma-CDF
identities, conditional-CDF limits, calibration-root residuals, agreement between
set membership and direct candidate tests, degenerate cases, and doubled-order
quadrature for the first 64 datasets in every scenario. The largest of those
CDF refinement differences was below \(6.7\times10^{-13}\).

The calibration root is bracketed by the proven \(\chi^2_1\) and \(\chi^2_K\)
quantile bounds. Quadrature resolves the chi-square transition and integrates
only up to an absolute standard-normal value of 8; the omitted probability is
at most \(1.25\times10^{-15}\). The integration scale follows \(c\), not \(c+s\).
Independent review caught a false-convergence failure at extreme scales and
the corrected implementation includes explicit regression checks for it.
The current implementation and quadratic roots remain research numerics, with
scoped validation rather than universal floating-point error bounds. The exact
coverage statement concerns the mathematical construction.

Machine outputs and independent supporting derivations are under
`.tools/mass_equation_20260915/conditional_lr/`. Outputs record hashes of this
script and its imported investigation script. No baseline estimator, formal
uncertainty policy, or managed document was changed.

VD-docs `main` was resolved again to
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
The [catalog](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/catalog.yaml),
[working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md),
and [readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md)
were read at that snapshot. Local scope was the preceding investigation and
the new derivations and checks. Moreira's primary paper was checked for the
specific methodological antecedent. This is not an exhaustive literature or
wider VD-context review; the new contribution remains local working material.

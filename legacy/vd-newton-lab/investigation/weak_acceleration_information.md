# Mass information when acceleration is hard to resolve

**Working contribution 8 — 15 September 2026.** A focused contribution to the
shared investigation, prompted by the user's correction that a resolved force
with an unresolved acceleration conveys information about the object's mass.
Newton II is assumed throughout. The original ratio-of-means equation and its
flat positive-magnitude weighting are retained.

## Findings

1. **There is a useful coverage theorem in a specified domain.** With known
   force and direction, a numerical signed Gaussian acceleration reading, known
   measurement variance, and flat positive acceleration weighting, the induced
   posterior 95% lower mass bound covers the true mass in at least 95% of
   repeated experiments. Its coverage rises toward 100% at weak acceleration.
2. **The finite point has a different limitation.** At fixed force and noise,
   as true acceleration vanishes, the point divided by the true mass tends to
   zero. Useful lower-bound information can coexist with poor relative point
   accuracy. The information must be described accordingly.
3. **The reporting protocol matters.** A numerical reading of zero and a report
   saying only that the reading is below a threshold have different likelihoods.
   Correctly modeling threshold reporting changes the equation's inputs and can
   remove the conservative coverage guarantee. An explicit example falls to
   about 89.51% coverage for a nominal posterior 95% lower bound.
4. **Uncertain force needs its own assessment.** Correlated Gaussian signed
   readings admit an exactly calibrated one-sided confidence-set benchmark.
   The original conditional bound showed empirical coverage consistent with
   or above 95% in the tested strong-force settings, but has no uniform 95%
   guarantee across all masses.
5. **The vector problem has a direction-free benchmark.** A full Gaussian
   residual gives an exactly calibrated confidence set without estimating a
   direction from the same data and then pretending it was known.

These results establish a restricted theoretical use and concrete boundaries.
They do not establish superiority of the complete estimator. The confidence
benchmarks below evaluate uncertainty; they do not replace its mass point or
claim to satisfy the entire composition-based brief.

## 1. Start with a fully specified experiment

For the first result let the true net-force magnitude \(F>0\) and its direction
be known exactly. Newton II implies a true acceleration \(a=F/m>0\) along that
direction. The instrument supplies its **signed projection**

$$
Y\sim N(a,\sigma^2),\qquad \sigma>0\text{ known}.
$$

The observed numerical value is \(y\). It may be negative or zero. It is not
the norm of a noisy vector, and the direction has not been inferred from this
acceleration reading. Exactly known force is a deliberate limiting experiment;
Sections 5 and 6 restore uncertainty in force and direction.

The original flat magnitude measure reduces to flat \(da\) on \(a>0\). With
\(t=y/\sigma\), standard normal density \(\phi\), and CDF \(\Phi\),

$$
p(a\mid y)=\frac{\phi((a-y)/\sigma)}{\sigma\Phi(t)},\quad a>0,
\qquad
E[a\mid y]=\sigma H(t),\quad
H(t)=t+\frac{\phi(t)}{\Phi(t)}.
$$

Consequently the existing point and a posterior lower mass bound are

$$
\widehat m=\frac{F}{\sigma H(t)},\qquad
L_p(y)=\frac{F}{\sigma Q_p(t)}, \tag{1}
$$

where

$$
Q_p(t)=t+\Phi^{-1}\!\left(1-(1-p)\Phi(t)\right)
\quad\text{and}\quad
P(M\ge L_p\mid y)=p. \tag{2}
$$

These are two different summaries of the same distribution. In this special
case \(1/\widehat m=E[1/M\mid y]\): the point corresponds to the posterior
mean of inverse mass. It is not the posterior mean of mass.

At the numerical reading \(y=0\),

$$
\widehat m=1.253314\,F/\sigma,\qquad
L_{.95}=0.510213\,F/\sigma. \tag{3}
$$

For \(F=1\,\mathrm N\), \(\sigma=0.01\,\mathrm{m\,s^{-2}}\), these are
125.33 kg and 51.02 kg. The second number expresses the one-sided information
the user identified. The uncertainty scale participates in an inference about
the object; its presence does not make the result merely an apparatus fact.

## 2. A coverage theorem for that lower bound

For \(p>.5\), let \(z_p=\Phi^{-1}(p)\). Equation (2) gives

$$
Q_p(t)>t+z_p
$$

for every finite \(t\). Thus the posterior upper acceleration endpoint is
larger than the ordinary one-sided Gaussian upper endpoint. The event
\(Y+z_p\sigma\ge a\), which has probability exactly \(p\), implies
\(\sigma Q_p(Y/\sigma)\ge a\). Equivalently,

$$
\boxed{P_{a}\{L_p(Y)\le F/a\}\ge p.} \tag{4}
$$

This is an unconditional repeated-experiment guarantee at every fixed
\(a>0\), not a guarantee conditional on observing a nondetection or selecting
particular runs. The posterior statement and the coverage statement both hold
here, for different reasons.

An exact coverage expression gives more information. Put \(\theta=a/\sigma\)
and define \(t_*\) by \(Q_p(t_*)=\theta\). The truncated normal family is
stochastically increasing in \(t\), so its quantile is increasing and this root
is unique. For example the ratio of two posterior densities is proportional
to \(\exp((t_2-t_1)a/\sigma)\), proving the ordering. Hence

$$
C_p(\theta)=\Phi(\theta-t_*)
=1-(1-p)\Phi(t_*)>p. \tag{5}
$$

It approaches 1 as \(\theta\downarrow0\) and \(p\) as \(\theta\to\infty\).

| True acceleration / noise SD | Coverage of posterior 95% lower mass bound |
|---:|---:|
| 0.25 | effectively 100% |
| 1 | 99.9282% |
| 2 | 97.3801% |
| 4 | 95.0469% |
| 8 | 95.0000% |

### What the comparison costs

The Gaussian candidate set based on \(a\le y+z_{.95}\sigma\) has exactly 95%
coverage. When that upper endpoint is positive its mass lower endpoint is
\(F/(y+z_{.95}\sigma)\), which is higher than (1). At \(y=0\) it is
\(0.607957F/\sigma\), or 60.80 kg in the example. It can therefore supply a
stronger lower restriction for a specified coverage target.

When \(y+z_{.95}\sigma\le0\), however, this confidence set is empty in the
positive-mass model. This is an outcome of the confidence procedure, not an
estimate of infinite physical mass or a rule for rejecting Newton II. The
original posterior remains proper with a finite point and bound. This is a
real difference in output and intended use; coverage alone does not select
the entire estimator.

## 3. What the finite point can and cannot say at weak acceleration

Hold \(F,\sigma\) fixed and let the true acceleration decrease. Writing
\(Y/\sigma=\theta+Z\), \(Z\sim N(0,1)\), gives

$$
\frac{\widehat m}{m_{\rm true}}
=\frac{\theta}{H(\theta+Z)}\longrightarrow0
\quad\text{almost surely as }\theta\downarrow0. \tag{6}
$$

For each finite \(Z\), \(H(Z)>0\) and is continuous. The true mass grows
without bound while the finite readout stays on the \(F/\sigma\) scale.
This proves a substantial relative point-error limitation, even in the
otherwise favorable known-force experiment.

For example at \(a/\sigma=.05\), the exact median of
\(\widehat m/m_{\rm true}\) is 0.06125; its central 90% sampling range is
[0.02787, 0.11769]. Its lower mass bound still has essentially 100% coverage.
These statements are compatible: the experiment excludes small masses much
more effectively than it distinguishes among large masses.

At \(y=0\), the likelihood as a function of mass is

$$
\mathcal L(m)\propto\exp\!\left[-\frac{F^2}{2\sigma^2m^2}\right].
$$

It increases toward a plateau at large mass. The induced posterior density
also contains the Jacobian \(F/m^2\) from the flat acceleration measure.
Finite posterior quantiles therefore include the contribution of that
reference weighting; they are not upper mass exclusions supplied by the
likelihood alone.

This limit varies the physical mass at fixed noise. It is not the limit of
repeated measurements at one fixed positive acceleration. For repeated
independent readings of the same acceleration, using their mean and
\(\sigma/\sqrt n\) makes the point converge to \(F/a\). The experimental
protocol determines which limit is relevant.

## 4. A nondetection report requires a different likelihood

Suppose the instrument supplies only the report
\(E_k=\{|Y|\le k\sigma\}\), with a known threshold \(k>0\). With the same flat
positive reference, the posterior of \(\theta=a/\sigma\) is

$$
p(\theta\mid E_k)
=\frac{\Phi(k-\theta)-\Phi(-k-\theta)}{k},\quad\theta>0. \tag{7}
$$

The denominator is exactly \(k\): integrate the Gaussian likelihood over
\(-k\le Y/\sigma\le k\), then use \(\Phi(v)+\Phi(-v)=1\).
Integration of the first moment yields

$$
E[\theta\mid E_k]
=\frac{(k^2+1)(\Phi(k)-1/2)+k\phi(k)}{k}. \tag{8}
$$

This tends to \(\sqrt{2/\pi}\) as \(k\downarrow0\), recovering the zero-reading
calculation. A finite detection bin does not justify inserting \(y=0\).

| Report | Point / \((F/\sigma)\) | Posterior 95% lower mass / \((F/\sigma)\) |
|---|---:|---:|
| Numerical \(y=0\) | 1.253314 | 0.510213 |
| Only \(|Y|\le\sigma\) | 1.081478 | 0.443629 |
| Only \(|Y|\le2\sigma\) | 0.801850 | 0.344549 |

### A checked loss of the coverage guarantee

Specify the complete reporting experiment: the instrument reports the bin
inside \([-k\sigma,k\sigma]\), and reports the exact signed reading outside.
Use the appropriate posterior for every report. Let \(B_k\) be the bin
posterior's 95% upper quantile of \(a/\sigma\).
The bin posterior is a positive mixture of exact-reading posteriors over
\(-k<t<k\), so \(Q_{.95}(-k)<B_k<Q_{.95}(k)\).

As the true \(\theta\) increases just beyond \(B_k\), every bin report stops
covering it. The outside-bin positive readings cover it, and negative
outside-bin readings do not. Therefore

$$
\lim_{\theta\downarrow B_k,\ \theta>B_k} C_{\rm report}(\theta)
=\Phi(B_k-k). \tag{9}
$$

For \(k=1\), \(B_k=2.254137\), giving **89.5104%**. For \(k=2\),
\(B_k=2.902346\), giving **81.6563%**. These are right-hand limits; at the
endpoint itself the bin still covers. They are derived from the reporting
likelihood, not from a simulation accident or an incorrectly substituted zero.

A simple guaranteed comparison is to use upper acceleration endpoint
\((k+z_{.95})\sigma\) whenever the bin is reported, retaining the usual
Gaussian endpoint outside it. This increases the upper endpoint relative to
having the exact reading, so coverage remains at least 95%. Its lower mass
coefficients inside the bins are 0.378093 and 0.274359 respectively.

**Practical requirement:** obtain the numerical reading and its uncertainty
when available. Otherwise the data interface must supply the actual threshold,
rounding, or censoring rule. The probability model must describe that report.
The particular bin example does not claim that every real instrument uses
this protocol, or that every posterior bound after censoring undercovers.

## 5. Restore uncertain force, including correlation

In a known positive direction let signed measurements have known covariance
and satisfy

$$
\begin{pmatrix}X\\Y\end{pmatrix}
\sim N\!\left(\begin{pmatrix}ma\\a\end{pmatrix},
\begin{pmatrix}s_F^2&c\\c&s_a^2\end{pmatrix}\right),
\quad a,m>0,\quad |c|<s_Fs_a.
$$

The original law has the bivariate Gaussian likelihood restricted to
\(f,\alpha>0\) with flat \(df\,d\alpha\). Correlation changes both its point
and its mass probabilities. Separately truncating the two marginals is valid
only when the channels are independent.

At the true mass,

$$
T_m=\frac{X-mY}{\sqrt{s_F^2+m^2s_a^2-2mc}}\sim N(0,1).
$$

Thus the one-sided set

$$
\mathcal C=\{m>0:T_m\le z_{.95}\} \tag{10}
$$

has exactly 95% coverage at every fixed positive mass and acceleration.
Its upward hull, starting at \(L=\inf\mathcal C\), has at least 95% coverage.
The set can be disconnected or empty. The reproduction code retains the
empty-set outcome separately from a zero lower restriction.

For a numerical zero acceleration and \(x>z_{.95}s_F\), the lower endpoint is

$$
L=\frac{c+\sqrt{c^2+s_a^2(x^2/z_{.95}^2-s_F^2)}}{s_a^2}. \tag{11}
$$

For independent channels this becomes
\(\sqrt{x^2-z_{.95}^2s_F^2}/(z_{.95}s_a)\). It explicitly shows why the force
must itself be resolved before this benchmark excludes small masses.
The exact boundary \(x=z_{.95}s_F\) requires separate handling when correlated;
it is implemented in the companion code.

This is Gaussian test inversion for a ratio, related to established Fieller
confidence-set methods; see [von Luxburg and Franz, 2009](https://arxiv.org/abs/0711.0198).
The one-sided specialization and its application here are directly derived
above. No historical novelty is claimed.

### A limited repeated-experiment comparison

The companion script uses 16,384 experiments per setting, the same random
errors for compared methods, and the declared known-direction law. Here
\(s_F=s_a=1\); entries describe **true** signal strengths, not selection on
observed readings. Correlation is zero in this table.

| True force / SD | True acceleration / SD | Posterior lower-bound coverage | Pivot lower-hull coverage | Median point / true mass |
|---:|---:|---:|---:|---:|
| 8 | 0.05 | 100% observed | 94.96% | 0.0607 |
| 8 | 0.25 | 100% observed | 94.99% | 0.2766 |
| 8 | 1 | 99.87% | 94.98% | 0.7719 |
| 8 | 4 | 95.03% | 94.97% | 1.0003 |
| 1 | 1 | 97.85% | 95.43% | 1.0054 |
| 0.25 | 1 | 89.68% | 95.01% | 2.8440 |

Monte Carlo standard error near 95% is about 0.17 percentage points. A row
with no observed failures has a 95% Wilson lower coverage endpoint of about
99.9766%, not zero statistical uncertainty. Full output includes these
intervals. The exact pivot-set guarantee is proved; fluctuations around 95%
in this finite check are expected.

At true strengths \((3,.5)\), correlations \(-.8,0,.8\) give median point/true
mass ratios of 0.4281, 0.4795 and 0.5493. All three posterior lower bounds
covered in every simulated run. That is scoped numerical evidence, not a
uniform theorem for correlated data.

The last table row illustrates a boundary that can also be proved. With
fixed positive acceleration and fixed positive-definite measurement covariance,
let true mass tend to zero, so force becomes unresolved. Each finite data pair
has a strictly positive posterior lower mass quantile, varying continuously
with the readings. Under coupled measurement errors the readings converge to
a finite random pair while the true mass tends to zero. Consequently coverage
of that positive posterior lower bound tends to zero. A global 95% claim is
therefore impossible for this rule. This argument does not refute the favorable
known-force theorem or establish failure throughout a strong-force domain.

## 6. Retain vector direction uncertainty when it is present

For full vector readings with known positive-definite joint Gaussian
covariance, write its blocks as
\(\Sigma_F,\Sigma_a,C=\operatorname{Cov}(\widetilde{\mathbf F},\widetilde{\mathbf a})\).
At a candidate mass use

$$
D_m=\widetilde{\mathbf F}-m\widetilde{\mathbf a},\qquad
V_m=\Sigma_F+m^2\Sigma_a-m(C+C^T).
$$

At the true mass, \(D_m\sim N(0,V_m)\), so

$$
\mathcal C_d=\{m>0:D_m^TV_m^{-1}D_m\le\chi^2_{d,.95}\} \tag{12}
$$

has exactly 95% coverage. This is the absolute residual construction already
available in the investigation, applied here as an uncertainty benchmark.
Its use does not suppress the original conditional point or test whether
Newton II should remain a premise.

For independent isotropic channels and the observed vector
\(\widetilde{\mathbf a}=0\), it reduces to all positive masses when
\(\|\widetilde{\mathbf F}\|^2\le q_d s_F^2\); otherwise it is

$$
\left[\sqrt{\frac{\|\widetilde{\mathbf F}\|^2-q_ds_F^2}{q_ds_a^2}},\infty\right),
\qquad q_d=\chi^2_{d,.95}. \tag{13}
$$

This is a full-vector version of the force-resolution requirement, without
inventing a known direction from the measured force. It is not necessarily the
most informative calibrated set; the earlier conditional-LR contribution
studies a different benchmark.

The original 3D unknown-direction integrator, with unit channel uncertainties,
gives the following for observed force length 8 and zero observed acceleration:

| Output | Mass in units \(s_F/s_a\) |
|---|---:|
| Original point | 9.86452 |
| Original posterior 95% lower bound | 3.88696 |
| Vector residual set lower endpoint | 2.68135 |

The differing bounds have differing justifications. This deterministic example
does not establish 95% sampling coverage for the vector posterior bound, nor
superiority from its larger endpoint. Establishing that domain is a next task.

## 7. What this contribution supports

There is now a precise positive result for one practical interpretation: the
original construction supplies a conservative lower mass bound in the
known-force, exact signed Gaussian reading model. It captures the information
that a small mass would have produced a more readily detectable acceleration.

There is also a precise limit: a finite central point on the apparatus scale
does not become relatively accurate for arbitrarily large true masses at fixed
force and acceleration resolution. Reporting the bound and its meaning is
essential to presenting what this experiment has established.

The next useful work is to characterize the force-SNR, acceleration-SNR and
correlation domain in which the original vector conditional bounds have
acceptable coverage, and to compare their strength against calibrated sets.
The reporting protocol, covariance calibration, known versus inferred direction,
and intended error criterion must be fixed first. This contribution does not
choose a replacement weighting or certify the complete method for deployment.

## Reproduction, verification, and provenance

[weak_acceleration_checks.py](weak_acceleration_checks.py) uses Python and
NumPy. From the workspace root:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/weak_acceleration_checks.py' --output '.tools/mass_weak_acceleration_20260915/results.json'
```

It checks known-force means and quantiles by independent positive-axis
quadrature; the censored likelihood's normalization, mean and coverage;
correlated scalar posteriors against the existing estimator; independent
Cartesian CDF integration in uncorrelated cases; reciprocal transformations;
confidence-bound roots and boundary cases; and vector lower quantiles under
refinement. Sampling checks comprise nine settings of 16,384 runs plus three
censoring settings of 300,000 runs. The largest checked scalar point difference
was \(1.92\times10^{-12}\) relatively; the CDF difference was
\(2.27\times10^{-14}\). Refining a 128-run subset per sampling setting changed
points by at most \(8.63\times10^{-14}\) relatively and CDFs by
\(2.96\times10^{-14}\). These are checked-case errors, not uniform integration
bounds for arbitrary observations, correlations or signal strengths.
An independent review also used conditional-Gaussian Cartesian integration
on 30 scalar cases, including correlations from -0.98 to 0.98. It agreed with
the radial implementation within \(4.70\times10^{-14}\) relatively for points
and \(2.04\times10^{-14}\) absolutely for mass CDFs.

The two main analytical strands received independent derivations. Review notes,
additional checks and machine output are in
`.tools/mass_weak_acceleration_20260915/`. Numerical checks supplement the
proofs; simulated coverage is not promoted to a theorem. The exact-known-force
limit is derived separately because the empirical integrator requires strictly
positive-definite covariance.
The standalone check script uses fixed root and integration settings for its
declared scenarios; it is not a general production API. Its output identifies
runtime versions and hashes the research helpers it imports.

Managed `main` was resolved to
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
unchanged from the reconstruction earlier in this conversation. Source
selection followed `catalog.yaml`, the discussion README and these canonical
files at that snapshot:

- [Working specification, especially Sections 9–12 and 17–18](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md).
- [Law-assumed-correct steering](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md).
- [Composition theorem and its scope](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md).

Active local sources include `../estimator.py`, `../practical_formulas.py`,
the [Gaussian-conditioning analysis](gaussian_conditioning.md),
[weighting and calibration](weighting_and_calibration.md),
[joint inference](joint_inference_and_uncertainty.md),
[conditional confidence construction](conditional_uncertainty_and_design.md),
and the [operational analysis](operational_warrant.md) with their relevant
integrators. This is an active local contribution, not a managed policy ruling
or a submission to VD-docs. No managed source or baseline estimator was edited.
External reading was confined to the cited ratio-confidence paper's abstract
and bibliographic information; this is not an exhaustive literature or priority
review.

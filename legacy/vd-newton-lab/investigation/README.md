# Mass estimator investigation

**Current synthesis — 21 September 2026; extended 23 September.** We have a defined conditional
estimator, mathematical limits on what it can claim, and a reusable simulation
tester. A fresh operating-range study now distinguishes its relative gains
from absolute accuracy: its clearly demonstrated advantages over the norm
ratio occur in conditions that still fail the chosen 90% accuracy benchmarks.
Strong signals make it accurate, with similar performance from simpler
alternatives. This is a scoped result, not a universal ordering of estimators.
The [23 September contribution](prior_factors_and_sequential_updating.md)
scopes it further and makes two choices in the law explicit; see
[readout regimes, prior factors and updating](#readout-regimes-prior-factors-and-updating).

Start here for the argument, evidence and next question. Use the
[comparison guide](comparison_guide.md) for equations, assumptions and running
instructions. The linked mathematical notes supply derivations for particular
claims; they are not successive replacements for the estimator. This opening
replaces the former chronological contribution list. The original uniqueness
and uncertainty proofs remain [below](#first-result).

## The problem we are assessing

One object has a fixed positive mass. Its true force and acceleration obey
\(\mathbf F^*=m\mathbf a^*\), while both measured vectors contain error.
Newton II is assumed true throughout this empirical investigation. A large
measurement discrepancy describes the error needed to reconcile the readings;
it is not a reason to reject the physical law or silently drop the trial.

The construction has three distinct steps:

1. **Specify the experiment and conditional law.** Gaussian measurement errors
   and a declared latent reference measure give a joint law \(P\) for compatible
   force magnitude \(f\), acceleration magnitude \(\alpha\), and direction.
   The original choice is flat physical magnitude/direction volume
   \(df\,d\alpha\,d\Omega\).
2. **Compute an output from that law.** The existing point is
   \(\widehat m=E_P[f]/E_P[\alpha]\). The induced distribution of
   \(M=f/\alpha\) supplies conditional intervals. These are different summaries;
   the point is not \(E_P[M]\). It averages magnitudes before division,
   without selecting a single minimum-\(Q\) compatible pair.
3. **Assess the procedure against external truth.** Repeated simulated
   experiments reveal point error, interval coverage, prediction performance
   and failures under a stated design. Numerical integration error, Monte Carlo
   error in those performance estimates, and uncertainty about mass are three
   different quantities.

That separation is the organizing principle of this investigation. Changing
the law, changing its point summary, and changing the performance criterion
answer different questions and must be compared separately.

## What is established, and under which assumptions

**The readout theorem starts with a law; it does not select one.** The
[proof below](#2-a-stronger-uniqueness-theorem-with-substantial-premises) forces
the ratio of means on its broad domain under the stated composition and
magnitude-law dependence premises. Newton II and Gaussian errors alone leave
the reference weighting open. [Gaussian conditioning](gaussian_conditioning.md)
exhibits that freedom; [flat conditioning](why_flat_conditioning.md) gives
additional assumptions that select the original law; the
[operational analysis](operational_warrant.md) explains why other error-area
definitions select other laws. Those additional premises still need an
experimental warrant. Agreement between mass points cannot validate intervals.

**Finite output and mass information are different properties.** In the flat
Gaussian model with positive definite known covariance, every finite dataset
has divergent ordinary means of \(M\) and \(1/M\), although the ratio-of-means
point, central intervals and log moments exist. At true zero excitation, the data distribution is
independent of mass: no procedure can learn mass from those observations.
This is distinct from the special input of exactly zero *measured* vectors.
With resolved force and weak acceleration, useful lower-bound information can
remain even when relative point accuracy is poor; the
[proved bounds](weak_acceleration_information.md) have explicit known-force
and observation assumptions. [Crossover analysis](crossover_regime.md) proves
smooth behavior and identifies scoped performance gains, without a universal
accuracy guarantee.

**An uncertainty output needs an explicit promise.** A central 95% conditional
interval contains 95% of the declared mass law. A 95% confidence guarantee
concerns coverage of a fixed true mass in at least 95% of repeated experiments
under the procedure's assumptions.
We have [exact confidence-set constructions](conditional_uncertainty_and_design.md)
for known independent isotropic noise, with a common force/acceleration SD
ratio across the stacked trials. Sets may be unbounded, disconnected or empty.
They are not the conditional intervals scored by the main comparison
engine. Plugging estimated variances into those constructions does not inherit
their exact calibration. The [joint-trial analysis](joint_inference_and_uncertainty.md)
also demonstrates undercoverage and null false precision for a particular
joint Jeffreys model with a new unknown excitation vector in every trial.
That result does not describe stable repeated readings or all joint models.

**The unit of repetition determines how evidence combines.** With known
Gaussian covariance, independent readings of one stable vector pair reduce to
precision-weighted means. Separately fitted moment pooling can be inconsistent:
see [weighting and calibration](weighting_and_calibration.md). Trials with
different excitations but one common mass require a declared joint model.
[Exact curve combination](exact_mass_curve_combination.md) explains which
likelihood/prior factors and dependence must be retained; its two-function
merger is exact for a specified aggregation extension. Neither averaging mass
points nor multiplying normalized mass curves is a generic update rule.

**Calibration and application mechanics are part of the experiment.** The
implemented [variance-integration model](comparison_guide.md#pooled-isotropic-calibration-uncertainty)
propagates external pooled isotropic calibration uncertainty under explicit
variance priors. It uses stable measurement means and external scatter, omits
within-experiment scatter, and redraws calibration for every simulated dataset.
Shared calibration across deployments remains a different experiment. The
[robot payload connection](robot_payload_literature_bridge.md) identifies
mechanical and shared-reference requirements; it is not completed robotic
performance validation.

## What the completed simulations add

These findings quantify behavior within specified experiments. They do not
rediscover or strengthen the proofs merely by using more trials. Here \(n\)
is the number of stable readings in one experiment, \(R\) is the number of
simulated experiments per condition, and acceleration SNR uses the true
**per-reading** noise SD. All displayed ± values are one Monte Carlo standard
error (MCSE); `pp` means percentage points. They are pointwise, not simultaneous
uncertainty bounds across the searched conditions.

These studies use independent isotropic Gaussian errors in 3D, true channel
SDs of 1, and an unknown common direction. The main study supplies the true
SDs; the calibration studies vary that information as stated. In the table,
"joint" is the original ratio of magnitude means, and "geometric" is
\(\exp(E_P[\log M])\) under the same flat law.

| Question | Evidence from the completed studies | Consequence |
|---|---|---|
| Does one point summary win? | [Main study, n1][main-n1], \(m=.25\), SNR .5, \(R=10{,}000\): joint-minus-geometric mean loss is **+.012843 ± .000497** when squared log error is capped at \(\log^2 2\), and **−.092715 ± .003129** when capped at \(\log^2 20\). | The ranking reverses even between bounded criteria applied to the same points and trials. Lower loss is better; the acceptable cost of large errors matters. |
| Do stable repeats help? | At \(m=.25\), per-reading SNR 2, known SDs, \(R=10{,}000\) each, ratio-of-means factor-two failure falls from **66.30% ± .473 pp** for [n1][main-n1] to **28.19% ± .450 pp** for [n4][main-n4]. | Repeats improve this condition substantially. These two experiments use independent seeds; this comparison across reading counts is not paired. |
| Does broad null coverage imply identification? | [Null reuse diagnostic][null-sweep], \(n=1,R=10{,}000\): flat-law nominal 95% coverage is **99.71% ± .054 pp** at \(m=1\), but **6.93% ± .254 pp** at \(m=1/64\) and **7.26% ± .259 pp** at \(m=64\). | Coverage near the noise-ratio scale conceals failure elsewhere. These truth values reuse the same null trials; they are not independent confirmations. |
| Does accurate-looking point output validate uncertainty? | \(n=4,m=1\), SNR 8, \(R=1{,}024\): 95% coverage is **94.14% ± .734 pp** with [known SDs][known-cal] and **84.77% ± 1.124 pp** when [both supplied SDs are 25% too small][under-cal]. Both observe zero factor-two failures. | Calibration information can materially affect intervals while a coarse point tolerance reveals no errors. Zero observed failures does not imply zero risk. |
| Does integrating calibration uncertainty help? | [Verified pooled-calibration pilot][integrated-k4], \(n=4,m=1\), SNR 8, \(R=512\), four calibration vectors/channel: plug-in→integrated 95% coverage is **92.77%→94.53%**, paired gain **1.758 ± .581 pp**; mean log-width increases **13.2%**. | Higher coverage here accompanies wider intervals. Point loss barely changes here, and improves or worsens in other cells. No uniform coverage or risk guarantee follows. |

### What the full error curves show

The [archive sensitivity analysis][loss-sensitivity] examines six point rules
over all 24 main-study conditions without fitting again. It compares eight
tolerances and loss caps from 1.1 to 20 and checks survival-curve extrema at
every empirical absolute-log-error breakpoint. There is a substantial
advantage over the simple norm ratio in a particular partially resolved
condition, and a disadvantage in another condition with balanced signals:

| Condition (known SDs; \(R=10{,}000\)) | Allowed factor \(K\) | Our failure rate | Norm-ratio failure rate | Paired difference ± MCSE (pp) |
|---|---|---|---|---|
| \(n=4,m=.25\), acceleration SNR 2 | 1.5 | 53.40% | 65.28% | **−11.88 ± .49** |
| \(n=4,m=1\), acceleration SNR 2 | 1.5 | 28.12% | 22.52% | **+5.60 ± .23** |
| \(n=1,m=.25\), acceleration SNR .5 | 2 | 90.08% | 84.73% | **+5.35 ± .34** |
| Same weak-signal condition | 5 | 29.05% | 34.03% | **−4.98 ± .29** |

Success means \(m/K\le\widehat m\le Km\); negative differences favor our
readout. All four signs survive a conservative **95% simultaneous Monte Carlo
bound of ±4.352 pp**, covering every tolerance for every archived main-study
method and condition. This is stronger protection against selecting a
favorable threshold than the pointwise MCSE alone. The
[testing manual](comparison_guide.md#comparing-whole-error-curves) gives the
bound, its assumptions and reproduction command; the
[curve figure][loss-figure] shows the three distinct conditions.

In the first condition, the benefit is also present at factors 1.25 and 2,
and under normalized squared-log caps 1.25, 1.5, 2 and 3, beyond that same
simultaneous bound. This advantage is not an isolated favorable score.
Nevertheless, more than half of our estimates still miss a factor-1.5 target
there. The median from the same flat law has 56.06% failures: much of the
gain over the norm ratio is shared by another summary of that law. The
additional 2.66 pp benefit of our readout over that median does not clear the
conservative simultaneous bound.

The weak-signal rows exhibit an actual crossing of error curves: fewer very
large mistakes coexist with more moderate mistakes. No dominance claim follows
from a curve that lacks a resolved crossing. These are comparisons within
specified Gaussian truth conditions, with no pooled winner. In particular,
changing mass at fixed acceleration SNR also changes force SNR. The first
two rows differ in both physical mass and the balance of signal strengths;
they do not isolate a causal effect of mass alone.

### Fresh operating-range evidence

The [fixed follow-up protocol](studies/operating_range_followup.json) crosses
effective force and acceleration SNRs `[1,2,3,4,6,8,12,16]` in 64 cells, with
10,000 fresh simulated experiments per cell and no change to the estimator.
All numerical refinement checks passed. The
[accuracy map][operating-accuracy] and [figure][operating-figure] retain every
cell, method and failure. Its single 95% simultaneous Monte Carlo event covers
704 marginal error distributions: success bounds have halfwidth 2.263 pp;
method-difference bounds have halfwidth 4.527 pp.

**The earlier relative contrasts reproduce.** At effective SNRs \((F,a)=(1,4)\),
factor-1.5 success is 46.49% for our readout versus 35.14% for the norm ratio;
at \((4,4)\), it is 71.91% versus 77.08%. The paired failure differences are
−11.35 ± .50 pp and +5.17 ± .23 pp respectively (one MCSE). Both signs survive
the simultaneous bound. Fresh draws test conditions matching the earlier
effective means.

**Relative gains and the stated accuracy targets occupy different cells.**
For each of factors 1.25, 1.5 and 2, every cell with a simultaneously resolved
advantage over the norm ratio is also demonstrably below 90% success at that
same tolerance. The largest observed success among those advantages is 80.24%
at factor 2 and \((1,16)\); even its simultaneous upper bound is only 82.50%.
Thus no sampled cell establishes both that accuracy benchmark and an advantage
over the norm ratio. This does not rule out smaller improvements hidden by
the conservative bound, other loss functions, or other conditions.

Examples of the absolute-accuracy assessment are:

| Effective SNRs \((F,a)\) | Allowed factor | Observed success | Simultaneous lower bound | At least 95% success? |
|---|---|---|---|---|
| (16,16) | 1.25 | 98.76% | 96.50% | Supported |
| (6,6) | 2 | 99.14% | 96.88% | Supported |
| (12,16) | 1.25 | 96.46% | 94.20% | Unresolved |

Across the tested grid, 90% success at factor 1.25 is supported exactly at the
four cells with both SNRs in `{12,16}`; 95% is supported only at `(16,16)`.
For factor 1.5, the supported 95% cells have the smaller SNR at least 8 and
the larger at least 12 **among the sampled values**. `(8,8)` remains unresolved
under the simultaneous rule despite 97.03% observed success. These statements
do not interpolate or extrapolate between cells. The 90%/95% levels and
multiplicative factors are study benchmarks, not application requirements.

An independent audit recomputed all 2,112 method/cell/tolerance probabilities
and their MCSEs from the archived trial estimates, as well as the
reported flat-versus-norm comparisons; they matched the reusable analysis.
All 64 trial archives and their recorded hashes were checked.

The [main protocol](studies/mass_excitation_main.json) covers 24 conditions
at \(R=10{,}000\); the [calibration-sensitivity protocol](studies/calibration_sensitivity_pilot.json)
covers five noise-information regimes × 12 conditions at \(R=1{,}024\);
the [integrated-calibration protocol](studies/calibration_integration_pilot.json)
covers two calibration sizes × 12 conditions at \(R=512\). Each has saved
protocols, source snapshots and trial archives. The verified integrated run is
the preferred evidence: its earlier run omitted paired coverage summaries;
the [rerun check](../../../.tools/calibration_uncertainty_20260920/reporting_rerun_identity.json)
confirmed all 2,088 archived arrays reproduced exactly. It is one
experiment reported correctly, not an independent replication.

The linked results are local scratch evidence and may be absent in a fresh
checkout. The protocols and [running instructions](comparison_guide.md#run-checkpoint-and-resume-a-study)
are the route to reproduction; a run also records the execution source and
runtime, since later code may differ.

### Readout regimes, prior factors and updating

[Working contribution 12](prior_factors_and_sequential_updating.md) reanalyzes
the operating-range archive and three fresh studies (R = 10,000 per cell) with
new reusable modules. It keeps the readout and changes no estimator or earlier
result. Instead it makes two choices in the law explicit: a prior factor on
mass, and the per-reading nuisance used when readings are combined.

- **An unbeaten band, with edges.** With a paired Bonferroni bound over all
  eleven archived rules, the ratio of means beats every rule except its tube
  sibling where the weaker channel's effective SNR is 1, or 2 with the other
  channel strong. At factor 1.5, 15 of 16 archive cells reappear on fresh
  draws. On grids with both SNRs ≥ 1 it also has the smallest worst-case
  regret of the eleven rules (4.7–6.8 pp), although the floored norm ratio is
  narrowly ahead in about half the cells. With one channel at SNR 0.5 and the
  other at 2.5 or more, the same law's reciprocal-root point leads the ratio
  of means by up to 23 pp. The flat-law median and geometric points then
  become the minimax rules, and the ratio of means drops to 5th of 11.
- **Switch rules lose somewhere at every threshold.** Choosing |F|/|a| above an
  observed-SNR threshold and a fallback below is clearly worse than the ratio
  of means in 12–45 of 64 cells at every tested threshold. The largest gain
  in any cell is 7.0 pp. Per-reading switching can beat both branches in a
  cell, so the true-cell choice between branches (4.7–6.8 pp) is a reference,
  not a ceiling.
- **The tilt λ is a fitted prior factor on mass,** centred at σ_F/σ_a. Fitted
  by a declared criterion on the archive (λ\* = 0.35), it cuts worst-case
  regret by a third to a half on fresh draws with both SNRs ≥ 1 and worsens
  it below the fitted range. The hindsight-optimal value moves with the test
  set. A tilt with a free centre moves that centre towards the masses when
  both SNRs are at least 1.
- **Weak-signal points shrink towards σ_F/σ_a.** The median m̂/m is 3.6 for
  m = s/4 at acceleration SNR 1.
- **Declared log-normal priors are priced.** An overconfident prior is
  catastrophic. An honest wide prior that tempers the flat law helps slightly.
  An exact wide log-normal prior is worse than the flat default.
- **As an update rule,** contribution 11's eq. (18) settles below the true
  mass and eq. (20) between the truth and σ_F/σ_a. A symmetric rule that
  counts the prior once converges to the truth, although its 80% intervals
  undercover. Exact angle identities relate the three rules.

## How the tester should support the next question

The generic assessment discipline is to declare **aims, data-generating
mechanisms, estimands, methods and performance measures** before running a
study, then report simulation precision and failures. This is the ADEMP
framework in [Morris, White and Crowther (2019)](https://doi.org/10.1002/sim.8086),
the supplied paper. [Gneiting](https://arxiv.org/abs/0912.0902) explains why a
point rule must be judged with its error criterion in view;
[Gneiting and Raftery](https://sites.stat.washington.edu/people/raftery/Research/PDF/Gneiting2007jasa.pdf)
provide proper distribution scores. Prior-predictive simulation-based
calibration, as in [Talts et al.](https://arxiv.org/abs/1804.06788), checks a
different promise from fixed-truth coverage and requires a specified proper
generative prior; it cannot be silently applied to the improper flat reference.

The existing infrastructure gives that discipline a reusable implementation:

| Responsibility | Existing home | Rule for future work |
|---|---|---|
| Declare the scientific comparison | [Study protocols](studies/) and [configuration](experiment_config.py) | State truth, supplied information, repetition unit, loss, failure handling and fixed replication count. Put a new parameter sweep in configuration. |
| Execute and compare | [Study runner](study_runner.py) and [comparison engine](compare_estimators.py) | Freeze inputs/source, retain failures, compare methods on matched trials, and report paired MCSE separately from quadrature checks. Label true-noise oracles as having extra information. |
| Preserve reusable evidence | [Trial archive](trial_archive.py) | Keep observations, supplied noise, estimates, selected interval endpoints and hashes. Full probability curves are not archived; some new distributional questions need new inference. |
| Analyze and explain | [Analysis/manual](comparison_guide.md#rescore-saved-points-without-running-inference), results browser, this README | Reuse saved trials when they answer the question. Update the relevant claim here; retain detailed assumptions and commands in the guide. Generated outputs stay in `.tools/`. |

**The point assessment now has an explicit, limited operating map.** It shows
where the readout reaches declared accuracy benchmarks and where it improves
on a simpler comparator; the demonstrated gains do not coincide with those
benchmarks in this study. The next distinct assessment is uncertainty at a
common stated guarantee: compare the existing conditional intervals and
known-covariance confidence constructions on the saved observations, retaining
unbounded sets and both miss directions. That addresses the estimator's
remaining uncertainty promise without changing the point-comparison verdict.

The larger open decision is the intended measurement protocol and reporting
promise: a point within an application tolerance, a useful one-sided mass
bound, a confidence set with repeated-experiment coverage, or a prediction.
The existing factor-two criterion is a study choice, not a user-supplied
application requirement. Truth-indexed SNR maps cannot become an observed-data
quality rule without further validation. Any rule suggested by exploratory
archive analysis should be fixed before evaluation on fresh simulations.
Shared calibration, changing excitation and other noise models should follow
when that intended experiment calls for them.

## Scope and source authority

This synthesis covers the active local investigation, its mathematical notes,
the tester and the completed studies identified above. Managed context was
checked against `catalog.yaml`, the discussion README, the working
specification and the later law-assumed-correct steering at
[VD-docs 8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
resolved from `main` on 21 September. This is a scoped synthesis, not an
exhaustive managed/inbox audit. The managed discussion README distinguishes
discussion members from historical unsolicited reports; those reports are
not elevated into governing requirements here. Local investigation work is
active material and is not claimed to have been ingested into VD-docs.
No managed source or formal document relationship is changed.

[main-n1]: ../../../.tools/mass_estimator_studies/assessment-main-20260920/experiments/stable_n1/attempt-001/results.json
[main-n4]: ../../../.tools/mass_estimator_studies/assessment-main-20260920/experiments/stable_n4/attempt-001/results.json
[null-sweep]: ../../../.tools/assessment_main_20260920/main-null/null_truth_sweep.json
[known-cal]: ../../../.tools/mass_estimator_studies/assessment-calibration-20260920/experiments/known/attempt-001/results.json
[under-cal]: ../../../.tools/mass_estimator_studies/assessment-calibration-20260920/experiments/both_under/attempt-001/results.json
[integrated-k4]: ../../../.tools/mass_estimator_studies/assessment-calibration-integrated-20260920-verified/experiments/calibration_k4/attempt-001/results.json
[loss-sensitivity]: ../../../.tools/loss_sensitivity_20260921/main-reviewed/analysis.json
[loss-figure]: ../../../.tools/loss_sensitivity_20260921/main-reviewed/tolerance-comparison.png
[operating-accuracy]: ../../../.tools/operating_range_20260921/analysis/accuracy.json
[operating-figure]: ../../../.tools/operating_range_20260921/figures/operating-range.png

## First result

**Uniqueness of the mass point does not determine its uncertainty.** The
readout proof is sound on its stated domain and can be strengthened. Its
unresolved commitments concern the allowed domain and operations, further
direction weighting, and the probability law supplied to the readout.

A sensitivity calculation below keeps the measurements, Gaussian likelihood,
and mass point unchanged while changing the conditional 95% mass interval
enormously. Agreement on point values cannot establish that the remaining
freedom is harmless.

The following sections retain the original 15 September derivation. Read
their foundational questions together with the current synthesis above;
later notes and studies address parts of that original agenda.

## 1. Separate the mathematical objects

For one object, interval, and inertial-frame-qualified account, write

$$
\mathbf F^*=f\mathbf u,\qquad
\mathbf a^*=\alpha\mathbf u,\qquad f,\alpha>0,\quad\|\mathbf u\|=1.
$$

The empirical construction supplies

$$
dP=\frac1Z e^{-Q(f,\alpha,\mathbf u)/2}
       \rho(f,\alpha,\mathbf u)\,df\,d\alpha\,d\Omega.
$$

Here \(Q\) comes from the measured vectors and their declared Gaussian
covariance; \(\rho\) is the latent reference weighting. The existing numerical
laboratory uses \(\rho=1\). There are then two outputs:

$$
\widehat m=\frac{E_P[f]}{E_P[\alpha]},
\qquad
M=\frac f\alpha\ \text{with its induced distribution under }P.
$$

The first is a scalar summary; the second is the current conditional
uncertainty law. Repeated-experiment variation of the estimator is a third
object, requiring a sampling experiment. Integration error is a fourth.

The composition theorem takes \(P\) as given. It does not construct \(P\)
from measurements. The September 7 steering also distinguishes this averaging
operation from selecting a single minimum-\(Q\) pair.

## 2. A stronger uniqueness theorem, with substantial premises

Choose reference units to make \(X=f/F_0\) and \(Y=\alpha/a_0\) dimensionless.
Let \(R(X,Y)\) be positive and finite on **every** joint law of positive
integrable variables. The composition requirements are

$$
R(X_1+X_2,Y)=R(X_1,Y)+R(X_2,Y), \tag{C1}
$$

$$
\frac1{R(X,Y_1+Y_2)}
=\frac1{R(X,Y_1)}+\frac1{R(X,Y_2)}. \tag{C2}
$$

C3 says the output depends only on the magnitude-pair probability law, with
no further direction information or probability-space labels.

**Checked result:** these requirements, the universal domain, and the single
calibration \(R(1,1)=1\) force \(R(X,Y)=E[X]/E[Y]\).
The source's full deterministic calibration C4 and constant-limit continuity
C5 follow from the other assumptions. Appendix A proves the strengthening.
The original sufficient theorem remains valid.

### What this leaves open

- C1/C2 concern **one propagated joint law**, with the indicated common channel
  retained. They do not equate that law with a fresh fit to another experiment.
- Scalar force addition presupposes the relevant common direction. Vector
  addition does not license arbitrary addition of force magnitudes.
- Statewise Newton II alone does not require a chosen scalar summary of
  uncertainty to obey both composition rules. Their operational warrant
  remains substantive.
- The full domain includes arbitrary couplings, deterministic channels, and
  integrable laws far beyond the posteriors of one Gaussian model. Its
  application requires that extension to be justified, or a theorem on the
  actual restricted domain.
- C3 and the construction of \(P\) remain separate questions.

### The direction freedom noticed in the earlier work

Extend the source's counterexample by putting

$$
\mathbf v=E_P[\mathbf u],\quad b=(\mathbf u\cdot\mathbf v)^2,\quad
h_t=1+t b,\quad
R_t=\frac{E_P[h_t f]}{E_P[h_t\alpha]},\qquad t\ge0.
$$

For finite \(t\), this retains positivity, finiteness, units covariance,
rotation invariance, channel reciprocity, exactness at a fixed mass, and
C1/C2 when the common direction law is retained. It violates C3.

Define \(dQ_\alpha=\alpha\,dP/E_P[\alpha]\). Direct subtraction gives

$$
\boxed{
R_t-R_0=
\frac{t\,\operatorname{Cov}_{Q_\alpha}(b,M)}
     {1+tE_{Q_\alpha}[b]}.
}
$$

The freedom affects the point precisely when this direction feature and mass
co-vary under the acceleration-weighted law. Fixed-mass cases, and cases
with zero such covariance, hide it.

For the source's states \((\mathbf u,f,\alpha)=(\mathbf e_1,1,1)\) with
probability \(3/4\) and \((\mathbf e_2,2,1)\) with probability \(1/4\),

$$
R_t=\frac{80+29t}{64+28t};
\qquad R_0=\frac54,\quad R_1=\frac{109}{92}.
$$

This supplies sensitivity, not a reason to prefer any \(t\). C3 permits
direction to affect the construction of the magnitude law; it excludes a
second direction-dependent weighting at the readout stage.

## 3. The mass point cannot carry uncertainty by itself

Compare equiprobable-state laws in consistent reference units:

| Joint law | Force marginal | Acceleration marginal | Point | Induced mass |
|---|---|---|---:|---|
| \((1,1)\) or \((9,9)\) | \(1\) or \(9\) | \(1\) or \(9\) | \(1\) | \(1\) with certainty |
| \((1,9)\) or \((9,1)\) | \(1\) or \(9\) | \(1\) or \(9\) | \(1\) | \(1/9\) or \(9\), equally likely |

Both complete magnitude marginals are identical. Their pairing changes the
mass uncertainty from zero to a log standard deviation of \(\log9\).

These are exact examples in the theorem's universal domain, not asserted to
be exact finite-noise Gaussian posteriors. The point and marginals cannot
reconstruct mass uncertainty. This does not establish that using coupling
would improve point accuracy under a specified loss.

For the later information-use question, the comparison must specify the
complete output and intended use, rather than just compare scalar equations.

## 4. Same measurements and point, adjustable uncertainty

Take both **measured** vectors to be zero, with finite independent isotropic
Gaussian standard uncertainties \(\sigma_F,\sigma_a>0\). Set
\(s=\sigma_F/\sigma_a\), and use standardized polar coordinates

$$
f=\sigma_F r\sin\theta,\qquad
\alpha=\sigma_a r\cos\theta,\quad
r>0,\quad0<\theta<\pi/2.
$$

The likelihood is \(e^{-r^2/2}\), and the magnitude Jacobian is proportional
to \(r\). Consider the fixed-experiment sensitivity family

$$
\rho_\lambda=\exp(\lambda\sin2\theta),\qquad\lambda\in\mathbb R.
$$

Every finite \(\lambda\) gives bounded, strictly positive weighting, invariant
under rotations, covariant under unit changes, and symmetric in the
standardized channels. The existing flat measure is \(\lambda=0\).

The radial law remains Rayleigh, independently of

$$
p_\theta(\theta)=\frac{e^{\lambda\sin2\theta}}{Z_\lambda},\qquad
Z_\lambda=\int_0^{\pi/2}e^{\lambda\sin2\theta}\,d\theta.
$$

Symmetry about \(\pi/4\) makes \(E[\sin\theta]=E[\cos\theta]\), so

$$
\boxed{\widehat m=s\quad\text{for every }\lambda.}
$$

Yet \(L=\log(M/s)\) has density

$$
\boxed{
p_L(\ell)=\frac{e^{\lambda\operatorname{sech}\ell}}
                 {2Z_\lambda\cosh\ell}.
}
$$

Its median is zero; its width depends on \(\lambda\).

| \(\lambda\) | Point / \(s\) | Conditional central 95% interval / \(s\) | SD of \(\log(M/s)\) |
|---:|---:|---|---:|
| -4 | 1 | [0.006897, 144.989662] | 2.852883 |
| **0: existing flat measure** | **1** | **[0.039290, 25.451700]** | **1.570796** |
| 4 | 1 | [0.271496, 3.683292] | 0.666264 |
| 16 | 1 | [0.593688, 1.684385] | 0.263564 |
| 64 | 1 | [0.779872, 1.282261] | 0.126513 |

Multiply the interval by \(s\) to restore mass units. These are conditional
probabilities, not established repeated-experiment coverage statements.

The intervals shrink monotonically: the ratio of densities at two values of
\(\lambda\) is proportional to
\(e^{(\lambda_2-\lambda_1)\operatorname{sech}\ell}\), decreasing in \(|\ell|\)
when \(\lambda_2>\lambda_1\). Expansion near zero gives
\(\operatorname{Var}(L)\sim1/\lambda\) as \(\lambda\to+\infty\).

**No new measurement was supplied.** The weighting increasingly favors
similar standardized magnitudes. The checked symmetries and unchanged point
therefore leave a substantial uncertainty choice open.

### Qualifications that matter

The family uses the instrument-scale ratio. Changing that ratio changes its
physical reference measure; cross-experiment coherence may impose further
restrictions. Simultaneously scaling both uncertainties leaves the entire null
mass law unchanged for each fixed \(\lambda\), so that existing null-scaling
test does not select its width.

The weighting is smooth inside the quadrant and at the axes away from the
origin. Its limit at the origin is direction-dependent for nonzero \(\lambda\).
That point has zero measure, and continuity there is not a current requirement.
Bounds are finite for each \(\lambda\), but not uniform over the whole family.

Measured zeros still constrain excitation magnitudes. Their constant
fixed-mass profile does not by itself establish a general non-identifiability
theorem. At a **true zero-excitation state**, the Gaussian sampling law is
identical for every mass; that separate claim includes the zero pair in the
physical model's boundary.

The instrument-independent family \(\rho_k=(f\alpha)^k\), \(k>-1\), exposes
another sensitivity direction. Reference-unit constants cancel. At measured
zeros, the standardized squared magnitudes divided by two are independent
Gamma variables with shape \((k+1)/2\). Thus the point again equals \(s\), but

$$
\operatorname{Var}\log(M/s)=\tfrac12\psi_1((k+1)/2),
$$

where \(\psi_1\) is trigamma. The variance tends to infinity as \(k\to-1\)
and zero as \(k\to\infty\). This family also changes boundary behavior and,
for large \(k\), favors increasingly large latent magnitudes despite zero
readings. Neither family is endorsed as an estimator policy.

## 5. Infinite moments and practical spread

For the flat Gaussian construction, integrate out direction and call the
unnormalized magnitude density \(H(f,\alpha)\). With finite data and positive
definite covariance it is positive at the axes, with Gaussian decay. Then

$$
p_M(m)=\frac1Z\int_0^\infty\alpha H(m\alpha,\alpha)\,d\alpha
      =\frac1{Zm^2}\int_0^\infty fH(f,f/m)\,df.
$$

Dominated convergence gives \(p_M(m)\to C_0>0\) as \(m\downarrow0\) and
\(p_M(m)\sim C_\infty/m^2\) as \(m\to\infty\). The constants are positive
axis integrals; a Gaussian envelope dominates them uniformly.

Therefore **mass and inverse-mass means and the ordinary second mass moment
diverge for every finite dataset in the default model**, including strongly
aligned ones. At informative data the tail coefficient can be tiny and
central intervals narrow. All fixed-order absolute log-mass moments remain
finite. Divergent ordinary moments do not determine central interval width;
variance about a divergent mean is formally undefined.

The bounded angular family retains these tail exponents for finite \(\lambda\).
At measured zero-zero its density is explicitly

$$
p_M(m)=\frac{s}{Z_\lambda(s^2+m^2)}
\exp\!\left(\frac{2\lambda sm}{s^2+m^2}\right).
$$

Thus arbitrarily narrow central uncertainty can coexist with infinite
ordinary moments. The point-mass limit as \(\lambda\to\infty\) is nonuniform
in the tails; moments do not pass through that limit.

### A mean identity does not select an uncertainty law

Although \(E_P[M]\) diverges, \(dQ_\alpha=\alpha\,dP/E_P[\alpha]\) satisfies
\(E_{Q_\alpha}[M]=\widehat m\). This is a different law. At default measured
zeros, for \(R=M/s\),

$$
p_R(r)=\frac2{\pi(1+r^2)},\qquad
q_{\alpha,R}(r)=\frac1{(1+r^2)^{3/2}}.
$$

The latter has mean 1, median \(1/\sqrt3\), and infinite second moment.
Its quantile is \(p/\sqrt{1-p^2}\), whereas the former's is
\(\tan(\pi p/2)\). The mean identity supplies no argument to replace
\(P\)-uncertainty by \(Q_\alpha\)-uncertainty.

When reusing the original acceleration, retain its dependence with mass:
\(M\alpha=f\) state by state. Drawing independently from their marginals is
a different calculation. Independence for a new acceleration requires a
stated experiment and reuse assumption.

## 6. Remaining foundational questions

**What requirements should the entire uncertainty law satisfy?** Point
calibration, reciprocity, and null scaling leave the freedoms above open.

1. What operational argument requires C3 and excludes a second direction
   weighting? The covariance identity identifies where it matters.
2. What experiment or coherence requirement selects or restricts the latent
   measure? Test proposed requirements against both null families, including
   changes of instruments and repeated observations of one latent pair.
3. Which uncertainty outputs support the intended reuse of mass? Establish
   their behaviour and appropriate repeated-experiment calibration before
   treating narrower intervals as evidence of informational advantage.

Subsequent work has developed common-mass models, non-null sensitivity,
confidence constructions, performance comparisons and literature connections;
their scope is summarized above. Those developments do not by themselves
select the reference law or the intended uncertainty promise. The immediate
empirical priority is stated in [the testing agenda](#how-the-tester-should-support-the-next-question).

## Appendix A. Continuity follows on the full domain

Put \(r=R(X,Y)\), \(x=R(X,1)\), \(y=1/R(1,Y)\). Expanding
\(R(X+1,Y+1)\) in the two composition orders gives

$$
\frac{xr}{x+r}+\frac1{y+1},\qquad
\frac{(r+1/y)(x+1)}{r+1/y+x+1}.
$$

Their difference is

$$
-\frac{(ry-x)^2}{(x+r)(y+1)(yr+1+xy+y)}.
$$

The denominator is positive, so \(ry=x\):
\(R(X,Y)=T(X)/S(Y)\), with \(T(X)=R(X,1)\), \(S(Y)=1/R(1,Y)\).
Both are additive, positive, finite, law-dependent and normalized at 1.

For \(T\), positivity makes its additive action on constants increasing;
rational squeezing gives \(T(c)=c\). If \(X\le Y\), then
\(T(Y)+\epsilon=T(X)+T(Y-X+\epsilon)>T(X)\) for every \(\epsilon>0\), hence
\(T(Y)\ge T(X)\). Extend to nonnegative integrable variables by
\(T_0(Z)=T(Z+1)-1\). Additivity makes the positive shift immaterial;
\(T_0\) is nonnegative, additive and monotone.

Suppose \(E[Z_n]\to0\) but \(T_0(Z_n)\not\to0\). Choose a subsequence with
\(E[Z_n]\le2^{-n}\) and \(T_0(Z_n)\ge\epsilon>0\).
C3 permits realizing these laws on one probability space. Their sum \(Z\)
has finite mean and is finite almost surely (redefine its null infinite set).
The universal domain includes \(Z+1\). Finite additivity and monotonicity imply
\(T_0(Z)\ge\sum_{n=1}^N T_0(Z_n)\ge N\epsilon\), contradicting finiteness.
So \(E[Z_n]\to0\) forces \(T_0(Z_n)\to0\).

Consequently \(|T(X_n)-c|\le T_0(|X_n-c|)\to0\) when \(X_n\to c\) in \(L^1\).
This is C5. Countable additivity of \(T\) was not assumed; monotone convergence
was used for the underlying ordinary probability measure.

For iid copies \(X_i\), additivity, rational homogeneity and C3 give
\(T(n^{-1}\sum_iX_i)=T(X)\). The sample mean converges in \(L^1\) to \(E[X]\):
truncate at \(K\), bound its bounded-sample deviation by
\(\sqrt{\operatorname{Var}(\min(X,K))/n}\), and the two tail terms by
\(2E[(X-K)_+]\); let \(n\to\infty\), then \(K\to\infty\).
Derived continuity gives \(T(X)=E[X]\), and likewise \(S(Y)=E[Y]\).
Full deterministic calibration follows.

The countable-sum argument uses the universal finite-mean domain. A restricted
posterior family cannot silently borrow that construction.

## Original derivation: reproduction and sources

[checks.py](checks.py) reproduces the main table and exact discrete examples.
It uses Python and NumPy; it imports the existing estimator solely for the
baseline check and writes no files unless an output path is supplied:

~~~powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/checks.py'
~~~

All assertions passed on 15 September. Repeating the main quadrature with
192 and 384 nodes changed reported quantities by at most \(1.21\times10^{-13}\)
relatively. An independent Cartesian integral at 512 and 1024 nodes agreed
within \(10^{-7}\) in normalization. Checks also cover half-Cauchy recovery,
reciprocal quantiles, transformed log density, tail coefficients, asymptotics,
and rational direction/coupling examples. These bounds describe the checked
cases. Moment divergence and the uniqueness strengthening rest on proofs.

Proof and uncertainty derivations received independent parallel reviews.
Development scripts and machine outputs remain in the workspace scratch
folder .tools/mass_equation_20260915/. This README maintains the current
argument; the focused notes retain the detailed derivations and their checks.

All managed reads used
[VD-docs commit 8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
resolved once from main for the original derivation pass. Selection began with
catalog.yaml and the discussion README, following uniqueness, composition, measure, uncertainty,
null, and steering references:

- [Uniqueness proof](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md): proof target and existing C3 counterexample.
- [Readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md): composition motivation and measure freedom.
- [Working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md): construction and uncertainty objects.
- [Additional derivations](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md): null scaling, log spread, loss and weighted-law identities.
- [Chain audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md): premises versus consequences, read in the proof review.
- [User steering](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md): assume Newton II; distinguish averaging from least-error selection.

Local scope: laboratory README, estimator, practical formulas, uncertainty
walkthrough, and relevant tests. The older rejection contract does not govern
this investigation. Separated exploratory reports and literature files were
not used as discussion authority. This was a targeted investigation, with
no new literature survey, priority claim, exhaustive inbox audit, managed edit,
or entry-policy revision.

# Understanding the mass estimator: null uncertainty, reuse, and changing the experiment

**Subsequent user clarification:** this investigation conditions on Newton II
being correct. The gate/withholding interpretations below describe the earlier
specification and are not requirements for that conditional purpose. The
numerical results and dependence findings stand. See [the steering record](mass_estimator_law_assumed_correct_steering_2026-09-07.md),
including the distinction between averaging reconciliations and selecting a
minimum-error reconciliation.

**Date:** 2026-09-07. **Status:** local explanatory and test record; awaiting
transfer to VD-docs, not submitted. No new estimator/admission policy is ratified.

This continues the September 6 behavioral suite in response to the user's
request to work through the meaning of the results. The point equation remains
$\widehat m=E[f]/E[\alpha]$. The purpose here is to explain what that summary
retains, what it loses, and which uncertainty belongs in a downstream operation.

New executable work:

- [Experiments](../../03_trace_and_evaluator/mass_estimator/uncertainty_walkthrough.py)
- [Interpretation and propagation tests](../../03_trace_and_evaluator/mass_estimator/test_uncertainty_walkthrough.py)
- [Saved results](../../03_trace_and_evaluator/mass_estimator/results/uncertainty_walkthrough.json)
- [Suite verification](../../03_trace_and_evaluator/mass_estimator/results/tests.json)

## 1. Begin with an experiment, not the expression 0/0

Suppose both vector instruments read zero. Each has independent isotropic
Gaussian error, with per-component standard uncertainty 1 N and 1 m/s². The
latent measure is flat in the two positive magnitudes and uniform in their
common direction. These assumptions define the example; a reported zero alone
does not define its uncertainty law.

Under this policy the possible positive magnitudes are independent half-normal
variables. Their means are both numerically $\sqrt{2/\pi}\approx0.797885$ in
their respective units, so the point is 1 kg. The possible mass in each state
is $M=f/\alpha$; its law is half-Cauchy with scale 1 kg.

There are two different senses of uncertainty already:

1. The actual pair is unresolved. The fitted law assigns probabilities to
   neighboring positive force and acceleration values and thus to their ratios.
2. The Newton-II account at the zero central readings selects no positive
   mass by best fit. For every fixed $m>0$, the profile discrepancy is zero.

Neither statement erases the other. The first is conditional on how the model
weights nearby pairs; the second says every mass can fit arbitrarily well.

## 2. Why the null law is not “all masses equally probable”

In the standardized $(f/\sigma_F,\alpha/\sigma_a)$ plane, the latent positive
quadrant at null has a radially symmetric Gaussian density. Its polar angle is
uniform, and $M=\tan\theta$. Uniformity in this angle is not uniformity in
mass. The induced law has

$$P(M\le m)=\frac{2}{\pi}\arctan m\quad(m>0),$$

with dimensionless numerical $m$ in units of the 1 kg apparatus scale.

For example, the model assigns about 40.97% probability to 0.5–2 kg, but about
0.0955% to 500–2000 kg. Yet 1000 kg fits the zero readings exactly as well as
1 kg by the profile test. The probability difference is supplied by integration
over the declared neighboring-pair measure; it is not a difference in their
best Newton-II fit.

The point equals the median at this symmetric null, but it is not the ordinary
mean of $M$, which does not exist. Nor should it be called “the most likely
mass”: the half-Cauchy density per kilogram has its maximum at the zero
boundary, while its density per logarithmic interval peaks at the scale.
Density modes depend on the coordinate/measure used to describe them.

This explains the original surprise more precisely. A positive mass point and
non-identification can coexist. Writing “1 kg ± infinity” does not explain
that coexistence and cannot substitute for either probability or compatibility.

## 3. What a better instrument changes

If both instruments become 100 times more precise while both central readings
stay zero, both possible magnitudes shrink by 100. Their ratios have exactly
the same law. The force itself can become very tightly localized near zero
without mass being determined.

If only force uncertainty is halved, the null point and ratio quantiles halve.
If only acceleration uncertainty is halved, they double. In both cases every
positive mass still has zero profile discrepancy. The scale sensitivity is a
property of the declared null probability law, not a new measurement of the
object's mass.

At a physically exact null state $F^*=a^*=0$, changing the true mass cannot
change the distribution of Gaussian measurement noise. The new suite feeds
identical noise realizations generated for masses 0.001, 1, and 1000 into
the estimator and gets identical results. No repeated procedure can distinguish
these masses from this unexcited data distribution alone.

This statement concerns a fixed true null state. It is distinct from
conditioning the fitted positive-magnitude law on the particular central
readings, and from assuming that all future sample means are exactly zero.

## 4. Why infinite mass-ratio moments need not give infinite reconstructed force

The mass uncertainty was obtained from pairs. If we reuse the original latent
acceleration, the possible mass and that acceleration remain related:

$$M=\frac f\alpha,\qquad M\alpha=f.$$

A large mass often occurs because its accompanying acceleration is very small.
Those two possibilities belong together. Multiplying each mass by its original
acceleration reconstructs the original force in every sampled state. The
original force magnitude has a half-normal law, with finite mean and variance,
despite the mass ratio having no finite mean.

The experiment draws one million pairs from the exact null magnitude law and
checks this reconstruction. Maximum relative numerical error is about
$2.22\times10^{-16}$. Reconstruction is a consistency identity of the fitted
model, not an independent test of Newton II.

Now discard that dependence: draw a mass from the same marginal law, and pair
it with an independent acceleration drawn from the correct acceleration
marginal. That is a **different joint law**. Its force product has a long tail
that the original force did not have.

| Operation at unit null scales | Central 95% force-magnitude interval | Interpretation |
|---|---:|---|
| Original fitted force | 0.031338–2.24140 N | Analytic half-normal interval |
| Each mass with its original acceleration | Exactly the same distribution | Preserves the original pairs |
| Mass independently paired with the same acceleration marginal | Approximately 0.00844–20.3 N | Alters the joint inference |
| Multiply separate 95% endpoint intervals | 0.001231–57.0475 N | A range product, not a joint central 95% interval |

The independent-product interval is checked against a separate one-dimensional
quadrature oracle, as well as Monte Carlo. Its sample probability of force
above 10 N was about 5.05%; the original half-normal probability is only
$1.52\times10^{-23}$. A zero count in the original million draws does not
mean an exactly zero probability.

The log variances make the dependence particularly clear. All logarithms here
are of numerical ratios to the chosen reference units. At this null,

$$\operatorname{Var}(\log f)=\frac{\pi^2}{8},\quad
\operatorname{Var}(\log M)=\frac{\pi^2}{4},\quad
\operatorname{Cov}(\log M,\log\alpha)=-\frac{\pi^2}{8}.$$

Retaining the covariance recovers $\operatorname{Var}(\log(M\alpha))=
\pi^2/8$. Omitting it gives $3\pi^2/8$, three times as large. These are finite
log moments, not an illicit calculation with infinite raw mass moments.

**Operational lesson:** mass uncertainty cannot generally be detached from
the evidence that produced it when that same evidence is reused downstream.
A scalar and a marginal mass interval do not encode the necessary dependence.

## 5. A single point can preserve the force marginal while changing the relation

There is a second, subtler test. At unit null scales, replace $M$ by
$\widehat m=1$ and calculate $\widehat m\alpha$.

The result has exactly the same marginal half-normal law as $f$. It has the
same mean, variance, and quantiles. A test checking only the reconstructed
force distribution would pass.

But under the null fitted law, the positive magnitudes $f$ and $\alpha$ are
independent. After point substitution, the reconstructed force magnitude and
$\alpha$ are perfectly correlated. The original pairings have been replaced.
The simulation finds original magnitude correlation about -0.00034 and point
substitution correlation exactly 1.

This statement concerns **magnitudes**. The full vectors share a common
direction and are not independent vector random variables. The test does not
mistake zero magnitude correlation for independent original vector channels
after fitting the common direction.

For any supplied law with finite positive first moments, the chosen point
guarantees

$$E[\widehat m\alpha]=E[f],\qquad E[f/\widehat m]=E[\alpha].$$

This is a useful exact property of the readout: it preserves these first
magnitude moments in both directions. It does not guarantee each pair, the
joint distribution, or arbitrary downstream nonlinear quantities. If every
supported pair already shares one mass, point substitution does preserve
every pair; the suite includes that control.

The null coincidence of entire marginal laws is stronger than this mean
identity and relies on the equal-shape scale relationship. It must not be
generalized to every fitted law.

## 6. A genuinely new prediction is different from reconstructing an old pair

Assume the same mass persists, and ask what force would be required to produce
a newly specified exact acceleration of 1 m/s². This is a counterfactual
input, supplied without additional mass evidence. It is not the unresolved
acceleration from the old null trial.

The induced force magnitude is then $M\times1$ and retains the broad mass
law: central 95% interval 0.03929–25.45170 N. The point alone gives 1 N and
does not carry that uncertainty. Under the chosen law, the force mean and
variance do not exist; finite quantiles remain meaningful model summaries.

For a newly specified exact acceleration of **zero**, every finite positive
mass instead predicts exactly zero net force. No mass identification is needed
for that particular prediction. Conversely, exact zero net force predicts
exact zero acceleration for every finite positive mass.

| Downstream question after the null determination | What the same inference says |
|---|---|
| Reconstruct original latent force using original latent acceleration | Preserve joint pairs; recover the original finite-spread force law |
| Required force for externally specified exact acceleration 1 | Broad half-Cauchy-scaled force law |
| Required force for externally specified exact acceleration 0 | Exactly zero force under Newton II |
| Acceleration under externally specified exact force 1 | Reciprocal half-Cauchy-scaled acceleration law |

This is why treating mass uncertainty as one infinite number leads to trouble:
the same uncertain mass can support a sharp zero prediction, a broad nonzero
prediction, or a finite-spread reconstruction, depending on what is being asked
and which dependence is retained.

For a new **measured** acceleration, independence is not automatic merely
because the instrument errors are independent. The physical experiment can
couple acceleration to mass. Its joint model and any additional evidence must
be supplied before choosing a propagation rule. These experiments use an
externally specified exact input to make that boundary explicit.

The global compatibility gate still precedes reporting a prediction. In the
angle experiments, propagated intervals from a rejected account are explicitly
marked as conditional diagnostics; the post-gate prediction field is empty.
A new regression test verifies this withholding.

The 95% intervals in this record are induced-law probability intervals, not
claims of calibrated 95% coverage under repeated physical experiments.

## 7. Change excitation continuously while keeping the displayed ratio fixed

The next path keeps the measured vectors aligned and their magnitude ratio
equal to 2, but reduces their strength at fixed unit uncertainty:

| Measured F magnitude | Measured a magnitude | Mass point | Induced 95% mass interval |
|---:|---:|---:|---:|
| 20 | 10 | 2.0000 | 1.626–2.535 |
| 6 | 3 | 1.9966 | 1.044–5.920 |
| 2 | 1 | 1.5204 | 0.1560–25.39 |
| 0.6 | 0.3 | 1.0650 | 0.04440–25.76 |
| 0.02 | 0.01 | 1.000075 | 0.03930–25.45 |

The displayed quotient remains 2 throughout, but it ceases to be a stable
guide to the latent magnitudes as both measurements become unresolved. The
readout moves smoothly toward the apparatus scale 1. The point's movement is
not a physical change in the object's mass; it is a change in what these
uncertain inputs support under the policy.

If the readings **and** their uncertainty scales shrink together, all these
signal-to-noise ratios remain fixed. The mass point stays 2 and its interval
stays 1.62555–2.53459. Alternatively, holding the nonzero readings fixed while
reducing uncertainty concentrates the central law toward the ratio 2.

These paths also give different nonzero future-force intervals. A new exact
acceleration of 3 multiplies the mass endpoints by 3. The strong (20,10)
case gives 4.877–7.604 N; the weak (.02,.01) case gives about 0.1179–76.36 N.
The numerical readout's totality does not erase that practical difference.

## 8. Rotate vectors: smooth calculations, qualitative changes of standing

Under the default isotropic model, after choosing dimension and units, the
input behavior is controlled by three quantities: force signal-to-noise ratio,
acceleration signal-to-noise ratio, and their angle. Let these strengths be
$p,q$, let the angle be $\theta$, and write $t=m/(\sigma_F/\sigma_a)$.
The fixed-mass profile is

$$Q(t)=\frac{p^2+q^2t^2-2pq\cos\theta\,t}{1+t^2}.$$

This provides a compact map for testing the default policy; extra covariance
parameters or a known-direction constraint add further controls.

For equal strengths $p=q=5$, the point is the apparatus scale at every angle.
At the illustrative threshold 9, the free-mass gate changes from pass to reject
at approximately **50.2082 degrees**. Tests on either side (45 and 55 degrees)
leave the mass point at 1 while changing the fit decision. The numerical law
varies smoothly; thresholding is what creates a discrete decision.

For equal strengths $p=q=\sqrt6$, the global fit passes at all angles. Beyond
**120 degrees**, the point 1 itself fails the same fixed-mass threshold. The
compatibility set splits into two ranges. At 121 degrees they are approximately
$(0,0.7830]$ and $[1.2772,\infty)$; at opposition they are
$(0,0.267949]$ and $[3.732051,\infty)$.

This follows the anti-aligned finding through a continuous path rather than
treating opposition as a special numerical switch. The point remains a valid
output of the readout formula. Interpreting that point as one adequately fitting
physical explanation would be an additional claim the fixed-mass check rejects.

## 9. What these tests establish about use

- At null, numerical positivity is compatible with non-identification, but
  the induced mass probabilities are not uniform across masses.
- Dependence can make a downstream product well behaved even when a marginal
  ratio has no finite mean. Uncertainty arithmetic cannot discard that dependence.
- The ratio-of-means point preserves first magnitude moments for same-pair
  reconstruction. That is a specific capability, not a guarantee about each state.
- Even agreement of reconstructed marginal distributions is weaker than
  preservation of the joint law.
- An undetermined mass can still suffice for a zero force/acceleration
  prediction. Usefulness depends on the downstream demand.
- Strong-signal agreement, weak-signal apparatus dependence, and angular
  rejection can be followed as continuous paths through the same construction.

The next semantic contract needing a ruling is what an estimator result must
retain for downstream reuse: a marginal law suffices for certain newly
specified inputs, while original-input reuse can require the joint law or an
equivalent dependency-preserving representation. The current result description
names uncertainty and provenance but does not by itself prescribe that complete
representation. This is recorded for entry/evaluator development, not installed
as a new admission gate here.

These conclusions describe the declared Gaussian/flat-magnitude policy and the
stated prediction setups. They do not settle non-Gaussian policies, hard bounds,
mass changes, common-mass multi-trial posterior measures, or universal calibrated
uncertainty coverage. Existing real-data trials retain their provisional
uncertainty and mechanical-account limitations.

## 10. Reproduction, verification and custody

Run `python uncertainty_walkthrough.py` and `python run_suite.py` from the
laboratory directory. Seventeen new tests supplement the original forty.
The saved verification record captures outcomes, runtime, versions, and hashes
of the core and new behavioral sources. The experiment uses one million exact
null-law draws with seeds 20260907 and 20260908; simulation figures should not
be mistaken for analytic constants. Analytic half-normal/half-Cauchy quantiles,
log-moment identities, a separate independent-product CDF quadrature, and
state-by-state reconstruction provide independent checks.

Numerical non-null paths use the existing 2D integrator at 256 directions and
2048 ratio intervals. They remain subject to its documented refinement limits,
particularly in severe anti-aligned tails. A test's initial arbitrary “tenfold
upper-quantile inflation” expectation was rejected: the detached upper quantile
is roughly ninefold the original. It was replaced by an independently integrated
distributional oracle, rather than treating ten as a physical requirement.

**Managed context:** resolved `main` through the GitHub connector on 2026-09-07
to `canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`, unchanged from
the preceding test pass. The inventory and fetched relevant passages all use
that coherent snapshot:

- `catalog.yaml`: inventory/status for the construction and null discussion.
- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`:
  sections 2–3 and 8–12, especially provenance, the latent law, and the readout.
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md`:
  null scaling and same-law aggregation, sections 2–3.

The managed September 5 research and zero-zero discussion had been read in the
preceding pass; their tail/non-identification conclusions are retained through
that established context, not represented as freshly reread in full today.

**Local scope:** read `AGENTS.md`, the active reference implementation, README,
and runner; used the September 6 behavior report and saved experiments. Scoped
uncertainty/correlation/reuse/prediction/zero searches also inspected excerpts
from the newer local air-track, robot/cart, and arrow-data audit records. These
remain active/unsubmitted material and were not assumed absorbed into VD-docs.
Their datasets and code were not changed. No exhaustive inbox or historical
collection search was performed; this is a scoped behavioral continuation.

The new script/tests/results remain in `03_trace_and_evaluator/mass_estimator/`.
This report is in `pass to VD-docs/` under the workspace staging policy. Managed
context, the estimator equation, entry drafts, and `VDfirst/` were not edited.
No submission, commit, push, or publication is implied by this local save.

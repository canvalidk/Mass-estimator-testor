# Mass estimator: behavioral tests and a walkthrough

**Later steering, 2026-09-07:** the user clarified that the present investigation
assumes Newton II is correct and interprets disagreement as measurement error.
The rejection framing below records the earlier working specification; a point
outside a threshold set is not by itself unwanted under the clarified purpose.
See [the correction and remaining averaging/optimization distinction](mass_estimator_law_assumed_correct_steering_2026-09-07.md).

**Date:** 2026-09-06. **Status:** new local test report and interpretation;
awaiting transfer to VD-docs, not submitted. This report does not ratify a new
measure, admission rule, or correction to managed sources.

**Request:** stop trying to prove the estimator and test its behavior at zeros,
aligned and anti-aligned vectors, and under uncertainty.

**Executable laboratory:** [README and commands](../../03_trace_and_evaluator/mass_estimator/README.md).
**Numerical results:** [cases](../../03_trace_and_evaluator/mass_estimator/results/cases.csv),
[complete experiment data](../../03_trace_and_evaluator/mass_estimator/results/experiments.json),
[automated verification](../../03_trace_and_evaluator/mass_estimator/results/tests.json).

## 1. What exactly is being tested?

The reported point remains

$$\widehat m=\frac{E_P[f]}{E_P[\alpha]},$$

where the latent vectors are $f u$ and $\alpha u$, with positive magnitudes
and a common unit direction. Their law $P$ weights compatible pairs by the
measurement likelihood. The default tested policy is independent isotropic
Gaussian errors, uniform common direction, and flat measure $df\,d\alpha\,d\Omega$.
The integrator also supports declared anisotropic and correlated Gaussian
covariances; these are tested separately.

The mass is a ratio of two means. It is neither the ratio of measured vector
magnitudes nor $E_P[f/\alpha]$. For example, a law equally weighting magnitude
pairs (1,1) and (100,10) gives a readout of $101/11=9.1818$, although the mean
of its two mass ratios is 5.5. This is an explicit regression test.

Every numerical case in the main table uses **two spatial dimensions and unknown
common direction**. Magnitudes are in N and m/s²; masses are in kg. Standard
uncertainties are per component. The illustrative absolute discrepancy threshold
is 9, matching the working specification's playground example.

## 2. Four outputs to inspect in every experiment

| Output | Question it answers | What it does not establish |
|---|---|---|
| Point $\widehat m$ | What scalar does the agreed readout summarize? | Identification, concentration, or model fit |
| Induced law of $M=f/\alpha$ | How are ratios distributed under this fitted latent-pair law? | A calibrated repeated-sampling confidence interval |
| Fixed-mass profile $Q(m)$ and its threshold set | Which fixed positive masses can explain the measurements with small enough corrections? | Probabilities assigned to those masses |
| Repeated-sampling law of $\widehat m$ | What happens if the stated experiment is repeated at a stated true latent pair? | Uncertainty from one dataset without supplying the sampling setup |

The **free-mass fit** is $Q_{\min}=\inf_{m>0}Q(m)$. It checks whether any mass
fits; it is distinct from $Q(\widehat m)$, which checks the reported point.

For independent isotropic errors,

$$Q(m)=\frac{|\widetilde F-m\widetilde a|^2}
{\sigma_F^2+m^2\sigma_a^2}.$$

The set here is $\{m>0:Q(m)\le9\}$, an **absolute discrepancy set**. It is not
a profile-likelihood-drop interval and is not labelled a 95% confidence set.
If $Q_{\min}>9$, the working specification rejects the account and withholds
the conditional mass. The calculations retain that conditional point for inspection.

## 3. Start with exact inputs

| Exact force and acceleration | Outcome |
|---|---|
| Both nonzero and codirectional | Recover the unique positive magnitude ratio |
| Both zero | Every positive mass fits; no mass is identified |
| Exactly one zero | Contradiction with finite positive-mass Newton II |
| Nonzero perpendicular, anti-aligned, or otherwise non-codirectional | Contradiction with positive-mass Newton II |

Uncertainty changes the input, not just the numerical precision of that table.
An instrument display of zero with a positive standard uncertainty enters the
empirical branch. Supplying zero empirical covariance is explicitly rejected
by this numerical implementation; it does not silently regularize an exact premise.

## 4. The measured zero-zero experiment

With both readings zero and the default isotropic independent model, put
$s=\sigma_F/\sigma_a$. The test oracle is

$$\widehat m=s,\qquad M\sim\operatorname{HalfCauchy}(s),\qquad Q(m)=0
\quad\text{for every }m>0.$$

With both uncertainties equal to 1 in their respective units, the point is
1 kg and the induced central 95% interval is **0.03929–25.45170 kg**. Every
positive mass is nevertheless compatible by the profile criterion. Those
statements coexist because they answer different questions.

The point's scale comes from the apparatus/model. Doubling force uncertainty
doubles the point and every ratio quantile without excluding any mass from
the compatibility set. Reducing **both** uncertainties by a common factor,
even a millionfold, leaves the entire ratio distribution unchanged. The
positive magnitudes themselves shrink; their ratio does not become better known.

Repeated independent observations of **the same latent pair**, whose sufficient
sample means are both exactly zero, shrink both effective uncertainties by
$1/\sqrt n$ and likewise leave the mass law unchanged. This is a conditional
statement about those sufficient statistics, not a claim that actual repeated
noise readings will always average to exactly zero.

The automated tests check null point, quantiles, scaling, and profile behavior
in one, two, and three dimensions. The half-Cauchy result uses the flat
**magnitude** measure; inserting Cartesian volume factors changes the experiment.

### What the original infinity surprise does and does not mean

Under this default model, the induced ratio has an infinite positive mean and
second raw moment, so there is no conventional finite standard deviation about
a finite mean. Its mean-square error about the finite reported point is also
infinite. But this also occurs at every finite empirical input under the same
nondegenerate Gaussian/flat-magnitude policy, including excellent aligned data.

The mechanism is an arbitrarily small, nonzero probability density near
$\alpha=0$. An extremely tiny coefficient still leaves the $1/M^2$ tail with
divergent moments. High signal makes that coefficient tiny; it does not
mathematically remove it. A finite Monte Carlo sample can completely miss that
tail and return a reassuring finite sample standard deviation.

Consequently, **infinity alone is not a null detector**. The stronger null facts
are the flat zero profile and the absence of narrowing when both apparatus
scales shrink together. The September 5 managed research already established
this qualification; these tests preserve it rather than rediscovering it as a
new policy decision.

At null, $\log(M/s)$ has density $1/(\pi\cosh l)$ and standard deviation
$\pi/2$. This finite multiplicative-spread description is tested. It does not
replace the unrestricted compatibility set. Away from null, neither a geometric
mean nor the ratio median must equal the ratio-of-means point.

## 5. One zero is different from two

| Measured vectors; uncertainties both 1 | Point | Induced central 95% interval | Compatible masses at threshold 9 |
|---|---:|---:|---|
| F=(0,0), a=(10,0) | 0.08020 | 0.003150–0.2313 | $0<m\le0.3145$ |
| F=(10,0), a=(0,0) | 12.4693 | 4.324–317.5 | $m\ge3.1798$ |
| F=(0,0), a=(2,0) | 0.46576 | 0.01830–7.443 | Every positive mass |
| F=(0,0), a=(0,0) | 1 | 0.03929–25.452 | Every positive mass |

The first trial says that appreciable acceleration needed very little force;
the second says that appreciable force produced very little acceleration.
They constrain opposite ends of the mass range. The third shows that a nonzero
display alone is insufficient: the acceleration is not sufficiently resolved
to exclude either end at the chosen threshold.

For a well-resolved nonzero channel, the known-direction limits use
$c=\sqrt{2/\pi}$:

$$\widehat m\simeq c\sigma_F/|a|\quad(F\text{ reads zero}),\qquad
\widehat m\simeq |F|/(c\sigma_a)\quad(a\text{ reads zero}).$$

The finite-direction-uncertainty values above need not equal these limiting
expressions exactly. The known-direction formula is tested against independent
truncated-normal means; the full vector construction approaches the limit as
the resolved channel strengthens.

Free-mass $Q_{\min}$ is zero for these one-zero trials: its infimum is reached
as mass approaches zero or infinity. This does not produce an attained finite
best-fit mass. A previously known mass can still be rejected: for F=(0,0),
a=(10,0), unit uncertainties, the fixed mass 1 kg has $Q=50$.

## 6. Alignment, signal strength, and approaching zero

Strong aligned data F=(20,0), a=(10,0), unit uncertainties give a point of
2 kg and an induced 95% interval of **1.62555–2.53459 kg**. Weak aligned
data F=a=(1,0), with the same uncertainties, give a point of 1 kg but an
interval of **0.05323–18.78785 kg**. Perfect visual alignment is not precision.

The excitation sweep holds the true ratio at 1 and raises both measured
signal-to-noise ratios from zero to 20. The point stays 1 by symmetry, while
the induced interval narrows. This is information supplied by stronger signals,
not by numerical quadrature refinement.

There is no single undifferentiated "limit at zero":

- Shrink the readings at fixed positive uncertainties: excitation disappears,
  and the result approaches the apparatus-controlled null law.
- Shrink readings and uncertainties together: the signal-to-noise ratios and
  induced mass law stay the same.
- Keep a nonzero aligned pair and shrink uncertainties: recovery approaches
  the exact ratio with a concentrated central distribution.
- Stipulate both inputs exactly zero: the exact branch remains non-identifying.

The suite checks the first two as separate paths, alongside strong-signal
recovery. It does not assign one mass to the exact origin by continuity.

## 7. Rotate the vectors without changing their strengths

Equal standardized strengths are a revealing symmetry test. If
$|F|/\sigma_F=|a|/\sigma_a$, the isotropic model's channel-exchange symmetry
gives $\widehat m=\sigma_F/\sigma_a$ at every angle. Therefore the point can
remain unchanged while the physical account deteriorates drastically.

For the angle sweep $|F|=|a|=5$ and both uncertainties 1:

- At 0 degrees the point is 1, the profile fits at that point, and the central
  ratio distribution is concentrated.
- At 90 degrees the point is still 1, but $Q(m)=25$ for every mass and the
  model is rejected at threshold 9.
- At 180 degrees the point is still 1, $Q_{\min}=25$ and $Q(1)=50$; again
  the model is rejected. The conditional ratio law is broad and split across
  small- and large-mass explanations.

### Perpendicular vectors can reproduce the null ratio law

The balanced perpendicular case F=(0,$\sqrt6$), a=($\sqrt6$,0), with unit
uncertainties, has the **same half-Cauchy ratio law as measured zero-zero**.
But its constant profile is 6, not 0. Stronger balanced perpendicular readings
can keep that same normalized ratio law while failing the compatibility gate.

This is visible in the independently computed curves. Geometrically, rotation
integration depends on the length of the combined standardized direction
evidence; equal orthogonal channel vectors make that length independent of
the ratio angle. This explanation is an interpretation of the tested model,
not a new assumption about all measurement policies.

It follows that even the **entire normalized ratio distribution** is not a
substitute for absolute model fit. Conditioning and normalization discard the
overall cost of forcing a compatible account.

## 8. Anti-alignment: two different failure patterns

The specification's stress case F=(-10,0), a=(10,0), both uncertainties 2,
gives $Q_{\min}=25>9$. The conditional point is 1 kg, but the result is
**rejection**. Its induced interval, approximately 0.002068–483.6 kg, is
diagnostic output under the rejected model, not an accepted mass determination.

A less extreme case exposes a separate issue:

$$F=(-\sqrt6,0),\quad a=(\sqrt6,0),\quad \sigma_F=\sigma_a=1.$$

Here the tested outputs are

$$\widehat m=1,\qquad Q_{\min}=6,\qquad Q(\widehat m)=12.$$

At threshold 9, the full compatibility set is

$$ (0,2-\sqrt3]\ \cup\ [2+\sqrt3,\infty)
\;\approx\;(0,0.267949]\cup[3.732051,\infty).$$

**The overall fit passes, but its reported point is outside the compatible
set at the same threshold.** One group of latent explanations nearly erases
force and favors small mass; another nearly erases acceleration and favors
large mass. The ratio of their mean magnitudes lies between the groups. It is
not required to minimize the discrepancy or lie inside a nonconvex compatibility
set. The ratio median is also 1 here by symmetry, so switching to the median
does not remove the issue in this case.

The induced central 95% interval is about 0.009374–106.677 kg. As an interval,
it crosses the central region between the two explanations. Inspecting the
shape conveys information that this pair of quantiles conceals.

**Handoff question:** should a globally admissible estimate whose central
readout fails a fixed-mass compatibility check carry an explicit warning or
different downstream standing? This laboratory records the mismatch; it does
not silently change the estimator, add a new rejection rule, or settle that
policy. Requiring the readout itself to belong to the threshold set would be
an additional contract beyond the current free-mass gate.

Unequal channel resolution removes the symmetric ambiguity. F=(-2,0),
a=(5,0), unit uncertainties, gives a small point near 0.07748 kg: moving
the less resolved force toward zero is cheaper. Anti-alignment alone is not
an empirical rejection rule.

## 9. How to think about uncertainty in general

1. **Start with the measurement law.** State units, dimension, covariance,
   whether direction is known, and the latent measure. Standard uncertainties
   without these details do not define the whole inference.
2. **Look at shape, not just a width.** Positive ratio distributions are
   asymmetric and may have separated regions of support. Equal-tail intervals
   can hide a low-density central region.
3. **Use multiplicative language when appropriate.** A lower/upper interval
   or log spread often describes ratios better than an additive ± number.
   At high signal and for a locally concentrated law, the delta approximation
   is $\operatorname{Var}(\log M)\approx
   \operatorname{Var}(f)/E[f]^2+\operatorname{Var}(\alpha)/E[\alpha]^2
   -2\operatorname{Cov}(f,\alpha)/(E[f]E[\alpha])$.
   This is a local description, not a claim that the exact raw ratio variance
   exists. Covariances here belong to the fitted magnitude law.
4. **Keep fit and identification visible.** Broad, narrow, or finite uncertainty
   is not itself a goodness-of-fit test. A pass says some account fits, not
   necessarily that the point is well supported or that the data determine mass.
5. **Transform inverse mass from the same inference.** The point becomes
   $1/\widehat m$, and interval endpoints become $[1/U,1/L]$. The tests check
   both transformations, including correlated channels.
6. **Separate numerical error from physical uncertainty.** More quadrature
   nodes reduce integration error. More Monte Carlo draws reduce simulation
   error. Neither changes the data's physical information.

Anisotropy is tested by rotating both readings and their covariance: this
preserves the answer. Rotating just the readings in a fixed anisotropic
apparatus can change it. Correlation is tested against an independent
positive-quadrant Monte Carlo calculation. Symmetric positively correlated
null channels can retain the point 1 while producing a narrower induced
ratio distribution than independent channels. Therefore even the half-Cauchy
null law is conditional on the declared independence/isotropy policy.

The source research's alternative null measure with chi-3 magnitudes keeps
the same symmetric point but yields finite ratio variance. That is a
documented sensitivity example, **not an implemented alternative policy in
this suite**. Non-Gaussian errors, hard bounds, and other latent measures
require their own experiments before generalizing today's results.

For general metrological background, the [JCGM publications](https://www.bipm.org/en/committees/jc/jcgm/publications)
include GUM guidance and Supplement 1 on propagating distributions with Monte
Carlo. They are background references, not sources endorsing this VD estimator.

## 10. Repeated trials and gate calibration

The saved experiment performs 400 full 2D reruns in each of three Gaussian
sampling regimes with true mass 2, unit channel uncertainties, and true
acceleration 0, 1, or 5. These are **conditional candidates including rejected
trials**, to avoid silently changing the sample by selection.

| True acceleration | Mean candidate | Median candidate | Candidate sample SD | Empirical coverage of induced 95% interval |
|---:|---:|---:|---:|---:|
| 0 | 1.1786 | 1.0064 | 0.6887 | 99.0% |
| 1 | 1.7013 | 1.3857 | 1.2093 | 98.5% |
| 5 | 2.1296 | 2.0512 | 0.5544 | 95.25% |

Monte Carlo standard errors of the mean candidates are 0.0344, 0.0605, and
0.0277 respectively. Coverage standard errors are about 0.50, 0.61, and 1.06
percentage points. These are a pilot sampling study, not a proof of unbiasedness
or universal coverage. The low-excitation results visibly depend on apparatus
scale, not just the true mass. Finite sample SD does not establish existence
of an infinite-population moment. At the null true pair the data distribution
is independent of the stipulated positive mass.

A separate 100,000-trial simulation per dimension/excitation combination checks
the free-mass gate (900,000 trials total). With strongly excited equal channels,
threshold 9 rejected approximately **0.250% in 2D** and **1.035% in 3D**.
At the null, those rates were about 0.009% and 0.052%. Their binomial Monte
Carlo standard errors are saved in the experiment data. A zero observed count
in some 1D regimes is not a zero true probability; the saved plug-in standard
error is then zero and should not be interpreted as a confidence bound.

Thus 9 is not a universal "three-sigma" gate with one false-rejection rate.
Dimension, excitation, and boundary geometry matter. These simulations explore
its calibration; they do not install a replacement threshold.

Two individually aligned trials can each fit perfectly but fail to admit one
common mass. The suite checks F/a=2 and F/a=8 with tight uncertainty. Their
joint profile is a sum at **one shared mass**, not a sum of separate minima.
This does not define a common-mass posterior or authorize averaging separately
fitted mass points. The working specification leaves that measure/model choice
open, and this laboratory preserves the boundary.

## 11. Verification, coverage, and what is ready

The suite contains 40 named tests, many with subcases and randomized checks.
They include independent truncated-normal formulas, direct radial quadrature,
direct full-vector integration in different coordinates, a correlated
Monte Carlo oracle, direct generalized least-squares verification of the
profile, and 3,750 randomized set-membership comparisons.

The report's 59 saved cases comprise 12 named cases, 37 angles, and 10
excitation levels. Refinement from 128/1024 to 256/2048 direction/ratio
resolution changed the point by at most **1.26e-9 relative** and the upper
97.5% quantile by at most **0.285% relative** in this set. Quantiles in severe
anti-alignment converge more slowly than the point. Reported precision is
limited accordingly. These are observed refinement checks, not rigorous
global error bounds for arbitrary inputs.

This is broad coverage of the declared Gaussian compatible-pair construction.
It is not exhaustive across measurement policies, arbitrary conditioning,
extreme numerical ranges, all possible physical experiments, or every VD
provenance/admission rule. Numerical stability at still higher SNR, extreme
covariance condition numbers, non-Gaussian models, hard bounds, and common-mass
posterior integration remain explicit extension targets.

The principal new entry-facing issue is the globally passing, locally
incompatible central point in section 8. Preserve the whole set and shape
when reviewing that case. The equation has been left intact.

## 12. Provenance and managed-context coverage

**Managed snapshot:** `canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`.
`main` was resolved through the installed GitHub connector on 2026-09-06.
`catalog.yaml` supplied inventory and record status. All fetched managed
task sources used this commit:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
  — construction, exact branch, direction integration, covariance, gate,
  anti-alignment fixture, uncertainty and output contract.
- `Newton-analysis/_dscn_mass_estimator/mass_estimation_zero_zero_two_cases.md`
  — competing null admission positions and the original uncertainty surprise.
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`
  — inspected especially sections 5–6: universal finite-input tail issue and
  probability/compatibility/sampling distinction.
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md`
  — null scaling, aggregation, finite log spread, and scope limits.

The construction is a working synthesis record, not an assertion that every
policy item is settled. Historical stationary-root and median proposals were
identified through the catalog and specification; they were not substituted
for the current readout or exhaustively reread.

**Local active/unsubmitted scope:** filename and text searches excluding
`VDfirst/`, using mass/estimator/readout, zero-zero, 0/0, uncertainty, and
anti-alignment terms. Inspected the full local
`pass to VD-docs/Mass estimator/mass_readout_composition_uniqueness_2026-09-06.md` and the
source/drafting record `mass_estimation_introduction_candidate_sources_2026-09-06.md`;
scoped search excerpts from the other September 6 introduction/development and
shared-quantity notes supplied active context. These are newer local authoring
material and are not assumed absorbed into the managed snapshot. The proof
was used only to fix which readout to test, not extended.

The local engine/evaluator surfaces were searched for an existing estimator;
none was found in that scoped search. The specification refers to an earlier
interactive playground, whose implementation was not recovered or assumed to
be this code. This is a new independent reference implementation. No separate
inbox/unclassified collection or historical local archive was exhaustively
audited, and no claim of exhaustive VD-document coverage is made.

**Custody:** the tests, reference code, and results remain on the local active
evaluator surface. This report is in the mandated `pass to VD-docs/` staging
folder. No managed sources, prior entry packets, existing local authoring
documents, or `VDfirst` files were edited. Saving the report does not submit it
to VD-docs, and no commit or push is implied.

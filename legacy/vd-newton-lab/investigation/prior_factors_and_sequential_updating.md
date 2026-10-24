# Readout regimes, prior factors and sequential updating

**Working contribution 12 — 23 September 2026.** User-requested record of a
conversation with Claude, held from the VD-docs project, which re-examined what
the estimator can claim and then built reusable code for it. Newton II is
assumed throughout. Every simulation uses known independent isotropic Gaussian
noise in 3D with one effective reading per channel. Not submitted to VD-docs.
No existing module, protocol, result or estimator was changed; the new code and
studies are additive.

Three of the user's questions steered the work:

> "When should I stop using the mean version and start using the median
> version. That is the question! … The issue is thinking that what is implied
> is 1 equation, when it's an infinite class. And that infinite class of
> piecewise functions needs to be calibrated using an added parameter. Which
> includes it's own uncertainties that people will then ignore."

> "As long as we don't derive lambda, it is a parameter that is being fitted.
> Lets say that more clearly: 'It is a parameter of our equation.'"

> "Why not interpret M as simply… you previous mass measurement best fit. Then
> the equation looks to be just an update rule. Previous estimate of mass + new
> info -> new estimate."

All numbers below come from the modules and saved outputs listed in
[Code, studies and reproduction](#code-studies-and-reproduction). Exploratory
scripts used earlier in the conversation were not reusable and are not the
evidence for anything here; two of their conclusions changed when rerun
properly, and [§9](#9-what-changed-from-the-exploratory-numbers) says which.

## Summary

1. **Among the eleven archived rules, the ratio of means is unbeaten in a
   band, and it is the minimax rule on grids with both SNRs ≥ 1.** With a
   paired Bonferroni bound, it beats every other archived rule except its
   tube-measure sibling where the weaker channel's effective SNR is 1, or 2
   with the other channel at 12–16. At factor 1.5, 15 of 16 archive cells
   reappear on fresh noise (all 12 at factor 2; 7 of 9 at factor 1.25). On
   these grids its worst-case regret, against the per-cell best of the eleven
   rules, is the smallest of any rule: 4.7–6.8 pp.
2. **The band has edges.** With one channel at SNR 0.5 and the other at 2.5 or
   more, the flat law's own reciprocal-root point is clearly ahead of every
   rule in 6 of these 8 cells at factor 1.5 and in all 8 at factor 2 (none at
   1.25), by 0.8–2.9 percentage points (`pp`) over the next best. The ratio of
   means trails it by up to 23 pp. Once such cells are included, the minimax
   rule becomes another summary of the same law: the median at factors 1.25
   and 1.5, the geometric point at factor 2, tied on average over factors. The
   ratio of means drops to 5th of 11 (15–23 pp), the worst of the flat-law
   summaries; the non-flat rules are worse still. Where both channels are moderate, the floored
   vector-length ratio is resolved ahead in about half of the fresh cells, by
   0.2–6.8 pp.
3. **Every switch threshold loses somewhere.** Two-branch pretest rules
   (|F|/|a| above an observed-SNR threshold, a fallback below) are clearly
   worse than the ratio of means in 12–45 of 64 cells at every threshold. The
   largest gain in any cell is 7.0 pp. A switch can beat both of its branches
   in a cell (one does at (1,1)), so choosing between the fixed branches with
   knowledge of the true cell (at most 4.7–6.8 pp) is a reference, not a
   ceiling.
4. **The tilt λ is a fitted parameter of the equation, and it carries its
   fitting range.** Fitted on the operating-range archive by a declared
   criterion, λ\* = 0.35 cuts worst-case regret by a third to a half on fresh
   draws with both SNRs ≥ 1 and makes it worse below the fitted range. The
   hindsight-optimal λ moves with the test set (0.35, 0.4, 0.7, and −2 or
   lower).
5. **λ is a prior factor on latent mass, centred at the instrument ratio
   σ_F/σ_a.** Its value can depend on quantities known before the experiment
   (the SDs and what is believed about the object) but must not depend on the
   readings, or it becomes a switch. A two-parameter version (λ, m₀) moves its
   centre towards the masses (2s for heavy objects, s/2 for light) when both
   channels have SNR ≥ 1.
6. **The flat law shrinks towards σ_F/σ_a at weak signal.** A light object seen
   through a weak force channel reads 3.6× too heavy (median).
7. **A declared log-normal prior is priced.** Overconfidence is catastrophic.
   An honest wide prior that tempers the flat law helps a little on average;
   its worst cells lose up to 14 pp when its centre is off by a factor of four
   or more. An exact *wide* log-normal prior is worse than the flat default on
   average.
8. **As an update rule, only one of three combination rules converges.** The
   symmetric rule (both channels' excitation as one nuisance, prior counted
   once) converges to the true mass, though its 80% intervals undercover.
   Contribution 11's eq. (18) settles below the true mass, and eq. (20)
   settles between the truth and σ_F/σ_a. Exact identities explain why.

## 1. Where the ratio of means is unbeaten

[`paired_dominance.py`](paired_dominance.py) compares one focal rule with
every other archived point rule on identical replicates. For each cell and
tolerance factor K it tests the paired failure difference against
z × MCSE, with z the two-sided Bonferroni normal quantile over cells ×
factors × compared rules (z = 4.18 over 1,728 comparisons without tube; 4.21
over 1,920 with it). This uses the pairing, unlike the simultaneous DKW band
of the operating-range protocol, which bounds each marginal separately and is
more conservative. It says nothing about unsampled cells.

Cells (force SNR, acceleration SNR) where `flat_joint` beats every rule except
`tube_joint`:

| K | Operating-range archive | Fresh validation grid |
|---|---|---|
| 1.25 | (1,1) (2,8) (2,12) (2,16) (8,1) (12,1) (12,2) (16,1) (16,2) | (1,1) (1,8) (1,12) (1,16) (2,12) (2,16) (8,2) (12,1) (12,2) (16,1) (16,2) |
| 1.5 | (1,1) (1,3)–(1,16) (2,8) (2,12) (2,16) (6,1) (8,1) (12,1) (12,2) (16,1) (16,2) | same, with (8,2) in place of (2,8) |
| 2 | (1,1) (1,4)–(1,16) (3,1) (4,1) (6,1) (8,1) (12,1) (16,1) | same plus (1,3) |

"(1,3)–(1,16)" means every sampled acceleration SNR from 3 to 16 at force
SNR 1. The band is nearly symmetric under exchanging the channels. The
remaining asymmetry, with (1,3) and (1,4) in the band at K = 1.5 but (3,1) and
(4,1) only at K = 2, is Monte Carlo noise at the resolution threshold. There
the unresolved rivals are `flat_median` and `flat_geometric`, which are
reciprocal-symmetric, and the margins are 0.6–1.9 pp against MCSEs of
0.3–0.5 pp. The swap of (2,8) and (8,2) between the archive and the fresh draws
is the same effect. Within the band at K = 1.5 on fresh draws:

- the closest rival is usually another summary of the same flat law
  (`flat_median` or `flat_geometric`), 1.4–4.9 pp behind;
- the vector-length ratio is 2.4–22.9 pp behind;
- absolute success is 43–66%. Unbeaten is not accurate.

Counting `tube_joint` as a rival leaves only (1,1) at every factor, plus
(2,8) in the archive at K = 1.25 and 1.5, and (8,2) on fresh draws at K = 1.5. The unbeaten object is the
average-then-divide readout rather than the particular flat measure.

**The edges**, from the off-grid study (effective SNRs 0.5, 1.5, 2.5, 5, 10, 20):

- With one channel at SNR 0.5 and the other at 2.5 or more,
  `flat_reciprocal_root` = E[√M]/E[1/√M] under the same flat law beats every
  rule in 6 of these 8 cells at K = 1.5 and in all 8 at K = 2. Its lead over
  the next best rule, `flat_geometric`, is only 0.8–2.9 pp. At K = 1.5 the
  ratio of means trails the best rule there by 8.2–22.6 pp. The cause is the
  shrinkage of [§6](#6-the-flat-law-shrinks-towards-the-instrument-ratio):
  median m̂/m is 2.0–2.6 for the light objects and 0.39–0.49 for the heavy
  ones in those cells.
- At K = 1.5, in 11 off-grid cells with both SNRs between 1.5 and 20, the
  floored vector-length ratio is resolved ahead by 1.4–5.5 pp, most at
  (2.5, 2.5).

**Minimax.** Among the eleven archived rules, the ratio of means has the
smallest worst-case regret on both SNR ≥ 1 grids:

- archive, K = 1.25 / 1.5 / 2: 4.8 / 6.1 / 4.8 pp;
- fresh draws: 5.2 / 6.8 / 4.7 pp.

In each case the worst cell is (3, 3) against `norm_floor`, and the next best
rule is `flat_median` (5.5–7.9 pp). This fits with `norm_floor` being resolved
ahead in 28–35 of the 64 fresh cells, by 0.2–6.8 pp: it wins many cells narrowly
and loses the weak-channel band heavily.

On the off-grid set that includes SNR 0.5 the order changes:

- `flat_median` has the smallest worst-case regret at K = 1.25 and 1.5
  (8.2 and 12.5 pp);
- `flat_geometric` has it at K = 2 (13.1 pp against the median's 15.0), and
  the two are tied on the mean over K;
- the ratio of means ranks 5th of 11 (14.8 / 22.6 / 22.9 pp), the worst of the
  flat-law summaries;
- the non-flat rules are worse still (24–72 pp).

This is the regime split the user recalled from the earlier median proposal. The reciprocal-root
point is itself a same-law readout that does not satisfy the composition
requirements. Whether a composition-respecting construction can recover its
subunit performance is open (see §4 for a partial answer).

## 2. Switch (pretest) rules

[`switch_rules.py`](switch_rules.py) rescores saved trials without refitting.
The rule reports |F|/|a| when min(|F|/σ_F, |a|/σ_a) of the *observed* readings
exceeds τ, and a fallback point otherwise. It is compared with the ratio of
means on identical replicates, with a Bonferroni family of 3,456 comparisons
(z = 4.34). Fresh validation grid, fallback `flat_median`, K = 1.5:

| τ | Cells switch clearly worse | Clearly better | Worst cell | Best cell | Mean over cells |
|---:|---:|---:|---:|---:|---:|
| 1 | 17 | 36 | −22.6 pp | +6.8 pp | −1.7 pp |
| 2 | 24 | 32 | −13.9 | +6.8 | −2.1 |
| 3 | 35 | 18 | −12.8 | +4.1 | −2.2 |
| 4 | 41 | 10 | −4.9 | +1.9 | −1.7 |
| 6 | 44 | 1 | −4.9 | +0.3 | −1.4 |

The archive gives the same table to within about a percentage point, and
factors 1.25 and 2 have the same shape. At K = 2 on fresh draws, no threshold
of 3.5 or more is clearly better in any cell.

An oracle that picks, in each cell, the better of the two fixed branches *with
knowledge of the true cell* gains at most 6.8 pp (K = 1.5), 5.2 pp (K = 1.25)
and 4.7 pp (K = 2) over the ratio of means. It is a reference point, not a
ceiling. A switch chooses per reading, and that choice can carry information:

- the best tested switch gains 7.0 pp, at (3, 3);
- at (1, 1) the τ = 1.5 switch is clearly better than the ratio of means
  (+2.1 pp), although each of its branches alone is clearly worse there.

So "unbeaten" in §1 means unbeaten among the archived rules. What no threshold
achieves is safety: each is clearly worse somewhere. This is the
preliminary-test estimator problem of statistics. The uncertainty of the branch
choice goes unreported when the chosen branch is treated as fixed; that is the
"own uncertainties" of the user's question.

The `noise_ratio` fallback (σ_F/σ_a) is also saved. With unit SDs it equals
the true mass wherever the two true SNRs are equal, so its apparent gains on
this grid are an artefact and are not interpreted.

## 3. Fitting λ: a parameter of the equation

[`weighted_posterior.py`](weighted_posterior.py) multiplies the flat law by a
declared factor of latent mass. The readout remains E[f]/E[α] under the
weighted law. The sech tilt is

$$
w_\lambda(M)=\exp\!\bigl(\lambda\,[\operatorname{sech}(\log(M/c))-1]\bigr),
$$

which for c = 1 is the README §4 family ρ_λ = exp(λ sin 2θ) up to a constant.
[`prior_weight_study.py`](prior_weight_study.py) scores a grid of λ on saved
trials. The regret in a cell is the best success among all eleven archived
rules minus the tilted law's success.

**Fit.** The declared criterion is the mean over K ∈ {1.25, 1.5, 2} of the
largest cell regret. On the operating-range archive (fit set, 64 cells) it
selects **λ\* = 0.35**. It was frozen before any λ was evaluated on
validation data; the λ-free dominance and switch analyses of the validation
grid ran first.

**Validation.** Largest cell regret at K = 1.25 / 1.5 / 2, in pp:

| Data | λ = 0 | λ\* = 0.35 | λ minimizing the criterion on this data |
|---|---|---|---:|
| Fit set: operating-range archive | 4.8 / 6.1 / 4.8 | 3.2 / 3.4 / 2.9 | 0.35 |
| Fresh grid, same cells, new seeds | 5.2 / 6.8 / 4.7 | 3.5 / 4.1 / 3.1 | 0.4 |
| Off-grid, both SNRs ≥ 1.5 (25 cells; SNR 20 exceeds the fitted grid's 16) | 4.0 / 5.5 / 3.7 | 2.5 / 2.5 / 1.6 | 0.7 |
| Off-grid below the range (11 cells with one SNR 0.5) | 14.8 / 22.6 / 22.9 | 15.8 / 24.1 / 25.9 | −2 (smallest tested) |

Mean success changes by only +0.2 to +0.7 pp in the first three rows. The gain
is the removal of the worst balanced-cell losses. Below the fitted range λ\*
makes things worse, and the best tilt there is *negative*: it pushes
estimates away from σ_F/σ_a and so counters the shrinkage of §6. −2 is the
smallest λ tested, so the optimum there is not bracketed.

The fitted value is therefore not a constant of nature. It encodes the
distribution of conditions it was fitted on. That is the user's "added
parameter" with its own uncertainty, measured.

## 4. What λ is: a prior factor on mass

In standardized coordinates the latent mass is M = s tan θ with
s = σ_F/σ_a, and

$$
\sin 2\theta=\frac{2\tan\theta}{1+\tan^2\theta}=\operatorname{sech}\bigl(\log(M/s)\bigr).
$$

So ρ_λ depends on M alone: it is a prior factor on the mass, centred at the
instrument ratio, with λ setting its concentration. The implementation
reproduces the README §4 null table exactly (points 1; log SDs 2.852883,
1.570796, 0.666264, 0.263564, 0.126513; the same intervals) through the
general weighted-law code path.

This answers the user's follow-up question: "does that imply lambda is
effectively a function of concepts knowable before the experiment…? Does that
mean it might depend on F, a their uncertainties somehow?"

- It may depend on the supplied SDs (through s), and on designed or expected
  force and acceleration, i.e. on what is believed about the object before
  weighing it.
- It must not depend on the measured readings. A λ chosen from the readings
  is a data-dependent switch (§2), with the same unreported selection
  uncertainty. With a single pair, a prior also cannot be estimated from the
  data.

**Two parameters: the centre as well.** The tilt
exp(λ[sech(log(M/m₀)) − 1]) was fitted on the first 5,000 replicates of the
shifted-mass study and scored on the other 5,000. "Heavy" cells have
m = 2s, 4s or 8s and "light" cells the reciprocals, with the weaker channel's
SNR 0.5–3. The values are held-out means over K of the largest regret, in pp.

| Cells | λ = 0 | Best tilt centred at s | Best (λ, m₀) |
|---|---:|---:|---:|
| Heavy, both SNRs ≥ 1 (12 cells) | 2.8 | 2.5 (λ = 0.25) | **1.8** (λ = 0.25, m₀ = 2s) |
| Light, both SNRs ≥ 1 (12 cells) | 2.6 | 2.3 (λ = 0.25) | **1.6** (λ = 0.25, m₀ = s/2) |
| Heavy, including SNR 0.5 cells (15) | 14.8 | **7.5** (λ = −1) | 7.9 (λ = −1, m₀ = s/2) |
| Light, including SNR 0.5 cells (15) | 15.7 | **7.4** (λ = −1) | 8.4 (λ = −1, m₀ = 2s) |

With both channels at SNR ≥ 1 the fitted centre moves towards the masses and
mirrors exactly between the heavy and light grids. It sits at the near edge of
each grid's mass range (masses 2s–8s and s/8–s/2), and there are only two
grids. Once subunit cells are included, those cells dominate the worst case:

- the fit turns negative (−1 is the smallest λ in this grid, so it is not
  bracketed);
- a centre adds nothing on held-out data;
- the tilt's job becomes undoing shrinkage rather than expressing knowledge of
  the object.

**A correction to the conversation.** Claude said the flat law "is a
half-Cauchy prior on M". That holds under the symmetric split of the flat
reference, df dα = σ_F σ_a r dr dθ, with uniform θ. Under the split
α dα dm the same measure is flat in m with weight α dα on the acceleration
magnitude. The two readings agree for one reading and differ when readings are
combined; §7 shows which one survives repetition.

## 5. Pricing a declared log-normal prior

Two ways to add a normal belief about log mass of width w (a factor e^w of
uncertainty):

- **A, tempering:** the flat law times exp(−(log(M/c))²/(2w²));
- **B, replacing:** the flat law's implied angle factor sech(z)/2 is replaced
  by an exact normal density on log mass.

The centre c is placed at a stated multiple of each cell's true mass. That
oracle placement prices a prior that is wrong by a known factor; it is not an
estimator anyone can run. The data are the 15 excited cells of the fresh mass
ladder (m ∈ {s/4, s, 4s}, acceleration SNR 1–6). Entries are the change in
success against the flat law at K = 1.5, mean over cells / worst cell, in pp:

| A, width | centre ×⅛ | ×¼ | ×½ | **correct** | ×2 | ×4 | ×8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ×1.3 (w = 0.25) | −56 / −90 | −52 / −89 | −34 / −57 | **+42 / +3** | −40 / −61 | −58 / −94 | −58 / −97 |
| ×1.6 (w = 0.5) | −40 / −66 | −32 / −56 | +10 / −32 | **+34 / +2** | −8 / −33 | −39 / −59 | −50 / −76 |
| ×2.7 (w = 1) | −14 / −50 | +8 / −33 | +17 / −8 | **+12 / +1** | +3 / −11 | −9 / −34 | −21 / −51 |
| ×7.4 (w = 2) | +6 / −12 | +6 / −5 | +5 / 0 | **+4 / 0** | +2 / −3 | 0 / −8 | −3 / −14 |

Variant B at w = 2 is worse than the flat law on average at every centre,
including the correct one (−2.4 mean, −15 worst). Per cell it still helps the
three light cells with force SNR below 1 (+10 to +28 pp when correctly
centred), where it counters the shrinkage of §6. Restricting to the 12 cells with
both SNRs ≥ 1, it is −8.2 / −15 when correctly centred. A wide exact log-normal
prior tends to a flat prior on log mass, the measure excluded on 2 September
for collapsing weak readings onto the boundaries.

What follows:

- **A prior claiming more certainty than it has is worse than none.** It gains
  42 pp when right and loses 40 pp when off by a factor of two.
- **An honest prior is fairly safe.** Variant A helps on average when its
  centre is within roughly its stated width of the truth. In this grid it is
  more tolerant of centres below the truth than above. At w = 2 the worst cell
  loses 5–14 pp only when the centre is off by a factor of four or more.
- **"Vague" is not "ignorant."** The flat law's implied half-Cauchy is a much
  better default than a vague log-mass prior, so a stated prior should *temper*
  the flat law (A), not replace its angle factor (B).

## 6. The flat law shrinks towards the instrument ratio

Median m̂/m of the unweighted ratio of means on the fresh mass ladder
(10,000 replicates per cell):

| True mass | acc SNR 1 | 2 | 3 | 4 | 6 |
|---|---:|---:|---:|---:|---:|
| s/4 (force SNR a/4) | 3.60 | 2.63 | 1.82 | 1.42 | 1.12 |
| s | 1.00 | 1.01 | 1.01 | 1.01 | 1.00 |
| 4s (force SNR 4a) | 0.72 | 0.97 | 1.00 | 1.00 | 1.00 |

At weak signal the null value σ_F/σ_a leaks into informative data. The
crossover analysis and the operating-range protocol already note that accuracy
near the noise-ratio mass can reflect this shrinkage. Here it is quantified
away from that mass, and it explains both the subunit losses of §1 and the
negative λ fitted in §4.

## 7. The estimator as an update rule

Reading π(M) as what was known before a measurement makes the single-reading
equation one Bayesian update from a default prior: the half-Cauchy at
σ_F/σ_a. Contribution 11 already fixes two requirements:

- **Carry curves, not points.** A point plus log spread cannot be combined
  exactly (its §6).
- **Carry a second function for the readout.** Preserving the
  average-then-divide readout across readings needs E[α | m] as well (its §7).

[`sequential_update.py`](sequential_update.py) implements three declared
rules for N independent readings of one mass. Each reading keeps its own
latent excitation.

- **symmetric**: prior(u) × ∏ L_i(u). Each reading's nuisance is its
  standardized radius r, which treats force and acceleration alike. The
  likelihood is L_i ∝ J₁(h_i(θ_i)), the radial kernel of
  `comparison_posterior`. The prior enters once; the default is the first
  reading's half-Cauchy.
- **mass_matching**: contribution 11 eq. (18), the product of single-reading
  mass densities (acceleration magnitude as nuisance).
- **log_matching**: eq. (20), the product of log-mass densities.

All three return the flat law at N = 1. With a common s they differ by exact
factors, verified in the tests:

$$
p_{(18)}=p_{\rm sym}\,\cos^{2(N-1)}\theta,\qquad
p_{(20)}=p_{\rm sym}\,\Bigl(\tfrac{\sin 2\theta}{2}\Bigr)^{N-1}\quad(\text{up to normalization}).
$$

Each single-reading law contains a half-Cauchy factor. The matching rules
count it once per reading, so the combined prior sharpens with N whatever the
data say. cos²θ pulls towards m → 0, and sin 2θ pulls towards m = s. The
2 September "m vs μ" asymmetry, a ratio m^(2N−2) between acceleration-nuisance
and force-nuisance joint references, is the same issue. Eq. (20) is exactly
the geometric mean of those two references; the symmetric rule carries no such
factor.

**Repetition study.** Each saved cell of the fresh mass ladder is split into
200 consecutive series of 50 readings, treated as 50 independent readings of
one mass. The readout is eq. (28) of contribution 11 applied to each rule's
joint law. Coverage MCSEs are 1.3–3.1 pp (95%) and 2.9–3.5 pp (80%).
At N = 50:

| Cell | symmetric: readout/m, 95% cov, 80% cov | eq. (18) | eq. (20) |
|---|---|---|---|
| m = s, acc 1, force 1 | 1.06, 80%, 60% | 0.34, 0%, 0% | 1.02, 100%, 88% |
| m = s, acc 3, force 3 | 1.01, 95%, 74% | 0.83, 15%, 3% | 1.01, 97%, 80% |
| m = 4s, acc 1, force 4 | 1.03, 94%, 78% | 0.48, 0%, 0% | 0.64, 0%, 0% |
| m = 4s, acc 3, force 12 | 1.01, 96%, 78% | 0.85, 1%, 0% | 0.92, 54%, 26% |
| m = s/4, acc 4, force 1 | 1.02, 92%, 74% | 0.91, 84%, 72% | 1.61, 0%, 0% |
| m = s/4, acc 1, force 0.25 | 1.18, 74%, 55% | 0.41, 52%, 26% | 3.42, 0%, 0% |

Across all 15 excited cells the symmetric readout at N = 50 is 1.00–1.06 of
the true mass when the force SNR is at least 1, and 1.03–1.18 in the three
cells with force SNR below 1. Its 95% coverage is 74–96% and its 80% coverage 55–79%. Observed 80%
coverage is below 80% in every excited cell; in the strongest cells the
shortfall is within about one MCSE. In the weakest cells coverage falls as N
grows: in the lightest cell, 95% coverage goes from 0.95 at N = 1 to 0.73 at
N = 50. Eq. (18) is below the true
mass in every cell at N = 50 (0.34–0.96 of truth). Eq. (20) is correct only when the
true mass is σ_F/σ_a.

**True zero excitation** (the data cannot depend on mass). After 50 readings:

- eq. (18) gives a posterior median of 0.07 s (readout 0.08 s) with a typical
  95% interval of [0.003, 0.23] s;
- eq. (20) reports m ≈ s, interval [0.76, 1.3] s;
- the symmetric posterior stays wide (typical interval [0.17, 10] s) but not
  unchanged. Its series medians scatter (interquartile range of log medians
  3.1), so random noise still produces some spurious concentration at random
  places. With *exactly zero* readings it is exactly unchanged, as tested.

Contribution 11 presented eq. (18) as exact for its declared joint reference
(17) and noted that (17) is not uniquely forced. That remains true. This
assessment adds that, *used as an update rule*, it counts a per-reading prior
factor N times and converges to the wrong mass. Under the symmetric rule, the
anchoring at σ_F/σ_a from §6 lives only in the first prior and washes out as
readings accumulate: the heavy and light cells above converge to the truth.

## 8. What this changes, and what is open

The readout E[f]/E[α] remains the composition-forced part. The law now has two
explicitly declared choices, each with measured consequences:

- **the prior factor on mass.** The default is the flat law's half-Cauchy at
  σ_F/σ_a; a stated prior may temper it. Fitted tilts are parameters whose
  values carry their fitting distribution.
- **the per-reading nuisance used when readings are combined.** The symmetric
  rule is the one that converges.

Open:

1. **A default prior that composes.** A prior anchored at the instrument ratio
   gives every object the same centre, whereas the parts of a lumped object and
   the whole should not share one. An object-anchored or scale-free default
   might be derivable from composition applied to the law rather than the
   readout. Not attempted here.
2. **The subunit regime.** The reciprocal-root summary leads the ratio of
   means there by up to 23 pp. The flat-law median and geometric points become
   the minimax rules. A negative tilt recovers part of the loss. Whether a
   composition-respecting construction recovers all of it is open.
3. **Symmetric-rule calibration.** 80% intervals undercover (55–79%). Some
   concentration at random locations remains under true null excitation.
4. **Scope.** Readings with different s, correlated or calibration-uncertain
   noise, non-Gaussian noise and real data are untested by these modules.

## 9. What changed from the exploratory numbers

Two conversational results were revised when rerun through the modules:

- The first off-grid check of λ\* compared against eight of the eleven rules.
  It omitted the flat law's median, geometric and reciprocal-root points, and
  reported below-range regrets of about 5–8 pp. Against all eleven they are
  15–26 pp, because the reciprocal-root point wins there.
- The first shifted-grid fit used weaker-channel SNRs from 1 to 8 only and
  found positive tilts centred towards the mass distribution. The declared
  study includes weaker-channel SNR 0.5, where the fit turns negative. Both
  are reported in §4.

Every other exploratory result reproduced: the dominance band, the switch
rules, λ\* = 0.35, the log-normal table, the shrinkage table and the
sequential comparison.

## Code, studies and reproduction

| File | Role |
|---|---|
| [`weighted_posterior.py`](weighted_posterior.py) | Flat law times a declared mass factor (`SechTilt`, `LogNormalFactor`, `LogNormalPrior`); points and full continuous summaries. It reuses `comparison_posterior`'s panels, kernels and summaries, and returns them unchanged when no factor is given. |
| [`sequential_update.py`](sequential_update.py) | Per-reading likelihoods; the symmetric, eq. (18) and eq. (20) rules; the eq. (28) readout; series analysis of saved cells. |
| [`archive_rescoring.py`](archive_rescoring.py) | Shared loading, failure indicators, paired differences and Bonferroni bounds, built on `comparison_scores` and `trial_archive`. |
| [`paired_dominance.py`](paired_dominance.py) | Cells where a focal rule beats every other saved rule. |
| [`switch_rules.py`](switch_rules.py) | Two-branch pretest rules and their oracle. |
| [`prior_weight_study.py`](prior_weight_study.py) | `sweep`, `fit` and `summarize` for prior-factor grids, with replicate halves and an SNR filter. |
| [`test_weighted_posterior.py`](test_weighted_posterior.py), [`test_sequential_update.py`](test_sequential_update.py), [`test_prior_update_analyses.py`](test_prior_update_analyses.py) | 18 tests (below). |
| [`studies/prior_update_validation_grid.json`](studies/prior_update_validation_grid.json), [`studies/prior_update_offgrid_shifted.json`](studies/prior_update_offgrid_shifted.json), [`studies/prior_update_mass_ladder.json`](studies/prior_update_mass_ladder.json) | ADEMP protocols for the three fresh studies (R = 10,000 per cell, new seeds, existing engine unchanged). |

**Independent checks built into the tests:**

- no factor reproduces `infer_batch`'s flat point and every distribution
  summary exactly;
- the sech tilt reproduces the README §4 closed-form null law
  (λ = −4, 0, 4, 16, 64);
- an exact log-normal prior with zero readings returns that normal law;
- the sequential rules reproduce contribution 11's §7 values (1.993026469631,
  2.373589604907, 2.008902161134 and the scaled 1.997267503687);
- with two zero readings they give its readouts s/2 and s;
- every rule reduces to the flat law at N = 1;
- the module's own rule combination satisfies the two angle identities;
- the symmetric rule is unchanged by exactly-zero readings;
- quadrature refinement agrees;
- switch limits (τ = 0 and τ = ∞) and the oracle agree with direct rescoring.

Saved outputs are in `.tools/` (scratch, not committed):

- fresh studies: `.tools/mass_estimator_studies/prior-update-validation-20260923`,
  `prior-update-offgrid-shifted-20260923` and `prior-update-mass-ladder-20260923`;
- analyses: `.tools/mass_prior_update_20260923/` (subfolders `dominance`,
  `switch`, `lambda`, `centre`, `lognormal`, `sequential`).

All three studies were run in scratch space on the linked computer, copied here
and re-verified with `study_runner.py --resume`, which checks every saved
file's hash. From the workspace root:

```powershell
$S = "03_trace_and_evaluator/mass_estimator/investigation"
$O = ".tools/mass_prior_update_20260923"
$OR = ".tools/mass_estimator_studies/operating-range-20260921/experiments/effective_snr_grid/attempt-001"
python -B -m unittest discover -s $S -p "test_*update*.py"
python -B -m unittest discover -s $S -p "test_weighted_posterior.py"
python -B $S/study_runner.py --study $S/studies/prior_update_validation_grid.json --output .tools/mass_estimator_studies/prior-update-validation-20260923
python -B $S/paired_dominance.py --run $OR --exclude tube_joint --output $O/dominance/operating_range_excluding_tube.json
python -B $S/switch_rules.py --run $OR --output $O/switch/operating_range.json
python -B $S/prior_weight_study.py sweep --run $OR --family sech --lams -2 -1 -0.5 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 1 1.5 2 3 4 --output $O/lambda/fit_sweep_operating_range.json
python -B $S/prior_weight_study.py fit --sweep $O/lambda/fit_sweep_operating_range.json --output $O/lambda/fit_selection.json
python -B $S/sequential_update.py --run <mass-ladder attempt> --cells 00_mass_excitation 01_mass_excitation --output $O/sequential/example.json
```

The validation, off-grid, shifted and log-normal analyses use the same
commands with the other study runs, `--cells`, `--half first|second`,
`--centres`, `--family lognormal_factor|lognormal_prior --widths …
--centre-factors …` and `summarize --min-snr 1`. Sweep, dominance, switch and
sequential outputs record their input runs, manifests, source hashes and
runtime. Fit and summarize outputs record their input sweep files; hashes of
those sweeps and of the sources were added to them after this study's own fit
and summary files were written. A full sweep over a 64-cell run takes
about a minute; the sequential analysis about 14 s per cell.

## Provenance and scope

- **Workspace.** Before any change, every file under
  `03_trace_and_evaluator/mass_estimator/` matched `canvalidk/VD-Newton@0525c5a`
  (the user's push of 23 September), apart from line endings.
- **Managed context.** It was read during the conversation from the user's
  local VD-docs working copy rather than a resolved `main` commit: the
  mass-estimator discussion README and steering record, the equation README,
  the introduction candidate, the originality assessment, and the Claude
  project notes on the forcing results. The protocols carry forward the
  laboratory's snapshot `8782a293297330f7a354d87ebe62792018fa8866`, as earlier
  studies did. No managed source was edited.
- **Scratch files.** Two small files in `.tools/mass_prior_update_20260923/`
  record tooling rather than results: `fs_probe.json`, a filesystem check, and
  `sequential/timing_probe_cell07.json`, which duplicates part of
  `series_cells_00_08.json`. Both can be deleted.
- **Code changed after some outputs.** After some outputs were written,
  `prior_weight_study.py` gained parameter grouping (for speed), the
  `--min-snr` filter and provenance in fit and summarize outputs. The
  docstring of `switch_rules.py` was also corrected. Outputs record the source
  hashes current when they were written. An independent re-run of `sweep` on
  four cells with the final code reproduced the saved values exactly.
- **Review.** An independent review of this record against the saved outputs
  found the tables consistent. It led to the corrections of interpretation now
  in the summary and §§1, 2, 5 and 7: the oracle is not a ceiling; eq. (20)
  settles between the truth and σ_F/σ_a; the reciprocal-root point's 23 pp is
  its lead over the ratio of means; and the minimax comparison was added.
- **Not committed.** Nothing here is committed or pushed; that is the user's
  or the VD-Newton workspace's step.

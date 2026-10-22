# The ours26 sweep: every flat record against the lemma list

**2026-09-27.** Record, written by Claude at the user's request: "Use the lemma
list as the premise table, with ours26 as the equation. Check each flat record
against it: confirmed, moved or broken." It is step 1 of the order of work in
this folder's `README.md`. It edits nothing and binds nothing. The corrections it
finds in the lemma record are listed in §4. They have not been made, because
that record changes by correction, and that is his call.

**What was checked.** The 25 records at the discussion's flat level: 24 files
and the 21 September record folder (`mass_estimator_equation/`). The lemma record
is the table, so it is not one of them. The premise table is
`../mass_estimator_lemmas_2026-09-27.md` (L1–L14, T1–T20). The equation is
ours26 (`ours26_per_axis_noise_2026-09-26.md`, with its §4 as corrected by
`ours26_starting_law_2026-09-27.md`).

**Method, and its limit.** Five readers took four or five records each and
classified every exact, algebraic, numeric or structural claim, about 500 rows
in all.
Each reader first wrote its own ours26 implementation and checked it against
numbers already on record. Between them the readers reproduced:

- ours26 §5's zero-reading readouts 1.759, 1.730 and 1.720;
- the 3D worked example, 4.9937539062;
- the audit's 1D/3D one-zero table;
- the prior table's 1.33213.

Every **broken** call and every new **[one pair]** call was checked, by algebra
or by a short computation. Claude read all five ledgers and rechecked the
findings against the lemma record (§4) and a sample of numbers (§8). Claude also
classified `type_theoretic_analysis_newton_II_mass_recovery_d2.md` directly: it
is hard-linked on disk, so it could not be staged to a reader.

The ledgers themselves are not filed. This record keeps every broken and moved
item, and the confirmed items that carry weight.

---

## 0. The verdicts

- **confirmed.** The claim holds under ours26, and what it rests on is in the
  table.
- **moved.** The claim still holds, but its standing changed in one of these ways:
  - *[one pair]*: it holds only when every component shares one noise pair,
    which is ours25b's ground;
  - *limit*: it holds only as a limit;
  - *re-based*: it now rests on a different L or T than the record says;
  - *older equation*: it is true of ours21 or ours24 as stated;
  - *open*: something settled became open, or the reverse;
  - *numbers*: the conclusion stands but a stated number changes.
- **broken.** The claim is false under ours26, or it rests on something the
  table rejects.

A record's verdict is the worst status among the claims its conclusion needs. A
broken aside that nothing rests on does not break the record. It is still listed.

**In one paragraph.** 2 records are confirmed, 13 moved and 10 broken.
Nothing found breaks ours26 or the chain. Every broken item is one of three
things:

- one of the three retired moves: the rejection gate, picking the best
  candidate, or combining separately computed answers;
- a sentence in an older record that later records already overtook;
- one public-claim wording.

What moves most is **[one pair]**. Almost every exact result that names $s$
holds only with one noise pair: zero readings give $s$, equal energies, the
four-number state, one von Mises vector, $\kappa_0$, Liseo's identity, the width
factor and frame invariance. ours26's replacements are per pair, or open.

The sweep also finds **six errors and four gaps in the lemma record** (§4). The
one that matters most: T14, "zero readings give $s$", does not follow from the
listed premises. L3 now lets a declared starting law carry its own mass scale.
Take one pair and a starting law centred at $c=3s$. It passes L2–L4, T5, T7 and
T20, and returns $2.0950\,s$ at zero readings. T14 needs back the premise the
old P2 carried: with no knowledge of the object, the starting law is built from
the instruments alone.

## 1. At a glance

| record | verdict | what breaks or moves | what survives |
|---|---|---|---|
| `vector_division_inertial_mass_measurement_d1` (26 Jul) | **broken** | the boxed argmin (picks the best candidate); the common-mass test (gate); averaging separate 1D answers; "$g(1)=1$" on measured inputs | the vector-division framing; OLS facts; "one joint fit" |
| `type_theoretic_analysis_newton_II_mass_recovery_d2` (29 Jul) | moved | §7's policy "selecting an aligned pair" → weighting all candidates; the policy slot shrinks to $\pi$ and L13 | §§1–6 whole (T1) |
| `why_reconciliation_cannot_be_upstream_of_E10` (1 Aug) | moved | "estimated codirectional pair" → the pair of mean sizes; non-uniqueness shrinks to $\pi$ + L13; **ruling owed** on "outside the VD's authority" | the ratified dilemma |
| `why_the_estimator_is_the_definition` (1 Aug) | **broken** | Result 1 (the positive root is the misfit *maximum* when $\tilde F\cdot\tilde a<0$, and absent when $=0$) | its conclusions, re-based on L2, T13, T14 |
| `median_based_estimator_inertial_mass` (27 Aug) | **broken** | the median as the point (fails lumping); profile vs marginal; the physical-units latent weight; "0/0 has no universal interpretation"; "mean ≈ $m_0$"; the compatibility verdict | integrate, don't profile; prior sensitivity at weak signal |
| `disagreement_rule_mass_estimator_trace` (29 Aug) | **broken** | §7.5's incompatibility return (gate); §8's dispatch order | the trichotomy and the repair rule (§§2–6, §7.1–7.4) |
| `inertial_mass_estimator_full_working_spec_d1` (30 Aug) | **broken** | the gate, throughout; §23 (median needs more inputs) | the readout, now forced (T4); the exact branch |
| `what_forces_the_chain_to_the_estimator` (2 Sep) | **broken** | the thesis clause "an absolute compatibility test"; S9's statistic; S11's $\chi^2$ reading | the necessity chain S1–S8, S12 |
| `what_forces_the_ratio_of_means_readout` (2 Sep) | **broken** | §6's pooled first moments; §7's gate calibration | §§1–4: reciprocity selects nothing, composition does |
| `mass_estimation_zero_zero_two_cases` (5 Sep) | **broken** | Case I ("mass undefined" on measured input); Case II §8 (gate) | Case II's thesis, value [one pair] |
| `mass_estimator_originality_research_2026-09-05` | moved | ours21's measure; zero values [one pair]; the composition "sketch" is now T3/T4; §10.2 (admission) broken | the literature |
| `mass_estimator_additional_derivations_2026-09-05` | moved | $\kappa_0$ [one pair]; the retained gate (broken aside) | the loss, the aggregation, Hellinger: exact for any law |
| `mass_readout_composition_uniqueness_2026-09-06` | **confirmed** | §6.5 "gate validation" (broken aside) | the theorem, noise-free: C1–C5 = L5, L14, L6, L7 |
| `mass_estimator_law_assumed_correct_steering_2026-09-07` | **confirmed** | "full declared covariance" is open (§3.4) | L1's source |
| `mass_estimator_equation/` (21 Sep, ours21) | moved | measure, $\Sigma$ type, worked-example interval (older equation) | the readout form; no threshold |
| `mass_estimator_regimes_priors_and_updating_2026-09-23` (contribution 12) | moved | all numbers ours21; λ and centre [one pair]; §7's per-reading $s_i$ default broken (T20) | T4's corroboration; each matching rule counts the prior once per reading |
| `mass_estimator_default_weighting_direction_2026-09-23` | moved | its opening argument (instrument-free ⇒ $(f\alpha)^k$) broken | its conclusion: reference each reading's own instrument |
| `old_information_in_the_prior_slot_2026-09-24` | **broken** | §2.2 additive weight (known); §2.3 "not a double-counted prior" | §2.1 old × new, for any noise structure |
| `vonmises_state_2026-09-24` | moved | one vector → one per pair; TLS peak [one pair]; "does not carry $A$" broken (reversed) | exact per pair, now in every $d$ |
| `excitation_weighting_premise_2026-09-24` | moved | re-based on L8–L10 + T7; §2 step 2 broken for $N\ge2$ | the flat-in-$\mathbf v$ weighting |
| `axis_splitting_and_the_24_september_change_2026-09-25` | moved | "axis splitting forces the mass factor" re-based on T7; width [one pair] | L8 and T6, per-axis included |
| `ours25b_the_many_reading_readout_2026-09-25` | moved | collapse, frame, zero readings, width [one pair] | T9's proof, unchanged per-axis |
| `what_ours24_and_ours25b_change_in_earlier_results_2026-09-26` | moved | nearly every exact result [one pair] | the audit's breaks and moves |
| `the_combination_equation_is_exact_2026-09-26` | moved | four numbers → four per pair, plus the pair's SDs | additivity (T12) |
| `credit_and_originality_of_the_mass_estimator_2026-09-26` | **broken** (as worded) | §5's public claim omits C3 and most premises | the attributions, with [one pair] on the prior |

## 2. Broken, by cause

### 2.1 The rejection gate (L1: Newton II assumed, not tested; ruling of 7 Sept)

- **Where it appears:**
  - `vector_division` §4 ("a test of whether one common scalar mass describes
    all trials") and §6 ("a failed fit says…");
  - `median` §15 and §19 ("model compatibility = poor");
  - `disagreement_rule` §7.5 and §8's `PoorModelCompatibility`;
  - the spec: §§1, 7 (the profile's use), 13–15, 17, 19–20, 22, 25–27, 29,
    31, 32.10–11, 32.13, 32.15, 34–35;
  - `chain` thesis, S9 and S11;
  - `ratio_of_means` §6 (the Level-3 test), §7 and §9;
  - `zero_zero` I.6 and II.8;
  - `originality` §10.2;
  - `additional_derivations` (scope and §4);
  - `composition` §6.5.
- **What ours26 does instead:** it returns a mass for every measured input
  (L2, T13). Precise anti-aligned $(5,0)/(-2,0)$ at $\sigma=0.01$ gives about
  $1.0\times10^5$, with the law piled at the heavy end. The equal-energy case
  gives exactly $s$.
- **What survives:** the misfit, as the law's exponent and as description kept
  with the result. Where that belongs in the table is gap G2 (§4).

### 2.2 Picking the best candidate (d2 step 4; T4 averages)

- **Where it appears:**
  - `vector_division` §4, §8 and §9 (the boxed argmin);
  - `why_the_estimator` Result 1;
  - `disagreement_rule`'s single `latent_pair`.
- **Failures, checked:**
  - $\tilde F=(5,0)$, $\tilde a=0$: the misfit falls $12.5\to0.248\to0.0025$
    at $m=1,10,100$, so no minimiser exists. ours26 gives 6.137 in 1D and
    6.130 in 2D.
  - Anti-aligned $(5,0)/(-2,0)$, $\sigma=1$: the stationary roots are $-2.5$
    and $0.4$, with $Q(0.4)=29>Q(0^+)=25>Q(\infty)=4$. So the positive root is
    the *worst* positive mass.
  - When $\tilde F\cdot\tilde a=0$ there is no positive root at all.

### 2.3 Combining separately computed answers instead of one law (T8–T10)

- **`vector_division` §2.2, orthogonal 1D trials averaged:**
  - $F=(3,0.4,-1)$, $a=(1.5,0.1,0.3)$, $\sigma=1$;
  - the separate answers 1.808, 1.049 and 1.330 average to 1.396;
  - the joint ours26 value is 1.983.
- **`ratio_of_means` §6, $\sum_i\bar f_i/\sum_i\bar\alpha_i$ from separately
  fitted readings:** each reading gets its own $m$, against the one shared
  mass. The previous audit already flagged it as superseded.

### 2.4 The rest, one by one

| record § | claim | why it breaks | check |
|---|---|---|---|
| `vector_division` §3 | on *measured* aligned inputs, all sensible methods return $\lvert F\rvert/\lvert a\rvert$ | ours26 shrinks: $(2,0)/(1,0)$ at $\sigma=1$ gives 1.5204, not 2; $(20,0)/(4,0)$ gives 4.99983. It holds only as T17's limit. True of Deming, false of ours26 | numeric (reader and Claude) |
| `median` §1, §7, §22 | the posterior median is the reported point | fails lumping. For three equally likely states $(X_1,X_2)\in\{(1,1),(2,5),(3,2)\}$ with $Y\equiv1$: the medians add to 4, the median of the sum is 5, and the means add exactly. The median survives as a law summary | exact |
| `median` §§2, 4–5, 20 | integrating differs from profiling | with the flat push in noise units (T6), marginal $=\sqrt{2\pi}\times$ profile at every $m$. The contrast is empty for ours26; what separates the August estimators is $\pi$ and the readout | sympy |
| `median` §5 | a latent weight on the physical $\mathbf a^\ast$, independent of $m$ | leaves a factor $\propto\cos\theta$ per axis in every reading, which T7 bans | sympy |
| `median` §9.1, §16 | with strong data, mean ≈ $m_0$ | under the half-Cauchy, $E_p[m]=\infty$ for every finite data set (T19) | algebra |
| `median` §12 | a 0/0 reading has no universal interpretation | the push weight is forced (T6), and a zero reading has likelihood exactly 1 (T20): it is "no evidence" | sympy |
| `median` §13 | $\mathbf a^\ast\sim N(0,\tau^2)$ | the formula is right, but the model is not: $\tau$ is a built-in scale (L3, L10), it is a per-reading mass factor (T7), and it is the physical-units branch of the units fork | sympy |
| `disagreement_rule` §8 | dispatch on the compatibility set first | a measured aligned pair returns $\lvert\tilde F\rvert/\lvert\tilde a\rvert$ (ours26: 1.520), and a measured (0,0) returns "undefined" (against L2). Fix: exact vs measured first, as d2 step 1 does | numeric |
| `spec` §23 | the median needs stronger inputs than the ratio of means | both are summaries of the same law, with the same inputs. The ratio of means wins by T4 | algebra |
| `zero_zero` I.1, I.6, I.7 | a measured input may return "mass not identified" | L2; T13 gives a finite positive answer for every proper $\pi$ | algebra |
| `regimes` (contribution 12) §7 | default prior = the first reading's half-Cauchy, with per-reading $s_i$ | this is T20's order-dependence clause. None of its numbers is affected (one-$s$ studies) | algebra (T20) |
| `default_weighting` origin | an instrument-free, unit-respecting, symmetric weighting must be $(f\alpha)^k$ | T6 references each component's own instrument; T7 excludes the family | algebra |
| `old_information` §2.2 | $A_{1:N+1}=A_{1:N}+E[\alpha_{N+1}\mid m]$ | already found by the audit (§2). Under ours26 the exact state is $\pi$ plus $(C_g,T_g,M_g)$ per pair. Check: 2.8946 against 4.0502 | numeric |
| `old_information` §2.3 | eq. 18's drift comes from its nuisance split, "not a double-counted prior" | eq. 18 $=p_{\rm sym}\,\pi(m)^{N-1}$ in every $d$: the prior counted $N$ times as a density in $m$. The nuisance-split identity holds only in 2D, where both descriptions are true. So the "not" fails. Contribution 12's wording ("counts it once per reading") is the one that survives | algebra |
| `vonmises` §2 | the vector does not carry $A(m)$ | reversed by the audit. Under ours26, $A$ is exact from $(C_g,T_g,M_g,s_g,\sigma_{a,g})$: quadrature against Monte Carlo 2.8946/2.8944 at $m=2$ | numeric |
| `excitation` §2 step 2 | the flat-in-$F$ vs flat-in-$\alpha$ factor "belongs to the mass prior", so the choice doesn't matter | true for one reading. For $N\ge2$ it is a per-reading mass factor, which T7 retires. The record's own §3 uses the right weighting | algebra, numeric |
| `credit` ¶ and §5 | the readout "is selected exactly by the composition requirements", from Newton II, $m>0$, declared uncertainty and two consistency requirements | $R_h=E[hf]/E[h\alpha]$ with $h=1+(\mathbf u\cdot\mathbf v)^2$ passes lumping, common force, calibration, units and the swap, and is not $E[f]/E[\alpha]$. Only C3 (L6) excludes it. The premise list also omits L4, L7 and L9–L12. Fix: name C3 and the premises | algebra |

## 3. Moved, by pattern

### 3.1 [one pair]

Every row holds exactly on ours25b ground (one pair). The ours26 column is what
replaces it. The numbers use axis ratios 5, 2 and 0.8 with declared centre 2,
unless a row says otherwise.

| result | where it is stated | under ours26 | check |
|---|---|---|---|
| zero readings give $s$ | d2 step 2; spec §21; `zero_zero`; `originality` §3–4; `why_the_estimator` §6; `ours25b` §4; `combination` §1; `credit` #7; audit §5 | the law is $\pi$, whose median is $c$. The readout is 1.759, 1.730 and 1.720 at $N=1,5,50$, with limit 1.7186 | numeric (three readers) |
| the padding limit $s\,E[\sin\theta]/E[\cos\theta]$ | `ours25b` §7 | $E_\pi\big[m\sqrt{\textstyle\sum_gn_g\sigma_{a,g}^2\cos^2\theta_g}\big]/E_\pi\big[\sqrt{\cdot}\big]$, with $n_g$ the share of components in pair $g$ | numeric |
| equal energies ($P=Q$) give readout = median = $s$ | audit §6.1; d2 §4; `why_the_estimator` Result 2 | even $P_g=Q_g$ in every pair gives neither $c$ nor any $s_g$: 2.8086 (anti-aligned) and 3.7047 (aligned), for 2D ratios 5 and 2 with $c=\sqrt{10}$ | numeric |
| the four-number state $(C,T,M)$ | audit §2; `combination`; d2 §4; `credit` #12 | $(C_g,T_g,M_g)$ per pair, plus $(\sigma_{F,g},\sigma_{a,g})$ and the declared $\pi$. The rebuilt readout matches the direct one to $10^{-15}$ | numeric |
| one von Mises vector carries the law | `vonmises` §2; audit §2 | $\log L=\sum_g\mathrm{Re}[C_ge^{-2i\theta_g}]$, one term per pair in its own angle. The best single-angle fit leaves a residual of 0.116 on a range of 2.68 | numeric |
| the likelihood's peak is the TLS (Deming) slope | `vonmises` §2; `credit` §1; d2 §4 | the heteroscedastic errors-in-variables ML slope (York-type). Two pairs: ML 2.644 against TLS 1.679, 2.112 and 1.873 in the three standardisations | numeric |
| $N$ readings in $d$ dims are one reading in $Nd$; the stacked space rotates freely | `ours25b` §5.3–5.4; d2 §4 | one block per pair; rotations only within a block ($10^{-10}$, ours26 §5) | tester |
| the recording frame does not matter | `ours25b` §5.2; `vector_division` §2.2; `chain` S8 | lost. A 2D reading rotated by 0, 0.5 and 1 rad gives 2.809, 2.911 and 3.202. It may return with the full covariance form (§3.4 of the lemma record: the law is frame-covariant there; the readout is unchecked) | numeric; algebra |
| $\kappa_0=\Gamma(3/4)^2/\sqrt\pi=0.847213$ | `additional_derivations` §4; audit §4.5 | 0.8352. It also depends on $\pi$'s shape | numeric (two readers) |
| the 1D law is Liseo's eq. (14) | audit §4.1; `credit` §5; d2 §5 | needs one ratio, $c=s$ and the current half-Cauchy shape. A 1D series from mixed instruments is outside it | algebra |
| the width factor $1+d/r^{\ast2}$, $1+M/R^{\ast2}$ (T18) | `axis_splitting` §6; `ours25b` §5.5; d2 §4 | $1+\sum_c\theta_c'^2\big/\sum_c\theta_c'^2r^{\ast2}_c$ with $\theta_c'=s_c/(s_c^2+m^2)$ (§5). The one-pair formula can be far off: 1.296 against a true 2.20 | Monte Carlo |
| the profile $\lvert\tilde F-m\tilde a\rvert^2/(\sigma_F^2+m^2\sigma_a^2)$, one denominator | spec §7; `zero_zero` I.3; audit §5 | a sum of per-component Beale terms (T8) | algebra |
| λ is a prior factor centred at $s$; $(m_0,\sigma_{\ln m})$ start at $(s,\pi/2)$; weak data shrink toward $s$ | `regimes` §4, §6; `old_information` §1, §2.4, §4; audit §2 | all three become the declared $c$ (T20). The zero-reading law's interval $[0.0393c,\,25.45c]$ and its log SD $\pi/2$ survive | numeric |
| one reading: every version's point lands on $s\tan\theta_{\rm mode}$ | audit §4.4 | with per-axis noise one reading has one angle per axis and no single mode angle | algebra |

### 3.2 Re-based

- **The swap is now two things:** L4 (a premise on $\pi$) and T5 (a theorem of
  the readout). Records that use it as one requirement:
  - spec §§11, 12, 32.9;
  - audit §5;
  - the steering ruling;
  - `ours25b` §1;
  - `credit` #9.

  "Any two of lumping, common force and swap" is lopsided (lemma §4.1). The
  claim that the reciprocal-invariance survey "selects nothing" holds for the
  survey's candidates. T5 does exclude $h_1\ne h_2$ direction weightings.
- **C3 left out of composition summaries.** In `axis_splitting` §1, `ours25b`
  §1, `credit` #9 and `chain` S10, the composition laws alone are said to force
  the readout. T4 needs L6 and L7 as well. In `credit` §5 this becomes a
  public claim (§2.4).
- **"$\rho_\Pi$ is policy" is obsolete.** Records: spec §§9 and 33, `chain` §4,
  `ratio_of_means` §5, `originality` §1 and the 21 Sept README. The push weight
  is a theorem (T6). What stays declared is $\pi$ and L13. L9's units and L10's
  $\omega$ stay open.
- **"No prior on $m$" was never true.** Spec §11 and `ratio_of_means` §2 say the
  magnitude measure carries no prior. In noise units,
  $df\,d\alpha\,d\Omega=\sigma_F\sigma_a\,r\,dr\,d\theta\,d\Omega$. Its
  $\theta$-marginal is uniform, so the half-Cauchy at $s$ was always there, and
  these claims rest on L4. In 2D this measure equals T6's, so every 2D
  one-reading number in the spec and in `ratio_of_means` is a live ours26
  number: all seven of its cases reproduce.
- **Totality comes from T13**, which needs a proper $\pi$, not from $\sigma>0$
  alone. Records: spec §10, `zero_zero` II.5, `additional_derivations` (scope),
  the 21 Sept README.
- **T6's noise units come from T7.** In `axis_splitting` §4(c), "axis splitting
  forces the mass factor" is really the combination rule. L8 alone allows
  $\prod_k\sigma_{a,k}\cos\theta_k$ per reading, and that factor splits over
  axes. The same applies to `excitation` §2.
- **T19's infinite mean is a property of $\pi$'s tail**, not of the model.
  Records: the 21 Sept README, `zero_zero`, `originality` §§0 and 5. A
  lighter-tailed $\pi$ would remove it.
- **`why_the_estimator` Result 3** now rests on L2, the type, not on
  detection-limit practice.
- **`type_theoretic` §7**: the "estimation policy" no longer selects an aligned
  pair. It weights all of them, and only $\pi$ and L13 remain declared.
- **`reconciliation`:**
  - "the positive scalar relating an estimated codirectional pair" now means
    the pair of mean sizes. $E[\mathbf F^\ast]$ and $E[\mathbf a^\ast]$ are
    not codirectional in general.
  - The non-uniqueness shrinks to $\pi$ (shape; centre with several pairs)
    plus L13.
  - **Ruling owed:** the 1 August ruling says estimation "lies outside the
    VD's authority", as a policy-neutral external competence. The readout is
    now forced by dictionary commitments (L5, and L14 via Newton III). Is the
    estimator still "named, not supplied" by E10?

### 3.3 Older equation, limits, numbers

- **ours21's measure $df\,d\alpha\,d\Omega$**, in the spec, `originality`, the
  21 Sept README, `ratio_of_means` §5 and `regimes` §4: it agrees with ours26
  in 2D for one reading. The 1D and 3D numbers in these records are ours21's.
- **The one-zero closed forms**, $\sigma_F\sqrt{2/\pi}/\lvert a\rvert$ and
  $\lvert F\rvert/(\sigma_a\sqrt{2/\pi})$ (spec §§17–18, `zero_zero`,
  `originality` §3):
  - they are limits (audit §4.3);
  - they also depend on $\pi$'s tails: they need $\pi(0^+)$ finite and nonzero,
    and $m^2\pi(m)\to$ const;
  - a lognormal $\pi$ (log-SD 1) gives $1.215/A$, $1.743/A$ and $2.241/A$ at
    $A=5,50,500$, with no limit near $\sqrt{2/\pi}$;
  - so they are conditional on the open shape.
- **"Independent half-normals at zero readings"**, in `zero_zero` I.4 and II.9,
  `originality` §4 and `additional_derivations` §2: 2D only.
- **Contribution 12's $m^{2N-2}$ ratio** is $m^{M-2}$ under ours26, with
  $M=Nd$. It equals $m^{2N-2}$ only in 2D, and the geometric-mean identity is
  also 2D only.
- **The 21 Sept worked example:** the point is unchanged at 4.993753906. The
  interval becomes $[4.344623,\,5.739732]$ and the log SD 0.070979. The
  README's own interval came from an unconverged grid.
- **The spec's §20 "near 1"** is exactly 1 (T15).
- **The previous audit missed one thing.** The 3D study tables in contribution
  12, `old_information` §3 and `vonmises` §4 ran on ours21's 3D law. Their law
  columns (medians, intervals, coverage) move too, not only the ratio-of-means
  columns the audit listed. On single weak 3D readings, medians move 0.4–8% and
  the 2.5% quantile up to 26%. None was rerun.
- **`why_the_estimator` §3's null spread:** ours26's is much tighter, with
  interquartile range 1.34–2.84 against the old root's 0.86–4.92 (3D,
  $s=2$). ours26's point also has a finite sampling mean there (2.34 over 600
  draws), unlike the old root. Three different "means" are in play in the
  corpus and should not be cited for one another:
  - the old root's sampling mean;
  - $E_p[m]$ (T19);
  - ours26's sampling mean.

### 3.4 Open ↔ settled

**Settled by ours26 or the table, though records list them as open:**

- several readings of one mass (21 Sept README, `originality` §10.3,
  `composition` §6.4);
- readings with different instruments (`regimes` §8.4, `vonmises` §6,
  `default_weighting`);
- carrying the readout compactly (`vonmises` §6);
- axis splitting under unequal noise (`axis_splitting` §9);
- a mass at measured 0/0 (`disagreement_rule` §10, now L2);
- the composition proof (`originality` §8: from sketch to T3/T4).

**Still open, though records state them as settled or supported:**

- **The full covariance**, in the steering ruling ("using the full declared
  covariance"), spec §6, the 21 Sept README's $\Sigma$ type and
  `vector_division` §4. ours26 is diagonal (L13), and the full form is lemma
  §3.4.
- **The half-Cauchy**, presented in several records as *the* null law. It is
  the current, declared shape.
- **One unchanged pair measured many times.** Lemma §3.4 lists it as open. The
  21 Sept README already gives the rule (vector mean, $\Sigma/n$), which is
  correct under L13. §3.4 could cite it.

## 4. What the sweep finds in the lemma record

Proposed, not made.

### 4.1 Errors

- **E1. T14 needs a premise the table does not have.** It came up from three
  readers independently. Two problems:
  - **Sources.** The argument needs the swap at the readout, which is T5, not
    only L4. And "for any estimator of this type" needs a readout-level swap
    premise, which the table no longer has, since T5 is proved for this
    readout only.
  - **Scale.** L3 now admits "the declared $\pi$" as a mass scale. So one
    pair with $\pi$ centred at $c\ne s$ satisfies:
    - L2 and L3;
    - L4 (symmetric about $c$);
    - T5 (the swap partner, centred at $1/c$, gives the reciprocal);
    - T7 and T20.

    Yet the zero-reading readout is neither $s$ nor $c$.
  - **Checked by Claude:**
    - $s=1$, $c=3$: readout 2.09496;
    - swap partner, $c=1/3$: 0.47734;
    - product 1.0000;
    - $c=s$ gives exactly $s$.
  - **Proposed fix:** split L3 in two, and scope T14 to this readout, from L2,
    L3′, L4 and T5.
    - **L3**, units covariance, as now.
    - **L3′, knowledge-free default:** with no knowledge of the object, $\pi$'s
      scale comes from the instruments' SDs alone. Knowledge of the object
      enters only as old information in the prior slot.
  - **What else rests on L3′:**
    - lemma §3.1's "with one noise pair, centred on $s$";
    - ours26 §4's derivation of the $c_p$ family;
    - `default_weighting`'s "borrow the instrument's scale".
- **E2. T6's sources lack T7.** Flat in the physical $\alpha_k$ is also a
  product over axes, isotropic at $\omega\to\infty$ and scale-free. It differs
  from flat in noise units by $\prod_k\sigma_{a,k}\cos\theta_k$ per reading,
  and only T7 (L11) excludes that factor. ours26 §2 already argues it this way.
  Numeric check: flat-$v$/Beale $=\sqrt{2\pi}$ at every $m$, while
  flat-$\alpha$/Beale $=\sigma_a\cos\theta\sqrt{2\pi}$.
- **E3. T16 is narrower than it needs to be.** A tilt $\lambda\sin2\theta_c$,
  centred at any $c$, is exactly the evidence of one reading with $P=Q$ and
  $D=2\lambda$ on an instrument of ratio $c$ (identity to $10^{-16}$). The
  equivalence holds at the level of the law only, since a real reading would
  also change $M$ and $H$ in the readout. It settles contribution 12's
  two-parameter tilt: $(\lambda,m_0)$ with $m_0\ne s$ is an unseen reading of
  an object, so it is old information for the prior slot and never the default.
- **E4. L7's "continuous at constants" needs its topology.** The 6 Sept theorem
  (C5) uses convergence in mean absolute error, and the choice matters. Take
  $X_n=c$ with probability $1-1/n$ and $c+n$ with probability $1/n$. Then
  $X_n\to c$ in probability, but $E[X_n]=c+1$ for every $n$. Read loosely, L7
  contradicts T4. Proposed wording: "continuous at constants in mean absolute
  error".
- **E5. L2's "so $\widehat m\in(0,\infty)$" belongs to T13.** With $\sigma>0$
  and an improper $\pi\propto1/m$, the law at zero readings is improper and
  there is no readout.
- **E6. T18 has a per-component form.** It is
  $1+\sum_c\theta_c'^2\big/\sum_c\theta_c'^2r^{\ast2}_c$ (§5). Only the
  shorthand $1+M/R^{\ast2}$ needs the [one pair] mark.

### 4.2 Gaps

- **G1. L1 has no addressing or persistence.** d2 P1 had "one object, one
  interval, one inertial frame", and d2 §1 had "all sharing one mass". T8 uses
  one $m$ across readings, and the breaks in §2.3 are breaks precisely because
  that mass is shared.
- **G2. No home for the kept misfit.** d2 P3's "the misfit is kept with the
  result" is needed by spec §§30 and 32.14, `chain` S9, `disagreement_rule`
  §7.2 and `reconciliation`. It is an output clause, not a premise the
  mathematics rests on.
- **G3. The inherited premises are not named.** Provenance independence, laws
  as predictions and the authority ruling come from `_note_why_mass_estimation`.
  The table should say it inherits them rather than leave them out silently.
- **G4. Two omissions in §3 and §2 of the lemma record:**
  - §3 lacks law-level composition of the default $\pi$ (contribution 12
    §8.1);
  - the prior-slot factorisation (old law × new likelihood) has no T number,
    though T20's proof uses it. It holds for any noise structure while the
    readings are independent (checked to $2\times10^{-14}$).

**Minor.** L13's "never argued" is too strong for its cross-channel clause.
`reconciliation` §§1–2 argue that the channels are determined independently.
That is not the same as independent errors (a shared timing error would
correlate them), so the clause stays declared, but it is not unargued.

## 5. New results from the sweep

All of these were checked, and none is filed elsewhere yet.

- **The width factor, per component.** The score and Hessian of
  $\ell(m)=\tfrac12\sum_cw_c^2$ at the truth give
  $J=\sum_c\theta_c'^2(1+r^{\ast2}_c)$ and $H=\sum_c\theta_c'^2r^{\ast2}_c$. So
  $\mathrm{Var}_{\rm true}/\mathrm{Var}_{\rm law}=1+\sum_c\theta_c'^2/\sum_c\theta_c'^2r^{\ast2}_c$.
  Monte Carlo (3D, fixed-length pushes):

  | $\sigma_a$ per axis, $m$ | $N$ | measured | per-component formula | $1+M/R^{\ast2}$ | 95% coverage, measured |
  |---|---|---|---|---|---|
  | 0.2, 0.5, 1.25; $m=0.5$ | 400 | 2.197 ± 0.069 | 2.210 | 1.296 | 81.4% |
  | 0.1, 1, 3; $m=2$ | 1600 | 1.605 ± 0.041 | 1.564 | 1.318 | 87.8% |

- **Tempering needs one $\omega$ per noise pair.** In the first case above, with
  oracle $\beta$:
  - one $\beta$ for all components leaves a ratio of 1.693;
  - one $\beta_g$ per pair gives 0.974, against a predicted 1.

  The algebra: with $\beta_g=R^{\ast2}_g/(M_g+R^{\ast2}_g)$,
  $J_\beta=\sum_g\beta_g^2\theta_g'^2(M_g+R^{\ast2}_g)=H_\beta$. This is a fact
  about the hint, not an adoption. It firms up the README's "one $\omega$ per
  noise pair" item.
- **Profile = marginal.** Under the flat noise-unit push, integrating out the
  push gives the profile likelihood times $\sqrt{2\pi}$, with no $\theta$. So
  ours26's law is $\pi\times$ profile likelihood, and d1's argmin is its
  interior mode.
- **The zero-reading padding limit per pair** (§3.1): 1.7186 for ratios 5, 2
  and 0.8 with $c=2$.
- **The several-pair likelihood peak** is the heteroscedastic
  errors-in-variables ML slope, not TLS in any one standardisation (§3.1).

## 6. Per record, what d3 needs

Only what §§2–3 do not already say.

- **`vector_division` d1.** Confirmed: §1, §2.1, §2.3–2.5, §3's form (as
  mathematics; ours26 is outside its hypothesis, since it uses $\sigma$), and
  §5's magnitude-only measurand.
  - The $\Sigma_a\to0$ limit is older equation and a limit: ours26 gives the
    truncated-normal mean. At $\sigma_F=1$, $\sigma_a=10^{-3}$, $F=0.5$, $a=1$
    it gives 1.0096, against OLS 0.5.
- **`type_theoretic` d2.**
  - §§1–6 are confirmed: the total operations follow from L1; the checked mass
    operation and its four result types are T1.
  - §7's `MassEstimate` becomes total (L2, T13). Its "policy selecting an
    aligned underlying pair" is re-based (§3.2).
- **`reconciliation`.** §§1–4 are confirmed, upstream of the table. §2's
  independence is L13.
- **`why_the_estimator`.**
  - The §2 quadratic is confirmed: it gives the stationary points of T8's
    exponent.
  - Result 4 (path-dependence as $\sigma\to0$) is confirmed. With several
    pairs the path runs through every $s_g$ and $c$.
  - §8.1 is confirmed. §8.2–8.3 are re-based on T13.
- **`median`.** Confirmed: §8's decision theory (it does not select the
  readout), §§10–11's direction of the one-zero shifts, §14's apparatus scale,
  and §15's large-$\sigma$ case.
  - §14's "a null answer of $\sigma_F/\sigma_a$ is an additional design
    requirement, not a theorem" is **right** for an arbitrary estimator: that
    is E1.
- **`disagreement_rule`.** Confirmed: §§1–6 and §7.1–7.4. d2 step 1 should cite
  only these.
- **The spec.** Confirmed: §§2–5, 8, 10, 16, 21 (the exact branch), 24 and
  32.2–32.7. The rest is in §2.1 above, or [one pair].
- **`chain`.** Confirmed: S1–S7, S12, §3, §4's "horn two needs *predictive*",
  §5 R1, R2 and R4, and §8.
  - §5 R3 is re-based: it needs only a finite, apparatus-scale null value,
    which is T13.
- **`ratio_of_means`.**
  - Confirmed:
    - §§1–4, with the uniqueness claim re-based on the 6 Sept theorem;
    - Table A, all seven cases reproduced;
    - §3's $\alpha^k$ analysis.
  - §3's "Fisher weight $\propto\alpha^2/\sigma_F^2$" is known-$\alpha$
    information. The efficient information is
    $\alpha^2/(\sigma_F^2+m^2\sigma_a^2)$, and the conclusion stands.
  - §6's $m^{2N-2}$ asymmetry disappears under T6. The swap fixes each
    reading's shared push and each Beale term for every $N$.
- **`zero_zero`.** Confirmed:
  - at zero readings, $Q\equiv0$ and the compatibility set is $(0,\infty)$, for
    any number of pairs;
  - the quantiles $s\tan(\pi p/2)$, [one pair] with $c$;
  - Case II's averaging.
- **`originality`.** The literature is confirmed. §5's weighted-moment example
  is correct arithmetic; as a policy it is the $(f\alpha)^k$ family, which T7
  rejects.
- **`additional_derivations`.** Confirmed:
  - §1's loss characterisation, for stacked sizes too;
  - §3's associative bookkeeping, for mixtures of one law, not reading by
    reading;
  - §4's Hellinger affinity;
  - §5's log-moments, with $c$ in place of $s$.
  - §2's "common scaling leaves the null unchanged" holds for any number of
    pairs.
- **`composition`.** Confirmed in full.
  - C1 = L5, C2 = L14, C3 = L6. C4, C5 and the domain are L7's three clauses.
  - §2 is T3 and §3 is T4, with exactly the table's sources.
  - The two-order identity and 109/92 are confirmed exactly.
- **Steering.** Confirmed. Its anti-aligned example gives exactly 1, by T15.
- **21 Sept README.** Confirmed: the readout form, "95% of the law, not of
  repeated experiments", and the exact branch.
- **Contribution 12.** Confirmed:
  - the matching rules each bring a per-reading factor, and only the symmetric
    rule converges;
  - zero readings leave the law unchanged, for any pairs;
  - a wide log-normal tends to flat in $\log m$, which is excluded (improper;
    the likelihood tends to positive constants at both ends).
  - §§1–3's switch rules are excluded by T4, and their numbers are older
    equation.
- **`old_information`.**
  - §2.1 is confirmed for any noise structure: two noise pairs, 3 old readings
    plus 1 new, spread $2\times10^{-14}$.
  - §2.4's tilt-family table is confirmed. The $1/\sqrt{1+\lambda}$ rule is
    within 2% from λ ≈ 64.
- **`vonmises`.** It is exact per pair in every $d$ (§3.1). With several pairs
  the conjugate family is products of von Mises terms in the $\theta_g$.
- **`excitation`.** The Jacobians and the Maxwell–Boltzmann reading are
  confirmed. Herschel–Maxwell is the table's own route: L8, L9, L10.
- **`axis_splitting`.**
  - §4(b), the product over axes, isotropy and no scale giving flat, is
    confirmed per-axis.
  - The §6 consistency argument is confirmed per-axis: each $\theta_c(m)$ is
    monotone, so every term peaks at the true $m$.
- **`ours25b`.** §3's proof of T9 is confirmed per-axis. It is noise-free, read
  in physical units. §4's composition checks are confirmed.
- **`combination`.**
  - Additivity is confirmed per pair.
  - Its note that states add only with shared $s$ **and** $\sigma_a$ is right:
    the readout needs each pair's $\sigma_{a,g}$.
  - Its 76.5% coverage cell is explained by T18.
- **`credit`.**
  - The fit measure (#1), the totality argument (#2) and the flat per-point
    prior (#4) are confirmed.
  - Beale's weighting (#3) extends component by component. Whether Lindley or
    Beale treat unequal ratios was not checked, because the papers are not in
    `papers/`.

## 7. What d3 must carry

- **The premise table.** It is the lemma list, with E1–E6 and G1–G4 if he
  accepts them.
- **Step 1.** Cite `disagreement_rule` §§2–6 and §7.1–7.4 only. Its §7.5 and §8
  are broken.
- **Step 2.** Scope it to this readout and to the knowledge-free default (L3′,
  E1).
- **Step 4.** Cite its own run-off argument plus T13, not `median` §2, whose
  mechanism is empty under ours26 (profile = marginal).
- **Step 5.** Name C3 and L7, with its topology, wherever the readout's forcing
  is summarised. Keep the lopsided "any two of three" table.
- **Step 7.** Add T7 as the source of "in noise units".
- **Step 8.** With several pairs the centre is declared (T20). "Centred on $s$"
  is [one pair] and needs L3′.
- **Step 10.** Use the several-pair form of the tilt identity (E3).
- **§4's results.** Mark them [one pair] and give the ours26 replacements from
  §3.1 here.
- **No gate anywhere.** The misfit is kept as description (G2).
- **§5's credit summary.** Fix it as in §2.4.

## 8. Checks and reproduction

- **Reader implementations.** Each reader's ours26 was independent: the law on
  a $\log m$ grid, the one-pair weight in closed form, and the several-pair
  weight by ours26 §6's Laplace-transform integral. Before use, each reproduced
  numbers already on record (see "Method" above).
- **Rechecked by Claude:**
  - E1's counterexample: 2.09496, 0.47734, product 1.0000, and $c=s$ giving
    exactly $s$;
  - E3's identity: exact at four masses;
  - E4's counterexample;
  - $(2,0)/(1,0)$ giving 1.52036;
  - $(5,0)/(0,0)$ in 2D giving 6.13036, the audit's 6.1304;
  - $(20,0)/(4,0)$ giving 4.99983;
  - the wording of `credit` §5, `old_information` §2.3 and `composition` C5,
    against the files.
- **Not re-derived:** rows the ledgers mark "not checked", mostly
  literature and confirmed structural claims. The previous audit's own *not
  re-derived* items stay as they were, except $m^{2N-2}$, now derived (§3.3).
- **Not run:**
  - the 3D study tables of contribution 12, `old_information` §3 and
    `vonmises` §4 under ours26;
  - contribution 11's two-instrument example (2.0089…) under ours26, which by
    T20 would need a declared $\pi$.
- **Scripts.** They were run in Claude's workspace and are not kept. They are
  numpy/scipy quadrature on the formulas of ours26 §1, with Monte Carlo for
  the width factor. The tester was not touched.

## 9. How we got here

1. On 27 September the lemma record was filed and ours26 was adopted as the
   working equation, "knowingly flawed", to be pushed through until its flaw
   shows. This folder's README set the order of work, with this sweep first.
2. The user asked for it in the words at the top.
3. Five readers each took four or five records, with the lemma record, ours26, its
   starting-law correction, the previous audit and d2 as context. Claude
   classified `type_theoretic` itself.
4. Three readers independently found the T14 problem. Claude settled it with
   the $c\ne s$ counterexample, and wrote this record.

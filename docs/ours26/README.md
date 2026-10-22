# ours26/ — pushing ours26 through

**2026-09-27.** Written by Claude on the user's instruction: "create new folder
in the discussion called ours26 and do our remaining investigation in there.
Move the existing relevant ours 26 files there. The lemma record would be
outside of that. You can create a markdown in ours26 folder that briefly
explains what we are doing here. How the next thing to do is push ours26
through. And the coverage factor and omega hypothesis points would go in here."

## What this folder is

The working folder for the rest of the ours26 investigation. It is a grouping
inside `_dscn_mass_estimator`, like `papers/` and `literature_review/`. It is
not a nested discussion.

- **Here now:**
  - `ours26_per_axis_noise_2026-09-26.md`, the record that names ours26 (moved
    in 27 Sept; it was the only file about ours26);
  - `ours26_starting_law_2026-09-27.md`, what is forced about the starting
    law, the mixed-instrument no-go, the arc-length candidate and the Jeffreys
    route (it corrects ours26 §4);
  - `ours26_full_noise_covariance_2026-09-27.md`, ours26 with correlated
    errors;
  - `ours26_learned_push_scale_2026-09-27.md`, the units fork decided by
    consistency, the tempering identity, and the coverage a learned $\omega$
    restores;
  - `ours26_learned_push_scale_in_the_tester_2026-09-27.md`, which corrects it
    (the tester's `hierarchical` rule already was the learned $\omega$), extends
    that rule to ours26, and reports the study;
  - `ours26_learned_omega_is_a_hint_2026-09-27.md`, the user's ruling that
    learned $\omega$ is a hint, not the solution;
  - `ours26_sweep_of_the_flat_records_2026-09-27.md`, step 1 below: every flat
    record checked against the lemma record (2 confirmed, 13 moved, 10
    broken), with six errors and four gaps it finds in the lemma record itself
    (proposed, not made);
  - (`ours27_full_covariance_2026-09-27.md` was filed here the evening of 27
    Sept and moved the same night to `../ours27/`, created for it.)
  - this README.
- **The spine, one level up:** `../mass_estimator_lemmas_2026-09-27.md`, the
  premises (L1–L14), what they give (T1–T20) and what is left. It stays outside
  this folder. It is an ordinary record: changes arrive as corrections or
  supersessions.

## What we are doing

**Since the evening of 27 Sept the working equation is ours27; its work is in
`../ours27/`.** Everything below written for ours26 applies to it unchanged,
since ours26 is ours27 with a diagonal covariance.

We adopt ours26 knowing it is flawed, and push it through: clean it up, rerun
it, and check every earlier result against it, until its flaw shows. ours26 is
ours25b exactly whenever the axes share one noise pair, so moving to it loses
nothing on old ground. Its known flaws at the start:

- the prior's centre is open, and by the mixed-instrument no-go (T20) it has to
  be declared before data;
- the zero-reading readout is not the centre;
- the per-pair $(C,T,M)$ carry is not built in the tester;
- one cell covered 91%.

**Order of work** (proposed 27 Sept, agreed):

1. **Sweep.** Check each flat record against the lemma record: confirmed, moved
   or broken. The output is d3 of the equation note, featuring ours26 with the
   lemma list as its premise table.
   Sweep done 27 Sept (`ours26_sweep_of_the_flat_records_2026-09-27.md`).
   Before d3: his ruling on the lemma-record corrections it proposes (§4
   there), since they change d3's premise table.
2. **A record for C3**, the "no second weighting" argument (lemma record §3.2).
3. **The learned-$\omega$ study** on ours26 (below). Done 27 Sept; ruled a
   hint, not the solution.
4. **The per-pair $(C,T,M)$ carry** in the tester.

## Coverage: two different failures

**(a) Fieller / Gleser–Hwang: accepted, not fixed.**

- Where: the weak-signal corner, as the true push goes to zero and the mass
  stops being identified.
- What: no interval that is always finite (the type's promise) can guarantee
  coverage there. Any scheme that moves the shortfall around is choosing one
  arbitrary shortfall over another. That is policy, not a forced move.
- Standing: a 95% interval holds 95% of the law; it is not a
  repeated-experiment guarantee (equation note d2).

**(b) Many readings: correctable.**

- Where: the strong-information region. It grows with $N$. At SNR 1 the cell
  average went 95.8% → 91.1% → 86.3% for $N$ = 2, 6, 21, and the worst cell was
  76.5% (`../the_combination_equation_is_exact_2026-09-26.md` §4).
- What: the law is too narrow by $\sqrt{1+M/R^{\ast2}}$ (T18), and per reading
  the factor $1+d/r^{\ast2}$ does not shrink with $N$.
- The formula predicts the number: 3D, $m=s$, SNR 1 gives $r^{\ast2}\approx2$,
  a factor 2.5 and a width 1.58 too small, so a nominal $\pm1.96\sigma$ covers
  $\pm1.24\sigma$, about 78.5%. Measured: 76.5%.
- Gleser–Hwang does not forbid good coverage here.
- **Not fixed by a multiplier.** The sandwich widening is bolted on from
  outside the chain, and the tester's version breaks axis splitting (it uses
  per-reading scores). The user's position: wait for the fix to drop out
  naturally rather than impose it. The $\omega$ hypothesis is the candidate for
  that.

**The confidence set $\{m:Q(m)\le\chi^2_{M,1-\alpha}\}$** (from the 27 Sept
completion message, not filed). Exact coverage under the declared model, since
$Q(m_{\rm true})\sim\chi^2_M$. Two problems:

- An empty set means no mass fits at level $\alpha$: rejecting Newton II, the
  gate the 7 Sept ruling removed. Descriptive at most, never a verdict.
- It shrinks like $N^{-1/4}$, not $N^{-1/2}$. $E[Q(m)]=M+\delta^2\sum_i\mathbf a^{\ast\top}_iV_m^{-1}\mathbf a^\ast_i$,
  which grows like $N\delta^2$, while $Q(m_{\rm true})$ fluctuates by
  $\sqrt{2M}\sim\sqrt N$. (Derived in conversation, not simulated.)

## The $\omega$ hypothesis

**Ruled 27 Sept (evening): a hint, not the solution**
(`ours26_learned_omega_is_a_hint_2026-09-27.md`). The hint: the many-reading
undercoverage is exactly the $\omega=\infty$ claim, so the width problem lives in
L10. Read the notes below as the investigation that led there, not as a plan to
adopt learned $\omega$.

- **L10 as used:** $\omega\to\infty$, the flat push weight. That is the
  infinite-temperature limit the width formula blames.
- **L10 as stated:** no preferred push size *before data*. The candidate is
  **$\omega$ unknown and learned from the readings**. Correction (Claude, 27
  Sept, later): a scale-free prior $d\omega/\omega$ does **not** work. The
  likelihood stays positive as $\omega\to0$, where the data say nothing about
  the mass, so the posterior is improper and the readings wash out (checked
  numerically). $\omega$ needs a proper prior, or a plug-in estimate.
- **First checks (Claude, 27 Sept, standalone simulation, not the tester):** in
  noise units the $\omega$ law is the flat law raised to the power
  $\omega^2/(1+\omega^2)$, so its peak is ours26's for every $\omega$, and a
  learned $\omega$ widens it by $1+M/R^{\ast2}$, the missing factor. Coverage
  in the 3D, $N=21$, SNR 1 cell went from 78% to 96%. Physical isotropy is
  inconsistent when the pushes are isotropic in noise units instead. Recorded
  in `ours26_learned_push_scale_2026-09-27.md`.
- **What survives:** axis splitting ($N(0,\omega^2I)$ is a product over axes and
  isotropic: Herschel–Maxwell's own answer, with the flat weight as its limit);
  L3 (see the fork below); probably the $(C,T,M)$ state, since everything stays
  Gaussian and isotropic (unchecked).
- **Prediction:** calibrated widths when the world's pushes match the model.
  The tester's worlds use pushes of fixed length in random directions, so a
  Gaussian $\omega$ is slightly misspecified there. Expect most of the coverage
  back, not all.
- It does nothing for failure (a).

## The fork ours26 forces: in which units is isotropy?

The flat weight is flat in any units, so ours26 could weight each axis in its
own noise units without choosing. A finite $\omega$ cannot: a Gaussian isotropic
in noise units is not isotropic in physical units once the axes have different
SDs.

- Isotropic in **physical** units: $\omega$ is a physical push scale with units,
  and each axis sees width $\omega/\sigma_{a,k}$ in its own noise units. Fine
  under L3 if $\omega$ is learned under a scale-free prior.
- Isotropic in **noise** units: $\omega$ is dimensionless, but the physical
  weighting then depends on the instruments.

The same fork decides the arc-length prior's Jeffreys justification (lemma
record §3.1), which needs equal pushes per axis in noise units. Settle it once,
for both.

## Studies owed here

- **Learned $\omega$ against ours26 (flat):** done 27 Sept
  (`presets/ours26_learned_omega.yaml`, fixed-length pushes; see the tester
  record). Still owed: the same grid with Gaussian pushes, and **one $\omega$ per
  noise pair**, which may calibrate the per-axis case and stay consistent in both
  worlds of the fork.
- **A repeated-measurement world:** one unchanged pair measured many times.
  The current law treats every reading as a new push.
- **The 91% cell** of the ours26 preset ($m=0.5$, SNR 1, $N=5$).

## Not in this folder

The completion message of 27 September (VD-Newton) is not filed. Its results
are recorded by thread: the starting law and the full covariance in the two
records above; the confidence set and the repeated-measurement caveat here.

## Moved out, 29 Sept (night)

- `ours26_learned_omega_is_a_hint_2026-09-27.md` moved to
  `../noise_vs_signal/` on his instruction, with the other records on the
  signal-versus-noise problem. The two records it corrects stay here.

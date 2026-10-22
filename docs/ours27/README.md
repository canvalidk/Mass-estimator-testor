# ours27/ — the ours27 investigation

**2026-09-27 (night).** Written by Claude on the user's instruction: "Make a
folder called ours27 put it there. Make a record of it."

## What this folder is

The working folder for ours27, the working equation since the evening of 27
September. Like `ours26/`, `papers/` and `literature_review/`, it is a grouping
inside `_dscn_mass_estimator`, not a nested discussion.

## Here now

- `ours27_full_covariance_2026-09-27.md` names ours27: ours26 with a declared
  full error covariance per reading. The change is forced by the property that
  mass cannot depend on the choice of axes. ours26 and ours25b are its special
  cases. The record was moved here from `ours26/` the night it was filed,
  because this folder was created for it.
- `ours27_perturbation_matrix_2026-09-27.md` bends each move of the chain and
  tests every property. Six moves are forced and three directions are free:
  the starting law's shape, C3, and the push scale ω. It also shows that
  Gaussian noise is forced by consistency, given the flat push.
- `ours27_candidate_size_and_null_readings_2026-09-28.md`: seven short
  reminders from the night of 27–28 Sept and the next morning. They cover
  candidate size (formerly "push"), repeated null readings, the flat measure as
  a delta at infinity, rotating a reading to 0/0 (known-empty against
  measured-empty), impossible data under correct σ, and per-reading vs uniform
  down-weighting. The width test of the signal-to-noise-angle rule is there
  too.
- `ours28_candidate_push_scale_as_nuisance_2026-09-28.md` (evening): **ours28,
  a candidate.** L10 ("no push scale") is not a symmetry but the assumption
  that every push is infinitely large, and T18's width factor is exactly its
  error. ours28 keeps the push scale as a per-reading nuisance averaged under
  the size-angle rule (uniform in $\tan\varphi=\omega$, the rule that gives
  the half-Cauchy). Coverage 92–98% where ours27 gives 66–87%; a zero reading
  leaves the law unchanged exactly; T11's four-number state is the price. Open
  after it: C3, and the angle rule (now one choice, with a Jeffreys check
  owed). Filed here because it is ours27's open item 3; an `ours28/` folder is
  the user's call.
- `ours28_size_measure_jeffreys_check_2026-09-28.md` (night): the Jeffreys
  check the ours28 record owed. Under ours28's own model the mass and the
  push scale are orthogonal and the reference prior on the mass is uniform in
  $\theta$ for any size prior (the half-Cauchy, without holding the push
  fixed). For the size, Jeffreys gives $\tan\varphi$, not uniform, and
  under-covers. What calibrates is the push–noise symmetry
  $\omega\leftrightarrow1/\omega$: every symmetric prior covers 94–98%,
  every asymmetric one fails. Shape open, as for the mass. Proposes L4′.
- `ours28_known_empty_vs_measured_empty_2026-09-28.md` (night): his rotation
  story. Under ours28 a 1D-declared reading is narrower than the same
  numbers embedded in 3D with two null axes (r = 4: 3.04 vs 3.65 log-width;
  equal by r = 16), because the nulls inform the reading's shared push scale.
  Under ours27 they were equal. With per-axis noise the same story finds the
  units fork: precise empty axes wrongly discount the signal axis. Candidate
  repair: one $\omega$ per noise pair per reading. This moves the several-pair
  form of $\omega_i$ to the top of the ours28 owed list.
- `ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md` (evening,
  written in a session that ran in parallel with the ours28 one): C3 derived
  from betweenness, own frames and sign conventions, grouping as bookkeeping
  and same-law-same-mass (all but grouping ruled 28 Sept); the per-axis
  survivor that needs same-law-same-mass; self-consistency coverage; the size
  measure forced to one Gaussian scale per series (Kingman) if grouping is
  bookkeeping. Then **the fork with ours28**: grouping as bookkeeping, the
  zero-reading property and a proper size measure cannot all hold; ours27,
  a shared scale and ours28 each drop one. The ruling owed is whether which
  components were measured together is bookkeeping or physics.
- **Four Codex records of 28 September** (VD-Newton side, written at the
  user's request against VD-docs@b4bfc20, so they do not cite the ours28
  records or the betweenness record; sorted here the same evening). Their
  closing "not submitted" banners are stale.
  - `ours27_axis_splitting_does_not_force_independent_candidates_2026-09-28.md`:
    representation invariance does not imply integrating each axis's
    candidate state separately; a shared size state can travel through the
    split (proper counterexample). Bears on the owed ruling of the
    betweenness record: is grouping bookkeeping or physics?
  - `ours27_full_covariance_width_and_candidate_geometry_2026-09-28.md`: T18
    in full-covariance form (score variance $J+G$ against curvature $J$), an
    exact conditional test pivot as a benchmark, identities linking candidate
    size and direction to curvature, and the Fieller limit kept separate from
    the avoidable width defect.
  - `ours27_limits_of_scalar_and_tensor_size_repairs_2026-09-28.md`: one shared
    scalar size does not calibrate widths with heterogeneous instruments
    (regular counterexample, SD ratio 1.61); a shared tensor carries through
    recodings, but its reference measure is improper for few readings and its
    rotation-invariant continuations disagree on the mass law.
  - `ours27_zero_readings_and_shared_candidate_information_2026-09-28.md`: his
    rotation story, reached independently of
    `ours28_known_empty_vs_measured_empty_2026-09-28.md`. An all-zero reading
    can leave the mass law equal to the starting law and still update the
    shared size; omission, precise conditioning and a fresh lower-dimensional
    model are different operations; the precise-null limit is not the 1D fit
    (2.0697 against 2.2064).
- **Moved to `../ours28b/`** (28 Sept, night, on his instruction, with the
  record naming ours28b): `ours27_one_scale_per_noise_pair_2026-09-28.md` (night, same session as the
  betweenness record): asks the law's width to equal the spread of its own
  centre over repeats, for any pushes (Bartlett at fixed candidates, regular
  limit). No per-reading rule can do this (an Ornstein–Uhlenbeck equation with
  only constant smooth radial solutions); shared scales do it exactly when
  they don't mix sensitivities; the precise-null limit picks classes by
  declared noise pair. Result: one Gaussian scale per noise pair, shared over
  the series; $(P,Q,D,M)$ per pair kept. On the repair-limits counterexample
  it covers 94–95% at 100–1600 repetitions where the other rules sit at
  71–88%. Price: the strong zero-reading property (the weak one holds). Open:
  full covariance, the two shapes.
- `ours28_full_covariance_form_and_declared_subspace_2026-09-28.md` (late
  night, after the five records above): ours28 in two lines, $H_i=T_i-Q_i$
  with $Q_i$ ours27's misfit and evidence ${}_1F_1(\tfrac12;\tfrac k2+1;H_i/2)$,
  $k$ the dimension of the reading's **declared candidate subspace** $N_i$, a
  new reading input. Codex's heterogeneous counterexample is a fixed
  direction left undeclared: declared, ours27 and ours28 both calibrate
  (1.01, 0.95). Withdraws the known-empty record's §3 verdict (the transverse
  widening is a solid-angle effect, correct when undeclared). Recommends on
  the fork: grouping is physics (a reading is one push), Z and properness
  kept, C3 narrow folded into L12. Calibration judged in the definition's own
  world; the two shapes are declared. Owed list reordered in its §6.
- **Pointer (28 Sept, late night):** `../ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` records
  the user's ruling that the definition borrows nothing between readings. Under
  it the ours28 line in this folder is the live one, and ours28b is withdrawn.
- `ours28_reference_implementation_and_readout_2026-09-28.md` (late night), with `ours28_reference.py`,
  `ours28_reference_checks.py` and `ours28_reference_sim.py`: a reference
  implementation of ours28 in the full-covariance form, with its **readout**
  computed for the first time (candidates $\alpha_i\mid m,\beta\sim N(\beta G_m^{-1}b,\beta G_m^{-1})$).
  Swap, units, rotation with cross-channel correlation, the zero-reading
  property and totality all pass. A declared $k=1$ subspace gives exactly the
  1D answer. Axis splitting holds only with the scale carried. Gaussian noise
  is still forced, by the zero-reading property. In ours28's own world the
  readout goes to the truth, and from 21 readings up its scatter matches the
  law's width, so the law's interval serves as the readout's.
- `ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md` (late
  night): C3's narrow form under the no-borrowing ruling. Without regrouping
  across readings, the only surviving readouts are ratios of means weighted by
  a function of how the candidate's size is shared among the readings. On
  ours28 mild weights move the point 0.1–0.8%, and the premises put no bound
  on the weight. The family is excluded if the reading structure counts, like
  $\Sigma_i$ and $N_i$, as a declared instrument fact that the score does not
  see (same-law-same-mass); C3 then follows in full. **Ruled the same night:**
  the user confirmed that reading, so C3 is derived (record §9).
- this README.

## Where things stand

The uniqueness statement the matrix supports: **given C3 and "no push scale",
ours27 is forced up to the starting law's shape.** The open list was three
items; the ours28 record (28 Sept, evening) proposes that item 3 was a mistake
in the premise table and folds item 1's shape question into it, leaving C3 and
one angle rule. As written on 27 Sept:

1. **The starting law's shape.** With several noise pairs, its centre too.
2. **C3's justification.**
3. **The push scale ω**, which is also where the many-reading width problem
   lives.

What is owed, from the two records:

- the tester build of ours27;
- a record defending axis splitting, since it is the sole support of two moves;
- the next lemma record, carrying the sweep's proposed corrections and this
  matrix's;
- a check on pairs of bends;
- widths under the bends.

## Not in this folder

- The earlier ours26 work, including the sweep of the flat records, stays in
  `../ours26/`.
- The lemma record stays at flat level.

## Moved out, 29 Sept (night)

- To `../noise_vs_signal/`, on his instruction, with the other records on
  the signal-versus-noise problem:
  `ours27_candidate_size_and_null_readings_2026-09-28.md`,
  `ours27_zero_readings_and_shared_candidate_information_2026-09-28.md` and
  `ours28_known_empty_vs_measured_empty_2026-09-28.md`. Links to them from
  records in this folder are relative and now one level off.
- `../ours28b/` is now `../noise_vs_signal/ours28b/`.

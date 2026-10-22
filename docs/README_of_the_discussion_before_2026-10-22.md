# _dscn_mass_estimator — read this before treating the folder as one argument

This is a discussion: a folder grouping records that address one another. The
FLAT files are its members — records the user chose, one by one, to file here.
They scope each other; nothing is bound (no `supersedes` / `corrects` asserted
anywhere in the folder). Content, provenance and relations live in
`catalog.yaml` — the `files:` rows and the `folders:` entry for this path. Read
the rows, not just the filenames.

## The two subfolders are not members in the same sense

`literature_review/` (22 files) and `unsanctioned/` (7 files) hold the output of
one VD-Newton Codex session, 2026-09-06 to 2026-09-08, received in the
2026-09-08 pass. That session produced records on its own initiative. In this
repo a record is something the user decides to make, and these were not decided
one by one. The user's ruling of 2026-09-15: keep them, do not formalise them,
keep them out of the flat view. That is a RECORDED EXCEPTION under
ORGANIZATION.md "Exceptions" — nothing more, no new object kind, no "log" model.

What that means if you are an AI working here:

- Do not read the subfolder files as part of the discussion's argument unless
  the user brings one up or moves it flat. Do not promote, bind, stamp a kind
  on, or reorganise them on your own initiative.
- Do not delete them either. A pass surfaces fossils, never deletes.
- Their in-file "Status: local, awaiting transfer to VD-docs, not submitted"
  banners are stale — this repo IS the transfer. Ignore the banners.
- Two files in `unsanctioned/` (behavior_test_report, uncertainty_and_reuse_
  walkthrough) were written under a gate/rejection spec that the flat steering
  record of 2026-09-07 retired the next day. They carry a banner saying so and
  were not rewritten. The banner links point at the steering record by bare
  filename; the file now sits one level up. Filenames are the identifiers;
  paths in citations decay (ORGANIZATION.md). VD-Newton's lab README links
  these files by path — that is theirs to repoint, not a reason to move
  anything back.

## Two of the session's files stayed flat, by decision

- `mass_estimator_law_assumed_correct_steering_2026-09-07.md` — the user's own
  ruling, quoted verbatim. His content, not the session's.
- `mass_readout_composition_uniqueness_2026-09-06.md` — a self-contained theorem
  (C1–C5 force E[f]/E[α]) answering a question the flat records themselves pose
  (`what_forces_the_ratio_of_means_readout`). Claude's call, 2026-09-15; cheap
  to push down if the user disagrees.

The three 2026-09-05 files (`mass_estimation_zero_zero_two_cases`,
`mass_estimator_originality_research_2026-09-05`,
`mass_estimator_additional_derivations_2026-09-05`) predate the session, were
filed on the user's instruction the day they were written, and are flat members.

## Added 2026-09-23: one record folder and one reference paper at flat level

- `mass_estimator_equation/` is a **member**: the 2026-09-21 specification of
  the current estimator (ratio of mean compatible magnitudes), received from
  VD-Newton. It is one record in its own folder. `README.md` is the record, and
  the PDF, equation images and standalone code are parts of it. It is not bound
  to `inertial_mass_estimator_full_working_spec_d1.md`; that call is the user's.
- `SIM-38-2074.pdf` is a **reference paper** (Morris, White & Crowther 2019, on
  simulation studies), placed here by the user. It is not a record of the
  argument, and not part of `literature_review/`.

## Added 2026-09-23 (evening): one flat record

- `mass_estimator_regimes_priors_and_updating_2026-09-23.md` is a **flat
  member**, filed as `kind: record` at the user's request. It is VD-Newton's
  working contribution 12: readout regimes (dominance, minimax, switch rules),
  the fitted prior weight λ and what it is (a prior factor on mass centred at
  the instrument scale s, not knowledge of the object), log-normal priors, and
  the estimator read as an update rule. Received from VD-Newton @`ecd1a1c`;
  receipt `outbox/VD-Newton/PASS_RECEIPT_2026-09-23_ecd1a1c.md`. It scopes and
  binds nothing. Its code and study links point into VD-Newton's
  `investigation/` folder, which is where the reproducible studies live.

## Added 2026-09-23 (night): `papers/`

- `papers/` holds **published papers**, collected on the user's instruction as
  standard citations. They are references, not records and not members:
  nothing in this discussion is bound to them. `SIM-38-2074.pdf` (listed above)
  moved into it the same night, alongside Wang et al. (2021) on penalized
  Fieller intervals. `papers/README.md` carries the reading list, including the
  citations not yet downloaded.
- Not to be confused with `literature_review/`, which is the unsanctioned Codex
  search *about* the literature.

## Added 2026-09-23 (night): one investigation direction

- `mass_estimator_default_weighting_direction_2026-09-23.md` is a short,
  deliberately non-conclusive **investigation direction**, filed flat at the
  user's request. It asks what the default weighting may reference when
  readings are combined. It scopes and binds nothing.

## Added 2026-09-24: one flat record, and where the studies now live

- `old_information_in_the_prior_slot_2026-09-24.md` is a **flat member**,
  written at the user's request. It shows algebraically that the exact
  combination splits into old × new, so the old readings' law sits exactly
  in the prior slot, and reads contribution 12's tilt parameters as carried
  state (m₀, σ of log m). Its evidence is a preset of the tester repo. It
  scopes and binds nothing.
- **Reproducible studies now live in the tester repo,
  `canvalidk/Mass-estimator-testor`**: one program where each study is a
  recorded preset file. Records here cite a preset and a commit. It replaces
  VD-Newton's `investigation/` scripts for new work.

## Added 2026-09-24 (later): a second flat record

- `vonmises_state_2026-09-24.md` is a **flat member**, written at the
  user's request. It shows that in 2D each reading's likelihood is exactly a
  von Mises law in the doubled angle, that contribution 12's λ is a von Mises
  concentration, and that readings sharing only the mass combine by adding one
  vector each (exact in 2D, approximate in 3D). Its peak is the
  total-least-squares axis. It carries a reading list for following this
  into the literature. It scopes the morning record and binds nothing.

## Added 2026-09-24 (night): a third flat record

- `excitation_weighting_premise_2026-09-24.md` is a **flat member**,
  written at the user's request. From his premise (force and acceleration
  are symmetric, and either can take any value), it argues that the weighting
  over the true pushes must be flat in the vector, as Maxwell–Boltzmann is
  flat in velocity. The 21 September equation is unchanged in 2D. The record
  proposes a change to it in 1D and 3D, where it currently favours large or
  small pushes. It names the argument's weak spot (accelerating frames) and a
  second route (Herschel–Maxwell). Tests are a tie. It scopes and binds
  nothing.

## Added 2026-09-25: a fourth flat record, and a premise

- `axis_splitting_and_the_24_september_change_2026-09-25.md` is a **flat
  member**, written at the user's request. It isolates the 24 September
  change as one weight on each candidate pair ($r^{d-2}$, identity in 2D),
  and shows the new law of mass is the 21 September 2D law in every
  dimension. It proposes **axis splitting** as a premise: a vector reading
  and its axes read separately must give the same law. The user judged it
  "too elementary to be wrong". It forces both the flat-in-vector weighting
  and the mass factor the previous record left open. It also proves the
  symmetric rule consistent and gives the exact width factor
  $1+d/r^{\ast2}$ for many readings. On accuracy the two equations tie; the
  gain is consistency. Evidence: four presets in the tester repo. It scopes
  and binds nothing.

## Added 2026-09-25 (evening): a fifth flat record, and a new equation

- `ours25b_the_many_reading_readout_2026-09-25.md` is a **flat member**,
  written at the user's request. It pays the debt the axis-splitting record
  left: with many readings, the readout's "size of the series" must be the
  length of all the readings' true vectors stacked into one, not the plain
  sum of their sizes. Axis splitting applied to the answer forces this,
  given one rule for joining sets of readings. It names the result
  **ours25b** (ours24's law with that readout) and presents the equation in
  full. Consequences: a series and its axis readings agree exactly, the
  recording frame does not matter, and the answer depends on the data only
  through three sums. On accuracy it ties with ours24. Evidence: two presets
  and a test file in the tester repo. It scopes the axis-splitting record and
  binds nothing.

## Added 2026-09-26: an audit record

- `what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` is a
  **flat member**, written at the user's request after he asked whether the
  new equation breaks earlier exact results. One proof broke: the additive
  readout weight of `old_information_in_the_prior_slot_2026-09-24.md` §2.2.
  Its exact replacement is a four-number state. Everything about the law
  stands, including $m_0$, $\sigma_{\ln m}$ and λ as carried state. It also
  lists what the 24 September law change moved (the 1D Liseo identity now
  matches the reference-prior posterior, not Jeffreys; some 1D/3D closed
  forms hold only as limits), a new exact result (equal channel energies
  give exactly $s$), and what the equation note's next draft must not copy.
  It scopes several records and binds nothing.

## Added 2026-09-26: a credit and originality record

- `credit_and_originality_of_the_mass_estimator_2026-09-26.md` is a **flat
  member**, written at the user's request after the 26 September literature
  check. It splits the estimator into its law and its readout, and credits
  each part. The law is Lindley and Beale's (1966 discussion of Sprent), with
  a prior shape that matches known defaults. The readout device (a ratio of
  posterior means) is known. What was not found is the selection principle:
  composition forcing the readout, and axis splitting forcing the weighting
  and the many-reading pooling. It also lists who in the project contributed
  what, and what can honestly be claimed in public. It extends the
  5 September originality record and binds nothing.

## Added 2026-09-26 (evening): ours26, per-axis noise

- `ours26_per_axis_noise_2026-09-26.md` is a **flat member**, written at the
  user's request. Every equation up to ours25b assumed one noise pair for all
  axes; ours26 drops that. Walked by what stays forced: the law becomes a sum
  of Beale terms, one per axis (axis splitting plus the combination rule); the
  readout keeps ours25b's definition, its weight now a weighted noncentral chi
  mean computed by one Laplace-transform integral; the prior's centre comes
  off the wall (the symmetries that fixed it leave a family), so it is
  declared. Also recorded: the answer depends only on per-noise-pair sums, so
  the exact combination state survives as one (C, T, M) per pair; zero
  readings no longer return the centre through the readout. Evidence: a
  test file and a preset in the tester repo (axis splitting to 5e-11; the open
  centre moves weak-data answers by a factor of up to about 1.5, strong-data
  ones by under 1%). It extends the ours25b record and binds nothing.

## Added 2026-09-27: the combination equation, built and rerun

- `the_combination_equation_is_exact_2026-09-26.md` is a **flat member**,
  written at the user's request after he asked what in the 26 September
  session was unrecorded. It delivers the audit record's owed item: ours25b's
  carried state, four real numbers per series $(C,T,M)$, built in the tester
  and tested. Combining sets of readings is adding their states; the law, the
  median, the intervals and ours25b's readout all follow from the state and
  match the rule on all readings to $10^{-14}$ over 21,600 series. The rerun of
  the old-information study under ours25b measures what the earlier carries
  lose: keeping the law but not the old weight moves the readout by up to about
  6% (cell median), the two-number fits move the median by up to about 19%.
  It also notes, incidentally, 95% intervals covering as little as 76.5% with
  21 weak readings. It extends the audit and the old-information records and
  binds nothing.

## Added 2026-09-27 (afternoon): the lemma record, and `ours26/`

- `mass_estimator_lemmas_2026-09-27.md` is a **flat member**, a record at the
  user's request ("It won't be a note but a record. It'll be dated, and we can
  imagine the changed versions go in as replacements"). It lists the premises
  the estimator rests on (L1–L14), what they give (T1–T20) and what is left,
  for ours26. An ordinary record: changes arrive as corrections or
  supersessions like any other. It is the spine of the ours26 sweep.
- `ours26/` is a **working folder** for the rest of the ours26 investigation,
  created on his instruction. It is a grouping inside this discussion, like
  `papers/`; it is not a nested discussion. It holds
  `ours26_per_axis_noise_2026-09-26.md` (moved in from flat level) and a
  README with the plan, the two coverage failures and the $\omega$ hypothesis.
  The lemma record stays outside it by his instruction.

## Added 2026-09-27 (evening): the sweep, and ours27

- `ours26/ours26_sweep_of_the_flat_records_2026-09-27.md`: every flat record
  checked against the lemma record under ours26 (2 confirmed, 13 moved, 10
  broken); it also proposes six corrections and four additions to the lemma
  record, not made.
- `ours27/ours27_full_covariance_2026-09-27.md` names **ours27**, the working
  equation from that evening: ours26 with a declared full error covariance per
  reading, forced by the property that mass cannot depend on the choice of
  axes. ours26 and ours25b are its special cases.

## Added 2026-09-27 (night): `ours27/`

- `ours27/` is a **working folder** for the ours27 investigation, created on
  his instruction, like `ours26/`. It holds the ours27 record (moved in from
  `ours26/` the same night), a README, and
  `ours27_perturbation_matrix_2026-09-27.md`: each move of the chain bent and
  tested against every property. Six moves forced, three free directions (the
  starting law's shape, C3, the push scale); Gaussian noise forced by
  consistency given the flat push.

## Added 2026-09-28 (night): `ours28b/`

- `ours28b/` is a **working folder** for ours28b, a candidate the user named
  on 28 September: ours27 with one Gaussian candidate-size scale per declared
  noise pair, shared over the series. It holds the naming record, its
  derivation (moved in from `ours27/`), a reference implementation with its
  checks, and a README. It is not the ours28 in `ours27/`, which has one
  scale per reading.

## Added 2026-09-29: the lemma record superseded

- `mass_estimator_lemmas_2026-09-29.md` **supersedes**
  `mass_estimator_lemmas_2026-09-27.md`, on his instruction ("You can just
  turn this into a superseding record"). The 27 Sept record said its changed
  versions would arrive as replacements, so this is the first `supersedes`
  asserted in this folder. It stays in place, unedited; the relation is in
  the catalog and in the new record's header.
- The new record is the current list of the properties of mass (L1–L20),
  what follows from them (T1–T29), the declared shapes, and the open items in
  the order the definition needs them. It carries the corrections the
  28 Sept records proposed to the old list, without correcting those records.
  Its first open item is new: the readout fails the zero-reading property
  (T29).

## Added 2026-09-29 (later): why mass can't be undefined

- `why_mass_cannot_be_undefined_2026-09-29.md` is a **flat member**,
  written at the user's request. It records his argument that mass is never
  "undefined" on measured inputs (the argument behind his O1 ruling on the
  lemma record), quoted verbatim, and sharpens it: an undefined output is
  contradictory, redundant or arbitrary under the definition's own
  combination rule. It flags that comparing the starting law's intervals
  would bring the contradiction back. It offers two consequences for O2:
  the all-zero point is forced to be $c$, and the best-guess repair's limit
  at zero depends on direction. It scopes the 5 Sept two-cases record and
  binds nothing.

## Added 2026-09-29 (evening): the noise-measure view

- `noise_vs_signal/the_noise_measure_view_2026-09-29.md`, written at the
  user's request, now in `noise_vs_signal/`, a **working folder** created the
  same evening on his instruction (like `ours26/`, not a nested discussion). A way of looking, not a result: the equation has no
  notion of how much of a reading is push and how much is noise, and every
  open item in the lemma record is that absence showing through. It records
  the afternoon's `==^` thread (why the law's width is not the result's
  uncertainty), the finding that a series of $N$ noise-sized readings tilts
  the law by $O(\sqrt N)$ while the data-only resolution stays the whole
  half-line, his property N, and a map of what a noise measure would touch.
  It parks five items, including a correction owed to the undefined record.
  It scopes the lemma record and binds nothing.
- `noise_vs_signal/null_series_tilt_and_property_N_2026-09-29.md` (night)
  develops the view's first finding: the tilt a null series gives the law
  has an exact form (the carried state $C$; Rayleigh under noise), ours28
  halves it and keeps the $\sqrt N$, the point moves by a factor 2.6 at
  $N=100$, the data-only resolution stays whole, and any fix at the law's
  level must share one scalar across readings. It binds nothing.
- `noise_vs_signal/c_is_the_noise_measure_at_zero_signal_2026-09-29.md`
  (night, short): $c$ is what the definition returns when the noise
  measure reads zero, but noise picks no scale and Z forbids $c$ from
  depending on the series' instruments, so $c$ is declared before data and
  L3′ is a convention for the declaration. It binds nothing.
- `noise_vs_signal/alignment_what_one_reading_can_say_about_noise_2026-09-29.md`
  (night): per reading the law sees $(T_i, C_i)$; in 1D every reading is
  aligned with some mass; in $d\ge2$ misalignment says almost nothing
  about weak signal and the residual is weak evidence *for* a push; the
  only per-reading signal evidence is size along the line; evidence that
  a series is noise is collective, $|C|/T$ against $1/\sqrt{dN}$. So the
  noise measure collides with L19 in every dimension. It binds nothing.

## Added 2026-09-29 (night, later): `noise_vs_signal/` filled

- On his instruction ("do A and B lists") the folder gained a pointer
  record, `older_hints_on_noise_and_signal_2026-09-29.md` (eleven earlier
  records, one quoted sentence each), and five moves: the learned-$\omega$
  hint ruling from `ours26/`, three candidate-size and null-reading records
  from `ours27/`, and **`ours28b/` whole**, which now sits at
  `noise_vs_signal/ours28b/`. The section above naming `ours28b/` is
  history. `noise_vs_signal/README.md` lists everything. Relative links
  to the moved files are one level off; filenames are the identifiers.

## Added 2026-09-29 (night, latest): readings are evidence, just less

- `noise_vs_signal/readings_are_evidence_just_less_2026-09-29.md`, written
  at the user's request. The user ruled that the true, errorless vector is
  called the **signal** (not a push). The record finds a defect separate
  from the null-series tilt: ours28 treats a reading misaligned far beyond
  its declared noise as strong signal, and the law collapses with unbounded
  confidence as $\sigma\to0$. His ruling ("The readings are evidence. Just
  less.") forces each reading's factor to be flattened by its off-line
  misfit, $E_i(m)^{g(Q_i)}$: $g(0)=1$ by Z, a $\kappa/Q$ tail by bounded
  confidence, the middle free. Six forms tested in 2D: all $\kappa/Q$ forms
  fix the impossible reading (coverage 8.5% to 91-100% in a rescue test),
  and forms flat inside the noise band cost 0-6% width. Inert in 1D; does
  not touch the $\sqrt N$ tilt. Axis splitting collides (parked). It binds
  nothing.

## Added 2026-09-29 (night, latest): one reading is one completeness claim

- `noise_vs_signal/readings_are_evidence_just_less_2026-09-29.md` gained a
  §10 the same night, and `mass_estimator_lemmas_2026-09-29.md` was **revised
  in place** on his instruction (a revision note sits under its Supersedes
  paragraph). His ruling: a $k$-dimensional reading declares that the signal
  has no component outside its declared subspace; that claim imposes $k-1$
  codirectionality conditions, and their violation is the reading's off-line
  misfit. A 1D reading is a projection and imposes none. So one 3D reading is
  not three 1D readings (same data, different declaration), and the strong
  form of axis splitting (L8) is false. In the lemma record: L8 rewritten,
  L21 (readings are evidence, just less) and T30 (the tempered law) added,
  O3 reopened (the flat weighting and the stacked length need their own
  support), O12 and O13 added. The 25 Sept axis-splitting record is not
  edited; its premise is withdrawn as a forcing argument in O3 (offered).

## Added 2026-09-30 (morning): scattered observations

- `noise_vs_signal/scattered_observations_2026-09-30.md`, written at the
  user's request as a list of the morning's observations. Its centre: after
  his correction ("you were wrong that the solution has this idyllic
  property"), the property identified is that **the law carries its own
  uncertainty**. Every impossibility about the $\sqrt N$ tilt is a theorem
  about the law's width and a wall only under that property, which the
  `==^` thread had already given up. So the noise measure lives beside the
  law (resolution, $|C|/(2\sqrt N)$, or $\eta=\min(1,\mathcal H/J)$, the
  Bartlett power, which is Z-exact and removes the $\sqrt N$ growth). Also:
  what the instrument declares ($k_i$ must be declared, not inferred), and a
  Herschel–Maxwell candidate for O3. It binds nothing; his call is owed on
  whether the tester's standard moves from the law's interval to the
  resolution.

## Added 2026-09-30 (afternoon): ours29

- `noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md`
  and three scripts beside it. After his rulings of 11:57 the working
  equation is **ours29**: ours27's flat factor per reading, weighted by the
  reading's own misfit (L21), with the product raised to the Bartlett power
  $\eta$ computed from the readings' own residuals. ours28's per-reading size
  law is dropped as a declared population of setups, which L19's reason
  excludes. L19's scope as ruled: nothing about any reading's signal is
  inferred from the others; the width may use each reading's own residual.
  In the tester the quoted width covers 93.7–97% wherever there is signal at
  $N=200$, in 1D, 2D and 3D. Lemma changes are offered in its §11, not
  placed. It binds nothing.

## Added 2026-09-30 (afternoon, later): `ours29/`

- A folder for the ours29 equation, created on his instruction at this
  level (not inside `noise_vs_signal/`). It is kept exclusive to the
  equation and mostly empty: only a README pointing to the defining record
  and reference implementation, which stay in `noise_vs_signal/`, where the
  work continues.

## Revised 2026-09-30 (afternoon, later): the lemma record placed on ours29

- `mass_estimator_lemmas_2026-09-29.md` revised in place on his instruction,
  carrying the ours29 record's §11: T8 is ours29's law, L10 restated (no size
  law is declared), L10′ struck and L4′ withdrawn, L19's reason and scope
  added, T31 (the Bartlett power) added, O14 and O15 opened, O3 now carries
  the equation. A second revision note sits under the first.

## Revised 2026-09-30 (late afternoon): O3 and O8 closed

- On his ruling ("Yes I like it. I agree.") the lemma record closes O3 and
  O8 from the O3 record in `noise_vs_signal/`: the flat weighting from L8's
  conditional form, L9 and L10 (Herschel–Maxwell), the stacked length from
  L12, noise units for any covariance from T7. Third revision note.

## Revised 2026-09-30 (16:26): O10 and O12 closed

- On his ruling, persistence is a premise of L1 (O10 closed, no argument
  owed) and the middle of the misfit weight is forced by the meaning that
  fixed $\kappa$ (O12 closed, from the O12 record in `noise_vs_signal/`).
  The declared-shapes table now holds $\pi$ alone. Fourth revision note.

## Added 2026-10-01: what can be rerun, and where tester work goes

- On his instruction ("Ok go ahead. Make safe changes", 30 Sept), the tester now holds this
  discussion's reference implementations in `reference/`, byte for byte:
  `ours27/ours28_reference*.py`, `noise_vs_signal/ours29_reference*.py`,
  `noise_vs_signal/ours28b/ours28b*.py` and the 21 September
  `mass_estimator_equation/mass_estimator.py` with `example.py`. Its tests
  check that each still prints what it printed, and that it reproduces the
  numbers these records state. Tester commit `eae695b`. The copies here are
  unchanged.
- A staging inventory of every computed result in this discussion, and
  whether it can be rerun now, is in
  `inbox/mass_estimator_reproducibility_inventory_2026-10-01.md`. Nothing was
  recomputed; what was not kept stays lost.
- For anyone adding work, the rule is the tester README's "How work gets
  into the tester": code lands as pushed tester commits (or a git bundle),
  every study is a seeded preset, reference implementations are copied into
  `reference/`, and a record cites the pushed commit and the preset, test or
  script with its seed. A sandbox commit hash is not a citation.

## Added 2026-10-10: one flat record, on units and axes

- `mass_as_force_paired_with_reciprocal_acceleration_2026-10-10.md` is a
  **flat member**, written at the user's request. His idea: every cancelling
  of units is a pairing with a reciprocal, not a division, so mass is force
  paired with the reciprocal of acceleration, $m=\mathbf F\cdot\mathbf a'$
  with $\mathbf a'=\mathbf a/|\mathbf a|^2$. It covers the **aligned case
  only**, on his instruction, and tells how the session got there: the plane
  wave as the $\mathrm m\cdot\mathrm m^{-1}$ example, what the literature
  already establishes, the two magnitudes, constant rates built one
  multiplication at a time ($m=\Delta\mathbf p\cdot(\Delta\mathbf v)'$), and
  the terminology **units** (M, L, T) versus **axes**: only L has axes, so
  only the $\mathrm L\cdot\mathrm L^{-1}$ pairing can need a choice. §8 parks
  the unaligned material and a units argument offered for O11. It scopes and
  binds nothing.
- Revised the same afternoon on his instruction: §9 corrects the comparison
  with the wave (the mass corner is the forward pairing, in two versions,
  $m=\mathbf a'(\mathbf F)$ and $1/m=\mathbf F'(\mathbf a)$, agreeing only on
  aligned pairs), reads L5 and L14 as the linearity of the two pairings, and
  notes the L14/swap trade, which he ruled stays unruled: L14 remains
  provisional.

## Added 2026-10-10 (later afternoon): the lemma record superseded again

- `mass_estimator_lemmas_2026-10-10.md` **supersedes**
  `mass_estimator_lemmas_2026-09-29.md`, on his instruction ("Just copy the
  old lemma list and edit the new version. That way you can write all of it
  as accepted and just put a provisional above all of it"). It is the 29
  September record copied, every row keeping its 29–30 September state, with
  the 10 October changes marked inside the rows and a **provisional** banner
  over all of them: the reading's window enters the type and L18; L5 and
  L14 are read as the linearity of the two pairings (L14's status unchanged,
  the swap trade unruled, by his ruling); T1 is restated in pairing form;
  T32 (the misfit is the pairings' disagreement in noise units) and T33
  (same window) are new; O9 gets the window criterion as a proposal; O11
  gets the units argument. The second `supersedes` in this discussion; the
  29 September record stays in place, unedited, as the 27 September one did.
- The units-and-axes record gained a §9.4 the same afternoon: the start on
  its four next moves, with the checks the lemma rows cite. The check script
  ran against `noise_vs_signal/ours29_reference.py` on his laptop and is
  **not yet kept**; where it goes (beside the record, or into the tester as
  a test) is still his call.

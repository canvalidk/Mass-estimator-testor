# noise_vs_signal/ — the signal-versus-noise investigation

A **working folder** inside `_dscn_mass_estimator`, created 2026-09-29 on the
user's instruction ("Can you make a folder called noise vs signal and put it
there"). Like `ours26/` and `ours27/`, it is a grouping, not a nested
discussion.

**The question.** The equation has no notion of how much of a reading is
push and how much is noise. A series of pure instrument noise tilts the law
by about $\sqrt N$, and every open item in the lemma record touches this.
Start with the view record, then the three developments.

## Written 29 Sept

- `the_noise_measure_view_2026-09-29.md`: the view, how it was reached (the
  `==^` thread), the map of what it touches, five parked items.
- `null_series_tilt_and_property_N_2026-09-29.md`: the tilt's exact form
  (the carried state $C$, Rayleigh under noise), ours27 against ours28,
  property N, and the dichotomy: keep L19 and keep the $\sqrt N$, or share
  one scalar.
- `alignment_what_one_reading_can_say_about_noise_2026-09-29.md`: per
  reading, only size along the best line is signal evidence; evidence that a
  series is noise is collective in every dimension.
- `c_is_the_noise_measure_at_zero_signal_2026-09-29.md`: $c$ is the answer
  at zero signal, but the noise supplies no scale; $c$ is declared.
- `older_hints_on_noise_and_signal_2026-09-29.md`: eleven earlier records
  that touched the problem, one quoted sentence each. Pointer only.

## Moved in 29 Sept (night), on his instruction ("do A and B lists")

- `ours26_learned_omega_is_a_hint_2026-09-27.md`, from `../ours26/`: his
  ruling that learned $\omega$ is a hint, and the tempering identity. The two
  records it corrects stay in `../ours26/`.
- `ours27_candidate_size_and_null_readings_2026-09-28.md`,
  `ours27_zero_readings_and_shared_candidate_information_2026-09-28.md` and
  `ours28_known_empty_vs_measured_empty_2026-09-28.md`, from `../ours27/`:
  candidate size, null readings, and known-empty against measured-empty.
- `ours28b/`, the whole folder, from `../ours28b/`: the one shared-scale
  candidate and the ruling that withdrew it (L19).

**Links.** Some records elsewhere link these files by relative path; those
links now point one level off. Filenames are the identifiers
(`ORGANIZATION.md`); the records were not edited. The same holds for
`../../ours27/` links inside `ours28b/` records.

## Written 29 Sept (night, later)

- `readings_are_evidence_just_less_2026-09-29.md`: a defect separate from
  the $\sqrt N$ tilt. A reading misaligned far beyond its noise is treated
  as strong signal and collapses the law as $\sigma\to0$. His ruling, "The
  readings are evidence. Just less.", forces tempering by the off-line
  misfit, $E_i^{g(Q_i)}$ with $g(0)=1$ and $g\sim\kappa/Q$; six forms
  tested. Also records the evening's thread: signal (his term, replacing
  "push"), the back-calculated-error check, the non-codirectional averaged
  pair, and the plug-in signal share. Parks the collision with axis
  splitting.
- Same record, §10 (21:13–21:25): the error-inflation side path closed
  (misfit at a mass vs misfit at the best mass); his ruling "one reading is
  one completeness claim"; $\nu=k_i-1$; strong axis splitting false; the
  lemma record revised in place (L8, L21, T30, O3, O12, O13).

## Written 30 Sept (morning)

- `scattered_observations_2026-09-30.md`: a list, not an argument. The
  argument under F1 with three open joints; $\eta=\min(1,\mathcal H/J)$,
  the Bartlett power (Z exact, closed form from the carried state, the
  $\sqrt N$ gone with a ~2 nat floor; ours25's sandwich by another route);
  rotation kept, axis splitting dropped; what is forced in that chain; null
  evidence does not accumulate under Z and L19; what the instrument declares
  (same space, $k_i$, $\Sigma_i$; $k_i$ must be declared); his
  $d$-intuition split; the shape sought is one radial $\rho_k(r)$; a
  Herschel–Maxwell candidate for O3; and the reframe after his correction:
  the property to drop is that the law carries its own uncertainty. The
  noise measure lives beside the law. His call owed on what the tester
  counts as a failure.

## Written 30 Sept (afternoon): ours29

- `ours29_misfit_weight_and_bartlett_power_2026-09-30.md`: the working
  equation after his three rulings of 11:57. Each reading gives the flat
  factor weighted by its own misfit, and the product is raised to the
  Bartlett power $\eta=\min(1,\mathcal H/J)$, computed from the readings'
  own residuals; ours28's size law is dropped. L19's reason ("it can't rely
  on things like a repeated experimental setup") and its ruled scope: nothing
  about a reading's signal is borrowed; the width may use each reading's own
  residual. Z forces the raw $J$ and the whole projective line. Also the
  eigenvalue dictionary, the combination case, the tester study (quoted width
  93.7–97% with signal at $N=200$; a bounded null floor of 81–85%) and the
  readout (the point stops wandering on pure noise).
- `ours29_reference.py`, `ours29_reference_checks.py`,
  `ours29_reference_worlds.py`: the reference implementation and its checks.
  The reproducible study is in the tester (`misfit`, `calibrate: bartlett`,
  `presets/ours29_noise_vs_signal.yaml`).

## Written 30 Sept (late afternoon): O3

- `o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md`:
  an offered closure of O3. The flat-in-vector weighting (ours29's whole
  per-reading factor) follows from L8's conditional form and L9 by
  Herschel–Maxwell, with L10 setting the scale to infinity and L7 removing
  the improper members; additive state (T12's exactness) is then a theorem,
  and the exponential is the only kernel with it. The stacked length follows
  from L12 applied to components as atoms. Rests on L8's conditional form
  (ruled) and L12 (medium). His ruling owed to close O3. Its §4b (added
  the same afternoon) shows O8 dissolves too: under L10 the weight is flat,
  and T7 fixes noise units for any covariance, checked by quadrature under a
  strongly anisotropic, correlated $\Sigma$.

## Written 30 Sept (late afternoon): O12

- `o12_the_middle_of_the_misfit_weight_2026-09-30.md`: an offered closure of
  O12. The meaning that fixed $\kappa=\nu$ (the share of a reading's
  self-contradiction its declared noise accounts for), read at every misfit,
  forces the kink $\min(1,\nu/Q)$: the same rule as $\eta$ one level down and
  the PDG scale-factor convention. Its cost on ordinary data is measured
  (about 7% width at strong signal in 2D/3D, coverage 94–96%). Retracts the
  earlier "threshold too low" reading. His ruling owed.

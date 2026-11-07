# VD mass estimator — the goal

**Written 29 Oct 2026 by Claude from VD-docs and the tester, for Can to correct.** This is what every session starts from. §0 is the short form (paste it into the project's custom instructions); the rest is the full statement.

---

## 0. Short form

**Target.** The body of the VD-Newton entry `mass := mass_estimator(F, a)`: the one function from readings (force and acceleration vectors with declared error covariances, a declared candidate subspace, a reading structure, a window) and a declared starting law π, to a law of mass on (0,∞) and a point — *and the proof that it is the only one* the version's commitments admit. The working equation is ours30 and the newest probe ours33 (4 Nov, on main 7 Nov); the proof is `definition/proof.md` in the tester, with its holes typed (G1–G27).

**How the work goes (Can, 29 Oct).** The candidates are probes. Each is built to the current premise list and is *supposed to break*; the break names a premise Newtonian practice has been running on unwritten. The growing premise list is the output, not a symptom. Sixteen probes so far (ours33, 4 Nov, is the latest, run beside ours30); each will fail too. Phase two, not yet started as a document: sort every find into three bins — a hidden commitment of Newton (persistence, additivity, positivity, one object), a hidden commitment of measuring anything (isotropy in noise units, the Gaussian model, no borrowing), or an input. `ladder.md` §2 draws the first line. The proof is relative to that sort. "Agrees with experiment" is a check on the survivor, not on a probe; no experiment is being run.

**Why it is hard.** There is no literature answer: the textbook has a type-level hole (m = F/a is typed on exact reals; no function from readings to mass exists in it), standard statistics assumes a loss function and the whole point is that the properties of mass (composition, the swap, Z, units, frame independence) replace the loss function, and the premises that do the work were not written anywhere before a probe broke on them.

**Standing rules.**
1. This is a long-term math problem, a sub-problem of the frontier. Every probe is presumptively dead; what it is for is the premise it surfaces. The only permitted report on a probe is *passes / fails / unknown*, per premise, plus the premise its failure names. Never "this resolves it."
2. At a fork, sort it: forced, input, or wrong premise (TASK.md §3). Never pick by test score; never hand Can a fork to pick a side on. The sort is Claude's deduction from the goal and provisional, like every premise, his or Can's (R-2026-11-07-1); Can corrects it where it misreads the goal, and a correction is itself a conjecture, made with more of the goal in view. Where the sort turns on what the goal means, ask that, as a question (`definition/bounds.md`).
3. The equation's existence is the one working assumption never dropped (TASK.md §4). It directs the search; it never discounts a proof. Bedrock premises (L1, L2, L3) are never dropped; if a strong premise has to go to save consistency, say so, and say which of two things happened: the deduction was wrong (the premise misread the goal: restate it, keep the old wording struck in its slot; L11, 7 Nov) or the goal was wrong (nothing can satisfy what it means: evidence against the assumption). The first is ordinary; the second is the finding.
4. VD-docs holds decisions and is the output; the tester holds the work. Rulings are Can's, verbatim, dated. Mark what is Claude's as Claude's.
5. Read `docs/existence_argument_audit_2026-10-29.md` before the textbook prior fires. Then `TASK.md`, then `definition/proof.md`.
6. **Bring the standard mathematics in first:** before working a fork, name what it is in standard mathematics or statistics and use what that theory already knows. The project's vocabulary is for the dictionary's own premises; it is not a licence to re-derive known results, and Can should not have to supply them (29 Oct; T36).

---

## 1. The target, precisely

The VD appends commitments in order, and an entry is complete when appended. The entry for inertial mass is written at Newton II (E10) as `mass := mass_estimator(F, a)`, and its body has to be the one function that the commitments *later* in this version (persistence and lumping at R3, the force sum at R4, Newton III at R5) leave admissible. So the estimator's form is a theorem of the whole version, worked out ahead of the entry (F4: commitment order suspended for the *work*, not in the log). Energy is not in this version; a version with energy before mass is a different version (R-2026-10-27-2).

Two halves to the deliverable, both in the tester's `definition/`:

- **The frame** (T36, kept central, R-2026-10-29-3): a reading is a vector in the stacked force–acceleration space, whose only metric is the noise's Σᵢ⁻¹ (its axes share no unit); a mass is a subspace of that space; Hᵢ is the projection onto it and Eᵢ − Hᵢ the residual, with kᵢ dimensions per reading. η's disagreement and the readings' misalignment are two components of one residual.
- **(A) The equation** — `equation.md`, currently ours30. Law: p(m) ∝ π(m) · exp(η Σᵢ gᵢ Hᵢ(m)/2), with Hᵢ the reading's explained energy in noise units on the whole projective line of slopes, gᵢ = min(1, κᵢ/Q₊ᵢ) its misfit weight at its best *positive* mass (κᵢ = kᵢ − ½), η = min(1, 𝓗/J) the Bartlett power from the series' own disagreement about m. π half-Cauchy with centre c = σ_F/σ_a by default (L3′, unruled). Point: the law's median (28 Oct, provisional, L23). Exact inputs (σ = 0) go to the three-way check; "undefined" lives only there, at exact 0/0.
- **(B) The proof that nothing else does** — `proof.md`, eighteen steps from the type and premises to ours30, every step marked ruled / forced / forced-in-archive / input / provisional / gap. What is proved today: *on isotropic uncorrelated channels, any function meeting the premises as ruled, with L22, L23 and η's finite-N form taken as given, is ours30 up to the declared inputs.* Three things taken as given, two of them rulings. **Since 30 Oct that statement is conditional** on reading R1 as not constraining the misfit weight: the rulings at steps 9–10 conflict (G18), left standing until something forces one side (R-2026-11-01-1); if R1 constrains g, steps 9–11 are open and ours30 is not the function (`proof.md` §4). **Since 31 Oct** ours30 also fails L11 at fixed per-reading signal (G22), so steps 8–10 are open on that ground too.

**The uniqueness principle** (TASK.md §1–2, R-2026-10-24-3): mass from readings is not a function of the readings alone, as log is not a function of w alone. The missing pieces — the instrument (Σᵢ), the structure, the subspace, the window, the prior (π) — are inputs, declared before data and reported with the result. The function is unique *given* them. The all-zero experiment returns the centre of π as log returns the branch.

**The guards on that move** (TASK.md §5): an input must be something real measurements need and practice supplies, declared before the data; moving a fork to the input because it is hard is a dodge. And the survivor must agree with experiment — detection limits on null measurements, combination with the density route and weighing, and real data. That is a check on the survivor, not on a probe, and it has not been made. The public robot dataset is not this check: it is a phase-one testbed, where probes break and name premises (O16; R-2026-10-31-1). No experiment of Can's is being run.

---

## 2. Why we think it exists, in one paragraph

See `docs/existence_argument_audit_2026-10-29.md` for the chain and its audit. Short version: the textbook's mass corner is overdetermined by d−1 and typed on exact inputs, so estimation is internal to the mass entry (forced); a definition cannot return "undefined" where its own other routes return 5 ± 0.5 kg, so the estimator is the definition and the law is its limit (ruled, with a VD-internal argument); a VD entry has one body, so every choice the properties leave open is an alternative mechanics and uniqueness is the stake (ruled, forced by what a dictionary is); and declared inputs are branch cuts, not choices. What that establishes is a well-posed uniqueness problem of the Cox/Shannon kind. That the premise list is *consistent* — has a model — is the conjecture, and the proof is the proof.

---

## 3. The project's shape: three repos

| Repo | Holds | Who reads it | Rule |
|---|---|---|---|
| **VD-docs** (private) | the decisions: founding arguments, theses, proofs accepted, rulings; the equation note and lemma note as draft lines; the Theory and Newton-analysis records. **This is the output of the project.** | Can, for human reading and absorbing ideas | only made decisions enter (`_dscn_mass_estimator/README.md`); work does not come here; a new statement of the lemmas or equation is a new draft of its note, not a record |
| **Mass-estimator-testor** (public) | the work: `definition/` (type, premises, theorems, equation, proof, ladder), `problems/` (one file per open item), `rulings.md` (Can's, verbatim, append-only), `briefs/` (Codex queue), `met/` + `tests/` + `presets/` (code), `docs/` (archive, moved from VD-docs 22 Oct) | Claude; Codex (Sol) for briefs | the living layer is the state; do not reconstruct it from records; a theorem is proved only when its file names its test; decisions go back to VD-docs through `outbox/VD-docs/` |
| **VD-Newton** (public) | the dictionary version itself: the ledger, the entry packets (E10 d4 is the mass entry, candidate, not promoted), the NML lesson drafts | Can | the rungs the proof rests on live here; R3 (persistence, lumping) is unwritten, R5 (Newton III) provisional, so the body cannot be written into E10 until both are settled (`ladder.md`) |

Division of labour (29 Oct): proof-first; derivations briefed to Sol as one sorted fork with a named deliverable each; sorts and write-ups Claude's; rulings Can's. Intake (provisional, 7 Nov; `briefs/README.md`): Claude merges brief branches as delivered and owns the file organisation, not a re-check; exploration branches are input, harvested by need.

---

## 4. The premise list so far (what a probe is checked against)

A probe is reported against every row; a row is a find from an earlier probe or a rung of the version. The phase-two sort (Newton's / measurement's / input) is not yet a column here — it should be. The right column is where the check lives today.

| premise | one line | grade / status | checked by |
|---|---|---|---|
| L1 | Newton II, one object, one m > 0 shared by the series; assumed, not tested (no gate) | bedrock; persistence a premise (30 Sept) | construction |
| L2 | every declared σ > 0; on measured inputs never "undefined", all-zero included | bedrock; ruled 25 Sept, 29 Sept | `test_every_measured_input_gives_a_mass` |
| L3 | units covariance: law and point scale by λ/μ | bedrock | `test_swap_and_units` |
| L4 | π symmetric under m/c ↔ c/m | strong | same |
| L5 | lumping: masses add at common acceleration | strong (rung R3, unwritten) | `test_composition_law` (law level since 28 Oct) |
| L14 | common force: inverse masses add | strong, provisional (rung R5) | same |
| L8 | one reading is one completeness claim; k−1 codirectionality conditions | ruled 29 Sept | T6 tests |
| L9 | isotropy in noise units | strong | T6/T7 tests |
| L10 | no size law declared | strong, ruled 30 Sept | same |
| L11 | consistency: m̂ → m when the readings identify m, $\Delta_N/\sqrt N\to\infty$ (restated 7 Nov; was "as total push → ∞", struck: T41) | strong; **ours30 fails it** at fixed signal as N grows (G22). The restatement is Claude's, provisional: the old wording took energy for information (T41, Codex, not checked by Claude); the rate is the literature's (B19; = T42, T53). Open: the full type, and whether the rate is Newton's or measurement's (G25, O19) | `test_t17_*`; `test_weight_consistency.py` |
| L13 | Gaussian noise model | stipulated (31 Oct: T23's forcing routes refuted, B05); conditionally forced by Herschel–Maxwell | T23's tests |
| L16 | frame independence: rotations, signs, Galilean boosts | ruled 28 Sept | `test_rotation_with_covariance_changes_nothing` |
| L17 | same candidate law, same mass | ruled 28 Sept | construction |
| L18 | declared inputs are correct (σ, structure, window; Nᵢ proposed) | ruled | construction |
| L19 | no borrowing between readings; joined only through m | ruled 28/30 Sept | `test_readout_weight` |
| L20 | Z: a reading with no information does not change the law (point-level Z *not* ruled, R-2026-10-29-2) | ruled 27 Sept; scope O7 | `test_a_zero_reading_leaves_the_law_alone…` |
| L21 | readings are evidence, just less: 0 < g ≤ 1, decreasing | ruled 29 Sept | `test_no_gate_*` |
| L22 | mass is a scalar: candidates share one direction | provisional (G10) | construction |
| L23 | the point is parametrisation-free | provisional (G6) | `test_t35_*` |
| N | the width tracks how far the point moves under repeats (η); a true statement, not nominal coverage (R-2026-09-29-4). Uniform coverage in the corner is impossible (Gleser–Hwang, O15, B13: infimum 0 at every finite N). ours30's 0.93–0.95 at N = 80–200 is finite-N only: its point is inconsistent at fixed signal (G22) | property | `test_ours29`; `tools/ours30_coverage_vs_n.py`; O15 |
| R-09-07-1 | smallest error under the law: Q₊ | ruled 7 Sept | `test_q_plus_*` |
| T30 | bounded confidence on a fixed misaligned pair as σ → 0 | forced by meaning | `test_t30_bounded_tilt_*` |
| guard | inputs cost: every input is something practice supplies, declared before data | TASK.md §5 | by sort |
| guard | agrees with experiment | TASK.md §5 | **not yet made**; a check on the survivor (the robot is a testbed, R-2026-10-31-1) |

Not yet checked for ours30: its coverage on its own world (only ours29's was measured).

---

## 5. The core samples (the probe history)

Every probe, what it was, what broke it. The breaker is the find: a premise that was not written down before. Read the right-hand column as the list of what Newtonian practice was silently relying on, in the order it surfaced.

| # | candidate | date | what it was | what killed it |
|---|---|---|---|---|
| 1 | Deming best point | 1 Aug | errors-in-variables profile minimiser; the positive root of the quadratic | boundary-seeking: argmin at 0 or ∞ for any zero reading, and the positive root is the profile's *maximum* for anti-aligned data (Sept 2 audit). Fails L5 and L14 by type (a mode) |
| 2 | median with excitation model | 27 Aug | median of the induced mass under an excitation model and a mass prior | needed an excitation model and a mass prior the type does not carry; superseded by integrating over the compatible-pair law (spec d1 §23) |
| 3 | spec d1 / ratio of means with gate | 30 Aug | integrate over compatible pairs under the flat measure; point E[f]/E[α]; absolute compatibility gate Q_crit = 9 | the gate: Newton II is a premise, not tested (R-2026-09-07-1); the readout form survived into ours21 |
| 4 | ours21 | 21 Sept | the equation spec with code; flat magnitudes, uniform direction; 1D = Liseo's flat posterior, positive part | the measure: axis splitting warrants a different joint law in 1D and 3D (2D identical) → ours24 |
| 5 | ours24 | 24 Sept | joint law via axis splitting (`cartesian2`) | the many-reading readout summed per-reading acceleration sizes → ours25b |
| 6 | ours25 / ours25b | 25 Sept | identification-limit study; ours25b: the many-reading readout on the stacked length (T9) | per-axis noise not handled → ours26 |
| 7 | ours26 | 26 Sept | per-axis noise, diagonal covariance; learned ω studied ("a hint, not the solution") | depends on the choice of axes: a rotated frame gives a different answer (L16) → full covariance |
| 8 | ours27 | 27 Sept | full error covariance per reading (forced by L16); flat push | under-covers (66–87%): the flat push is "every push is infinitely large" and T18's width factor is exactly its error |
| 9 | ours28 | 28 Sept | push scale as a per-reading nuisance under the arcsine size law; coverage 92–98% | a declared size law is a declared population of setups, which L19's reason forbids; withdrawn 30 Sept with L10′ and G (R-2026-09-30-3). **4 Nov (O14):** the cheat made precise: it declared the law of the one quantity η estimates, the signal share Δ/(Δ+M); η → that share on long series, so the fudge back to ours28b exists |
| 10 | ours28b | 28 Sept | one candidate-size scale per noise pair, shared over the series | borrows between readings (L19, ruled 28 Sept night) |
| 11 | ours29 | 30 Sept | ours27's flat factor × misfit weight g(Q) on the whole-line misfit × Bartlett power η | the whole-line misfit cannot see the sign of F·a: an anti-parallel pair has Q = 0 ("we know ours29 is wrong"); and its point (the cloud readout) fails Z (T29) |
| 12 | the cloud readout (T10) | 28 Oct | ratio of means over the candidate cloud on the stacked length | fails Z at the point and reads the instrument's s; replaced by the median (provisional). Took T3/T4, T9, T21, L12 off the path |
| 13 | ours30 | 29 Oct | ours29 with the misfit at the best *positive* mass, κ = k − ½; point = median | **working; a break found 30 Oct (provisional, T37):** with the σ's believed as the covariance (R-2026-10-30-1), its κ = k − ½ needs a noise size with no mean and its cap g ≤ 1 holds only for the fixed-size Gaussian. Names the premise: what a declared σ is. At strong signal the cost is efficiency (10–17% width in 3D), not coverage. **Second break, 31 Oct (B05):** the sign of F·a in the weight is correlated with the score, so at fixed per-reading signal the point is inconsistent (truth 2, signal 1, 3D: converges to 1.81). Names the conflict: the sign in the weight (T30) against consistency (L11); G22 |
| 14 | **ours31** (probe, beside ours30) | 30 Oct | T37's exact cut factor per reading, ν_eff declared (GUM), no η, median; `met/ours31.py`, O18 | **breaks at weak signal in 3D** with many readings (coverage 0.57–0.71 at N = 20; ours30 0.81–0.88): re-attributed the same day: it is ours27's flat-latent-measure problem (L10's flat endpoint, T18), which η had patched; ours31's size posterior slightly helps. η corrected the latent measure, not the noise. Passes L2 (B11's inputs), Z, T30, units, swap, L24 |
| 15 | **ours32** (probe, beside ours30) | 1 Nov | ours31's reading level with η put back (R-2026-11-01-1); `met/ours32.py`, O18 | **not yet broken:** ours31's weak-signal break goes (0.82–0.87 in 3D at N = 20, ours30's level or better); better than ours30 where the noise really varies (ν_eff = 10, N = 5). Fails L21's power form by construction (G19, standing). Its reading level passes the consistency check ours30 fails (G22; 4 Nov). What it does not fix: many-reading weak-signal coverage, η's finite-N form and the latent measure (O14, O15) |
| 16 | **ours33** (probe, beside ours30) | 4 Nov (main 7 Nov) | T48's exponent (ν+k)/2 on ours32's residual, η as T47's normal-score warp of the raw law instead of a power; point = raw median; `met/ours33.py`, O18 | **broken at the report:** at small η the warped width is π's Cauchy tail (T39's pair: 78 decades around a raw interval of 1.5%); names the output type (the width is not a quantile of π·L). And no width exists without π (B16). Passes L2, Z, units, swap, T30 and both T39 breaks. Its latent measure is not forced (G19: B17, E5); η's estimate is still open (G3: B15); its warp does not commute with composition (B18, T51) |

Struck premises, kept in their slots: L4′ (push–noise symmetry), L6 (C3, now T21), L10′ (push scale as nuisance), G (the size law). Each went with the candidate it served.

The one recorded case where three things could not all stand: O1 (Z at the point, a readout at the all-zero experiment, the cloud readout). Resolved by dropping the lowest-graded — and then, on 29 Oct, by finding that point-level Z had never been ruled at all.

---

## 6. Open holes

**In the proof** (`proof.md` §5), by type:

- *standing, until forced* (R-2026-11-01-1): **G18** R1 against O12's cap and κ's meaning; **G19** L21's power form, T37's cut factor, or T48's normalized residual, as the body (4 Nov: L11 with T30 is a candidate for what forces it, O18; T48's latent measure is forced neither by T36 nor by the instrument's moments, B17, E5).
- *sort* (Claude), formerly typed "rule (Can)" (R-2026-11-07-1): G25 L11's rate on the full type and its bin (restated 7 Nov); G10 L22 scalar mass as a clause of L1 (O17); G11 Nᵢ as an input (TASK §5's test); G9 π's centre (O11); G6 L23 (O2; a draft sort on `claude/b21-brief`). Where one of these turns on what the goal means and the page does not say, it is a *goal* question to Can (`definition/bounds.md`).
- *sort* (Claude): **G22 the sign in the weight against L11**; G17 the 1D factor; G13 why the weight is a power on Eᵢ and which premise makes bounded confidence a requirement; G3 η's cap and power form at finite N (O14).
- *from the 3–4 Nov exploration* (harvested 7 Nov, not checked by Claude): G23 tempered tail index against amplitude consistency; G24 local η does not select a global correction; G26 the fitted J against sampling variance under leverage (T43); G27 the joint composition protocol (T51).
- *derive* (Codex briefs): G1 a Gaussian-forcing argument under named premises (B05 refuted both stated routes; L13 stipulated); G5 κ under general Σᵢ (O13, only if G18 keeps O12's meaning); G14 agreement with experiment (later: a check on the survivor; B01 merged 31 Oct as a robot input check).
- *write* (Claude): G15 19 theorem files still "carried from the records", ten on the path (T1, T2, T7, T8, T20, T26, T27, T30, T31, T33); G7 done (B10). G2's remainder, the coefficient, rides on G18 (B07); G8, G12 done (B08, B09).

**Problem files** (`problems/`): O2 readout, O4 (Gaussian forcing, B05), O5 (parked, B09), O6 π's shape (family derived, B08) (needs a real idea, not more tests), O7 Z's scope, O9 window criterion, O11 centre, O13 contradiction between readings, O14 η at small N, O15 the null floor (corner: B13; cure table closed: B14), O16 shared errors (the robot's miss, B01), O17 anisotropy, O18 ours31/32/33, O19 what growth makes L11 attainable.

**Where the argument is most exposed** (7 Nov): L11's old wording had no model on the reading type (T41); it is restated with the literature's rate (provisional, Claude), and what is open is whether that rate is a property of mass or of measurement, and its form on the full type (G25, O19). Then (4 Nov): ours30 fails L11 because the sign of F·a entered its power weight (G22, a sort owed), and that also makes its pure-noise law concentrate on c (O15); ours32's cut factor does not share it. And (30 Oct) the rulings at the misfit weight conflict (G18), left standing (R-2026-11-01-1); what proof.md proves is conditional on reading R1 as not constraining g. T37's tail index α must pass the inputs-cost guard (G20). Earlier (from the audit): the point readout is the least-forced step since 28 Oct (L23 is a *rule* hole); and π's shape and the zeros are one question, decided by the input-vs-choice criterion of TASK.md §5.

---

## 7. What Claude does with this

- Start every session here, then `TASK.md`, then `proof.md` §4–5 and `docs/STATE.md`. The living layer is the state.
- Deduce the day's work from what Can says that day, against this document. If the day's ask conflicts with a standing rule, say which.
- Turn ideas into contributions: a sort written into a problem file; a brief for Codex with one deliverable; a theorem file with its test; a decision through the outbox. Not a record at the end of a conversation.
- When a probe breaks, the first question is not "how to fix it" but "what did it find" — name the premise, then ask which bin it goes in (Newton's / measurement's / input).
- Keep this document and the audit true: when a hole closes or a probe breaks, the row moves.

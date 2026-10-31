# The ladder

What determines the body of `mass := mass_estimator(F, a)`: the commitments of this version,
in their order, and what each removes from the family of admissible bodies (`../TASK.md` §7,
R-2026-10-27-2). Every freedom left in the definition gets its status from here: a rung of
the version removes it, or it is an input, or a rung is mis-stated, or "unknown yet"
(R-2026-10-27-1).

Built 27 Oct by Claude from VD-Newton at commit `61afb9e` (2026-10-24): the Design 2.1 ledger
(`07_design_2_1/DESIGN_2_1_ENTRY_STRUCTURE.md`), the E10 packet d4
(`07_design_2_1/entry_packets/ENTRY_PACKET_E10_INERTIAL_MASS_MATERIAL_OBJECT_d4.md`) and the
NML3–NML6 forward map (`08_nm_lesson_drafts/NML3_TO_NML6_PROVISIONAL_FORWARD_MAP_2026-08-08.md`).
Statuses below are VD-Newton's own. Everything in the "removes" and "findings" columns is
Claude's reading, provisional.

## 1. Rungs of the version (mechanics)

| Rung | Where in VD-Newton | Status there | What it commits, for mass | Tester premises | What it removes from the body | Theorems |
| --- | --- | --- | --- | --- | --- | --- |
| R0 | E-2..E0, K1–K4 (raw kinds, time, reference-frame, displacement, straight-line) | in ledger | coordinates and units are conventions | L3, L16 (rotation) | any body that depends on frame, sign convention or units beyond m's own scaling | T22 (full covariance, from L16) |
| R1 | E1–E6, Newton I | in ledger | the acceleration channel is inertial acceleration; Galilean frames | L16 (boosts) | bodies that read the frame's motion | — |
| R2 | E7–E11, Newton II; **E10 is the mass entry** | E10 d4 candidate, not promoted | one material object; mass is a **positive scalar**; $F = m\,a$; an exact pair recovers $m$ by check; a measured pair takes the joint branch, **selected by provenance, not by alignment** | L1, L22 (L1's "scalar"), L2's measured branch | non-scalar and non-positive bodies; per-reading alignment gates; the exact branch's three outcomes | T1, T2, T34 |
| R3 | E12 wall / NML2.2: address, persistence, trackability, lumping | **deferred; not written** | one mass per addressed object over time; masses of lumped objects add | L1's persistence clause, L5 (mass additivity) | per-reading masses (one $m$ per series, so readings combine); bodies that do not add under lumping | T3, T4 (with R5), T12 |
| R4 | E13–E17, NML3 force sum | active, under correction (aggregate-vs-member bug, 12 Aug) | net force is the vector sum of the impressed forces of a closed interaction-set | — directly | with R5: the lump's force is the sum of its parts' forces (finding 1) | — |
| R5 | E18–E23, NML4 action-reaction (Newton III) | provisional | reaction force is the negative of the action force | L14 (common force, reduced mass), F4 | bodies whose inverse does not add under a common force | T3, T4, T5 (swap), T21 with L15–L18 |
| R6 | E24–E41: NML5 acting objects, NML6 closure, meta-typing | provisional | acting-object identity; mechanical closure | none found | **unknown yet**: no use for the estimator found | — |
| — | energy, kinetic-energy form | **not in this version** (R-2026-10-27-2) | — | — | nothing; "mass is the metric" (O17) is a lens on R2's scalar, not a rung | — |

Two rungs the body rests on are not fixed in VD itself: R3 is unwritten and R5 is
provisional. So anything forced through them is provisional by inheritance, and the body
cannot be written into E10 until both are settled. That is the ladder run to the end, which
`../TASK.md` §7 requires.

## 2. The measurement constitution (not rungs)

Most of the premises are not commitments of Newton's version at all. They are about how any
measured quantity is read: they hold for every entry with measured inputs, and they come from
the constitution and Can's rulings (F1, F2, F3; R-2026-09-27 to R-2026-10-10), not from the
version's order.

| Premise | One line | Grade or status | Bin (29 Oct sort) |
| --- | --- | --- | --- |
| L2 | every declared σ > 0 (the noisy type) | bedrock, ruled 25 Sept | **M** for σ > 0 (JCGM: a measurement is a triple). The "never undefined on measured inputs" clause is a theorem of **N** (lumping) and **D** (F1): `why_mass_cannot_be_undefined` |
| L3′ | the default π's scale comes from the instruments | proposed | **I**: which member of the input family is the principal branch |
| L4 | π symmetric under m/c ↔ c/m | strong | **N** constraining **I**: the swap is Newton II's own symmetry (F and a on equal footing, R-2026-09-07-1), applied to the admissible inputs |
| L7 | calibration: R = f/α on sure inputs | mild | **N** (the exact branch, R2); its continuity clause is technical |
| L8 | one reading is one completeness claim (count mis-stated, 27 Oct sort) | ruled 29 Sept | **M**: what a reading declares about its own dimensionality. Written nowhere before 29 Sept |
| L9 | isotropy of the noise, in noise units | strong | **M**: about instruments, not bodies. That the law *cannot be read without it* is a fact about Newton's text (the pattern this project exists to surface) |
| L10 | no size law is declared | strong, ruled 30 Sept | **D**, via L19's reason: a definition cannot rely on a population of setups |
| L11 | consistency | strong | **N**: the exact law is the estimator's limit (F1's "the law is its limit") |
| L12 | one joining rule for series sizes | **medium** | **M**; off the path since 28 Oct |
| L13 | the Gaussian noise model | forced (T23) | derived from L11 (**N**) and L20 (**M**): no bin of its own |
| L15 | betweenness | ruled provisionally 28 Sept | **M**: a property of combining evidence, not of mass |
| L17 | same candidate law, same mass | ruled 28 Sept | **D**: the readout is a function of the law; it cannot read instrument facts |
| L18 | declared inputs are correct | ruled 28 Sept (σ, structure), 10 Oct (window) | **I** by construction |
| L19 | no borrowing between readings | ruled 28 Sept, 30 Sept | **D**: "the equation must act as the definition of mass, so it can't rely on things like a repeated experimental setup" (his reason, 30 Sept) |
| L20 | Z | ruled 27 Sept; scope O7 | **M**: what "no information" means is fixed by the declared inputs (`TASK.md` §6) |
| L21 | readings are evidence, just less | ruled 29 Sept | split: no gate is **D** (Newton II is a premise of the definition, R-2026-09-07-1); decreasing in misfit is **M** |
| L23 | the point is parametrisation-free (listed with the rungs in proof.md; sorted here) | provisional | **D**, provisional: a requirement on what a definition's *point* is, not on mass or on instruments. The probe that tests it varies the point functional (median vs geometric mean), and that is a legitimate probe |

**The bins** (Claude's sort, 29 Oct, provisional; the question is Can's of the same day: "have
they always been in Newton or no?"). Three were asked for; the sort needs four:

- **N** — a hidden commitment of Newton's mechanics. Everything in §1 (L1, L3, L5, L14, L16,
  L22) and, from this table, L4, L7, L11. All rungs or readings of rungs. *A probe may not
  vary these.*
- **M** — a hidden commitment of measuring *anything*: it would hold for Hooke's k or for G
  exactly as for mass. L8, L9, L12, L15, L20, L2's type, half of L21. JCGM-shaped, but L8, L9
  and Z are not in the GUM; they are new. *A probe may vary these, and what it finds is a
  result about measurement.*
- **D** — a commitment of the dictionary: what makes the body a *definition* rather than a
  procedure. L10, L17, L19, L23, R-2026-09-07-1, F1, F2, the no-gate half of L21. Neither
  Newton's nor metrology's. **This bin is the project's own contribution**, and it is exactly
  the set a statistician would call arbitrary axioms. *A probe may vary these only to test
  whether the VD's notion of definition is the right one.*
- **I** — the caller's. L3′, L18, and everything in `type.md`. *A probe that varies an input
  learns nothing about the definition* — which is why O6 (π's shape) needs an idea and not
  more tests.

Two things the sort shows. First, the smuggling check this section exists for passes by
inspection: no M or D premise makes a mechanical claim, and the N bin is exactly the rungs.
Second, the bins say which probe is worth building next: the open D items (L23; whether
L10's reading of L19 is the only one) and the open M items (O7, Z's scope; O14, η at finite
N) are where a break would name something new. A probe aimed at an N premise is testing
Newton, and one aimed at an I premise is testing the caller.

The body is a theorem of all four columns. The check this section exists for: **a measurement
premise must never force a mechanical choice.** If one does, it is doing a rung's job, which
is smuggling. Where a measurement premise and a rung both bear on one freedom, the rung
decides.

## 3. Every open freedom, by the ladder

| Freedom | Rung that bears on it | Status from the ladder |
| --- | --- | --- |
| The readout's form (O2) | R3 + R5, but as exact-value commitments, delivered by T1 and L11 for any consistent point | **revised 28 Oct:** composition is a law-level commitment (R-2026-10-28-1), checked by coverage; the point is a law functional (L19), and L23 + T5 make it the median (T35), provisionally |
| The series' size inside that ratio (O2) | none: R4 joins sources at one time, not readings; composition passes every 1-homogeneous join | **unknown yet, low stakes** (27 Oct): not smuggling, under 0.2% effect (O2, evening) |
| The weight's noise part (O2's blocker) | none; L17's reason (a measurement premise) bears on it | **provisional, outcome 3:** T21's domain mis-stated against L17's reason; repair unknown yet (O2, evening) |
| The all-zero value (O2) | none; the constitution (Z, units, swap) | derived, not a rung: c |
| π's shape (O6) and centre (O11) | none | input (`../TASK.md` §2) |
| The sign weight's numerator (O13) | R2 forces positivity to enter the weight | that it enters: forced from R2. **Ruled 29 Oct (R-2026-10-29-1, ours30):** the misfit at the best positive mass, numerator k − ½ (the value proved for isotropic channels: T38) |
| Z's scope for noise (O7) | none | input, provisional (O7) |
| η at small N (O14) | none | measurement side; a proof |
| Shared errors (O16) | none | input (a series covariance), Claude's sort |
| Which mean under anisotropy (O17) | R2's scalar | does not arise in this version |
| L14 | R5 | provisional, matching R5's own status in VD |

## 4. Findings from building it

1. **L5's readout form may sit at R5, not R3.** Mass additivity, $m_{1+2}=m_1+m_2$, is R3. The
   tester's L5 is the readout law $R(X_1+X_2,Y)=R(X_1,Y)+R(X_2,Y)$, which also needs the lump's
   force to be $X_1+X_2$. If the $X_i$ are each part's net force, that holds only when internal
   forces cancel in pairs: R4 (the sum) and R5 (action-reaction). If the $X_i$ are external
   forces only, it needs R6 to say which are external. Either way L5 as used is later than the
   additivity it expresses. **Unknown yet** which, until NML2.2 is written. It does not change
   T4's conclusion; it changes which rung T4 is provisional on.
2. **E10 d4 has a slot this work fills.** The packet routes a measured pair to "a declared,
   licensed estimation policy" and leaves "how an estimation policy is selected or certified"
   to "later joint-measurement/estimation machinery" (d4, lines 365, 373). The tester is that
   machinery, and under F2 its answer is not a policy to declare: it is the one body the ladder
   leaves, with the declared inputs of `type.md` (Σ, structure, window, π) in the slot where a
   policy would go. That is a decision about E10's wording. It belongs in VD-docs when Can
   checks it, through `outbox/VD-docs/`.
3. **Persistence.** Can, 27 Oct: "clearly the definition of mass can't actually assume that the
   mass of something will remain the same with every measurement." On the ladder that is R3's
   commitment, after E10, so E10's body does not assume it: the one-reading body is R2's, and
   joining readings into one series is licensed by R3 through F4. Consistent with L1 as ruled
   on 30 Sept.

## History

- 2026-10-27: created by Claude from VD-Newton `61afb9e`.
- 2026-10-27 (evening): the series join and the weight's noise part, from the L12 check.
- 2026-10-28: the readout row revised (O2, 28 Oct).
- 2026-10-29: the sign numerator row, ruled (ours30).
- 2026-10-29 (night): §2 gets a bin column — the sort of every measurement-constitution premise into Newton's / measurement's / the dictionary's / input (Claude, provisional, on Can's question "have they always been in Newton?").

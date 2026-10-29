# Why we think the mass equation exists — assembled and audited

**Written 29 Oct 2026 by Claude**, from VD-docs (`_note_why_mass_estimation/`, `_dscn_mass_estimator/`, `_note_mass_estimator_lemmas/`, `Theory/theories_and_belief/`) and the tester's `definition/proof.md` as of commit "R-2026-10-29-2". Can's argument; Claude's assembly and audit. **For future sessions: read this before the textbook prior fires.** Can corrects it; the corrections are the content.

**Verdict in four lines.** There are three existence claims travelling under one name. Two of them are established — the problem is well-posed as a uniqueness problem, and the right research program is the Cox/Shannon/Jaynes genre (declare the properties, prove the function is forced). The third — that the property list is *consistent*, i.e. has a model at all — is the conjecture, and it is exactly what the proof would show. Nobody needs to believe the third to work on it; it is the working assumption, and the method for working under it is sound provided premises are dropped by grade and the graveyard is kept.

---

## 1. Three claims under one name

| | claim | status | where |
|---|---|---|---|
| **X1** | *Some* total function fills the textbook's gap: NM as practised is total at 0/0 and in noise; the textbook is a partial transcription. | argued, non-constructive; one weak joint (§3, P1) | `theories_and_belief.md` §5 |
| **X2** | Given the premises as ruled and the declared inputs ($\pi$, $N_i$, structure, window), the function is unique — and it is ours30. | partially proved; holes typed G1–G17 | tester `proof.md` §§2, 4, 5 |
| **X3** | The premise list is consistent: a function satisfying all of it exists. | **conjecture** — the licensed working assumption | nowhere, correctly |

X1 motivates. X2 is the theorem being written. X3 is the bet. Treating them as one claim is what makes a fresh reader either over-believe (X1 is "proved", so the equation "exists") or under-believe (X3 is unproved, so the whole thing is "just a framing"). Keep them apart.

---

## 2. The chain, with load-bearing steps marked

Each step: the claim; its status; **what it forbids** — because that is what makes it a selection principle when two networks of argument contradict.

**E1. Belief flows to generators; a theory is a total function $(\text{question},\text{situation})\to\text{string}$, believed at untried inputs.** *Philosophical premise.* Forbids: scoring a theory on its track record alone; treating a hypothesis as the thing confirmed.

**E2. NM has persisted at inputs where the written text fails ($0/0$, signal in noise), and nobody bankrupts physicists there. So practitioners execute something total there; it exists because it is run.** *Argued; the weak joint of X1* (§3, P1). Forbids: "the text is fine, experimentalists just know when not to use it" — that is text plus unwritten judgement, which is a total function with half the code off the books.

**E3. The textbook has a type-level hole.** $m=F/a$ is typed on exact reals with $a\neq0$; experiments produce $(\tilde{\mathbf x},\Sigma,P)$ triples; the textbook contains no function from readings to mass. *Forced — this is simply true of the text.* Forbids: treating the hole as a philosophical mood rather than a type error; treating "mass must be estimated" as a generic reminder that experiments are noisy.

**E4. Estimation is internal to the mass corner.** From three premises a VD reader already holds (addressed positive mass; measurements as triples; laws runnable as predictions), the audit axiom, and the authority ruling: the mass corner is overdetermined by $d-1$, exact recovery is a checked sum-type, independent readings are generically off the law, and reconciliation cannot move upstream without either leaving the exact definition inert or contaminating the other two corners (horn two bites for *predictive* traces). A compatible pair exists only as a distinct latent object built where both inputs are in scope, and building it *is* mass determination. *Forced* (`what_forces_the_chain_to_the_estimator.md`, S1–S12). Forbids: a short mass entry; a preprocessing step that "cleans" force and acceleration toward each other; dividing magnitudes automatically; "not codirectional → no mass".

**E5. The estimator is the definition; the law is its limit (F1).** Two arguments. (a) The law-rearranged candidate and the estimator candidate disagree in truth value on a reachable input (a null trial); an outside standard — detection-limit practice — adjudicates for the estimator. Premise (a) rests on, undefended: that a sensitivity bound is a legitimate determination of mass. (b) Stronger and VD-internal: a definition cannot disagree with another of its own routes. The density route (lumping over volume elements + homogeneity + persistence) gives a stationary object $5\pm0.5$ kg; an "undefined" from the estimator would be the definition contradicting itself. And "undefined" combined with other information must be absorbing (erases the density result — contradiction), identity (the constant factor a null series already contributes — redundant), or neither (arbitrary). *Ruled 27 Sept; argued in `why_the_estimator_is_the_definition.md` and `why_mass_cannot_be_undefined_2026-09-29.md`.* Forbids: "there is a super meta mass the routes fail to reach"; an output type that depends on the data; `undefined` on any measured input.

**E6. Uniqueness is the stake (F2).** Any choice the properties leave open is an alternative definition of mass, i.e. an alternative Newtonian mechanics; "the choices only differ where nobody measures" is not an answer. *Ruled; forced by what a VD entry is* — the body of an entry is one string, two bodies are two dictionaries. Forbids: picking a readout by test score; handing Can a fork to pick a side on; "good enough" estimators.

**E7. Declared shapes are inputs, not choices.** $\pi$'s shape and centre, $N_i$, the reading structure, the window: once declared the function is unique, as the complex logarithm is unique once its branch cut is. *Can's position, Oct 2026; its test case is open (O6).* Forbids: counting a declared input as non-uniqueness. **But see P3** — this step and E6 pull against each other, and the corpus has not yet said where the line falls.

**E8. Given the premises as ruled, the function is ours30.** Proof.md's honest statement today: *on isotropic uncorrelated channels, any function satisfying the premises as ruled, with L22, L23 and the finite-$N$ form of $\eta$ taken as given, is ours30 up to the declared inputs.* Three things taken as given; two of them are rulings, not mathematics. *Partially proved; 16 holes typed (write / derive / sort / rule / VD-Newton).* Forbids: calling ours30 "the" function before G6, G3, G10 close; citing a passing test as a proof.

**E9. The premise list has a model.** *Conjecture.* Nothing in E1–E8 establishes it and nothing could short of the proof. The evidence for it is indirect: thirteen candidate equations (the `oursNN` numbers are dates, not a count; the lineage is in `GOAL.md` §5), each killed by a *named* premise, never, in the records read, by two premises that could not both be satisfied — with one recorded exception (O1: Z at the point, a readout at the all-zero experiment, and the cloud readout could not all stand), resolved by dropping the lowest-graded of the three. That is what a consistent-but-unfound list looks like from the inside. It is also what an inconsistent list looks like until the last premise.

---

## 3. Audit

**On E1–E4: no hole found.** E4 is the strongest piece of the corpus and the Sept 2 audit already did the work: three premises, one ruling, theorems between. I re-read it looking for a smuggled step and did not find one. The one thing I would add is already in that audit's §8 and should be owned loudly: **nothing in E4 is specific to mass.** Every scalar-parameter corner of every law-triplet — Hooke's $k$, $G$, resistivity, any constant of proportionality — has two total forward corners and one overdetermined inverse corner, the same dilemma, the same retention requirement. Mass is where the house first noticed. This is a contribution to pass on, and it changes what the proof is *of*.

**P1. E2's weak joint: abstention.** The pressure list already has it — at $0/0$ practitioners abstain, abstention is unexploitable, so survival supports "unidentified here" as well as a number. The corpus has the reply but has not connected it: practitioners do *not* abstain at $0/0$; they go and get the mass by another route (density, weighing) and combine. Combination with a null series is the identity (T24), which is a *value*, not a gap — it is $\pi$ unchanged. So the thing practitioners run at $0/0$ is "the starting law, updated by nothing," and that is a total function. E5(b) is the argument; it closes P1 and the two records should cite each other. What P1 still wins: it says nothing about *which* $\pi$. See P3.

**P2. The descriptive framing overclaims; the explication framing is what the work does.** "What have physicists believed all along" cannot pick a unique total function from a finite record of practice (Kripke; the pressure list's #2, which the corpus concedes: "what does not survive: that finding it is pure description"). What the project actually does is *explication*: fix the properties of mass from NM's own other entries (lumping, Newton III, frame independence, units), then prove the function those properties force. That is a stronger position than description, not a weaker one — it is why there is a theorem to prove at all — and it costs nothing in the mathematics. It changes how the goal is stated to outsiders: not "we discovered what mass always was" but "here are the properties NM already commits to; here is the unique function they admit; here is that it agrees with experiment." Recommend stating it that way.

**P3. E6 and E7 are in tension, and O6 is where it is decided.** F2 says "differs only where nobody measures" is not an answer. The branch-cut position says a declared input is fine precisely because its influence is confined to where the readings are silent (consistency, T17: $\widehat m\to m$ whatever $\pi$; the all-zero law *is* $\pi$). These are the same sentence with opposite verdicts. The resolution is a *criterion* for when a freedom is an argument and when it is a choice — and the tester's `TASK.md` §5–6 has one (correcting an earlier draft of this paragraph, which said the corpus lacked it): F2 is true of choices *the function* makes, not the caller; an input must be something real measurements need and practice supplies, declared before the data and reported with the result; and the function must agree with experiment. What I would add as a sharper, checkable form: **a freedom is an input iff (i) the readings cannot decide it, (ii) its influence vanishes in the limit of total push, and (iii) every physical output in the identifying regime is invariant under it** — the three things that make a branch cut harmless. TASK.md already concedes (iii) strains for π ("π does change answers at weak signal"), and answers it by making the output type honest: a mass is (m̂, law, π, Σ). $N_i$ and the window pass trivially (they are instrument facts). $\pi$'s *centre* passes if O11's units argument holds ($c=s$ is then not free). $\pi$'s *shape* is the hard case: the family $m^{-1}g(\ln m/c)$, $g$ even, is forced by L3, L4 and properness (G8), but $g$ itself is not, and $g$ is visible in every weak-signal output. If $g$ cannot be forced by a property of mass, then either the test above admits it (and F2 is narrowed to "choices in the function, not in its arguments" — a real narrowing that should be ruled, not slid into) or the definition has a parameter no physicist ever declared, and the descriptive thesis (P2) takes the cost. This is the single point where the existence argument is most exposed — and it is the zeros question seen from the other side: the all-zero output *is* $\pi$ (T24), so "what should happen at the zeros" and "what is $\pi$'s shape" are one question, and the test above is how it gets decided.

**P4. The uniqueness leg was weakened on 28 Oct and this should be visible.** The strongest *forced* result in the corpus — the Sept 6 composition-uniqueness proof, where C1–C5 force the ratio-of-means readout with no loss function chosen — is now "kept for comparison" because the point switched to the median. Composition (L5, L14) now applies at the law level, and the point rests on L23 (parametrisation-free), which is *provisional*. In exchange the median passes Z at the point exactly and is a law functional (L17). That may be the right trade. But the chain audit's verdict of 2 Sept ("the form was being selected by the composition laws all along") no longer describes the point, and nothing of equal forcedness has replaced it: G6 (L23) is a *rule* hole, not a *derive* hole. So as of today, the point readout is the least-forced step in the chain. Either L23 gets a forcing argument — the swap (L4/T5) narrows the point to functionals that commute with $m\leftrightarrow1/m$ (the median, the geometric mean, any symmetric functional of $\ln m$), and L23's invariance under *all* monotone reparametrisations is what then singles out the median; so the question is whether L23 follows from a property of mass or is a declaration — or the ratio of means comes back at the point with a Z repair, or the point is declared an input. One of these should be the next sort.

**P5. The method — assume existence, drop a premise — has a known failure mode, and the guard is already in place but should be named.** Under the working assumption, any contradiction between networks can be resolved by dropping a premise, and if dropping continues the list becomes trivially consistent. The guard: premises are graded (bedrock / strong / medium / mild), the drop goes by grade, and struck items keep their slot. The rule to state: **bedrock is never dropped; if a strong premise has to go to save consistency, that is evidence against E9, not a move licensed by it.** The graveyard (thirteen candidates, each with the premise that killed it, plus every struck L; `GOAL.md` §5) is what lets a reader check that the assumption has not been used to launder an inconsistency. Keep it as a first-class document, not as history.

**P6. One thing the corpus should stop doing.** `why_the_estimator_is_the_definition.md` Results 1–2 computed the positive root of the Deming quadratic, which is the profile's maximum on half of all null trials; the Sept 2 audit caught it, the record is unchanged. It is a *decision* record now. Either it carries a correction line or a fresh reader re-derives the wrong object. Small, but it is the kind of thing that costs a session.

---

## 4. What Claude is convinced of, and what it does not need to be

**Convinced** (meaning: read the argument, looked for the hole, did not find one):
- E3, E4: the text has a type-level hole and estimation is internal. Forced.
- E5: the estimator is the definition. Ruled, with a VD-internal argument (b) that is better than the detection-limit argument (a).
- E6: uniqueness is the stake, because a VD entry has one body.
- That the program — declare the properties of mass from NM's other entries, prove the function is forced — is a well-posed research problem of the Cox/Shannon/Jaynes kind, and that nothing in standard statistics ("estimators are choices") bears on it, because that objection presupposes a loss function and the whole point is that composition and the swap replace the loss function.

**Not needed:** belief in E9. It is the conjecture. The honest stance is the one Can named: work under it, drop by grade, and the proof is the proof.

**The selection principle, as Claude will apply it when two networks contradict:**
1. Never drop bedrock (L1, L2, L3). If a resolution requires it, stop and say E9 is under threat.
2. Drop by grade; among equal grades, drop the one whose *forcing argument* is weakest, not the one whose test score is worst.
3. Before counting a freedom as non-uniqueness, run the P3 test. Before typing it as an input, run the P3 test. If it passes neither cleanly, it is a *sort*, and the sort is written down before anything else.
4. Prefer "forced by a property" over "forced by meaning" (T30, T31 are the latter). "By meaning" is a premise that has not been named yet.
5. A test is a check, not a proof. A passing candidate is "not yet killed".
6. The default report on any candidate: passes / fails / unknown, per premise. Never "this resolves it".

---

## 5. Where things live

- The chain, steps 0–18, with statuses and the 16 holes: tester `definition/proof.md`. **Read this second.**
- The premises as ruled, graded: tester `definition/premises.md`; the 10 Oct lemma list in VD-docs is its ancestor.
- The existence argument's parts: `theories_and_belief.md` (E1, E2, pressure list); `why_inertial_mass_must_be_estimated.md` and `what_forces_the_chain_to_the_estimator.md` (E3, E4); `why_the_estimator_is_the_definition.md` and `why_mass_cannot_be_undefined_2026-09-29.md` (E5); `mass_estimator_lemmas_2026-10-10.md` §0 (F1–F4, E6); `mass_readout_composition_uniqueness_2026-09-06.md` and `what_forces_the_ratio_of_means_readout.md` (the demoted composition result, P4).
- Rulings, verbatim and dated: tester `rulings.md`.
- The graveyard and the acceptance checklist: `GOAL.md` §§4–5 (project doc, 29 Oct), extracted from the tester's `docs/` folders and `proof.md` §3.
- The input-vs-choice criterion: tester `TASK.md` §§3–6.

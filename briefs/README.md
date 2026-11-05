# Briefs: work queued for other agents

One file per brief. Each brief is **one sorted fork and a named deliverable** (CLAUDE.md):
a test, a counterexample, a derivation or code, written into a named theorem or problem
file. Claude writes a brief only after sorting the fork; a question that is not sorted stays
in its problem file and does not come here. Can may add, reorder or strike any brief.

## Intake: what Claude owns (provisional, 7 Nov 2026; revised as it gets used)

Can, 7 Nov 2026: "I mean responsibility of file organization and file sprawl. If a briefs
returns with lines and lines, the understanding is that that's what you asked for. So you
merge it into main regardless and have to find ways to successfully coexist with it in main.
This is not the case for file sprawl by exploration." Also: "you don't need to check the
briefs. Codex is trustworthy to an extent and adding incorrect ideas to main isn't really a
problem since we have other ways of checking."

- **A brief is merged as delivered.** Claude does not re-verify it before merging; main's other
  checks (probes that break, strict xfails, later briefs, Can's reading) catch what is wrong.
  Claude owns where the result sits, what links to it, and that main stays navigable.
- **The size of a return is controlled when the brief is written:** the brief names the file,
  the slot, and roughly how long the write-up should be.
- **What lands unchecked says so** (Claude's addition): a theorem or result merged from a brief
  carries "Codex, Bnn; not checked by Claude" in its status until someone checks it.
- **Briefs branch from `main`.** If a brief needs material that exists only on an exploration
  branch, that material is harvested to `main` first. The 4 Nov exception (B15–B20 branched from
  `codex/ours32-consistency`) is withdrawn as of 7 Nov.
- **Explorations are input** (`../explore/README.md`): never merged as branches; left behind,
  scraped for parts, revised or iterated on before anything goes in. Harvested parts are taken
  by need (a brief, a hole in `proof.md`, a probe on main cites them) and say where they came
  from.

## The queue

| Brief | Item | Deliverable, in one line | Status |
| --- | --- | --- | --- |
| [B01](B01-robot-budget.md) | O16 | The robot's noise budget declared from calibration-side data, then the rerun | merged 2026-10-31: the calibration-side budget does not close the gap (all four answers no or mixed); the robot is a phase-one testbed (R-2026-10-31-1) |
| [B02](B02-positivity-supremum.md) | O13 | Prove S_k = k − ½ and sup E[Q] = k − 1 (the sorted numerator rests on both) | merged 2026-10-31: both proved (not attained at finite strength), T38; load-bearing for ours30 only if G18 keeps O12's meaning |
| [B03](B03-composition-preset.md) | O2, Gaps | `tools/composition_coverage.py` as a `met run` preset kind, same numbers | merged 2026-10-31: the preset reproduces the 28 Oct ours29 study exactly (same JSON as the original script); locked to that historical run |
| [B04](B04-eta-effective-readings.md) | O14 | Is J's finite-N bias a function of effective readings under the law? | merged 2026-10-31: the n_eff identity holds only for fixed nonnegative curvatures; ours30's are signed (target mean −∞); η′'s coverage reported, not adopted |
| [B05](B05-gaussian-consistency-route.md) | O4, T23 | The consistency route to the Gaussian model under the tempered evidence, written out | merged 2026-10-31: neither route forces the Gaussian as stated; the consistency route found ours30 fails L11 (G22) |
| [B06](B06-write-out-the-law-theorems.md) | T16–T19, T24, T25 | The proofs of six theorems about the law, written into their files | merged 2026-10-31: all six written; T16, T17, T25's statements need conditions (scope lines added); T17 holds for fixed N only |
| [B07](B07-t6-argument-into-the-file.md) | T6 | The per-reading factor's argument moved from the O3 record into T06.md, with its two stated limits | merged 2026-10-31: written; the brief's L7 citation for the coefficient was mis-sorted (L7 calibrates a point); the remainder attaches to G18 via T36's homogeneity |
| [B08](B08-starting-law-family.md) | O6, T28 | The family of the starting law derived in one place (L3, L4, properness) | merged 2026-11-01: the family derived (density and general-law forms); T28 written; no member selected |
| [B09](B09-demoted-theorems-bookkeeping.md) | audit, O2, O5, ladder | The nine theorems and L12 that the median took off the path, marked as such; the audit reports it | merged 2026-11-04: the nine comparison-readout theorems, L12 and O5 marked, not removed; the audit tags T21 as not load-bearing |
| [B10](B10-short-proofs.md) | T5, T12, T13, T22 | Four short proofs written into their files; T5 re-derived for the median | merged 2026-11-04: T5 (for the median), T12, T13, T22 written; statements unchanged |
| [B11](B11-scale-branch-law.md) | O14, O7 | The raw law on the "Σ up to a scale" branch, π·K^−(n/2+a−1): derived, tested, coverage against ours30; null readings | merged 2026-10-30: the law is improper on valid inputs (false premise, stopped correctly) |
| [B12](B12-admissible-tempering.md) | L24, O14 | Which temperings are cut posteriors of a mean-1 noise-size law (Bernstein test): ours30's g, the smooth form, the floor, T37 controls; and η | merged 2026-10-30 |
| [B13](B13-gleser-hwang-hypotheses.md) | O15, property N | Gleser–Hwang's hypotheses checked against the type; the corner's corollary proved; one test | merged 2026-11-04: infimum coverage 0 at every finite N; the planted ce^{±3} floors are positive (not 0), and under ours30 they erode with N (G22) |
| [B14](B14-modified-profile-likelihood.md) | O15 | Barndorff-Nielsen's modified profile likelihood for the latent pairs: the one non-hierarchical cure not yet checked; m-free or not | merged 2026-11-04: not flat (Jacobian 1/|cos(θ − θ̂)|) but not a law: poles at positive masses on anti-parallel readings; O15's table closed |
| [B20](B20-anchor-line-bookkeeping.md) | T54–T65, G24, G27 | The score-law and anchor constructions marked as benchmarks | struck 2026-11-04: moot; exploration branches are never merged (`../explore/README.md`) |
| [B15](B15-eta-with-the-conditional-score-variance.md) | O14, O18, G3, G26 | η with T45's exact conditional score variance in place of the fitted J: does it reach the target Δ/(Δ+M) at small N and under leverage (T43), and on T39's pair? | open |
| [B16](B16-the-width-without-the-starting-laws-tails.md) | O18, step 18, R-2026-09-29-3 | The quoted width computed on the evidence L alone: when it is proper, its coverage, and the fraction "unresolved"; beside ours33's π·L width | open |
| [B17](B17-what-fixes-t48s-latent-measure.md) | T48, T36, G19 | What singles out T48's latent measure (b = 0, exponent (ν+k)/2): T36 and L13, an input, or nothing | open |
| [B18](B18-composition-coverage-under-ours32-and-ours33.md) | O2, L5, L14 | Composition at the level of the law under ours30, ours32, ours33 (clause (c) as coverage), truth at c and 4c; warp before or after the pushforward | open |
| [B19](B19-l11-against-the-tls-consistency-literature.md) | O19, L11, G25 | L11's growth condition translated from Gleser, Gallo, Kukush–Van Huffel and Shklyar; which is T42/T53's threshold | merged 2026-11-07, as delivered (not checked by Claude): the literature's scalar rate is Δ_N/√N → ∞ (Gallo, Kukush–Van Huffel, Shklyar), equal to T42's and T53's threshold; T41's Δ_N ~ log N fails it; L11 unchanged (a rule, G25) |

Order is priority. Since 29 Oct briefs are cut from the holes table in `../definition/proof.md` §5 (the *derive* rows); a brief names its hole. Exploration directions, which are not briefs and are never merged, are on `../explore/README.md`. Several can run at once: each branches from `main` and writes to
different files, except the one log line each adds to `docs/STATE.md`. (The 4 Nov exception, B15–B20 branched from `codex/ours32-consistency`, is withdrawn: see Intake above.)

## For the agent running a brief

1. Read `CLAUDE.md`, `TASK.md`, the brief, and the problem file it names.
2. Branch from the current `origin/main` as `codex/<brief id>` (for example `codex/b02`).
   One brief per branch. Never branch from another agent's branch.
3. Do what the brief says, nothing outside it. The brief's "Do not" list is binding. If the
   brief turns out to be mis-sorted (its premise is false), stop, write why into the problem
   file's attempts, and finish as below: that is a result.
4. Run `python tools/audit.py` and the full suite. Both must pass; audit findings unchanged
   unless the brief says otherwise.
5. Last commit: set this brief's status line to `done, awaiting check: <one line of what
   was delivered>` on your branch. Push the branch. **Never push to `main`, never merge.**

## For Claude, at each catch-up

`git fetch`; for every brief whose `codex/<id>` branch exists: merge it to `main` (Intake
above: as delivered), resolving conflicts so it coexists with what is there; run the audit and
the full suite, which check coexistence, not the brief's mathematics; read it enough to write its
status line and one STATE log line saying what it found and which premise or hole it names; set
the status to `merged <date>`. `returned <date>: <reason>` is for a branch that cannot coexist
(it breaks main's suite, or writes outside the files its brief names). Then sort what has come up
since and add briefs.

## Status values

`open` · `done, awaiting check: …` (agent, on its branch) · `merged <date>` (Claude) ·
`returned <date>: <reason>` (Claude) · `struck <date>: <reason>` (Can or Claude)

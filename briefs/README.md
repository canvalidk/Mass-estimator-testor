# Briefs: work queued for other agents

One file per brief. Each brief is **one sorted fork and a named deliverable** (CLAUDE.md):
a test, a counterexample, a derivation or code, written into a named theorem or problem
file. Claude writes a brief only after sorting the fork; a question that is not sorted stays
in its problem file and does not come here. Can may add, reorder or strike any brief.

## The queue

| Brief | Item | Deliverable, in one line | Status |
| --- | --- | --- | --- |
| [B01](B01-robot-budget.md) | O16 | The robot's noise budget declared from calibration-side data, then the rerun | merged 2026-10-31: the calibration-side budget does not close the gap (all four answers no or mixed); the robot is a phase-one testbed (R-2026-10-31-1) |
| [B02](B02-positivity-supremum.md) | O13 | Prove S_k = k − ½ and sup E[Q] = k − 1 (the sorted numerator rests on both) | merged 2026-10-31: both proved (not attained at finite strength), T38; load-bearing for ours30 only if G18 keeps O12's meaning |
| [B03](B03-composition-preset.md) | O2, Gaps | `tools/composition_coverage.py` as a `met run` preset kind, same numbers | merged 2026-10-31: the preset reproduces the 28 Oct ours29 study exactly (same JSON as the original script); locked to that historical run |
| [B04](B04-eta-effective-readings.md) | O14 | Is J's finite-N bias a function of effective readings under the law? | merged 2026-10-31: the n_eff identity holds only for fixed nonnegative curvatures; ours30's are signed (target mean −∞); η′'s coverage reported, not adopted |
| [B05](B05-gaussian-consistency-route.md) | O4, T23 | The consistency route to the Gaussian model under the tempered evidence, written out | merged 2026-10-31: neither route forces the Gaussian as stated; the consistency route found ours30 fails L11 (G22) |
| [B06](B06-write-out-the-law-theorems.md) | T16–T19, T24, T25 | The proofs of six theorems about the law, written into their files | open |
| [B07](B07-t6-argument-into-the-file.md) | T6 | The per-reading factor's argument moved from the O3 record into T06.md, with its two stated limits | open |
| [B08](B08-starting-law-family.md) | O6, T28 | The family of the starting law derived in one place (L3, L4, properness) | open |
| [B09](B09-demoted-theorems-bookkeeping.md) | audit, O2, O5, ladder | The nine theorems and L12 that the median took off the path, marked as such; the audit reports it | open |
| [B10](B10-short-proofs.md) | T5, T12, T13, T22 | Four short proofs written into their files; T5 re-derived for the median | open |
| [B11](B11-scale-branch-law.md) | O14, O7 | The raw law on the "Σ up to a scale" branch, π·K^−(n/2+a−1): derived, tested, coverage against ours30; null readings | merged 2026-10-30: the law is improper on valid inputs (false premise, stopped correctly) |
| [B12](B12-admissible-tempering.md) | L24, O14 | Which temperings are cut posteriors of a mean-1 noise-size law (Bernstein test): ours30's g, the smooth form, the floor, T37 controls; and η | merged 2026-10-30 |
| [B13](B13-gleser-hwang-hypotheses.md) | O15, property N | Gleser–Hwang's hypotheses checked against the type; the corner's corollary proved; one test | open |
| [B14](B14-modified-profile-likelihood.md) | O15 | Barndorff-Nielsen's modified profile likelihood for the latent pairs: the one non-hierarchical cure not yet checked; m-free or not | open |

Order is priority. Since 29 Oct briefs are cut from the holes table in `../definition/proof.md` §5 (the *derive* rows); a brief names its hole. Several can run at once: each branches from `main` and writes to
different files, except the one log line each adds to `docs/STATE.md`.

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

`git fetch`; for every brief whose `codex/<id>` branch exists: check it out in its own
worktree, install from `requirements-lock.txt`, run the audit and the full suite, read what it
wrote into the named file against the brief's questions and "Do not" list. Then either merge
it to `main` with the status set to `merged <date>`, or set `returned <date>: <reason>` on
`main` and leave the branch. Then sort what has come up since and add briefs.

## Status values

`open` · `done, awaiting check: …` (agent, on its branch) · `merged <date>` (Claude) ·
`returned <date>: <reason>` (Claude) · `struck <date>: <reason>` (Can or Claude)

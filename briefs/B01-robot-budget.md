# B01: the robot's noise budget, declared from calibration-side data, then the rerun

**Status:** done, awaiting check: calibration-only budget and bootstrap declared; all 48 rows and four answers delivered, with the gap remaining
**Item:** `problems/O16.md`. **Branch:** `codex/b01`. **Written:** 2026-10-28, Claude.

Read `problems/O16.md` first, especially "28 Oct: what the robot's shared error is, and is not".

## The sorted fork (Claude's, provisional; do not re-sort it)

The robot readings' error covariance is an input (`TASK.md` §2), declared by Claude in
`met/robot.py`. Its white part, rest variance / samples in the window, assumes independent
samples; the rest recordings are autocorrelated (tau_int up to about 19 samples), so it is
understated 2–19 times per axis. Separately, the dataset's `R_FA` is nominal: the static rows
give a 1.153° force–IMU rotation. Both are calibration facts available before any dynamic
reading. So: declare the budget from calibration-side data only, then rerun. This fixes an
input, not the definition. (The data are Skrede's public dataset, Zenodo 10.5281/zenodo.11096791.)

## Deliverable

1. `met/robot.py`: a second declaration beside the existing one. The defaults of `budget`,
   `windows` and `static_heldout` stay exactly as they are; no recorded number may move.
   - White part per channel and axis: rest variance × tau_int / samples in the window, tau_int
     from `0-steady-state_wrench.csv` and `0-steady-state_accel.csv` by the
     initial-positive-sequence estimator (`_tau` in `tests/test_robot_budget_checks.py`). State
     the method in the docstring; do not change it after seeing dynamic results. Report tau per
     axis with a block-bootstrap spread (the rest files are one second long).
   - Optional rotation: rotate the IMU channel by R from F_fts ≈ s R A + b on the static rows
     (as in that test file). Report R and its bootstrap covariance; do not propagate it (the
     type has no shared slot yet).
   - Floor: keep the option unchanged. Correct the docstring sentence calling it "the same error
     in every reading": it is the static fit's zero-mean residual (O16, 28 Oct, item 1).
2. `tools/robot_run.py`: a flag selecting the declaration, and one table over runs
   {1-baseline, 2-vibrations} × windows {0.05, 0.2, 1 s} × white {as now, autocorrelated} ×
   {floor, no floor} × rotation {nominal R_FA, static R}, under the current equation, ours30
   (`tools/robot_run.py`'s default since 29 Oct). Per row: N, median, 95% interval, covers 0.9323,
   eta, min and median g. The numbers quoted in the questions below are ours29's; give ours29's
   median beside each row so they stay comparable.
3. Tests: the tau estimator on a synthetic AR(1) with known tau = (1+ρ)/(1−ρ); the new white
   part equals the old when tau = 1; the rotation fit recovers an injected rotation on synthetic
   static rows. Existing robot tests unchanged and passing; strict xfails left as they are.
4. Into `problems/O16.md`'s attempts, and one line in `docs/STATE.md`'s log (nothing on its
   main page): the table, and answers to these four questions, fixed before running:
   1. Does the floor/no-floor gap in the baseline median (1.4% now) fall below 0.5%?
   2. Do the medians at 0.05, 0.2 and 1 s agree within their intervals?
   3. Does the 95% interval cover 0.9323 at all three widths, with and without the floor?
   4. How far do the misfit weights rise (median g about 0.012 without the floor now)?

   A "no" anywhere is the finding. Report it as it comes.

## Do not

- Change `met/full.py`, the equation, anything in `definition/`, or `rulings.md`.
- Fit or tune anything on runs 1–3: every declared quantity comes from the rest files or the
  static rows.
- Choose among declarations by which one covers. Every row goes in the table.
- Add a cross-reading (series) covariance: that is O16's next step, and how g and eta apply
  under it is unsorted.
- Open a new problem file, or write study narration on STATE.md's main page.

# Exploration board

Directions for open exploration by other agents (Codex and the like). The brief board is `../briefs/README.md`:
briefs are sorted forks with a named deliverable and are merged; explorations are open, and are **never merged**.
Claude reads an exploration branch when there is a reason to and takes specific results across by hand.

## Rules

- Branch as `explore/<id>-<slug>` (for example `explore/e2-width-class`). Add the branch to the board's
  Branches column.
- Put a `SUMMARY.md` at the branch root, at most a screen long: what broke, which premise of the definition
  each break names, and which files hold each result.
- Do not redo what a brief is already doing: see the "Not to redo" column.
- Before working any fork, name its standard mathematics and cite it (`../GOAL.md`, standing rule 6); if the
  result exists in the literature, cite and translate it rather than re-deriving it.

## The board

| Id | Direction | What an answer would decide | Not to redo | Branches | Status |
| --- | --- | --- | --- | --- | --- |
| E1 | Break ours33 (`met/ours33.py`; `problems/O18.md`, 4 Nov) beyond its listed holes | which premise its next break names (G3, G19, the D bin) | holes 1–4 are B15–B18 | — | open |
| E2 | Is the finite-N width unique? The class of quoted-width rules that are Z-safe, free of π where the evidence is proper, calibrated in the limit and continuous in the data: one member or a family; if a family, what extra premise selects a member, and is it a property of mass or a reporting convention | what the proof is a proof of (`definition/proof.md` step 18, G3) | B16 measures one candidate; E2 asks about the class | — | open |
| E3 | Per-reading signal weighting without declaring a size law: ours28 weighted each reading by $\hat\beta_i=1-k/h_i^2$ from a declared law, and η is the series-level version recovered without one (`problems/O14.md`, 4 Nov). A per-reading statistic that is Z-safe, borrows nothing (L19) and declares no size law (L10), or a proof that none exists | the shared-vs-per-reading fork; whether η's form can be per reading | B15 tests one series-level estimate | — | open |
| E4 | The weak-signal centre: is any Z-safe, π-free point better than the raw median, which carries the flat latent measure's tilt? | ours33 hole 6; O15's corner | B13 proved the coverage limit (Gleser–Hwang) | — | open |

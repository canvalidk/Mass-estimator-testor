# B18: composition at the level of the law under ours32 and ours33

**Status:** done, awaiting check: later-equation preset arms and 400-replicate composition tables, centred worlds, three ours33 warp procedures and tail-safe propagation; weak sum coverage misses and order sensitivity recorded, 38 new checks.
**Item:** `problems/O2.md` (R-2026-10-28-1, clause (c)); L5, L14; `presets/composition_coverage.yaml` (locked to ours29, B03); T51; `problems/O18.md` (4 Nov, hole 4).
**Branch:** `codex/b18`, from the head of `codex/ours32-consistency` (needs `met/ours32.py`, `met/ours33.py`).
**Written:** 2026-11-04, Claude.

Read O2's 28 Oct sections and R-2026-10-28-1 in `rulings.md`, B03, `met/composition_coverage.py`, T51's statement.

## The sorted fork (Claude's, provisional; do not re-sort it)

R-2026-10-28-1 restated L5 and L14 at the level of the law: (a) exact inputs add exactly; (b) the lump's
law is the law of $m_1+m_2$ computed from the parts; (c) a new measurement of the lump is consistent with
(b), checked as coverage. The check has been run once, for ours29 (28 Oct; B03 locked the preset to that
run). It has not been run for ours30, ours32 or ours33, and ours33 changes what "the law" is: the quoted
law is the warp of the raw law, and T51 proves that a universal law-only warp commuting with
independent addition is the identity. So under ours33 clause (c) is a sharper question: **does the
pushforward of the quoted laws cover the composite's own quoted law at the stated rate, and does it matter
whether the warp is applied before or after the pushforward?** Clause (c) asks for consistency, not
equality (O2, 28 Oct, the coarsening bridge), so this is a coverage study with a stated procedure, not a
theorem. **Standard name:** propagation of distributions (GUM Supplement 1) against direct measurement
of a derived quantity; T51 is the statement that no fixed recalibration commutes with convolution.

## Deliverable

1. **Extend the preset kind** `composition_coverage` to take `equation: ours29 | ours30 | ours32 | ours33`
   (with `nu_eff` for the last two), keeping the ours29 run byte-identical (B03's test must still pass).
   For ours33, two procedures for the composed law, both reported: warp the parts then push forward;
   push the raw laws forward then warp with a stated η (say which: the parts' η's combined how, and
   why; if no combination is principled, say so and report both candidates).
2. **Run** the 28 Oct cells (seed 2810, 400 replicates per size) for ours30, ours32 and ours33, lumping and
   common force, strong and weak signal: coverage of the computed law against the lump's new
   measurement, and the overlap of the two laws, as O2's 28 Oct table. Truth at $c$ *and* at $4c$ (O18, 4 Nov:
   the truth-at-$c$ design rewards a pull toward $c$).
3. **Write** the table into O2's attempts with the two-procedure comparison; state in one paragraph what
   T51 does and does not say about the result (it is about a universal operator; clause (c) is a check on
   one supplied joint law).
4. Tests in `tests/test_composition_coverage.py` (the existing file; add, do not change the locked check);
   one line in `docs/STATE.md`'s log.

## Do not

- Change clause (c), L5, L14 or R-2026-10-28-1's wording.
- Adopt either warp order as the equation's; report both.
- Add a dependence model or a joint prior between the parts (T54's construction is a benchmark, not an
  input of this check).
- Extend the anchor/score-law constructions (T56–T65).

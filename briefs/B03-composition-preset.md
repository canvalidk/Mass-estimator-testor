# B03: the composition coverage study as a `met run` preset kind

**Status:** open
**Item:** `problems/O2.md` (28 Oct, later) and `docs/STATE.md`, Gaps. **Branch:** `codex/b03`.
**Written:** 2026-10-29, Claude.

## The sorted fork

None to sort: this is code. `tools/composition_coverage.py` (seed 2810) produced the coverage
table behind R-2026-10-28-1 (O2, "28 Oct (later): composition at the level of the law,
checked"). The repo's rule is that every study is a preset with its seed; this one is a script.

## Deliverable

1. A preset kind `composition_coverage` in `met/` (follow how `pooled_scale` and
   `readout_candidates` were added: `met/preset.py`, `met/runner.py`, `met/__main__.py`), and
   `presets/composition_coverage.yaml` with seed 2810, 400 replicates, the three signal sizes
   and the parts, lump and reduced mass of the script. Output under `results/` like the other
   kinds. `--workers` allowed if the output is identical for any count; otherwise serial only.
2. The preset reproduces the script's numbers: every entry of O2's 28 Oct table (coverage,
   overlap, the median point defects) to its printed digits. A test of that under `MET_SLOW`,
   and a fast test on a few replicates that the preset and the script agree exactly.
3. `tools/composition_coverage.py` becomes a thin wrapper of the preset kind, or is removed
   with `tests/test_composition_law.py` pointed at the new module; either way that test keeps
   passing unchanged in what it asserts.
4. Remove the Gaps line "(28 Oct) `tools/composition_coverage.py` is a seeded script, not a
   preset kind" from `docs/STATE.md`, add one log line, and change O2's pointer to the preset.

## Do not

- Change the study's computation or any number. If a number does not reproduce, that is a
  finding: a strict xfail naming both numbers, and stop there.
- Change `met/full.py` or anything in `definition/`.

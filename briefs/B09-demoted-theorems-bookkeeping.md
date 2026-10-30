# B09: the audit and the problem files after the median

**Status:** done, awaiting check: comparison-readout scope marked, O5 parked, and the audit retains and labels T21
**Item:** `tools/audit.py`, `problems/O2.md`, `problems/O5.md`, `definition/ladder.md`,
`definition/theorems/README.md`; `definition/proof.md` §4 "What is no longer load-bearing", hole G12.
**Branch:** `codex/b09`. **Written:** 2026-10-29, Claude.

## The sorted fork

None: bookkeeping. Since the point became the law's median (28 Oct) these are no longer on the
path from the premises to ours30: T3, T4, T9, T10, T11, T14, T15, T21, T29, the premise L12, and
the open item O5. They stay true of the cloud readout, kept for comparison. The files still
present them as the definition's.

## Deliverable

1. `definition/theorems/README.md`: a column or marker ("comparison readout only, since 28 Oct")
   on those nine rows; T8's From line and the others' are not changed.
2. Each of the nine theorem files: one line at the top of its Status saying it concerns the
   cloud readout, which is kept for comparison and is not the point since 28 Oct
   (R-2026-10-28-1, T35), with a pointer to `proof.md` §4. Statements and proofs untouched.
3. `tools/audit.py`: the three "untested theorem" findings (T6, T21, T27) are unchanged in
   kind, but T21 should be reported as "not load-bearing" rather than owed; add that state (read
   from the theorem file's Status line) and make `--table` show it. Do not drop any check.
4. `problems/O5.md`: status "parked: not needed by ours30's proof; needed only if the ratio of
   means returns". `problems/O2.md`'s requirements table: a note that R-a and R-d constrain the
   cloud readout only. `definition/ladder.md` §3: the two rows "the series' size inside that
   ratio" and "the weight's noise part" marked not load-bearing for the same reason.
5. One line in `docs/STATE.md`'s log.

## Do not

- Delete or rewrite any theorem, premise or problem file; mark, do not remove.
- Change `met/`.

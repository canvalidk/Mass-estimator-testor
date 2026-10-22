# docs/: the mass-estimator work surface

The documents of the mass-estimator investigation live here: equation states
(ours25b to ours29), the signal-versus-noise work, the literature review and
papers, the unsanctioned Codex records of 6 to 8 September, the 21 September
equation specification with its code, and the September work-state records.
They moved here from VD-docs (`Newton-analysis/_dscn_mass_estimator/`) on
**2026-10-22**, on his ruling that VD-docs holds made decisions only and that
doing this work there had been a mistake.

## Layout

`docs/` mirrors the old discussion root, so the relative paths the records use
between themselves (`ours27/...`, `noise_vs_signal/...`) still resolve. Paths
that start `_dscn_mass_estimator/` or `Newton-analysis/` point at VD-docs and
are stale for these files; the filename is the identifier.

- `README_of_the_discussion_before_2026-10-22.md`: the map of this material as
  it stood in VD-docs, with every relation recorded there. Start here.
- `SHA256SUMS`: hash of every file as it left VD-docs (`sha256sum -c SHA256SUMS`
  from this folder).
- `vd-docs-context/`: read-only copies of what stayed in VD-docs: the 12
  decision records and the current drafts of the lemma and equation notes.
- `mass_estimator_reproducibility_inventory_2026-10-01.md`: what in these
  records can still be rerun, and with what.

The research scripts here (`ours27/ours28_reference*.py`,
`noise_vs_signal/ours29_reference*.py`, `noise_vs_signal/ours28b/*.py`,
`mass_estimator_equation/*.py`) are the originals of the byte copies in
`../reference/`. The copies there remain the ones the tests check.

## What goes back

Only decisions go back to VD-docs, through `../outbox/VD-docs/` (its README
says what qualifies). The VD-docs record of this move, with every file and its
hash, is `outbox/Mass-estimator-testor/HANDOVER_2026-10-22.md` in VD-docs.

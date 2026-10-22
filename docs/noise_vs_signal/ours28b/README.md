# ours28b/ — the ours28b candidate

**2026-09-28 (night).** Written by Claude on the user's instruction: "make a
ours28b folder and put the appropriate files that should be there there."

## What this folder is

The working folder for ours28b, a candidate equation named by the user on 28
September. Like `ours26/` and `ours27/`, it is a grouping inside
`_dscn_mass_estimator`, not a nested discussion.

**ours28b in one line:** ours27 with the flat candidate measure replaced by one
Gaussian candidate-size scale per declared noise pair, shared over the series,
under the arcsine scale law.

It is not the ours28 of `../ours27/ours28_candidate_push_scale_as_nuisance_2026-09-28.md`,
which has one scale per reading. Those ours28 records stay in `../ours27/`.

## Here now

- `ours28b_one_scale_per_noise_pair_2026-09-28.md`: names ours28b. It covers
  the equation (types, classes, law, candidates, readout, carried state,
  reductions), the checks, the comparison with ours27, what is kept and given
  up, the premises not yet ruled, and why it is still a candidate.
- `ours27_one_scale_per_noise_pair_2026-09-28.md`: the derivation, moved here
  from `../ours27/` when this folder was made, because ours28b is its
  equation. It covers the calibration requirement, the per-reading
  impossibility, the exact per-class result, and pairs over ratios.
- `ours28b.py`: reference implementation (numpy, scipy). `ours28b(readings,
  c, flat=False)` returns the readout, the median, the 95% interval, the law,
  $A(m)$ and the per-pair classes. `flat=True` gives ours27 on the same input
  type.
- `ours28b_checks.py`: the checks behind the record's tables. It runs in
  about five minutes.
- `ours28b_no_borrowing_between_readings_2026-09-28.md` (late night): **ours28b withdrawn.** The user ruled that
  the definition borrows nothing between readings (readings may come from any
  regime in any combination; nothing learned about pushes carries over). That
  excludes ours28b, a learned population of per-reading scales and a shared
  direction tensor, and keeps ours27 and ours28. The record adds a result: no
  per-reading rule is never overconfident, because a pure-noise reading adds
  score noise and no curvature. So the per-reading scale law's shape is a
  declared trade, not an honesty-forced choice.
- this README.

The reference implementation is a working copy for the record. The tester
repo (`canvalidk/Mass-estimator-testor`) is where reproducible studies live;
ours28b is not built there yet.

## Where things stand

**Withdrawn, 28 Sept (late night):** ours28b borrows push scale across
readings, which the user's ruling forbids (see the no-borrowing record). The
live line is ours28 in `../ours27/`. The notes below describe ours28b as
it stood before the ruling.

- **Forced, with one noise pair,** given the premises of the derivation
  record: the structure of ours28b, up to the shapes of the starting law and
  the scale law.
- **Open:**
  - the jump when two noise pairs become equal (record §5.1);
  - general covariances (§5.2);
  - the two shapes (§5.3);
  - the readout's own interval, the matrix rerun and the tester build (§5.4).
- **His rulings owed:** calibration as a premise; dropping the strong
  zero-reading property.

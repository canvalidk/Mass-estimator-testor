# E10: null limits, prior selection, and consistency

Branch: `explore/e10-null-limit-prior`; exploration only, never merged.

- The averaged T42 null update has no full-support fixed prior on positive mass.
  Whole-circle uniformity does not survive per-result positive restriction.
  [Prior-selection proof](explore/e10-null-limit-prior.md).
- Stronger **null washout** (each returned law tends to pi) conflicts with
  consistency for every arbitrarily slowly identifying triangular array.
- With **permutation invariance explicitly assumed**, the conflict holds on
  one fixed infinite signed history too, despite Delta_N/sqrt(N) -> infinity.
  It reaches the current median's point consistency, not only law concentration.
- An ordered-history control achieves washout and fixed-history consistency,
  but fails permutation and exact-zero invariance. A smooth attenuation preserves
  those invariances but needs an extra growth margin and an unforced schedule.
  [Compatibility proof and controls](explore/e10-null-washout-consistency.md).
- Neither averaged neutrality nor strong washout is already established by Z.
  No prior, production correction or premise is adopted; full L19 compliance
  of the controls is not established. The null requirement's meaning is the fork.
- Checks: `tests/test_null_limit.py`, `tests/test_null_washout.py`.
  Latest focused run: 117 passed; audit unchanged. Independent reviews recorded.

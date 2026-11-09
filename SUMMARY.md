# E10: the null limit does not select a full-support prior

Branch: `explore/e10-null-limit-prior`, from main; exploration only, never merged.

- T42's pure-noise limit is a random angular tilt of the supplied prior.
  Testing whether its average equals that prior gives **no full-support
  fixed point on positive mass**, for every nonzero finite strength multiplier.
- Uniform angle is the unique full-support fixed point on the whole circle.
  Restricting each returned law to positive mass before averaging breaks it:
  random normalization and averaging do not commute.
- Half-Cauchy at c=s: the average law shifts about 0.0955 in total variation.
  Reciprocal symmetry preserves the average's centre, not its shape.
- Raw Gram tempering and raw Bartlett tempering share this scalar null limit.
  If instead null tilt vanishes, every prior survives; no member is selected.
- The failed combination is an **added averaged-null requirement**, positivity,
  full support and persistent isotropic evidence. Existing Z only establishes
  exact observed-zero neutrality. No premise, prior or production rule changed.
- Derivation and sources: [E10](explore/e10-null-limit-prior.md).
  Deterministic checks: `met/null_limit.py`, `tests/test_null_limit.py`.
  Validation: 14 new checks; 101 focused regression checks passed; audit unchanged.
  O6/O15 remain open; the next issue is the meaning of null neutrality.

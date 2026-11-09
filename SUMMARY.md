# E6: declared structure between readings

Exploration by Codex, 9 November 2026, from main `80cb9c1`.

- Ours28b prompted explicit L19: it pooled a signal-size scale. Known
  per-reading direction was explicitly allowed; ours28's population law
  was rejected separately on 30 September.
- A common frame supplies no signal relation. A common unknown line does:
  two datasets with identical per-reading Grams have different constrained
  mass profiles. This uses no signal-size population, but fails strict L19.
- That shared-line profile passes Z. Conversely, an observed zero changes
  mass evidence when the true vector is declared equal across acquisitions.
  Thus Z alone cannot classify permitted sharing or its neutrality.
- For block-independent Gaussian errors and fixed linear signal constraints,
  own-reading profile factorization holds iff the constraint is a Cartesian
  product. This includes relative evidence between distinct masses.
- Equal-vector acquisitions reduce exactly to separate channel means with
  their covariance; unrestricted averaging loses evidence. Reusing one
  acceleration observation requires retaining observation identity.

Proofs, history, sources and unresolved input-permission question:
[explore/e6-declared-structure.md](explore/e6-declared-structure.md).
Research helpers: `met/declared_structure.py`; independent checks:
`tests/test_declared_structure.py`. Production equations/premises unchanged.

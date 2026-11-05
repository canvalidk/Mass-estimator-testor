# E5 — B17's latent-scale separation

User-authorized exploration, 5 Nov 2026. Branch `explore/e5-latent-scale`,
based on completed B17. Findings are Codex's proposed measurement-bin sort.

- Normalization plus declared mean precision forces b=0 within power tilts.
  A positive cubic tilt preserves both covariance and mean precision while
  changing relative mass evidence. Instrument moments leave this freedom.
- Proper diffuse signal models can retain the original unconditional noise marginal
  and have different integrated-likelihood limits. Their residual sampling
  distribution retains the original noise law at every cutoff.
- A route to the Student exponent requires calibration conditional on every
  true signal, preserved in a common diffuse approximation. Another specifies
  the normalized residual experiment as evidence. L10/L18/L19 and R1/R2 leave
  that bridge open in G19.

Proofs and source references are in `definition/theorems/T48.md`, E5 section.
`problems/O18.md` records the premise diagnosis.
`tests/test_latent_scale_exploration.py` has 27 new checks. All 63 tests across
that file and B17's latent-measure/residual-mixture files pass. Audit findings
match the B17 baseline. Production code has no edits.

Per `explore/README.md`, this branch is for extracting specific results.

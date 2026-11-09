# E6: declared structure and the state retained between operations

Exploration by Codex, 9 November 2026, from main `80cb9c1`.

- Ours28b pooled signal size. A common frame supplies no signal relation.
  Gaussian factorization and conditional GLS results name what L19 excludes.
- Joint calibration passes Z; shared geometry changes P(sum=3) from 0.4733
  to 0.5. The Gaussian joint still fails T30 and literal L19.
- Born form admits arbitrary kernels; positive mass rays do not resolve
  a POVM, and real probes leave the pure complex completion undetermined.
- Allowed future operations fix the required state. Mass-only operations
  need only the mass marginal; addressed recalibration can require the full
  joint association. Even both marginals and the next-event law can agree
  while calibration gives mass posteriors (.61488,.38512) and the reverse.
- A normalized raw event distinguishes the old mass-neutral histories.
  This preserves an amplitude-indexed family, not an unconditional forecast.
- At a common null, one fixed normalized reduction removing all amplitude
  dependence also removes mass information. Candidate-wise pivots escape it.
- Fixture input permission remains open; no production/premise change.

Proofs: [structure](explore/e6-declared-structure.md), [Born bridge](explore/e6-born-bridge.md),
[predictive state](explore/e6-predictive-state.md). Research helpers: `met/declared_structure.py`.
Checks: structure (53), Born bridge (18), predictive state (24); 188 relevant checks pass.

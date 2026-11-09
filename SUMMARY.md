# E6: declared structure and the Born bridge

Exploration by Codex, 9 November 2026, from main `80cb9c1`.

- Ours28b pooled signal size. A frame alone supplies no signal relation;
  fixed linear constraints factor iff they form a product (L19).
- GLS gives an exact conditional reduction; an unknown shared line has
  no lossy sufficient reduction for the full signal experiment.
- Joint calibration passes Z; mass-only state loses future information.
  Shared geometry changes P(sum=3) from 0.4733 to 0.5. The Gaussian law
  still fails T30 and literal L19; addressed input permission remains open.
- Born bridge: positive mass rays cannot resolve a POVM, and real mass
  probes cannot identify T70's pure complex completion uniquely.
- E6 admits diagonal Born/Lüders form for any bounded positive kernel.
  The form cannot select the estimator. Binary normalization and covariance
  also admit non-Born rules; richer operational assumptions need justification.
- Next fork: which reductions preserve the declared experiment's outcome
  probabilities without supplying a signal population? No premise is adopted.

Proofs: [declared structure](explore/e6-declared-structure.md),
[Born bridge](explore/e6-born-bridge.md). Research helpers: `met/declared_structure.py`.
Checks: `tests/test_declared_structure.py` (53), `tests/test_born_operational_bridge.py` (18).
Production equations and premises unchanged.

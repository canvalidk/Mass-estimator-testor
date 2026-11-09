# E6: declared structure between readings

Exploration by Codex, 9 November 2026, from main `80cb9c1`.

- Ours28b prompted L19 by pooling signal size. Known per-reading direction
  was allowed. A common frame alone supplies no relation between signals.
- Fixed linear constraints factor iff they form a product and admit exact
  GLS reduction; its constraint provenance can still rely on Newton II.
- An unknown common line has no lossy sufficient reduction for the full
  signal experiment. A noisy direction installed as exact can give 0.8 at truth 2.
- Retaining calibration in joint (mass,direction) gives a proper, Z-safe
  research law. A mass-neutral history still changes later mass updates.
- Shared geometry couples inferred masses despite independent errors:
  P(sum=3) is 0.5 versus 0.4733 from multiplying even the correct marginals.
- Resetting a direction or counting calibration twice changes the experiment.
  Fixture identity, persistence and evidence identity must travel together.
- This Gaussian construction still fails T30 and literal L19. Permission for
  addressed calibration evidence remains open; no new estimator is adopted.

Proofs, history, sources and the unresolved input-permission question:
[explore/e6-declared-structure.md](explore/e6-declared-structure.md).
Research helpers: `met/declared_structure.py`; independent checks:
`tests/test_declared_structure.py` (53 checks). Production equations/premises unchanged.

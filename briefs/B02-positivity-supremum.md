# B02: the exact supremum S_k of the positive-mass misfit

**Status:** done, awaiting check: exact isotropic suprema derived, finite attainment classified, and deterministic quadrature checked
**Item:** `problems/O13.md`. **Branch:** `codex/b02`. **Written:** 2026-10-29, Claude; updated the same day after the sort.

Read `problems/O13.md`: "27 Oct: the sign fork sorted" and "27 Oct pressure check: exact null
moments and a finite counterexample to the numerator". T32 has Q_+.

## The sorted fork

That the sign of F·a enters the misfit weight is forced (L1's inequality, R2 on the ladder;
O13, 27 Oct). Since 29 Oct the numerator is sorted (provisional, O13 "29 Oct: the numerator
sorted"): it is the supremum over sizes of the expected misfit, and a seeded grid puts that at
**k − ½ for Q₊** and **k − 1 for the sign-blind Q**. Known by proof (O13, Codex): k − ½ ≤ S_k ≤ k.
This brief turns the numerical values into proofs or refutations.

Setup, as in O13: whitened x = R e₁ + ξ, y = b e₁ + ζ, with ξ, ζ ~ N(0, I_k) independent and
R, b ≥ 0. A = |x|², B = |y|², C = x·y; Q_+ = λ₋ when C > 0, min(A, B) when C ≤ 0.
f_k(R, b) = E Q_+, S_k = sup over R, b ≥ 0 of f_k (limits included).

## Deliverable

1. f_k(R, b) for k = 1, 2, 3 on a grid fine enough to locate its maximum, by deterministic
   quadrature with stated error bounds. Monte Carlo only as a cross-check, seeded.
2. S_k exactly, or the sharpest bounds you can prove, and where it is reached: an interior
   point, or a limit (R → ∞ with b → 0, as O13's lower bound uses, or another). In particular:
   is S_k = k − ½, or is there a point with f_k > k − ½?
   Also, for the sign-blind Q: is sup over R, b of E Q equal to k − 1 (reached only as both
   channels grow)? The sort rests on both.
3. Use and check the swap symmetry f_k(R, b) = f_k(b, R) (Q_+ is symmetric in x and y).
4. Into `problems/O13.md`'s attempts: the derivation, the grid's maximum and the answer to 2;
   tests in a new `tests/test_positivity_supremum.py`; one line in `docs/STATE.md`'s log.

## Do not

- Change `met/full.py`, g, or the equation; choose, propose or adopt a numerator.
- Extend to correlated channels or partial subspaces; isotropic, full subspace only.
- Write study narration on STATE.md's main page.

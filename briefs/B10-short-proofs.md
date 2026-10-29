# B10: the short proofs: T5 (for the median), T12, T13, T22

**Status:** open
**Item:** `definition/theorems/T05.md, T12.md, T13.md, T22.md`; `definition/proof.md` holes G7 and G15.
**Branch:** `codex/b10`. **Written:** 2026-10-29, Claude.

## The sorted fork

None: four theorems whose proofs are a few lines each from the law in `definition/equation.md`
(ours30), not yet written in their files.

## Deliverable

1. **T5, the swap, for the median.** The old derivation went through T4 (the ratio of means),
   which is no longer the point. New proof: the law is invariant under swapping force and
   acceleration with their covariances together with $\phi\mapsto\pi/2-\phi$ (show it for $H_i$,
   $Q_+$, $\kappa$, $\eta$ and $\pi$ by L4), so the swapped law is the image of the law under
   $m\mapsto c^2/m$, and the median of the image is $c^2/$median. State that the premise is
   R-2026-09-24-1 / R-2026-09-07-1 (force and acceleration are equal information) and that T5 no
   longer rests on T4; update the From line.
2. **T22, full covariance.** $H_i$ is invariant when a reading's data, covariance and subspace
   are rotated together (two lines), and a diagonal model is not (one counterexample with
   numbers, from the existing test).
3. **T12, combination.** The law's state adds across readings (the Gram matrices, weighted by
   $g_i$), and the tempered law of a union is not the product of the tempered laws of its parts
   because $\eta$ is recomputed on the union; one worked two-reading example with numbers.
4. **T13, totality.** With a proper $\pi$ and every $H_i$ bounded, the law is a law and its
   median is in $(0,\infty)$ for every measured input, the all-zero experiment included; say why
   $H_i$ is bounded (the explained energy is at most $E_i$).
5. Existing tests keep passing; one line in `docs/STATE.md`'s log; each file's History.

## Do not

- Change any Statement, `met/`, or `equation.md`.
- Touch T6, T30, T31 (holes G2, G13, G3: Claude's).

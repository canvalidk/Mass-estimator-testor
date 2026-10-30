# B07: T6's argument, from the O3 record into T06.md

**Status:** done, awaiting check: T06 proof and coefficient tests written; L7 lacks the archive's width condition, so independent coefficient forcing remains conditional (G2); G17/G13 preserved.
**Item:** `definition/theorems/T06.md`; `definition/proof.md` step 3, hole G2. **Branch:** `codex/b07`.
**Written:** 2026-10-29, Claude.

**Frame (added 29 Oct):** `definition/theorems/T36.md` (2) is what $e^{H_i/2}$ is geometrically; cite it where it
helps. Item 3 (the coefficient) is unaffected.

## The sorted fork

None: the argument exists and is not in the file. T6 (the per-reading factor is the flat Gaussian
endpoint, $E_i(m)\propto e^{H_i(m)/2}$) is forced from L8's conditional form, L9, L10, L7 with L13,
and T7. The argument is written in
`docs/vd-docs-context/o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md`
§2 (Herschel–Maxwell for the family; L10 for the endpoint; L7/L13 for the coefficient and the
improper members; T7 for noise units). T06.md's "Proof or support" says "carried from the records".

## Deliverable

1. T06.md's "Proof or support": the argument in full, in the file's own words, citing each
   premise by its row. Keep the Statement unchanged. State Herschel–Maxwell as used (a density on
   $\mathbb R^k$, $k\ge2$, that is a product of one-dimensional factors and rotation-invariant is
   Gaussian) with a reference, and give the one-page proof or a precise citation.
2. The two things `proof.md` step 3 says the argument does not cover, written as explicit
   conditions in the proof, not resolved: (i) $k=1$ is carried by T25 and rests on "a reading
   declared $k=1$ is the 1D reading", which no premise row states; (ii) L8's conditional form
   conditions on $Q_i$, so whether the weight may depend on $Q_i$ is settled later, by T30.
3. The coefficient: a test that the factor is $e^{H_i/2}$ and not $e^{H_i/(2(1-2c))}$ for $c\ne0$,
   by checking that the law's width at strong signal is the Gaussian propagation of $\Sigma_i$
   (L7, L13), in `tests/test_theorems.py` beside T6's test.
4. One line in `docs/STATE.md`'s log; T06.md's History.

## Do not

- Change the Statement, any premise, or `met/`.
- Resolve (i) or (ii): they are hole G17 and G13 in `proof.md`, sorts, not derivations.

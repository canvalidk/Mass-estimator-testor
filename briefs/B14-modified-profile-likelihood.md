# B14: Barndorff-Nielsen's modified profile likelihood for the latent pairs

**Status:** done, awaiting check: Jacobian is nonconstant (exact scalar control; general expected-score approximation), but the global law fails on allowed inputs; calibration branch stopped, with eight obstruction/invariance checks.
**Item:** `problems/O15.md` (30 Oct, later: the table of non-hierarchical cures). **Branch:** `codex/b14`.
**Written:** 2026-10-30, Claude.

Read O15's two 30 Oct sections, T36, and `definition/proof.md` step 3 (the flat candidate weight).

## The sorted fork (Claude's, provisional; do not re-sort it)

The law is the profile likelihood of $m$ over $N$ incidental push parameters $t_i$ (arc-length along the
line $u_m$, whitened). Every standard non-hierarchical repair of a profile likelihood was found to coincide
with the flat law here (Cox–Reid: $j_{tt}=1$ is $m$-free) or not to exist (conditional, marginal), **except
one, unchecked:** Barndorff-Nielsen's modified profile likelihood (1983, *Biometrika* 70:343–365),
$$L_{mp}(m)=L_p(m)\,\lvert j_{tt}(m)\rvert^{-1/2}\,\Big\lvert\frac{\partial\hat t}{\partial\hat t_m}\Big\rvert,$$
whose Jacobian term can depend on $m$ even when $j_{tt}$ does not. If it does here, it is a candidate
latent measure that does η's job in principle without a size law; if it reduces to the flat law too, the
table in O15 is complete and η's job is forced.

## Deliverable

1. **The derivation**, into O15's attempts: $L_{mp}$ for one reading in the whitened model
   $y=tu_m+\varepsilon$, $\varepsilon\sim N(0,I_{2k})$, with the ancillary chosen as the standard theory
   requires (say which, and why; if no exact ancillary exists, use the Severini / Cox–Reid approximation
   and say so). State whether the Jacobian is $m$-free. Then the series: the product over readings.
2. **If it is not $m$-free:** the modified law in closed form; its Bartlett identity at the true $m$,
   checked by simulation against the flat law and ours30 on O15's cells (3D, $N=20$, signal 0, 0.5, 1, 2;
   seeded); Z (a zero reading leaves it alone?), L3 (units), the swap.
3. **If it is $m$-free:** say so, with the two-line reason, and the table in O15 is closed.
4. Tests in `tests/test_modified_profile.py`; one line in `docs/STATE.md`'s log.

## Do not

- Change `met/`, any premise row, or η. A candidate measure is reported, not adopted.
- Use a size law or anything estimated across readings (L10, L19).
- Touch the positive-mass restriction or ν_eff (ours31's reading level); the flat Gaussian per reading only.

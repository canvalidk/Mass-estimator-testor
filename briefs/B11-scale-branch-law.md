# B11: the raw law on the "Σ declared up to a scale" branch

**Status:** done, awaiting check: scale derivation and 1D impropriety counterexamples; null-padding median finding; coverage study stopped under the false-premise protocol.
**Item:** `problems/O14.md` (primary), `problems/O7.md` (one finding). **Branch:** `codex/b11`.
**Written:** 2026-10-29, Claude.

Read `definition/theorems/T36.md` first (the reading space, its noise metric, one residual), then
`problems/O7.md` ("26–27 Oct: Z splits in two") and `problems/O14.md` ("29 Oct: the restart in the
reading space"). Background in standard terms: this is the integrated likelihood of a Gaussian
errors-in-variables model with an unknown common variance scale (Berger, Liseo & Wolpert 1999,
`docs/papers/`); integrating a Gaussian scale against a power-law prior gives a Student-type law.

## The sorted fork (Claude's, provisional; do not re-sort it)

Whether $\Sigma_i$ is declared exact or up to a common scale $v$ is an input (O7, 27 Oct). On the
"up to a scale" branch, with $g_i=\eta=1$, the law of $m$ is forced from the declared Gaussian
model (L13), the flat candidate weight (T6) and integrating the unknown scale out, as every
nuisance is (the output is a law), **up to one freedom**: the exponent $a$ of the prior
$v^{-a}$ on the scale. Rescaling the declaration ($\Sigma\to\lambda\Sigma$) cannot fix $a$ (every
power law is relatively invariant; Hartigan 1964), so $a$ is sorted as an **input**, with
$a=1$ (Jeffreys' scale prior) as its principal branch. The claimed result, to be derived and
checked, not assumed:

$$p(m)\ \propto\ \pi(m)\,K(m)^{-(n/2+a-1)},\qquad K(m)=\sum_iK_i(m),\quad n=\sum_i(2d_i-k_i),$$

with $K_i(m)=\lVert\mathbf y_i\rVert^2_{\Sigma_i^{-1}}-H_i(m)$, the residual against $L_i(m)$ in the
**whole** $2d_i$-dimensional reading space. With the full subspace ($k_i=d_i$) this is T36's
$E_i-H_i$ and $n=\sum_ik_i$. With a declared subspace ($k_i<d_i$) the data outside the span are
pure noise: ours30 drops them because they carry nothing about $m$, but here they carry
$2(d_i-k_i)$ dimensions of information about $v$. The exponent comes from integrating the
latent pair: $v^{-d_i}$ from the data density, $v^{+k_i/2}$ from the latent pair.

## Deliverable

1. **The derivation**, written into O14's attempts: from the declared model with $v\Sigma_i$, the
   flat candidate weight, and the prior $v^{-a}$, to the law above, every step named. Include
   the full covariance and declared subspaces ($H_i$ from `met.full`, not the isotropic closed form),
   and check the count $n=\sum_i(2d_i-k_i)$ for $k_i<d_i$; if it is wrong, say so and why.
2. **Tests** in `tests/test_scale_branch.py`:
   - the closed form against numerical integration over $v$, for random readings, $k=1,2,3$,
     correlated $\Sigma_i$;
   - invariance under $\Sigma_i\to\lambda\Sigma_i$ (all $i$ at once), units (L3) and the swap;
   - in 1D the law is defined and proper wherever ours30's is.
3. **A study** (a preset kind with its seed, following `pooled_scale`): a scalar world whose
   true noise is $v_0\Sigma_i$ with $v_0\in\{0.25,1,4\}$ and the shape declared correctly;
   $k=1,2,3$; $N=5,20,100$; signal from none to strong. Report the 95% interval's coverage of the
   true $m$ for three laws on the same draws: this branch's law ($a=1$); ours30's raw law
   ($g=\eta=1$, $\Sigma$ taken as exact); ours30 itself. Report widths beside coverage.
4. **Null readings (one finding for O7).** Append 1, 5, 20 exact-zero readings to a few fixed
   series: does this branch's law change, and does its median move? O7's table predicted
   "nothing to the point" on this branch. Report what happens with the mechanism (a null reading
   adds $k_i$ to $n$ and nothing to $K$). One attempt line in `problems/O7.md`.
5. One line in `docs/STATE.md`'s log.

## Do not

- Change `met/full.py`, `ours30`, `bartlett_eta` or the batch pipeline; the law lives in a study
  module (for example `met/scale_branch.py`), imported by nothing in production.
- Combine this law with $g_i$ or $\eta$, or say whether they are still needed on this branch.
  That is the open fork (O14, 29 Oct) and it is Claude's sort, not this brief's.
- Adopt the law, choose $a$, or state the coverage table as a verdict on any of the three.
- Judge the null-reading result against Z at the point: point-level Z was never ruled
  (R-2026-10-29-2). Report it.
- Touch π's shape (O6) or T36's parked item on the measure of $L_i(m)$.

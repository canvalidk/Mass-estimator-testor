# B15: η with the exact conditional score variance in place of the fitted J

**Status:** merged 2026-11-07, as delivered (not checked by Claude): T45 gives the working-model conditional score variance; η_C (opt-in in ours33, default unchanged) improves the flat leverage benchmark but misses the small-N target Δ/(Δ+M), and on T39's pair η_C = 0.036, not 1; G3 and G26 stay open
**Item:** `problems/O14.md` (4 Nov: η's target), `problems/O18.md` (4 Nov: ours33, hole 1); holes G3, G26; theorems T31, T43, T45.
**Branch:** `codex/b15`, **from the head of `codex/ours32-consistency`** (it needs `met/ours33.py` and `met/conditional_score.py`, which are not on `main` yet; recorded exception to the branch-from-main rule).
**Written:** 2026-11-04, Claude.

Read O14's 4 Nov section first, then T43, T45, O18's 4 Nov section, `met/ours33.py`.

## The sorted fork (Claude's, provisional; do not re-sort it)

η's target is now known: on long series $\eta\to\Delta/(\Delta+M)$, the share of the fitted explained
energy that is signal rather than the $M$ units of latent noise (O14, 4 Nov; the shrinkage ours28b got by
declaration). η estimates it as $\mathcal H/J$ with the **fitted** $J=\sum_i\ell_i'(\hat\phi)^2$, and that estimate
fails twice: at small $N$ it sits near 1 whatever the target ($N=5$, snr 0.5: 0.90 against 0.19–0.34), and
under leverage it under-covers in the limit (T43: one strong reading among nulls, 0.927). **Standard
name:** observed versus expected information in a sandwich; the fitted-residual middle term needs a
leverage correction (MacKinnon–White HC2; B04's pointer). T45 supplies the exact alternative: the score's
conditional variance $C$ given each reading's own candidate projection, which fits nothing, borrows
nothing across readings (L19) and is Z-safe. The fork: **is $\eta_C=\min(1,\mathcal H/C)$ the estimate of
$\Delta/(\Delta+M)$ that the form needs, or does it fail where $\mathcal H/J$ fails?** Not a replacement
equation: a measurement of one candidate estimate of a known target.

## Deliverable

1. **Derive $C$ for ours33's score.** T45 gives it for the flat Gaussian score. ours33's per-reading score
   is $a_i h_i'/(\nu_i+K_i)$ (met/ours33.py), nonlinear in the residual. Write the conditional variance of
   the summed score at the best line given each reading's own projection, under the declared noise with
   $v_i\sim\mathrm{IG}(\nu_i/2,\nu_i/2)$; if it has no closed form, give it by one-dimensional quadrature per
   reading and say so. Put the derivation in T45 as a section, with its test.
2. **Code:** `met/eta_conditional.py` (study module, imported by nothing in production): `eta_c(readings, c,
   nu_eff, factors="flat"|"ours33")` returning $(\eta_C,\mathcal H,C,\hat\phi)$; and a switch on `ours33` so the
   warp can be run with either η (default unchanged).
3. **Run and report** (seeded; a tool script beside `tools/eta_vs_size_law.py`):
   - O14's 4 Nov grid ($k=1,3$; $N=5,20,100$; snr 0.5–5): mean $\eta_C$ beside mean η and the target
     $\Delta/(\Delta+M)$. The question is the small-$N$ rows.
   - T43's sequence (one squared amplitude $4N$, $N-1$ nulls; $N=50,200,1000$): coverage of the warped
     ours33 interval and of the flat law's quoted interval with η and with $\eta_C$; T43's prediction is 0.927
     for η in the limit.
   - T39's opposed pair at $R=10^3$: $\eta$ is 0.0101; report $\eta_C$. If $\eta_C\approx1$, the pair's failure
     was the estimate, not the form, and say so in O18.
   - Z: append exact-zero readings; $\eta_C$ unchanged to the bit.
4. **Write the findings** into O14's attempts (one section) and one line each in O18 and `docs/STATE.md`'s
   log. Tests in `tests/test_eta_conditional.py`.

## Do not

- Adopt $\eta_C$ into ours30, ours32 or ours33's defaults, or change `met/full.py`.
- Touch the point (the raw median) or the cap at 1; report where the cap bites, do not remove it.
- Re-sort the fork: if $\eta_C$ fails somewhere, that is the result; write where and why.
- Extend the anchor/score-law constructions (T56–T65); this brief uses T45's variance only.

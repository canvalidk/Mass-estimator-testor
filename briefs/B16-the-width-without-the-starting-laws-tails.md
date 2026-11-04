# B16: the quoted width as a statement about the evidence, not a quantile of π·L

**Status:** merged 2026-11-07, as delivered (not checked by Claude): false-premise stop; ours33's evidence alone is improper in mass and in log mass on every finite valid series, disagreeing readings included, so no width exists without π; no study run
**Item:** `problems/O18.md` (4 Nov: ours33's break at the report); `definition/proof.md` step 18; R-2026-09-29-3; `docs/vd-docs-context/why_mass_cannot_be_undefined_2026-09-29.md` §3; T47.
**Branch:** `codex/b16`, from the head of `codex/ours32-consistency` (needs `met/ours33.py`).
**Written:** 2026-11-04, Claude.

Read O18's 4 Nov section, T47, R-2026-09-29-3 in `rulings.md`, the undefined-mass record §3.

## The sorted fork (Claude's, provisional; do not re-sort it)

ours33 applies η as T47's warp to the raw law $\pi\cdot L$. At small η the warped ranks
($\Phi(-z/\sqrt\eta)$: $10^{-7}$ at $\eta=0.16$, $10^{-85}$ at $\eta=1/90$) lie in π's Cauchy tails, which the data
never touch, so the quoted interval is π's tail: 78 decades on T39's pair while both readings sit at 2.
The power form (ours29–32) hid the same dependence by pulling the law toward $c$ instead. Can ruled on
29 Sept that the all-zero result's uncertainty is "unresolved", not the width of the law around $c$, and
the undefined-mass record says the result's relation to other information "must never be a comparison
against the starting law's width". **The sort:** the quoted width is a statement about the evidence
$L$; it is not a quantile of a law whose tails are the declared π's. **Standard name:** likelihood-based
(evidence) intervals against posterior credible intervals; at weak information the credible interval
reports the prior. The fork for this brief: **what does the width look like when it is computed on $L$
alone, and does property N hold for it where it is finite?** The point is untouched: it stays the median
of π·L.

## Deliverable

1. **When is $L$ alone proper on $m>0$?** For ours33's raw factors $\prod_ie^{\ell_i}$ (no π): derive the
   condition for finite central quantiles (each factor tends to a constant at $m\to0$ and $m\to\infty$, so the
   product is flat in the tails unless the readings disagree; state exactly when it is integrable in
   $\ln m$, and in $m$). The all-zero series and any series of nulls give a flat $L$: infinite width, which is
   R-2026-09-29-3's "unresolved". Write this into O18 as a short derivation with a test.
2. **The $L$-only interval**, in a study module (`met/evidence_width.py`): the raw law's central 95% and T47's
   warp of it, computed on $L$ with π removed (log-space tails as `met.ours33._tail_quantiles`), returning
   "unresolved" where $L$ is not proper. Report beside ours33's π·L interval on: T39's pair ($R=10^2,10^3$),
   the six-reading series of O18 (η = 0.16), and the cells of `tools/ours33_coverage.py` at truth $c$ and
   $4c$ (100 series): coverage, median log-width, fraction unresolved.
3. **Property N on the evidence interval.** On a null series ($m_0$ anything, signal 0) the evidence
   interval is unresolved by construction; on weak signal it is wide. Report whether its coverage, where
   resolved, is nominal across the signal grid, and how the fraction unresolved depends on $N$ and
   signal. No cure, no new form: the table.
4. **Compare with the power form on $L$ alone** ($L^\eta$ without π): which of the two is π-free *and*
   meets T47's Gaussian benchmark; if neither, say so.
5. Findings into O18's attempts; one line in `docs/STATE.md`'s log; tests in `tests/test_evidence_width.py`.

## Do not

- Change the point, π, or any production function; the median of π·L stays the point.
- Propose a new output type as adopted: this brief measures what an evidence-only width does; whether the
  output carries it is a sort (Claude's) and a ruling (Can's).
- Treat T47's law-only class as forced, or re-sort the fork above.
- Extend the anchor/score-law constructions (T56–T65).

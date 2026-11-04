# B17: what singles out T48's latent measure (the exponent (ν + k)/2)

**Status:** done, awaiting check: refuted T36's proposed forcing; b=-1 is also normalized, fixed-v normalization is additional; 20 counterexample checks, stopped per brief.
**Item:** `definition/theorems/T48.md` §2 ("a remaining scale-measure choice"), T36, T49; `problems/O18.md` (4 Nov, hole 3); G19.
**Branch:** `codex/b17`, from the head of `codex/ours32-consistency` (T48 and T49 live there).
**Written:** 2026-11-04, Claude.

Read T48 §§1–2, T36, T49 §5, O18's 4 Nov section (hole 3), L13 and L18 in `definition/premises.md`.

## The sorted fork (Claude's, provisional; do not re-sort it)

T48 integrates the normalised Gaussian over the candidate plane with arclength flat in the **declared**
metric and independent of the noise size $v$ ($b=0$), and gets the Student residual density with exponent
$(\nu+k)/2$. T48 §2 shows that multiplying the arclength by $v^b$ gives exponent $(\nu+k-2b)/2$ instead, and
Codex leaves $b$ unforced; T49 shows graph geometry does not pick it either. Every member of the family
passes Z, units, the swap and L2 (the $m$-dependence is the same $(1+K/\nu)^{-a}$; only $a$ moves), so the
usual checks cannot decide it. **The sort:** T36 decides it. A reading carries one metric, $\Sigma_i^{-1}$,
declared (L18); $v$ is latent, so a measure in the realised metric $v\Sigma_i$ builds an undeclared quantity
into the volume of the hypothesis, and L13 says the model uses the declared second moments and nothing
else. If that argument holds, $b=0$ is forced and the exponent is a theorem; if it does not, $b$ is an input
(then TASK §5's cost guard asks what practice supplies it) or a premise is mis-stated. **Standard name:**
the choice of dominating measure for an integrated (marginal) likelihood over incidental parameters;
Lindley (1966) and Berger–Liseo–Wolpert (1999) on integrated likelihoods for the functional
errors-in-variables model.

## Deliverable

1. **The argument, written out or refuted**, as a section of T48: for each $b$, which premise its measure
   uses or violates (T36's declared metric, L13, L18, T7's "nothing attached to a reading may depend on
   $m$", L3). Say exactly what $b=0$ rests on. If it rests on a reading of L13/T36 that is not in their rows
   as written, say so: that is the result, and it goes to the holes table as a *rule* or *sort* row, not into
   the premise rows.
2. **The consequences, tested** (`tests/test_latent_measure_exponent.py`): for $b\in\{-1,-\tfrac12,0,\tfrac12,1\}$ and
   $k=1,3$, $\nu=10$: the raw law's totality on B11's inputs; T30 on the three pairs of `tests/test_ours33.py`;
   T39 part 4's single reading ($\nu=\tfrac12$: which $b$ keep the raw law concentrating, i.e. $\nu+k-2b>1$);
   and whether the per-reading factor is a normalised density of the $k$-dimensional residual (only $b=0$,
   per T48 §2; check it).
3. **The forward model check**: in `met/forward.py`'s world with $v_i\sim\mathrm{IG}(5,5)$, the calibration of the
   reading weight (O18, 2 Nov: ours31's weight is the true precision at strong signal) for each $b$. Report
   which $b$ gives $E[w]/E[1/v]=1$; this is a measurement, not the forcing argument.
4. Findings into T48 and O18 (hole 3's status); one line in `docs/STATE.md`'s log.

## Do not

- Change `met/ours33.py`'s exponent, or any premise row.
- Pick $b$ by the forward-model result or by coverage: item 3 checks a derivation, it does not make one.
- Re-sort the fork; if T36 does not force $b=0$, write that and stop.
- Extend the anchor/score-law constructions (T56–T65).

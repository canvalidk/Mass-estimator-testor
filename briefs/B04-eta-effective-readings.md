# B04: is J's finite-N bias a function of effective readings under ours29's law?

**Status:** open
**Item:** `problems/O14.md`. **Branch:** `codex/b04`. **Written:** 2026-10-29, Claude.

Read `problems/O14.md`, T31 (`definition/theorems/T31.md`) and `bartlett_eta` in `met/full.py`.
The current equation is ours30 (since 29 Oct: `ours30()` in `met/full.py`, `positivity: true` in the
batch pipeline); "ours29's law" below means the tempered law as ours30 has it.

**Frame (added 29 Oct):** `definition/theorems/T36.md`. The linear-Gaussian bias below is the count of T36 (4):
fitting the one best line takes one dimension. On the exact-Σ branch η stays and this question is
unchanged. Do not use T36's residual as a replacement for J: that is the other branch (B11) and an open sort.

## The sorted fork (Claude's, provisional)

eta = min(1, curv / J), with curv = Σ cᵢ (cᵢ = −lᵢ'' at the series' best line) and
J = Σ sᵢ² (sᵢ = lᵢ' there). J is biased low at finite N because the scores are taken at the
**fitted** best line, where they sum to zero, as a residual sum of squares is after fitting a
mean. In the linear-Gaussian case (readings xᵢ ~ N(θ, 1/cᵢ), sᵢ = cᵢ(xᵢ − θ̂)) this is exact:

  E[J] = C − Σcᵢ²/C = C (1 − 1/n_eff),  C = Σcᵢ,  n_eff = C² / Σcᵢ².

So the bias is forced by the fit and is a function of effective readings, not of N: lemma list
d4's conjecture for O14, here with its mechanism. A null reading has cᵢ = sᵢ = 0 and leaves C,
Σcᵢ² and J unchanged, so anything built on n_eff is Z-exact. **Not known:** whether this
survives ours29's law away from the linear limit (weak signal, the projective line, gᵢ < 1).
That is the question.

## Deliverable

1. The linear-Gaussian identity as a test, by exact moments, not sampling.
2. Under ours30's per-reading terms (`met/full.py`; isotropic noise; k = 1, 2, 3; c = s;
   signal sizes from pure noise to strong; N = 5, 20, 100, 400), seeded: E[J] / C against
   1 − 1/n_eff, with Monte Carlo standard errors. Where does the identity hold, and how does it
   fail at weak signal? Report the spread of J / C as well as its mean (J is heavy-tailed, O14).
3. Alongside eta, the candidate eta′ = min(1, curv / (J · n_eff / (n_eff − 1))) (Claude's,
   **not adopted**; eta′ = eta when n_eff ≤ 1), and the coverage of the quoted 95% interval
   under each, in the same cells. Lemma list d4 records 88–89% coverage at weak signal, N = 20.
4. A preset with its seed for 2–3 (a new kind, following `pooled_scale`); tests; into
   `problems/O14.md`'s attempts: the identity's range of validity and the coverage table;
   one line in `docs/STATE.md`'s log.

## Do not

- Change `bartlett_eta`, `ours29`, `ours30` or anything in `met/full.py`: eta′ lives in the study module.
- Adopt eta′, or state the coverage table as a verdict on it. The question is whether the
  identity holds; coverage is reported beside it.
- Touch the pooled-scale branch (O14's 27 Oct section) or O7.

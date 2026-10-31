# B13: does Gleser–Hwang apply to the type? (the corner of O15)

**Status:** done, awaiting check: zero uniform coverage proved directly for the type; original 1987 PDF unavailable, verified restatement cited; fixed-mass zero-limit clause disproved and positive-floor regression added.
**Item:** `problems/O15.md` (30 Oct sort, regime (a)); property N (R-2026-09-29-4). **Branch:** `codex/b13`.
**Written:** 2026-10-30, Claude.

Read O15's 30 Oct section, `definition/type.md`, L2, L10 and G8 in `definition/premises.md`, T24.

## The sorted fork (Claude's, provisional; do not re-sort it)

O15's corner is sorted as a case of Gleser & Hwang (1987), *Ann. Statist.* 15(4):1351–1362,
doi:10.1214/aos/1176350597 (and Dufour 1997, *Econometrica* 65(6), for the general form): near parameters
where $m$ is not identified, an always-finite confidence procedure has no positive guaranteed coverage.
The sort stands or falls on whether the paper's hypotheses hold for **our** model, which differs from the
textbook one: a line through the origin, $\mathbf F^\ast_i=m\,\mathbf a^\ast_i$, $m>0$ (a restricted
space); declared, known $\Sigma_i$; latent pairs free (functional, L10); $d\ge1$, $N$ readings, reading
structure $N_i$; and an interval that is always finite (L2 with π proper, G8). The PDF is not in this
public repo (`../legacy/README.md`); fetch it by its doi if you can.

## Deliverable

Written into O15's attempts as a new section "B13":

1. **The theorem as stated in the paper** (paraphrased, with its section and equation numbers) and its
   hypotheses, numbered.
2. **Each hypothesis checked against the type:** holds / fails / unknown, with the reason.
3. **The corollary for the type**, stated and proved, whether or not the paper applies verbatim. Claude's
   sketch, to be checked, not assumed: as every true push → 0 the readings' law tends to one free of $m$
   (continuity in the latent pairs), so the interval's law tends to one law for every $m$; it is a.s.
   finite, so for $m$ far enough from $c$ the coverage → 0, and the infimum coverage over $(m,\text{pushes})$
   of any always-finite procedure is 0. Say whether the infimum is 0 or only below $1-\alpha$, and whether
   the result needs $N$ fixed or holds uniformly in $N$.
4. **One test**, `tests/test_o15_corner.py`: ours30's coverage at planted $m=c\,e^{\pm3}$ falls toward 0 as
   the push shrinks (3D, a few hundred series, seeded; mark it slow if it is over 20 s).

If the paper cannot be fetched, do 2–4 from the sketch and say that 1 is missing.

## Do not

- Do not change `met/`, π, η or any premise. Do not re-sort regime (b) (Neyman–Scott); it is not this brief.
- Do not propose a fix for the corner. If the corollary holds, the shortfall there is a theorem, not a defect.
- Do not edit `docs/papers/README.md` beyond adding the PDF's line if you add the PDF (you may not: public repo).

# B08: the family of the starting law, derived in one place

**Status:** open
**Item:** `problems/O6.md`, `definition/theorems/T28.md`; `definition/proof.md` step 6, hole G8.
**Branch:** `codex/b08`. **Written:** 2026-10-29, Claude.

**Not this brief's (added 29 Oct):** T36's parked item, the measure the noise metric induces on the mass
subspace $L_i(m)$, which is where π's shape meets the metric. Leave it in T36.

## The sorted fork

None: $\pi$'s shape is an input (R-2026-10-24-3), and the open question O6 (which member is the
principal branch) is not this brief's. What is owed is the derivation of the *family* the input
ranges over, which is spread across O6 ("What is forced"), T28 and the Jeffreys-check record.

## Deliverable

1. Into `problems/O6.md`, section "What is forced", a derivation: from L3 (units covariance: a
   rescaling of force or acceleration units rescales $m$, and $\pi$ must transform with it), L4
   (symmetry under $m/c\leftrightarrow c/m$) and properness (T27's reason: on null input the law
   is $\pi$ itself and must be a law), every admissible $\pi$ is $m^{-1}g(\ln(m/c))$ with $g$ even
   and integrable, and conversely every such $\pi$ is admissible. State what "counted once" (T7)
   adds and what it does not.
2. The half-Cauchy as the member uniform in the angle $\theta=\arctan(m/c)$ (T28): shown, with
   what it does and does not follow from (T28 under ours28's model; under ours30 only that the
   reference form does not depend on a size law).
3. A test that the law on the all-zero experiment is $\pi$ and that its median is $c$ for three
   members of the family (half-Cauchy, log-normal in $\ln m$, log-logistic), using `met/full.py`
   only where it takes a general $\pi$; if it does not, say so and test the family algebra alone.
4. One line in `docs/STATE.md`'s log; O6's log.

## Do not

- Pick a member, argue for one, or change `equation.md`'s default. O6 stays open.
- Change `met/`.

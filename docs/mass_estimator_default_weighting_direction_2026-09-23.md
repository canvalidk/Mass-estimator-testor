# Investigation direction: what the default weighting can reference

**2026-09-23.** Filed at the user's request as an investigation direction, not a
result. Nothing here is concluded. It scopes and binds nothing.

## Where it came from

Earlier the same evening, an argument was proposed. If the default weighting may
not reference any instrument (so that combining readings from different
instruments does not depend on their order), respects unit changes, and treats
force and acceleration symmetrically, then it must be $(f\alpha)^k$, leaving a
single free exponent $k$. That argument holds only for repeated readings of one
unchanging force–acceleration pair.

## What a quick check suggested

For N readings of one mass, each with its own excitation, the same three
requirements with one factor per reading force the weighting
$m^{N(k+1)-1}\prod_i\alpha_i^{2k+1}$. For $k=0$ this is contribution 12's
eq. (20).

A large-N check (1D, $\sigma_F=\sigma_a=1$, fixed excitation per reading) found
that every tested $k$ settles away from the true mass whenever the mass is not
$\sigma_F/\sigma_a$. Each estimate is pulled towards $\sigma_F/\sigma_a$: for
example 0.55–0.84 × truth at $m=4s$ and 2.2–3.1 × truth at $m=s/4$, both at
acceleration SNR 1. The symmetric rule, whose per-reading weight uses that
reading's own noise ratio, settled within about 1% of the truth in every tested
cell. The $k=0$ values agree with contribution 12's 3D results.

## Possible reading, to be tested

- **Per-reading weighting.** Referencing each reading's own instrument may be
  needed for convergence. It does not break order-independence, because each
  reading brings its own factor.
- **The one-time prior on mass.** If it is scale-free and channel-symmetric, it
  is $1/m$, which gives a posterior that cannot be normalised. So a proper,
  knowledge-free prior would have to borrow a scale, most plainly the
  instrument's.
- **What would remain free.** Only the shape of that one prior (where λ lives),
  which washes out with enough readings.

## Open questions

- Prove or refute convergence, rather than rely on the numbers.
- Test weightings that are not one factor per reading (e.g. functions of
  excitation ratios).
- Test 3D, varying excitations, and readings from different instruments.
- Check whether the conclusion about the one-time prior survives those.

## Reproduction

VD-Newton `03_trace_and_evaluator/mass_estimator/investigation/default_weighting_consistency_check.py`
(exploratory script; `python default_weighting_consistency_check.py [R] [seed]`).
It was run with R = 4000, seed 1, and the symmetric rule was rerun with
R = 20000, seeds 2–3.

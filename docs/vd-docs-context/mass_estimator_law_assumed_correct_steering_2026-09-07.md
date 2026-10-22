# Mass estimator investigation: condition on Newton II being correct

**Date:** 2026-09-07. **Status:** user steering for the present investigation;
local record awaiting transfer to VD-docs, not submitted.

## User correction

The user clarified the purpose after the assistant called an anti-aligned
central mass outside a threshold compatibility set potentially unwanted:

> no the idea is to see what must be true if the law is correct. Given anti-aligned, the only way for the measured values to occur is with error. That means we are looking for the 'smallest' possible error. it treats force and accelerations as equal information and gives you the most likely reconciliation. It has to assume that no matter how unaligned, the relationship didn't fail, so it must simply be an unlikely thing to happen.

For this investigation, Newton II is a premise. For empirical inputs under the
declared nondegenerate Gaussian model, incompatibility of the central readings
is explained by measurement error. Large discrepancy describes the unusual
errors required by the law-constrained account; crossing an illustrative
discrepancy threshold is not a reason to stop conditional reconciliation.

The assistant's preceding objection used the wrong criterion: it asked whether
the output passed the older specification's empirical rejection/acceptance
contract. The user is asking what follows conditional on the law being correct.
The fact that the central point fails an arbitrary fixed-mass threshold does
not by itself demonstrate unwanted behavior for that purpose.

Force and acceleration should be treated symmetrically as sources of information.
Under the current Gaussian policy their corrections are measured relative to
their respective uncertainties, using the full declared covariance. Symmetry
does not identify their physical units or require equally precise instruments.

## Mathematical distinction still requiring care

The user's phrase “smallest possible error” must not silently change the
already agreed equation. Under Gaussian likelihood, minimizing the discrepancy
Q selects a maximum-likelihood reconciliation (or a boundary infimum). The
current estimator instead weights the full compatible-pair law and reports
E[f]/E[alpha]. This averaging rule need not select a minimum-Q pair or a mode
of the mass-ratio distribution.

The already tested symmetric anti-aligned example makes the distinction without
using an acceptance threshold. With F=(-sqrt(6),0), a=(sqrt(6),0), and unit
channel uncertainties, the infimum discrepancy is 6, approached by keeping one
channel and moving the other toward zero. The two alternatives are equally
favored by symmetry. The ratio-of-means point is 1; every fixed-mass-1
reconciliation has discrepancy at least 12. That is not the minimum error.

Consequently, two claims must be kept separate:

- The current construction assumes the law and gives less erroneous pairs
  greater likelihood weight.
- Its reported scalar is the single least-error or most-likely reconciliation.

The first describes the construction; the second does not hold in general.
The correction does not authorize replacing the ratio-of-means equation with
an optimizer. If literal least-error selection is a requirement, its relationship
to the existing readout is a substantive target for the next behavioral analysis.
No new point equation is selected in this record.

## Consequences for the existing laboratory

The September 6–7 calculations of Q, ratio laws, dependence, and central
readouts remain useful. Earlier gate/withholding tests reproduce the previously
inspected working specification. Their passing status is not evidence that
rejection is appropriate to the user's clarified conditional investigation.
Likewise, earlier “unwanted” interpretations based only on crossing threshold 9
are withdrawn for this purpose.

Preserve Q and the original measurement residuals to describe the error needed.
Do not use their size to imply that conditional law-based inference has failed.
This statement is about empirical Gaussian inputs. It does not redefine exact
contradictory premises as noisy data or guarantee support under a future bounded
error model.

The README and the two local behavior reports now point to this steering record
so their earlier gate framing is not silently treated as the current purpose.
Numerical reference and gate functions are retained as inspectable historical
behavior; no production entry or managed source is rewritten here.

## Provenance and custody

Authority: the user's clarification in this task, reproduced above. It governs
the present investigation over the earlier managed-specification framing.

Managed basis already read in this session:
`canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6`, resolved on
2026-09-07; `catalog.yaml` and
`Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`.
That specification separates integration/readout from an absolute rejection gate.
The conflict is recorded locally for handoff rather than changing its authority
or contents. No new managed discovery was required to record a user correction.

Local evidence: the September 6 behavior report, September 7 uncertainty/reuse
walkthrough, and their saved anti-aligned experiments. The numbers above were
already verified there. No new numerical or literature claim is made.

Saved in `pass to VD-docs/` under the workspace handoff policy. Awaiting transfer;
no submission, commit, push, or managed policy revision is implied.

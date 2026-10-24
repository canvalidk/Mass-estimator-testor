# Independent audit of weak-acceleration contribution

2026-09-15. Reviewed the proposed checks, saved results, mathematical review,
and one-sided benchmark. No production files edited. The lead manuscript was
not yet available when this initial audit was written.

## Findings

No mathematical defect found in the main results inspected:

- The exact-known-force positive-truncated Gaussian posterior and transformed
  mass density are correct. Its mass mean diverges at every finite reading;
  the ratio-of-means readout remains finite.
- The mass lower-bound coverage event is oriented correctly: posterior mass
  CDF at the true mass must be at least .05, whereas posterior acceleration
  CDF at the true acceleration must be at most .95.
- The exact-known-force continuous-reading coverage identity is correct and
  proves conservative sampling coverage. It is scoped to known force,
  direction, Gaussian standard uncertainty, and the stated reference measure.
- The high-mass boundary limit of the point relative to true mass is zero
  almost surely under a common-noise coupling. This is a weak-signal limit,
  distinct from inconsistency at a fixed finite mass as measurements accumulate.
- Censoring uses the correct bin likelihood, normalized CDF, and mean. The
  hybrid reporting rule's coverage jump is real and derives analytically.
- The scalar correlated-error one-sided pivot, concavity, possible set
  topology, and lower-hull implementation agree. The exact boundary at a force
  reading of the Gaussian critical value is handled correctly.
- The noisy-force lower-credible-bound failure as true force tends to zero
  is valid. It does not establish failure within a strong-force operating domain.

## Independent numerical verification

`check_audit.py` integrates the joint positive-quadrant law by conditioning
force on acceleration and performing Cartesian quadrature. This differs from
the radial/ratio-angle implementation under review. Thirty cases cover
correlations -.98, -.8, 0, .8, .98; negative/zero/positive readings; and mass
ratios through 160. Sharp conditional-CDF transitions are explicitly split
into integration intervals.

Maximum discrepancy from the radial implementation:

- Relative point: 4.70e-14.
- Absolute mass CDF: 2.04e-14.

The hybrid bin coverage was checked separately on both sides of its jump:

| Bin half-width k | Coverage at theta=B | Right-hand limit |
|---:|---:|---:|
| .1 | .98045446 | .96878483 |
| 1 | .99943131 | .89510391 |
| 2 | .99999953 | .81656341 |

`audit_results.json` contains the machine-readable results. These are numerical
consistency checks, not proofs of integration error over every possible input.

## Required presentation qualifications

1. A Monte Carlo coverage of 1 and plug-in standard error of 0 do not establish
   exact coverage or zero simulation uncertainty. With 16,384 successes, the
   one-sided 95% exact binomial lower confidence endpoint is
   `.05**(1/16384)`, approximately .999817. State that there were no observed
   failures or provide an appropriate interval.
2. Label `Phi(B-k)` as the right-hand coverage limit at theta=B. It is not the
   coverage at B, where inclusion of the reported bin gives `Phi(B+k)`.
3. Keep the point limitation next to the lower-bound success. A sound one-sided
   inference is compatible with severe relative underestimation by the finite
   central readout in weak signal.
4. Unknown-direction 3D tables are observed-data comparisons under independent
   isotropic known errors. They are not demonstrations of sampling coverage for
   the vector posterior bounds. The scalar theorem cannot supply that guarantee.
5. The classical comparator addresses one-sided information in the same scalar
   experiment. It neither supplies a competing compositional point readout nor
   establishes general practical superiority over the user's full method.
6. Treat the safe bin Gaussian endpoint as part of a specified reporting rule:
   use k+z on the bin and the ordinary Gaussian endpoint outside, or another
   explicitly stated conservative extension. A percentage attached solely to
   the observed bin is not conditional frequentist coverage given that bin.

## Optional implementation limits

The numerical root bracket in `known_force_results`/`coarsened_coverage` is
fixed to a lower observation of -100 and therefore supports the supplied
scenarios, not arbitrarily small positive true acceleration. This does not
affect the analytic theorem or current results. If exposed as a general-purpose
API later, grow the negative bracket adaptively.

The script's source hashes cover itself and the baseline estimator, but its
vector calculations also import existing research helpers. For a standalone
reproduction bundle, identify those helper versions or hash them too.

## Manuscript follow-up review

The completed `weak_acceleration_information.md` was subsequently reviewed.
Its numerical values, scalar/vector separation, confidence versus posterior
language, bin reporting protocol, and right-limit qualifications agree with
the results. Wilson intervals now address zero observed failures in simulations.

One required scope clarification was sent to the lead: Section 6 must require
positive-definite joint Gaussian covariance (or positive-definite residual
covariance for each candidate), since a merely known covariance can be singular
and then the displayed inverse and chi-square-with-d-degrees claim do not apply.
Section 5 should explicitly say that its covariance is known.

Two minor wording/proof improvements were also suggested: describe the noisy-
force Monte Carlo result as empirical coverage consistent with or above 95%,
and state `Q_.95(-k)<B_k<Q_.95(k)` before the censoring jump calculation. The
latter follows because the bin posterior mixes exact-reading posteriors and
makes explicit why all sufficiently nearby positive outside-bin readings cover
and all negative outside-bin readings fail to cover at that boundary.

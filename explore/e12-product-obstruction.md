# When bounded conflict evidence can compose

This is a necessary-condition argument for the scalar proof. Its assumptions concern the **returned law**. The conclusion changes if a separate report alone must stay wide. Keep the prior fixed while changing both channel standard deviations by a common factor \(\sigma\downarrow0\). Every finite \(\sigma\) is a measured experiment.

The standard setting is sequential updating by additive losses. Independent raw factors can be multiplied with the prior counted once. In this repo T12 explicitly permits a union-level correction to be recomputed from retained state. The stronger hypothesis below multiplies the **final** factors. [Bissiri, Holmes and Walker (2016)](https://doi.org/10.1111/rssb.12158) · [T12](../definition/theorems/T12.md)

## A duplicate-reading obstruction

Write \(P_{D,\sigma}(dm)\propto \pi(dm)e^{\ell_{D,\sigma}(m)}\) and \({\rm osc}(f)=\sup_{m,m'>0}|f(m)-f(m')|\). Suppose all factors are positive, and for disjoint lists \(D,E\),
\[
\ell_{D\cup E,\sigma}=\ell_{D,\sigma}+\ell_{E,\sigma}+c(D,E,\sigma)
\]
with \(c\) independent of mass. Duplicate numerical observations can be separately acquired and declared independent. Suppose for each *fixed* list admitting no common positive mass, its final log-evidence oscillation stays bounded as the declared errors shrink. The bound may depend on the list.

Take any nonzero compatible reading \(x\) and a reading \(y\) incompatible with it. The lists \((x,y)\) and \((x,x,y)\) are both incompatible. Cancellation gives
\[
 {\rm osc}(\ell_{x,\sigma})
 \le {\rm osc}(\ell_{x,x,y,\sigma})+{\rm osc}(\ell_{x,y,\sigma}).
\]
Thus the singleton oscillation is eventually bounded too. For any measurable set \(A\), \(P_{x,\sigma}(A)\le e^{{\rm osc}(\ell_{x,\sigma})}\pi(A)\). A continuous prior then precludes singleton law concentration at a single mass along this fixed observation. More generally, any fixed block can replace \(x\) if a suitable conflicting block can be appended. The obstruction uses two finite experiments and requires neither Gaussian factors nor a particular readout.

**Strength of the hypothesis.** Uniform boundedness here concerns *all candidate-mass contrasts for each observed list*. O14's words about the law staying spread can also be read more weakly. For example, multiplying a fixed positive angular prior by \(\exp[-t\,{\rm dist}(\theta,[a,b])^2]\) leaves a spread-out informative limit on \([a,b]\) and has unbounded oscillation. Nor does the cancellation argument restrict multiplication of raw factors when a series operation changes the final law. These are genuine differences in the theorem's premises. [Whole-series reading of T30](../problems/O14.md#17-nov-later-the-two-sorts-owed-on-the-discount-claude-provisional)

## From fixed observations to Gaussian experiments

The preceding proof addresses a deterministic noise-shrink path. A singleton observed *exactly* at its mean has probability zero. The following argument bridges that distinction. Assume measurable update factors and use angle coordinates on \((0,\pi/2)\).

Take the mean \(\mu\) of a single raw reading uniformly over a compact positive rectangle \(R\) bounded away from the axes, then observe \(X=\mu+\sigma Z\) with \(Z\sim N(0,I_2)\). The uniform distribution is a device for the proof. The resulting contradiction concerns at least one *fixed* true mean. Each such mean has a positive slope and is identifying as measurement error shrinks. Along \(\sigma_j=1/j\), the conditional law of \(X\) mixed over \(\mu\) tends in total variation to the uniform distribution on \(R\). Indeed, translating a uniform rectangle by a deterministic shift tending to zero moves a vanishing fraction of its area. Averaging this bound over \(Z\) gives the claim.

Fix the uniform-angle prior. For every raw \(x\in R\), cancellation gives an eventual finite bound \(M(x)\) for singleton oscillation. Any interval holding probability \(p\) then has angular length at least \(p(\pi/2)e^{-M(x)}\). Its bounded angular length therefore has positive liminf in expectation under the uniform distribution on \(R\), by Fatou's lemma. Total variation transfers this conclusion to the mixed observations. If every true \(\mu\in R\) gave singleton law concentration in probability at its slope, central \(p\)-interval lengths would converge to zero in mean for each \(\mu\), and bounded convergence over \(R\) would give zero for the same mixture. This proves the incompatibility in the Gaussian experiment, even without locally uniform bounds in \(x\).

**A median version.** Point consistency alone permits a spread-out law with the right median. There is a stronger corollary when the point is the median, the log factors at fixed prior centre are independent of the *shape* of the prior, and every positive reciprocal-symmetric density about this centre is admissible. Fix centre 1 and select \(R\) with all mean slopes above 2. Some positive-area subset has \(\limsup_j{\rm osc}(\ell_{x,\sigma_j})\le M\) for one finite \(M\), since the union of those subsets for integer \(M\) covers \(R\). Choose a sufficiently narrow lognormal prior centred at 1 so that \(\pi(m>1.5)<e^{-(M+1)}/2\). Shape independence keeps the same bound and subset. The median is eventually at most 1.5 for every observed point in that subset. Fatou and the same total variation transfer give positive liminf probability of a median at most 1.5 in the mixture. Consistency for every fixed mean in \(R\) would give probability zero in the limit. Thus even point-only L11 over these prior shapes conflicts with the two hypotheses about the final law. This corollary does not cover every possible point functional or a stipulated single prior shape.

## Two further boundaries of the requirement

* **The axes.** Suppose the final data factor for a perfectly compatible singleton at *each* finite \(\sigma\) is exactly its Gaussian factor, up to normalization. If the returned law is continuous in the measured row, then for the row \((F,0)\), \(F>0\), it must equal the Gaussian law by approaching with positive accelerations. Its angular kernel is proportional to \(\exp[-F^2\cos^2\theta/(2\sigma^2)]\). With a full-support positive prior the law tends to the infinite-mass boundary as \(\sigma\downarrow0\). No finite positive mass fits that row. Thus literal whole-series spread on this boundary, exact face-value Gaussian treatment of all perfectly compatible singletons at finite noise, and continuity in the data cannot all be required. Interchanging the channels gives the zero-mass boundary.
* **A balanced scatter.** The scalar rows \((1,0),(0,1)\) have Gram matrix \(I_2\). Their Gaussian data factor has identical value at every mass while the total residual is strictly positive at every mass. Any rule allowed *only to lower the magnitude of each Gaussian pairwise log contrast* must also have zero contrast on these observations at every noise scale. Literal informative evidence on every incompatible shape requires other reading information or a narrower meaning of “informative”. A pooled quadratic power cannot supply it. This argument treats the Gaussian contrast constraint as an additional hypothesis, separate from the ruling that the instrument errors are believed.

The [scalar construction](e12-scalar-proof.md) shows where each restriction is spent. The tests exercise the duplicate identity, the posterior mass bound and both boundaries. [Finite checks](../tests/test_e12_scalar_proof.py)

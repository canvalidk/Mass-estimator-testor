# Independent review of conditional LR implementation

2026-09-15. Scratch review; no production edits. Reviewed
`conditional_checks.py` with SHA-256
`1068b084d73243b96e91fbb5bc17877af2117447c18b9a00fce7e2c30ef66929`.
This records the pre-correction implementation; later corrections must be
checked separately.

## Mathematical and simulation checks

The implemented conditional CDF, stable LR algebra, decreasing-CDF inversion
to one fitted-energy threshold, and positive quadratic inversion agree with
the direct derivation in `derivation.md`. The simulated standardized group
averages have the correct means and identity channel covariance for the
stated balanced independent Gaussian repeated-readings model.

An independent conditional Monte Carlo used one million draws of independent
U~chi-square(1), W~chi-square(K-1), at K=1,3,24,1536,12288. For each K, the
same U,W sample was reused at conditioning strengths s=0,K,10K. Every
calibrated acceptance rate lay within 1.36 binomial standard errors of .95.
This validates the tested operating range, not arbitrary finite inputs.

The independent K=3 expression is

    F3(c;s)=erf(sqrt(c/2))
      -sqrt(2c/pi) exp(-c/2) integral_0^1 exp[-s(1-u^2)/2]du.

For s>200, the review substitutes v=s(1-u), integrates v from 0 to 100,
and bounds the omitted integral by 2 exp(-50)/s. The 1,722-point K=3
grid revealed an extreme-scale failure that self-refinement did not
always detect. The maximum error in the subgrid c,s<=1e5 was 7.2e-10,
at values whose true CDF is effectively one.

## Actionable numerical correction

The first Gaussian split in `cdf_weights` uses total L=c+s, but the
Gaussian factor has scale c. At K=3,c=1e8,s=1e12 the 64-node result is
0.0637620071 although the independent exact expression rounds to 1.

More seriously, `fitted_threshold(L=1e12,K=3)` returns

    c=L-b=60022794073.34546

although the correct c is approximately 3.841458821. Evaluating the wrong
root with both 64 and 256 nodes returns .95 to roundoff. Both quadratures
miss mass in the broad interval after the initial narrow Gaussian split.
The first retained interval has effective Gaussian cutoff 8 sqrt(c/L),
so the resulting false root is

    c/L = chi-square(1,.95)/64 = 0.0600227940733.

The same pathology occurs at other K because the chi-square factor is then
essentially one in the relevant Gaussian region.

Two concrete repairs follow from the mathematics:

1. Along the fitted-threshold root, bracket c in
   [chi-square(1,.95), chi-square(K,.95)], rather than [0,L]. The critical
   curve always lies in this interval, and the root is needed only when
   L exceeds the upper endpoint.
2. For direct CDF quadrature, make the Gaussian split depend on c. A useful
   bounded-error version integrates only z<=8 (equivalently beta<=
   asin(min(1,8/sqrt(c)))); the omitted probability is at most
   erfc(8/sqrt(2))=1.2442e-15. Clip the chi-square transition splits to this
   interval. Rebuilding c-dependent weights during the bounded root avoids
   the problematic cache assumption.

The existing 4096-by-11 reported simulation values remain supported by the
independent conditional Monte Carlo. The error is a numerical scope defect,
not a counterexample to the exact confidence construction.

Reproduction: `review_conditional.py`; numerical output: `review_results.json`.

## Corrected implementation recheck

The revised production script with SHA-256
`e4290180c800c0ee0de64e30e23459912ce8828f60fb947a7410cfa41f7dd22d`
was independently checked after both repairs were applied. The original
review, probe, and numerical evidence above are retained unchanged.

The same 1,722-point K=3 grid, spanning c from 1e-8 through 1e12 and s from
zero through 1e12, now agrees with the independent closed expression to a
maximum absolute difference of 6.31e-13. The former catastrophic examples
are repaired. All 15 previously simulated conditional critical values give
CDF .95 within 1.05e-13 under the corrected implementation, retaining the
independent million-draw Monte Carlo check without rerunning the same draws.

A threshold grid covers K=1,2,3,24,192,1536,12288, values just above the
whole-domain cutoff, ordinary signal scales, and L up to 1e12. At L=1e12,
the reconstructed c=L-b rounds to 3.8414306640625 because the float64
spacing of L is 0.0001220703125. Its CDF error of roughly 8.4e-7 is
consistent with that subtraction limitation, rather than the earlier
quadrature failure. No universal floating-point guarantee is inferred.

The Gaussian truncation tail argument and analytic critical-value bracket
in the revised implementation are correct. The fifth contribution,
`conditional_uncertainty_and_design.md`, was also reviewed for the conditional
test, positive-domain qualification, learned-strength curvature, null
interpretation, and repeated-group scope. No actionable mathematical or
scope correction was found. The full-simulation reproduction summary should
reflect the final rerun and can state the explicit Gaussian tail bound.

Recheck: `review_corrected.py`; new output: `review_corrected_results.json`.

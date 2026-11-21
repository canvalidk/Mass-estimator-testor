# The uncertainty correction for two labelled masses is a matrix

On a scalar push with one shared acceleration the true three-channel mean is \(a_i b(m)\), where \(m=(m_1,m_2)>0\) and \(b(m)=(m_1,m_2,1)^\top\). Whiten the observation by a factor \(W_i\) with \(W_i\Sigma_iW_i^\top=I\). Write \(y_i\) for the whitened data, \(u_i=W_i b/\|W_i b\|\), and \(D_i=\partial u_i/\partial(m_1,m_2)\), a matrix with two columns normal to \(u_i\). The pullback of the spatial noise metric is \(M_i=D_i^\top D_i>0\). A different whitening rotates all these vectors together. Every push may have its own signed strength and declared noise covariance. [Physical two-body experiment](e12-mass-observables.md) · [Reciprocal tangent](e12-joining-bodies.md)

## An exact covariance calculation

The log profile factor is \(\ell_i=(u_i^\top y_i)^2/2\) up to a mass-independent term. Its mass score is the covector \(s_i=(u_i^\top y_i)D_i^\top y_i\). At the true pair let \(r_i=a_i\|W_i b\|\). The longitudinal reading \(U_i=u_i^\top y_i\sim N(r_i,1)\) is independent of \(D_i^\top y_i\sim N(0,M_i)\). For independent pushes, writing \(S=\sum_i s_i\), we have
\[
 S\mid (U_i)_i\sim N(0,C),\quad C=\sum_i U_i^2 M_i,\qquad
 B=\operatorname{Cov}S=\sum_i(1+r_i^2)M_i,\qquad
 A=-E\,\partial_m^2\sum_i\ell_i=\sum_i r_i^2 M_i .
\]
For the Hessian differentiate \(E_{m_0}\ell_i(m)=[1+r_i^2(u_i(m)^\top u_i(m_0))^2]/2\). The second derivative of \(u^\top u_0\) at truth is \(-D^\top D\) by the unit-length identity. In particular **\(B-A=\sum_i M_i\) is known at each candidate from the instrument alone**. Both \(C\) and the observed score are available from the measurements. \(S^\top C^{-1} S\) has exactly the conditional \(\chi^2_2\) law when \(C>0\), almost surely for any fixed true pair and a nonempty block. An appended exactly zero observed row contributes zero to \(S,C\). This extends the scalar [T45 pivot](../definition/theorems/T45.md).

## Which width is forced near a regular fit?

Consider \(N\) independent repetitions of a finite design with an interior true pair. Hold both the physical drives and covariances fixed per design, and assume at least one nonzero drive. \(A,B\) above refer to one design. A consistent local maximum of the profile has asymptotic sampling covariance \(A^{-1}BA^{-1}/N\). Indeed the independent score sums obey the central limit theorem, and the empirical Hessian divided by \(N\) converges locally to its expectation by the law of large numbers applied to the channel scatter matrices. The expected criterion has a strict maximum at the true pair: equality for a nonzero push requires equality of the two lines, which inverts to the same positive masses. The score expansion about truth gives the formula.

A regular locally normal report centered there can have correct coverage for all its linear contrasts only with that covariance. Thus its local precision must be \(N A B^{-1} A\). The local quadratic form of the observable conditional pivot has exactly this limit, since \(C/N\to B\) and the derivative of the score divided by \(N\) tends to \(-A\). Equivalently, evaluate the observed negative Hessian and \(C\) at the regular fit and form the covariance \(\widehat A^{-1} C\widehat A^{-1}\) using block **sums**. Both ingredients ignore appended exactly zero observations. This is an estimator of local width. Small or indefinite observed curvature can occur on arbitrary finite data. No inverse-curvature readout is implied in that case. The sandwich formula is standard in estimating-equation theory. The calculation here identifies its two tensors directly in the mass geometry. [White, 1982](https://doi.org/10.2307/1912526) · [Scalar curvature argument](../definition/theorems/T31.md)

A scalar multiplication of the block log factor by \(\eta\) gives local covariance \(A^{-1}/(N\eta)\) under a smooth positive starting density. It matches every contrast exactly to first order iff \(A=\eta B\). The same requirement applies to any scalar chosen from the observed block and held fixed over candidate masses, if it has a positive probability limit. If all the \(M_i\) are proportional, such a scalar exists in the regular experiment. A single unlabelled mass always has this one-dimensional simplification.

**An experiment requiring two different discounts.** Use mass units in which the true labelled masses are both 1. Record \((F_1-A_{\rm obs},F_2-A_{\rm obs},A_{\rm obs})\) using a known unit mass in the force–acceleration subtraction. This is a fixed invertible readout of the three physical channels. For the first push the three errors in these coordinates are independent with standard deviations \((1/2,1,1)\), and for the second they are \((1,1/2,1)\). Take true accelerations 1 and 2 respectively. At truth the two whitened lines both point along the last coordinate. Their derivative metrics are \(\operatorname{diag}(4,1)\) and \(\operatorname{diag}(1,4)\). Thus
\[
 A=\operatorname{diag}(8,17),\qquad B=\operatorname{diag}(13,22).
\]
The first marginal requires \(\eta=8/13\) and the second \(17/22\). If the first width is calibrated, the second is too wide by a factor \(\sqrt{(17/22)/(8/13)}\), about 1.12. The discrepancy remains under arbitrarily many repeats of this same valid design. Declaring the two bodies' allocation and total in another smooth coordinate system changes the components of the tensors and preserves the impossibility of making them proportional.

## Which weighted covectors represent odds?

At a candidate write \(K_i=\|y_i\|^2-(u_i^\top y_i)^2\). The undamped score is \(-\tfrac12\sum_i dK_i\). Consider local discount rules of the explicit form
\[
 d\log E=-\tfrac12\sum_i w_i(K_1,\ldots,K_n)\,dK_i .
\]
The weights are continuously differentiable on the open positive cost orthant, continuous to any faces where used. They depend on candidate residuals alone for this theorem. Inputs that are fixed across candidates may be held as external parameters if the same set of experiments is available for each value. The question is whether these local instructions produce an unambiguous odds ratio between any two labelled mass pairs on every possible reading list.

**Classification.** This holds iff \(\partial_jw_i=\partial_iw_j\) throughout the cost orthant. Equivalently, \(w=\nabla\Psi\) for a single scalar cost \(\Psi\) and \(E\propto\exp[-\Psi(K)/2]\) wherever normalization is finite.

For necessity it suffices to consider the common three-channel unit-noise instrument at any interior mass pair. The derivative \(du\) identifies the two mass directions with the tangent plane of its unit vector. Choose \(y_i=t_i u+v_i\), \(v_i\perp u\), \(\|v_i\|^2=k_i>0\). Then \(dK_i=-2t_i v_i^\top du\). Any prescribed positive cost list is realized, and for any chosen \(i,j\) the corresponding differentials can be made independent, with every other \(t_l=0\). The exterior derivative of \(\sum w_i dK_i\) at that point therefore vanishes only if the two cross derivatives agree. A nonzero exterior derivative yields nonzero circulation around arbitrarily small rectangles of actual mass pairs. Conversely, symmetry gives a potential on the convex cost orthant, for example by integrating along straight segments from an interior base point. Pulling its differential back to mass space proves sufficiency. Positive finite factors times the supplied prior then determine a law. [Numerical loop and local tensor checks](../tests/test_e12_joint_score_geometry.py)

**Consequence for one common discount.** If \(w_i=d(K_1,\ldots,K_n)\) for all \(i\), consistency of odds requires every first partial derivative of \(d\) to agree. Hence \(d\) can depend only on the sum of the candidate costs (at a fixed block length). For example the seemingly reasonable symmetric downweight \(d=1/(1+K_1 K_2)\) produces different odds depending on the path whenever \(K_1\ne K_2\) and their differentials are independent. On the cost rectangle \([1,2]\times[2,3]\) the circulation of \(d\,(dK_1+dK_2)\) with positive orientation is strictly negative in the interior, since \(\partial_1d-\partial_2d=(K_1-K_2)/(1+K_1K_2)^2\). Physical data can realize the same defect on a local patch by the construction above. If the same residual-only prescription is used for a single body, is nonnegative, and neutral to appending an observed zero, the [total-cost obstruction](e12-reciprocal-evidence.md) applies. A block-wide weight computed once from the data and *held fixed as mass varies* is outside this derivative argument.

Even different conservative row weights have a restriction. Suppose \(0\le w_i\le1\) and write \(\bar w=n^{-1}\sum_i w_i\). For any positive vectors \(k,l\) and \(t>0\),
\[
 n\left|\int_0^t[\bar w(k+s{\bf1})-\bar w(l+s{\bf1})]\,ds\right|\le\|k-l\|_1.
\]
Integrating the gradient first along the two parallel segments rewrites the left-hand difference as \([\Psi(k+t{\bf1})-\Psi(l+t{\bf1})]-[\Psi(k)-\Psi(l)]\). Each bracket lies between \(-\sum_i(l_i-k_i)_+\) and \(\sum_i(k_i-l_i)_+\), by the component bounds on the gradient. Their difference obeys the claim. If the average weight depends only on the deviations \(k_i-\bar k\), it stays constant along each segment. Letting \(t\) grow forces the same average weight for every deviation pattern. The individual weights may still redistribute influence. For example \(\Psi=\sum_i\log(1+K_i)\) gives bounded, unequal and integrable weights. Neither this example nor local matrix calibration selects a global discount.

The two tests on a proposed equation are independent. Correct local coverage constrains its curvature to a matrix on regular data. Path independence constrains its first derivatives on every finite data set. A proposed probability equation for several masses must satisfy both wherever it claims both interpretations.

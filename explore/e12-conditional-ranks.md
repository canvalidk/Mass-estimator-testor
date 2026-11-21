# How much mass can the signed reciprocal score assign?

Under a candidate mass, [T45](../definition/theorems/T45.md) gives a conditionally standard normal signed score. Its Gaussian rank is therefore exactly uniform at the truth. Could it already be the CDF of the mass? Here the derivative of that rank determines the answer, an exact conditional selection probability, and a positive law constructed from the part of the rank that increases. The construction uses one isotropic scalar channel pair in noise units and any number of independent signed pushes. No push population is supplied.

## Where the exact ranks run forward

Stack the observed force and acceleration rows into \(Y\). For \(u=(\sin\theta,\cos\theta)\), \(v=u'\), \(\theta=\arctan(m/s)\), put \(U=Yu,\ V=Yv\) and \(g(\theta)=\|U\|\). At truth \(U\sim N(r,I_n)\), \(V\sim N(0,I_n)\) independently, with any fixed drive vector \(r\). T45's signed covector score and its candidate CDF rank are
\[
 Z=g'(\theta)=U^\top V/g,\qquad C(\theta)=\Phi(-Z).
\]
At every fixed true angle \(g>0\) almost surely. Conditional on \(U\), \(Z\) is exactly standard normal. Writing \(G=Y^\top Y\), \(L\ge l\ge0\) for its eigenvalues, a second derivative gives
\[
 g''=(\|V\|^2-Z^2)/g-g=\det(G)/g^3-g .
\]
Thus \(C\) increases where \(u^\top Gu>\sqrt{Ll}\) and decreases where the inequality reverses. For \(L>l\), in each projective period it increases on the arc about the leading eigenline with half width \(\arctan((L/l)^{1/4})\), interpreted as \(\pi/2\) at rank one. This follows by putting \(g^2=L\cos^2h+l\sin^2h\) in the eigenbasis with angle offset \(h\). At the two arc ends the signed score has opposite signs and absolute value \(\sqrt L-\sqrt l\) (use the limits from within the arc at rank one). The entire increasing variation over one period is therefore
\[
             W=2\Phi(\sqrt L-\sqrt l)-1 .
\]
The identical amount runs backwards. At rank one the entire decrease is one downward jump at the perpendicular line, where the score otherwise has no definition. At an eigenvalue tie the rank is constant \(1/2\). Only part of an increasing arc need lie within the positive physical interval \(I=(0,\pi/2)\). For instance \(G=\operatorname{diag}(9,1)\) gives \(C=1/2\) at both ends of \(I\) and \(C<1/2\) throughout its interior. These pointwise exact tests cannot all supply the ordinates of a single CDF.

**An exact selective statement.** At the true angle write \(R=\|U\|^2,\ J=\|V\|^2-Z^2\). For \(n\ge2\), \(R,J,Z\) are mutually independent with laws \(\chi'^2_n(\|r\|^2),\chi^2_{n-1},N(0,1)\). Conditional on \(U\), the component of \(V\) along \(U\) and the orthogonal projection are independent; neither the distribution of its signed component nor the projected norm depends on \(U\). The true candidate is on the increasing branch precisely when \(R>J\), since \(\det G=RJ\). Its probability is exactly \(\Pr\{F'_{n,n-1}(\|r\|^2)>(n-1)/n\}\). Even conditional on this event, \(C(\theta_0)\) is uniform. At \(n=1\) the branch condition holds almost surely at truth and \(J=0\). These are candidate-specific selection events. They supply no probabilities that one of several candidates is true.

## A positive law given the increasing variation

For \(0<d\le1\) constant in \(\theta\) on the given data, write \(C_d=\Phi(-d g')\). Take its positive variation on \(I\), the measure with density
\[
       H_d(\theta)=d\phi(dg')\,[g-Ll/g^3]_+ .
\]
Use zero at a zero of \(g\). The only possible jump of \(C_d\) is downward. The total \(M_d=\int_I H_d\) lies between zero and \(2\Phi(d(\sqrt L-\sqrt l))-1<1\) on finite data. Given the project's declared prior \(\pi\) on this angle, define
\[
            P_d(d\theta)=H_d(\theta)\,d\theta+(1-M_d)\pi(d\theta).
\]
There is a concrete sampling interpretation. Draw one uniform rank \(p\) and look for an upward crossing \(C_d(\theta)=p\) on the positive mass arc. In a full projective period almost every rank has either one such crossing or none. Hence the probability that an allowable crossing exists is exactly \(M_d\), by change of variables on the increasing arc. This is the **largest possible success probability** of any procedure required to realize the sampled rank with an upward crossing in the physical domain. Accepting every available crossing fixes its mass distribution to \(H_d\,d\theta\). Drawing independently from the declared prior on failure gives the displayed law. Other fallback rules or decisions to use decreasing crossings change the inference premise.

At \(d=1\) the procedure retains the underlying exact pivot levels. Its resulting CDF has different ordinates and has no claim to exact probability-rank calibration. Its failure on fixed conflict is measurable: multiply rows with Gram \(\begin{pmatrix}5&4\\4&5\end{pmatrix}\) by \(\sqrt t\). The original data have incompatible directions, but \(P_1\) concentrates at \(\pi/4\) as \(t\to\infty\).

Here is one choice using only the data that deals with that failure. Let \(Q=\inf_{\theta\in I}\|Yv_\theta\|^2\), \(D=\min(G_{11},G_{22})\) and \(B=4\sum_i(1-e^{-\|Y_i\|^2/2})\). Choose
\[
 d=\min\{1,\sqrt{BL}/Q,\sqrt{1+D^2}/\sqrt L\}.
\]
Omit the second term at \(Q=0\), and return \(\pi\) at \(L=0\). The fourth scalar \(B\) retains reading-level information; the previously established scatter obstruction applies if it is omitted. Other positive factors in the definition of \(B\) above 2 have the same limiting arguments. Neither the multiplier nor the conflict discount has been derived uniquely. The point of this witness is that the **measure being discounted** is determined by differentiating the signed reciprocal rank. [Four-coordinate result](e12-joining-gram.md) · [Earlier discount construction](e12-scalar-proof.md)

**Finite data and shrinking errors.** The law is proper, invariant under row permutations and individual reversals, equivariant under channel exchange, and exactly neutral on adding a measured zero. For each finite block it is continuous in total variation: the density is bounded above by \(d\sqrt L/\sqrt{2\pi}\), and the displayed formula is continuous pointwise except possibly at the isolated perpendicular to a limiting rank-one matrix. At a positive eigenvalue tie it converges to zero, and at the zero matrix the bound suffices. In the case \(Q\downarrow0,L>0\) the middle term eventually exceeds one. Ordinary dominated convergence now applies.

For any fixed nonzero shape incompatible on the *closed* positive arc, scaling \(Y=\sqrt t\,Y_0\) gives \(Q=tQ_0>0,\ L=tL_0,\ B\to4k\) for \(k\) nonzero rows. Then \(d\sqrt t\to\sqrt{4kL_0}/Q_0\) (both diagonals are positive). Consequently the law converges in total variation to the formula above for the unscaled shape with that finite multiplier. This follows also directly by homogeneity of \(g\) and both derivatives. The missing weight is strictly positive, providing a fixed prior contribution and preventing concentration, including for a full-support half-Cauchy prior. Isotropic ties return the entire prior. The remaining incompatible cases on the open mass domain are rank-one shapes exactly at an axis. Here \(Q=D=0\) and \(d\sqrt t\to1/\sqrt{L_0}\), yielding the same finite-multiplier conclusion. Other cases can also return the prior when their rising arc misses the interior. Informativeness on every fixed conflict is a stronger condition than retaining spread.

**Identification.** Under the [scalar proof's sampling experiment](e12-scalar-proof.md#experiment-and-quantifiers), fix an interior truth and let \(\Delta=\|r\|^2,\ I_n=\Delta^2/(n+\Delta)\to\infty\). Then the selected \(d=1\) with probability tending to one. For \(n\to\infty\), the independent bounded terms in \(B/4\) each have mean at least \(1/2\) (Gaussian integration gives \(1-e^{-r_i^2/4}/2\)), and \(Q\le\|V_{\theta_0}\|^2=n+O_p(\sqrt n)\). Since \(L\ge Q\), \(B>Q\) suffices for the middle cap to be inactive. For bounded \(n\), \(Q=O_p(1), L\to\infty\) and \(B\) stays bounded away from zero in probability because at least one signal diverges. The diagonal cap is inactive in both regimes: each diagonal is of order \(n+\Delta\) with probability tending to one at the fixed interior angle, while \(L=O_p(n+\Delta)\).

In the truth-aligned frame, the sample scatter differs from \(\operatorname{diag}(n+\Delta,n)\) by \(O_p(\sqrt{n+\Delta})\) in operator norm by the variances of its three entries. Thus the leading angle \(\hat\theta\) tends in probability to truth, \(L/(n+\Delta)\to1\) and \((L-l)/\Delta\to1\). Put \(A=(L-l)/\sqrt L\sim\sqrt{I_n}\) in probability. In leading-angle coordinates \(h\),
\[
 -g'(\hat\theta+h)=\frac{(L-l)\sin h\cos h}{\sqrt{L\cos^2 h+l\sin^2 h}} .
\]
At \(h=z/A\) the right side tends to \(z\), uniformly for bounded \(z\), and lies on the increasing branch whose half width is at least \(\pi/4\). Its derivative with respect to \(z\) tends uniformly to 1 there: dividing numerator and denominator by \(\sqrt L\) gives \(A\sin h\cos h/\sqrt{1-(1-l/L)\sin^2h}\), with the coefficient \(1-l/L\) uniformly in \([0,1]\). Both physical endpoints have fixed positive angular distance from truth. Consequently the rescaled increasing density tends uniformly to \(\phi(z)\) on each compact interval. Its mass between \(\hat\theta\pm K/A\) tends to \(2\Phi(K)-1\) for each fixed \(K\). Letting \(K\) grow gives \(M_d\to1\) and convergence in total variation, in probability, of the output law of \(A(\theta-\hat\theta)\) to \(N(0,1)\).

The sampling center has the matching fluctuation. In the truth-aligned frame \((L-l)\sin(2(\hat\theta-\theta_0))=2U^\top V\). Since \(\hat\theta\to\theta_0\), \(\|U\|/\sqrt L\to1\) in probability, and \(U^\top V/\|U\|\) is exactly standard normal at truth, \(A(\hat\theta-\theta_0)\) tends in law to a standard normal. This proves median consistency and nominal asymptotic coverage of all fixed central quantile intervals on the specified identifying arrays. [Exact CDF calculator](e12_conditional_ranks.py) · [Checks](../tests/test_e12_conditional_ranks.py)

The strictly decreasing arcs and the missing rank weight isolate two choices that exact covector calibration leaves open. The increase itself is fixed once the original signed quantile is specified. Its completion into a law and the amount of damping on conflicting readings still require a premise.

For the incompatible records \((4,1),(1,4)\), with equal channel errors and a unit half-Cauchy prior, the discounted construction gives median 1. With channel errors 1 its central 95% mass interval is \([0.177,5.664]\). With errors \(0.1\) the interval is \([0.391,2.556]\). With errors \(0.001\) the endpoints still equal those values to three decimals. These are probability intervals of the specified construction. Their finite-sample sampling coverage is a separate property.

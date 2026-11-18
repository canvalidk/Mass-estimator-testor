# A conservative scalar interval and the influence of a noise block

**Scope.** One scalar instrument, independent Gaussian errors of the declared size in both channels, fixed signed signals, and positive mass. The object here is a companion *report* at an explicit confidence level \(0<p<1\). Its promise is coverage **at least** \(p\), and its width has an operational consequence when pooling with an independent point. This is a possible statement for G28. It does not identify the intended content of property N or determine the raw law. The standard construction is conditional Gaussian score-test inversion, as derived in [T45, §§1–3](../definition/theorems/T45.md), followed by enlargement of a confidence set. The existing E2 results on score-state addition, whole-domain probability and the impossibility of vanishing null influence are [on their exploration branch](https://github.com/canvalidk/Mass-estimator-testor/blob/d450111c4716ff13690933a96f95e25356170a58/explore/e2-null-combination.md#L255-L375).

## Construction and coverage

Write \(s=\sigma_F/\sigma_a\), \(\theta=\arctan(m/s)\), \(u_\theta=(\sin\theta,\cos\theta)\), \(v_\theta=(\cos\theta,-\sin\theta)\). The finite positive mass domain corresponds to \(0<\theta<\pi/2\). Whiten each reading to \(Y_i=(F_i/\sigma_F,a_i/\sigma_a)\), so at truth \(\theta_0\) it has distribution \(r_i u_{\theta_0}+N(0,I_2)\). Write \(G=\sum_iY_iY_i^\top\), \(S_\theta=v_\theta^\top G u_\theta\), \(V_\theta=u_\theta^\top G u_\theta\), and \(z=\Phi^{-1}((1+p)/2)\). Define the accepted masses

\[
C_p=\{m>0:S_\theta^2\le z^2 V_\theta\}.
\]

At \(V_\theta=0\), \(S_\theta=0\) and the candidate is accepted. This is T45's explicit convention, including the empty or all-observed-zero experiment. At a fixed true candidate in a nonempty experiment the variance is positive almost surely. Conditional on all its parallel projections \(U_i=u_{\theta_0}^\top Y_i\), the perpendicular projections are independent standard normals, so \(S_{\theta_0}\mid (U_i)\sim N(0,\sum_iU_i^2)\). This proves exact pointwise conditional coverage for \(C_p\). The nominal probability concerns repetition of the declared experiment.

Take any proper raw mass law \(P_D\) with a positive continuous density on \((0,\infty)\). Let \(J_q\) be its central interval of probability \(q\), with \(p\le q<1\). A simple choice is \(q=p\). The finite coverage and null calculations also allow \(q\) to be a function of the observed block. Report

\[
H_{p,q}(D)=[\inf(C_p\cup J_q),\sup(C_p\cup J_q)]\cap(0,\infty).
\]

The infimum may be zero and the supremum infinite; neither is a possible mass. Every \(J_q\) has positive finite endpoints. Including endpoints of the hull inside the mass domain defines the interval even when \(C_p\) is disconnected. Computing just the smallest and largest accepted points on a finite mass grid would miss its unbounded ends.

* \(H_{p,q}\) covers with probability at least \(p\), also conditionally on the true \(U_i\), by containment. The empty experiment gives the whole domain and coverage one.
* If \(P_D\) respects exact-zero padding, permutations, mass units and channel swap, so does \(H_{p,q}\) for fixed \(q\), or for a choice of \(q(D)\) invariant under the same transformations. The Gram is additive with zero rows neutral, and taking the hull commutes with monotone rescaling and reciprocation.
* On a fixed conflicting series any nonvanishing raw interval width is retained by inclusion. This does not itself require the raw law to discount conflicts or prove that the reported interval measures the repeated variability of the point.

**Strict finite-sample overcoverage of every connected enlargement.** For every nonempty finite Gaussian series, every fixed signal vector and every interior truth \(\theta_0\), any interval containing all of \(C_p\) has unconditional coverage strictly greater than \(p\). Choose a first standardized row \(y=(f,a)\) just inside \((0,z)\times(-z,0)\), near \((z,-z)\), with \(u_{\theta_0}^{\top}y\ne0\). It can be chosen so \(v_{\theta_0}^{\top}y>z\), since at that corner the projection is \(z(\cos\theta_0+\sin\theta_0)>z\). Take all other rows zero initially. Wherever the parallel projection is nonzero the rank-one absolute score is \(|v_\theta^\top y|\). Its boundary values \(|f|,|a|\) are strictly less than \(z\), whereas its value at truth exceeds \(z\). The variances at all three angles are positive, so these strict comparisons persist on an open neighborhood of the entire observation vector. Both tails of the mass domain contain accepted points there. Every containing interval covers the truth there although the score set rejects it. The joint Gaussian law has positive density throughout that neighborhood for any fixed signals. Adding its positive probability to the exact score-set coverage proves the strict inequality. Thus finite-sample equality with the nominal level already excludes this entire class of connected enlargements.

## Contraction at the full scalar identifying rate

Let \(m_0\) be fixed in the interior and \(\Delta=\sum_i r_i^2\). Under \(\Delta^2/(N+\Delta)\to\infty\) the entire score acceptance set lies in any fixed angular neighborhood of \(\theta_0\) with probability tending to one.

Indeed \(G=NI+\Delta u_0u_0^\top+E\), with \(\mathbb E\|E\|_F^2=6N+6\Delta\) by independent Gaussian second and fourth moments ([O19](../problems/O19.md)). Uniformly on the closed positive-angle quadrant outside that neighborhood, \(|(v_\theta^\top u_0)(u_\theta^\top u_0)|\) is bounded below by a positive constant: the perpendicular line to \(u_0\) is strictly outside the quadrant. Thus \(|S_\theta|\ge c\Delta-\|E\|\) and \(V_\theta\le N+\Delta+\|E\|\). With probability tending to one the first quantity exceeds \(z\sqrt{N+\Delta+\|E\|}\) everywhere on this complement. This proves the assertion including fixed \(N\) with growing signal and arbitrary allocations of the energy.

If the raw law also concentrates in probability at \(m_0\), every fixed \(J_q\) contracts there, as does \(H_{p,q}\). Its log-width tends to zero in probability and it eventually has two finite endpoints with probability tending to one. The same holds for a variable \(q\to p\) in probability by sandwiching its quantiles between fixed ranks. This argument requires raw-law concentration; point consistency alone does not imply it.

**Nominal asymptotic coverage when the two intervals agree locally.** Write \(L_g\ge l_g\) for the eigenvalues of \(G\), \(h=L_g-l_g\), and \(A=h^2/L_g\), with \(A=0\) when \(L_g=0\). On identifying arrays \(A/I\to1\), where \(I=\Delta^2/(N+\Delta)\). With probability tending to one the principal line \(\hat\theta\) is near the interior truth. At offset \(\delta=\theta-\hat\theta\) the absolute score is
\[
\frac{|S_\theta|}{\sqrt{V_\theta}}
=\sqrt A\,\frac{|\sin\delta\,\cos\delta|}{\sqrt{1-(h/L_g)\sin^2\delta}}.
\]
In a small neighborhood of zero the fraction divided by \(|\delta|\) tends to one uniformly over \(0\le h/L_g\le1\). Outside a small neighborhood of \(\hat\theta\) the preceding exclusion argument applies with probability tending to one. Inside it the score is bounded below by a fixed positive constant times \(\sqrt A|\delta|\). On bounded \(x=\sqrt A\delta\) the statistic converges uniformly in absolute difference to \(|x|\). Thus the two extremes of \(C_p\) in the angle coordinate are \(\hat\theta+(\mp z+o_p(1))/\sqrt I\): for each \(\epsilon>0\) the rescaled set eventually includes \([-z+\epsilon,z-\epsilon]\) and lies in \([-z-\epsilon,z+\epsilon]\) with probability tending to one.

For [the pooled witness](e12-pooled-witness.md) the local law of \(\sqrt I(\theta-\hat\theta)\) tends to the standard normal, and \(\sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1)\). Any raw procedure with these same two limits has \(J_q\) endpoints \(\hat\theta+(\mp z+o_p(1))/\sqrt I\) whenever \(q\to p\) in probability. This follows at random ranks by trapping them between fixed nearby ranks and using continuity of the normal quantile there. Hence the hull has those endpoints too and sampling coverage tending to \(p\). For fixed \(q>p\) the raw endpoints instead dominate to first order, and under these local normality assumptions coverage tends to \(q\).

For a *deterministic incompatible* series the same score can be much too narrow to represent the whole-series T30 requirement. In the O14 witness with forces \((0.2,1,1.4)\), accelerations all one and common channel SD \(\sigma\),

\[
G=\sigma^{-2}
\begin{pmatrix}3&2.6\\2.6&3\end{pmatrix}.
\]

The positive principal line is \(m=1\). Away from any neighborhood of it, the score magnitude grows as \(\sigma^{-2}\), uniformly, while the standard deviation grows only as \(\sigma^{-1}\). Consequently \(C_p\) alone contracts to \(1\) as \(\sigma\downarrow0\). This explains the separate role of \(J_q\) for a raw law that stays spread on this input. It is a claim about conditioning on a fixed conflicting dataset as the declaration changes, distinct from sampling coverage at fixed signal.

## What the hull reports on random true-null data

When \(\Delta=0\), write \(f,a\) for the \(N\)-vectors of independent standardized force and acceleration noise. For \(N\ge1\) both diagonal Gram entries are positive almost surely. Except for endpoint equalities of probability zero, the score set reaches arbitrarily close to the zero-mass boundary iff \(G_{12}^2<z^2G_{22}\), and arbitrarily high masses iff \(G_{12}^2<z^2G_{11}\). This follows by continuity and evaluating the statistic at the limiting angles. The finite interior interval \(J_q\) changes neither fact. Hence the hull is the entire domain iff both boundary tests pass, even if interior candidates were excluded from \(C_p\).

Each boundary passes with probability \(p\): conditional on \(a\), \(f^\top a/\|a\|\) is standard normal, and interchange gives the other boundary. Let
\[
T=\frac{(f^\top a)^2}{\|f\|^2\|a\|^2},\qquad h_N(t)=F_{\chi^2_N}(z^2/t)
\]
with \(h_N(0)=1\). The two squared lengths are independent \(\chi^2_N\), independent of the two directions. Conditional on \(T\), therefore, both boundary tests pass with probability \(h_N(T)^2\). With \(b_N=\Pr_0\{H_{p,q}=(0,\infty)\}\),
\[
b_N=\mathbb E[h_N(T)^2],\quad \mathbb E[h_N(T)]=p,\quad p^2\le b_N\le p .
\]
The bounds use Jensen's inequality and \(0\le h_N\le1\). For \(N=1\), \(T=1\) and \(b_1=p^2\). For \(N\ge2\), \(T\) has the \(\operatorname{Beta}(1/2,(N-1)/2)\) law: fix one direction by rotational symmetry and square the first coordinate of an independent uniform unit vector. This gives a one-dimensional quadrature with no fitted signals.

As \(N\to\infty\), \(G_{12}/\sqrt N\) tends in law to a standard normal and \(G_{jj}/N\to1\) in probability. Therefore \(b_N\to p\). At 95% the **interval hull** tends to a 95% probability of reporting the entire mass domain. E2 gives \(1-\exp(-z^2/2)=85.35\%\) for the **original score set**: its acceptance must hold at *every* angle. The two probabilities differ because the hull fills gaps. For example \(G=\operatorname{diag}(16,1)\) accepts both boundaries but rejects \(m=s\) at 95%.

| True-null readings | Whole hull at \(p=.95\) | Two finite endpoints |
| --- | ---: | ---: |
| 1 | 0.90250 | 0.00250 |
| 5 | 0.91402 | 0.01402 |
| 20 | 0.92996 | 0.02996 |
| Limit | 0.95000 | 0.05000 |

## A finite bound on influence over an independent point

Here is one explicitly limited operation using the report. Let an independent route give a finite positive point \(B\) and a fixed finite positive external weight \(w_0\). Let the block give a finite positive point \(M\) and \(H_{p,q}\) with endpoints \(L,U\). If both endpoints are positive and finite give it weight \(w=1/\log^2(U/L)\), otherwise \(w=0\). The proper-law interval has strictly positive width for every finite dataset. Form the weighted mean of the two *log* points and exponentiate. This operation respects mass units and reciprocation. A reported \(w\) is a specified numerical weight, without a claim that it estimates inverse variance.

On random true-null Gaussian data the external answer can change only if neither boundary passes. By inclusion-exclusion,
\[
\Pr_0\{w>0\}=1-2p+b_N\le1-p.
\]
The probability tends to \(1-p\) and already has this upper bound at **every** sample size, for any independent external point. Under the full scalar identifying condition, if the raw law concentrates and \(M\to m_0\), the preceding contraction gives \(w\to\infty\) in probability. The combined log point differs from \(\log M\) by \(w_0(\log B-\log M)/(w_0+w)\to0\) in probability. The external error is tight because its law is fixed. A valid raw law and a consistent point thus suffice for this weaker, finite-error influence promise without a stronger signal rate. This bounds the probability of *any* change on pure noise. It does not bound the size of a change on that event. E2's obstruction is to the stronger vanishing-in-probability requirement for all identifying sequences.

## What remains unselected

Even fixing raw law, point, \(p\), asymptotically nominal coverage on identifying arrays and the finite influence bound leaves distinct reports. For any fixed integer \(k\ge0\) use
\[
q_k(D)=p+(1-p)\frac{k}{1+k+A(D)} .
\]
Here \(A=h^2/L_g\), set to zero at \(L_g=0\). Each choice respects the same invariances, has \(p\le q_k<1\), and \(q_k\to p\) along every identifying array. Thus every report has the finite coverage bound, whole-hull probability, influence bound and, under the preceding local normality assumptions, the same asymptotically nominal \(p\) coverage. On a scalar reading \(Y=R(1,1)/\sqrt2\) with \(R>\sqrt2 z\), the score is accepted exactly when \(R|\sin(\theta-\pi/4)|\le z\), so \(C_p\) is bounded strictly inside the positive domain and \(A=R^2\). As \(k\) increases on this fixed sample, \(q_k\uparrow1\). Because the same raw law has full support and positive continuous density, the raw interval eventually contains \(C_p\) and then strictly expands with each \(k\). This gives infinitely many distinct reports on reachable data. Finite-sample *equality* of coverage or a separate requirement comparing interval widths could change this conclusion.

The hull also loses state required for recomputing it after pooling. Both \(G=I_2\) and \(G=3I_2\) report the whole domain for every level. They have the same normalized Gaussian angular evidence curve (constant) and can give the same raw law under the same prior. Append \(Y=(3,-1)\). The new \(G_{12}^2\) is 9 in both blocks and the diagonal entries become \((10,2)\) and \((12,4)\), respectively. At 95% the first loses its zero-mass end, while the second still reaches both ends. Any positive finite \(J_q\) preserves this distinction. Hence even all original hulls do not determine the updated hull. Same-instrument score pooling can retain \(G\) and add the new outer products; the raw procedure may need its own additional state. The external point operation above supplies neither an uncertainty law for its output nor a general solution for joint composition or independent routes.

## The point readout also remains a choice

Fix the raw mass law \(P\) with a positive continuous density and let \(Q_a(P)\) be its mass quantile of rank \(a\). For fixed \(0<a\le1/2\), take
\[
 M_a(P)=\sqrt{Q_a(P)Q_{1-a}(P)}.
\]
Every member is finite and positive without moment assumptions. Rescaling mass by \(t>0\) multiplies both quantiles by \(t\). Channel interchange sends each quantile to the reciprocal of the complementary quantile. Thus \(M_a\) scales by \(t\) and reciprocates on the swap. On all-zero observations suppose the raw law is the prior and the prior is invariant under \(m\mapsto c^2/m\). Then \(Q_a Q_{1-a}=c^2\) and every readout returns \(c\). Exact measured-zero neutrality is inherited from the raw law on every dataset.

On the interior identifying arrays above, the angular quantiles at paired ranks are \(\hat\theta+(\pm\Phi^{-1}(1-a)+o_p(1))/\sqrt I\). Apply a first-order Taylor expansion of \(g(\theta)=\log(s\tan\theta)\) to both quantiles. Their average log mass is \(g(\hat\theta)+o_p(I^{-1/2})\) because the two leading offsets cancel and \(\hat\theta\to\theta_0\) in probability. The same expansion holds for the median \(M_{1/2}\). All these points agree to first order and are consistent.

They differ on reachable data. Take \(s=c=1\), a uniform angular prior, and one observed row \(\sqrt8\,(\sin 0.63,\cos 0.63)\). Its face-value Gaussian angular density is proportional to \(e^{-4\sin^2(\theta-0.63)}\), also the pooled witness law here. Its mass quartiles are approximately \(0.446553,0.763997,1.249750\), giving \(M_{1/4}=0.747047\) and \(M_{1/2}=0.763997\). The angular factor has absolute derivative at most four, so a composite midpoint integral on \([0,x]\) with \(n\) equal cells has absolute error at most \(x^2/n\). With \(n=16384\), the signs of the cumulative integral minus the target fraction of total mass certify \(Q_{1/4}<0.45,\ Q_{1/2}>0.760,\ Q_{3/4}<1.254\), and \(0.45\times1.254<0.760^2\). Consequently the units, swap, zero-block center and regular first-order behavior in this section leave [G6's readout choice](../definition/proof.md) open even for a single fixed raw law.

Finite identities, edge cases and the numerical entries are checked in [the report tests](../tests/test_e12_report.py).

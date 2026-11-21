# What physical composition determines about the joint mass law

L5 concerns different bodies that share a physical acceleration; L14 has the dual experiment at shared force. Repeated observations of one body instead join through one common mass. Here one physical experiment fixes a joint tangent law and gives an exact criterion for losing information about a lump when forces are summed. Outside the tangent experiment, even both part laws and both physical composite laws leave distinct *exchangeable* joint laws possible. [L5 and L14](../definition/premises.md) · [O2's Gaussian coarsening result](../problems/O2.md) · [T51](../definition/theorems/T51.md)

## The reciprocal of one shared measured acceleration

Use one spatial dimension, \(b\ge 2\) labelled bodies, and independent pushes indexed by \(i\). The source observation is \((F_i,A_i)\in\mathbb R^b\times\mathbb R\), with declared positive-definite joint covariance \(\Sigma_i\). All bodies share the same unknown physical acceleration \(a_i\) for that push. At a vector candidate \(m\in(0,\infty)^b\), the candidate means are \((a_i m,a_i)\). The same measured acceleration is present only once. Define
\[
 D=[I_b,-m],\quad d_i=F_i-mA_i,\quad
 V_i=D\Sigma_iD^\top,\quad C_i=\operatorname{Cov}(A_i,d_i),\quad
 x_i=A_i-C_iV_i^{-1}d_i .
\]
The shortest noise-metric correction consistent with the graph subtracts \(\Sigma_iD^\top V_i^{-1}d_i\) from the data. Its acceleration component leaves \(x_i\) as the fitted common drive. At truth \(d_i\sim N(0,V_i)\), independent of \(x_i\), since \(\operatorname{Cov}(x_i,d_i)=0\). This is E12's reciprocal tangent with one unknown common drive and multiple physical masses. [Single-mass reciprocal tangent](e12-tangent-pairings.md)

Conditional on the \(x_i\), consider the **tangent shift** \(\delta\in\mathbb R^b\) with \(d_i\sim N(x_i\delta,V_i)\) independently. Zero shift is the exact conditional model at the true candidate. A finite shift extrapolates the fixed tangents. For at least one nonzero \(x_i\) put
\[
 H=\sum_i x_i^2 V_i^{-1},\qquad u=H^{-1}\sum_i x_i V_i^{-1}d_i .
\]
Then \(u\sim N(\delta,H^{-1})\), independent of all orthogonal residuals by Gaussian projection. All-zero design has no reciprocal displacement. For one active push \(u=d_i/x_i\) and its cross-covariance is \(V_i/x_i^2\); the noisy acceleration usually couples the two inferred displacements even when the force errors are independent.

**Conditional uniqueness.** Fix the design and orthogonal residuals. Require translation equivariance of the returned tangent law and exact one-sided coverage at every rank in every linear contrast. For a direction \(c\ne0\), a rank-\(p\) endpoint must equal \(c^\top u + b_p(c)\). Its probability of covering \(c^\top\delta\) is \(\Phi(b_p(c)/\sqrt{c^\top H^{-1}c})\). Thus every projection has the corresponding normal law, which determines uniquely the joint report \(N(u,H^{-1})\) by its characteristic function. The Gaussian meets these requirements. Independent groups of pushes add their precision matrices and precision-weighted displacements. Exact-rank coverage for every contrast is an explicit condition of this result. [Scalar tangent result](e12-tangent-diagnostic.md) · [Cramér and Wold](https://doi.org/10.1112/jlms/s1-11.4.290)

The same identification follows if exact ranks are prescribed only for the individual bodies and every positively integer-weighted lump of the *same* bodies. The individual tangent marginals must then be Gaussian. Hölder's inequality applied to their Gaussian exponential moments makes the moment generating function of any joint linear combination finite everywhere. For each positive integer degree \(n\), the difference between the reported and Gaussian moments of \((c^\top\delta)^n\), with \(\delta\) drawn from each report, is a homogeneous polynomial in \(c\). Coverage makes this difference zero on all positive integer rays, hence on the positive rational orthant by homogeneity, then everywhere by continuity and the polynomial identity theorem. Moment determinacy of every linear combination and the projection argument give the same joint normal. Ordinary part and unit-weight lump coverage alone are weaker hypotheses.

## Does summing force observations preserve the reciprocal law of the lump?

The physical lump has displacement \(c^\top\delta\) with \(c=\mathbf1\). For any fixed nonzero \(c\) the pushforward variance is \(v_{\rm fine}=c^\top H^{-1}c\). Now *in the same conditional tangent experiment*, retain the same fitted drives \(x_i\) and only \(y_i=c^\top d_i\). The scalar experiment has
\[
 v_{\rm sum}=\left(\sum_i x_i^2/(c^\top V_i c)\right)^{-1}.
\]
Then \(v_{\rm fine}\le v_{\rm sum}\), with equality **if and only if**
\[
 \frac{V_i c}{c^\top V_i c}
 \quad\text{is identical for every active push }i .
\]
Under equality the two centers agree on every dataset as well. The condition concerns the direction of instrument covariance applied to the specific physical combination.

To see this, the fine reciprocal coefficient on \(d_i\) is \(w_i=x_i V_i^{-1}H^{-1}c\), whereas the summed experiment uses \(z_i=x_i v_{\rm sum}c/(c^\top V_i c)\). Each obeys \(\sum x_i w_i=\sum x_i z_i=c\). Orthogonality in the joint noise metric gives
\[
 v_{\rm sum}-v_{\rm fine}=\sum_i(z_i-w_i)^\top V_i(z_i-w_i).
\]
Equality forces every coefficient to coincide and hence the displayed matching condition. Conversely if all normalized \(V_i c\) equal \(v\), then \(Hv=c/v_{\rm sum}\) and the same coefficients follow. This is generalized least squares with a shared nuisance drive fixed by reciprocal graph projection. Complementary sensors are an elementary strict example. For \(m=(1,1)\), \(x_1=x_2=1\), let the two physical force covariance blocks be \(\operatorname{diag}(1,9)\), \(\operatorname{diag}(9,1)\), and acceleration noise have independent variance \(1/10\) in both. The coarse reading suppresses which body had the precise sensor, which the fine reciprocal estimator uses. [Aitken's least-squares theorem](https://doi.org/10.1017/S0370164600014346)

**One push before conditioning.** Coarsening the *raw* experiment to \((c^\top F_i,A_i)\) refits its drive,
\[
 x_{\rm sum}=A_i-\frac{C_i c}{c^\top V_i c}\,c^\top d_i .
\]
Fine and raw-coarse drives agree for every observation exactly if \(C_iV_i^{-1}=(C_ic)c^\top/(c^\top V_ic)\). With force covariance \(B>0\), independent acceleration error of variance \(v>0\), this reduces to \(m\propto Bc\), by \(C=-v m^\top,\; V=B+vmm^\top\). For a single shared push with both fits nonzero, the two squared standardized tests at the physical candidate total nevertheless agree on **every** such dataset:
\[
 \frac{(c^\top u)^2}{c^\top H^{-1}c}=\frac{(c^\top d)^2}{c^\top Vc}
 =\frac{(c^\top F-(c^\top m)A)^2}{\operatorname{Var}(c^\top F-(c^\top m)A)}.
\]
Their tangent laws need not agree, since their centers and scales use different drives. On an independent new measurement of a lump the analogous statistic also has the exact \(\chi^2_1\) benchmark at truth under the declared model. This is a coverage statement for the scalar candidate test. It does not give exact finite coverage for the pushforward of a purported global joint mass law. The common-force case uses acceleration numerators and one shared force denominator with coefficients \(\nu_j=1/m_j\); the quantity that adds is \(\sum\nu_j\). Inverting the resulting positive physical coefficient gives the reduced mass and its tangent uses the derivative of reciprocation.

## A spatial push and its reciprocal tensor

For a common spatial drive \(x_i\in\mathbb R^k\), the stacked force defect has \(bk\) components. Use \(D=[I_{bk},-m\otimes I_k]\) in the same graph correction. The resulting fitted drive \(x_i=A_i-C_iV_i^{-1}d_i\) and defect \(d_i\) are again independent at the true candidate. In the tangent model the operator taking a mass displacement to its expected defect is \(T_i=I_b\otimes x_i\), viewing \(x_i\) as a column. Its noise-metric left reciprocal across pushes gives
\[
 H=\sum_i T_i^\top V_i^{-1}T_i,\qquad u=H^{-1}\sum_i T_i^\top V_i^{-1}d_i .
\]
Each row of this reciprocal is a covector on the stacked defects. Its contraction with the measured vector gives one component of the mass displacement. With at least one nonzero drive \(H\) is positive definite. Under independent changes of orthonormal spatial frame \(R_i\), the defects and tangent maps transform by \(I_b\otimes R_i\), and their covariances transform by congruence. These factors cancel in each contraction. The inferred mass displacement and covariance are frame invariant even for anisotropic, correlated instruments.

Write \(L=c^\top\otimes I_k\), \(B_i=LV_iL^\top\) and \(h=\sum_i x_i^\top B_i^{-1}x_i\). Retaining the fitted drives and only the summed vector defects gives variance \(1/h\), compared with \(c^\top H^{-1}c\) on the full stack. The scalar-direction equality criterion now says that there must be one body-weight vector \(\beta\) satisfying
\[
             V_i L^\top B_i^{-1} x_i=T_i\beta
               \quad\hbox{for every active push }i.
\]
Indeed the full and summed coefficient covectors, written as columns, are \(V_i^{-1}T_iH^{-1}c\) and \(L^\top B_i^{-1}x_i/h\). The difference of variances equals the sum of their squared differences in the \(V_i\) metrics by the same unbiasedness identity. Equality forces the displayed condition with \(\beta=hH^{-1}c\). Conversely the condition implies \(H\beta=hc\) on multiplying by \(T_i^\top V_i^{-1}\) and summing. For \(k=1\) this recovers the normalized covariance condition above.

**A single genuinely spatial push can lose information on summing forces.** Take \(b=k=2\), \(c=m=(1,1)\), \(x=(1,1)\), independent force error covariances \(\operatorname{diag}(1,9)\), \(\operatorname{diag}(9,1)\) across the two bodies, and independent common acceleration error covariance \(I_2/10\). The fine variance is 2 and the summed variance is \(26/5\). Each sensor's precise spatial direction helps determine its body's mass. Summing forces obscures this information in either spatial frame. These are conditional variances at the same fitted drive.

## What if both physical composition laws were known globally?

Even supplying the entire individual marginals, the sum law, and the reduced-mass law fails to specify a two-body positive joint law. Consider equal weights on \((1,2),(2,3),(3,1)\) and on the three transposed pairs. Every univariate marginal and every symmetric composite is identical. The laws of \(2m_1+m_2\) differ. This strengthens the marginals-only ambiguity in T51. The ambiguity can persist even if both joint laws are required to be exchangeable.

**An exchangeable counterexample.** For each \(r=1,2,3,4\), use the unordered shape \((r,1/r)\) with a uniformly random labelling and a positive common scale \(t\). The three observables (one labelled part, lump, and reduced mass) multiply \(t\) by random factors whose laws are the column
\[
 T_r=\left(\tfrac12(\delta_r+\delta_{1/r}),\ \delta_{r+1/r},\ \delta_{1/(r+1/r)}\right)^\top .
\]
In the commutative ring of finite signed measures with multiplicative convolution as multiplication, form the \(3\times4\) matrix \(T\) with these columns. Set \(q_j=(-1)^{j-1}\det T_{\widehat{\jmath}}\) for \(j=1,\ldots,4\). These are signed radial measures. Expanding a determinant with a repeated row gives \(Tq=0\). Use \(q_j\) as the scale measure for shape \(j\); the resulting signed joint measure has zero pushforward for the three observables. It is nonzero. Applying the homomorphism \(\eta\mapsto\int t^2\,d\eta\) to the three rows gives \(((s^2-2)/2,s^2,s^{-2})\) for \(s=r+1/r\). Each three-column determinant is nonzero by the Vandermonde determinant at distinct \(s^2\), and the four unordered shapes have disjoint scale rays. Dividing its positive and negative parts by their common total mass therefore gives two distinct exchangeable joint laws with the required equalities. Direct fraction arithmetic gives common total mass 12; their probabilities for \(2m_1+m_2=80/17\) differ by \(1/24\). All masses are positive. This compares possible joint laws at a fixed declaration; reachability from a particular production estimator is a separate issue.

Even smooth positive full-support exchangeable laws can differ this way. In each column use a smooth shape density on a small disjoint open interval of \(r>1\) close to its original \(r_j\). Define the three column measures by averaging the same observables over that density. The determinant and zero-pushforward proof persist. Second-moment minors stay nonzero by continuity, and the four shape bands remain disjoint. Convolve all signed radial measures with one compactly supported smooth probability density in log scale, so their nonzero second moments persist. Independent drawing of that shape and signed scale produces a nonzero smooth compactly supported signed joint density, symmetric under interchange. Its listed pushforwards are zero. Adding sufficiently small positive or negative multiples to a smooth strictly positive exchangeable base gives the asserted laws.

**What extra composition information would suffice?** For positive masses, knowing the laws of \(p m_1+q m_2\) for all positive integer \(p,q\), together with the two marginals, identifies the full joint law. Indeed put \(X=e^{-m_1},Y=e^{-m_2}\) in fixed mass units. The mixed moments \(E X^pY^q\) are known from these lump laws, and moments on the axes come from the marginals. They determine integration of every polynomial, hence every continuous function on \([0,1]^2\) by polynomial approximation. The inverse coordinate map then recovers the joint. This assumes the same unknown mass is reused exactly when a coefficient exceeds one. Physical copies with independently uncertain masses constitute a different experiment. No supply of all these lump laws is in the stated L5/L14 clauses.

The two regimes identify what a composition proof can establish. A specified Gaussian shared-drive experiment fixes the local joint normal and supplies a sharp information comparison. The two physical pushforwards alone leave different global positive couplings even with their marginals fixed. Neither assertion licenses identifying a new measurement of a lump with coarsening data already taken. [Finite checks](../tests/test_e12_joining_bodies.py)

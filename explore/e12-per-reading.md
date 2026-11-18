# Shared-fit per-reading weights on one instrument

This branch studies the scalar rule in [O14, 17 November, evening](../problems/O14.md#17-nov-evening-the-sign-problem-is-about-where-the-leftover-is-read-claude-scratch-one-number-readings). A preliminary fit uses all readings. Each reading then receives a weight from its residual at that fit, and the weights stay fixed while fitting the law. Two claims follow under independent Gaussian instrument errors.

1. The raw law concentrates at the true mass on every identifying array at the [T53 energy threshold](../definition/theorems/T53.md), including sparse signals and fixed sample size with increasing signal. This proves the scalar sign-consistency claim for this rule.
2. At zero signal the raw law also concentrates, this time at a random angle. Central intervals have limiting coverage zero for every fixed interior mass. This extends to the whole critical weak-signal range and to two explicit ways of using the unweighted Gram signal share in the report.

The first claim concerns consistency of the point. The second concerns the reported uncertainty. The three-reading shrinking-noise example elsewhere on this branch concerns fixed incompatible observed data, a third experiment. None of these claims selects this weighting rule for the definition.

## Experiment and fit

Use noise coordinates and write the common instrument's positive mass as \(m=s\tan\theta\), \(0<\theta<\pi/2\). Fix an interior true direction \(u_0=(\sin\theta_0,\cos\theta_0)\). In an experiment with \(N\ge1\) rows observe \(Y_i=r_i u_0+Z_i\), with arbitrary fixed signed amplitudes and independent \(Z_i\sim N(0,I_2)\). Write \(\Delta=\sum_i r_i^2\). Arrays may change both \(N\) and the amplitudes. Fix a proper prior with a positive continuous density near the truth for the consistency result. For the coverage result assume the prior gives positive mass to every open subinterval of the positive arc. A fixed half-Cauchy prior meets both assumptions.

Let \(B=\sum_i Y_iY_i^T\). The pilot \(p\) is a unit principal eigenvector of \(B\), regarded as an unoriented line. We also consider the fit that maximizes \(p^TBp\) on the **closed positive quarter-circle**, including the endpoints as limiting slopes. This is the continuous version of the positive search in O14's scalar interval study. Its sign study uses the whole projective line. A fixed tie rule suffices. Both fits agree whenever the whole-line maximizer lies in the open positive arc.

For a unit normal \(v_p\), fix \(\kappa>0\) and let
\[
 w(z)=\min(1,\kappa/z^2),\quad w(0)=1,\qquad
 S=\sum_i w(v_p^TY_i)Y_iY_i^T.
\]
The raw law is proportional to \(\exp(u_\theta^TSu_\theta/2)\,\pi(d\theta)\) on the positive arc. Write \(G\) for the eigenvalue gap of \(S\) and \(q\) for its whole-line principal angle. Up to a constant its log evidence is exactly \(-G\sin^2(\theta-q)/2\). These statements concern the continuous law. A fixed numerical mass grid introduces separate resolution and boundary effects.

The standard tools are orthogonal regression, principal-subspace perturbation, and concentration of a one-dimensional quadratic likelihood. The pilot bound and its energy-class converse are given in [T42](../definition/theorems/T42.md) and [T53](../definition/theorems/T53.md). The perturbation and minimax background is given by [Yu, Wang and Samworth](https://arxiv.org/abs/1405.0680) and [Cai, Li and Ma](https://jmlr.org/papers/v22/20-610.html). The data-dependent second fit requires the argument below.

## 1. A population gap created by weighting

In coordinates along and perpendicular to truth, \(Y_i=(r_i+\xi_i,\zeta_i)\) with independent standard Gaussian coordinates. If the pilot is evaluated at truth, reflection of \(\zeta_i\) kills the weighted off-diagonal moment. With \(Z\sim N(0,1)\), write \(A=Ew(Z)\), \(B_\perp=E[w(Z)Z^2]\), and \(g=A-B_\perp\). The expected weighted scatter of row \(i\) is diagonal with entries \(A(1+r_i^2)\) and \(B_\perp\).

Put \(a=\sqrt{\kappa}\), with \(\varphi,\Phi\) the standard normal density and CDF. Splitting each integral at \(a\), using \(\int_a^\infty \varphi(z)z^{-2}dz=\varphi(a)/a-\Phi(-a)\), gives
\[
 g=4[a\varphi(a)-a^2\Phi(-a)]>0.
\]
Positivity follows from \(a\int_a^\infty\varphi(z)dz<\int_a^\infty z\varphi(z)dz=\varphi(a)\). At \(\kappa=1\), \(g\) is about 0.33326. Weighting creates a preferred axis even in a population of noise, once an axis is supplied. The data-dependent pilot means that this population calculation alone cannot prove either claim.

## 2. Consistency under arbitrary signal allocation

**Claim.** If \(\Delta^2/(N+\Delta)\to\infty\), the weighted principal line converges in probability to the true interior line, \(G\to\infty\) in probability, and the raw law and its median concentrate at truth. This is the sharp *growth condition* over energy classes. No optimal rate for the weighted angle error or coverage at truth is asserted.

By T42's scatter bound the pilot error relative to truth is \(O_p(e)\), \(e=\sqrt{N+\Delta}/\Delta\to0\). For the constrained pilot the positive and whole-line fits consequently coincide with probability tending to one. In the true basis write \(D=S_{11}-S_{22}\), \(C=2S_{12}\). The angular deviation is \(\tfrac12\operatorname{atan2}(C,D)\) modulo \(\pi\), and \(G=\sqrt{D^2+C^2}\).

**Growing \(N\).** Set \(M=\min(N^{1/16},e^{-1/6})\to\infty\). This threshold is used only to bound the error. Split the rows by \(|r_i|\le M\) and \(|r_i|>M\).

- For the low rows, replacing the pilot by truth changes each entry of their total weighted scatter by \(o_p(N)\). Indeed \(w\) is globally Lipschitz with constant \(2/\sqrt{\kappa}\), and the difference is bounded by that constant times the pilot angular error times \(\sum_{\rm low}\|Y_i\|^3=O_p(N(1+M)^3)\). At truth the rows are independent. Centering each entry costs \(O_p(\sqrt N(1+M)^2)=o_p(N)\) by the fourth-moment bound. Thus \(D_{\rm low}=A\sum_{\rm low}r_i^2+g n_{\rm low}+o_p(N)\) and \(C_{\rm low}=o_p(N)\).
- Call a high row good when \(\|Z_i\|\le |r_i|/\sqrt M\). Its projective angle with truth is bounded by a deterministic quantity tending to zero, since its perpendicular component is at most \(|r_i|/\sqrt M\) and its parallel component in absolute value at least \(|r_i|(1-M^{-1/2})\). For any pilot \(w(v_p^TY_i)\|Y_i\|^2\ge\min(\|Y_i\|^2,\kappa)\). Each such row therefore has \(D_i\ge c>0\) eventually and \(|C_i|\le o(1)D_i\), by the double-angle identities, uniformly over the good rows.
- A bad high row has \(\|Z_i\|>\sqrt M\) and \(\|Y_i\|\le(1+\sqrt M)\|Z_i\|\). Per row, its expected squared norm times its bad indicator is bounded by \((1+\sqrt M)^2 E[\|Z\|^2{\bf1}_{\{\|Z\|>\sqrt M\}}]=o(1)\) for a standard two-dimensional Gaussian \(Z\). Their total scatter and their count are both \(o_p(N)\), by Markov's inequality and the Gaussian tail bound.

Together \(D\ge c'N\) with probability tending to one, for a fixed \(c'>0\). Also \(D_{\rm good}\le D+o_p(N)\), because the low rows have nonnegative leading gap. Thus \(|C|/D\to0\) in probability even when a few good high rows dominate the entire gap.

**Bounded \(N\).** Here \(\Delta\to\infty\), the pilot error is \(O_p(1/\sqrt\Delta)\), and \(\max_i|r_i|\le\sqrt\Delta\). All residuals at the pilot are \(O_p(1)\), so the reciprocal of the smallest weight is \(O_p(1)\). The leading parallel sum \(\sum_i w_i r_i^2\) is at least \((\min_iw_i)\Delta\). The error in \(D\) and the entire \(C\) are \(O_p(\sqrt\Delta+1)\) by Cauchy–Schwarz and the bounded number of Gaussian errors. This proves the same line and gap limits. Every sequence of sample sizes admits subsequences in one of these two cases.

For any fixed angular neighborhood of the interior truth, a smaller interval around truth with positive prior mass eventually beats its complement outside the larger neighborhood by a fixed multiple of \(G\) in log evidence, with probability tending to one. This proves concentration of the raw law and consistency of every fixed interior quantile.

## 3. The raw law and two explicit reports under weak signal

Let \(N\to\infty\) and \(\Delta/\sqrt N=O(1)\). At a fixed interior true mass and any fixed nominal central level in \((0,1)\), the following intervals have sampling coverage tending to **zero** for either pilot.

- A central interval of the raw law.
- A central interval after multiplying the raw log evidence by the **unweighted** Gram share \(h=(\lambda_1(B)-\lambda_2(B))/\lambda_1(B)\).
- The central interval formed by applying the normal-score reporting warp \(F\mapsto\Phi(\sqrt h\,\Phi^{-1}(F))\) to the raw CDF, using that same share.

For the zero scatter or an exact tie these report comparisons may use any fallback. Both events have probability zero in the growing Gaussian experiments. These are the specified raw or share-only reports on the frozen per-reading scheme. Other rules for reporting uncertainty or combining results require their own analysis.

**Uniform comparison before choosing the pilot.** For pure-noise \(Z_i\), the empirical matrix \(N^{-1}\sum_i w(v^TZ_i)Z_iZ_i^T\) differs uniformly over unit normals \(v\) by \(o_p(1)\) from \(A(I-vv^T)+B_\perp vv^T\). A finite net of normals gives entrywise laws of large numbers. Interpolation costs at most the net mesh times \(\operatorname{Lip}(w)N^{-1}\sum_i\|Z_i\|^3\) for the empirical matrix, with the same finite expectation bound for its population version. Send the mesh to zero. Adding deterministic signal \(r_i u_0\) changes the summand by at most
\[
 \operatorname{Lip}(w)|r_i|\|Z_i\|^2+2|r_i|\|Z_i\|+r_i^2
\]
in operator norm, for every \(v\). The average of these bounds has expectation \(O(\sqrt{\Delta/N}+\Delta/N)=o(1)\). We can now substitute the fitted \(v_p\), which depends on **all** rows. Therefore
\[
 S/N=A pp^T+B_\perp(I-pp^T)+o_p(1),\qquad G/N\to g,\quad
 \operatorname{dist}_{\rm proj}(q,p)\to0 .
\]

**The pilot remains random.** In the true basis the traceless unweighted scatter vector \((B_{11}-B_{22},2B_{12})/\sqrt N\) has deterministic mean \((\Delta/\sqrt N,0)\). Its pure-noise part converges to \(N(0,4I_2)\) by the bivariate CLT. Its centered terms linear in \(r_i\) have variance \(4\Delta/N=o(1)\) in each component. Along a subsequence where \(\Delta/\sqrt N\to d\), the limiting vector is nonsingular Gaussian with mean \((d,0)\). Half its argument gives the limiting whole-line pilot, which has a continuous distribution on the projective circle. At pure noise the whole-line pilot is exactly uniform for each \(N\), by rotational symmetry. Projection onto the closed positive arc can add atoms at its endpoints. All interior masses retain probability zero as singleton maximizers. For the whole-line pilot the negative-arc angle equidistant from both endpoints has limit probability zero.

Since \(\lambda_1(B)/N\to1\), \(\sqrt N\,h\) converges along the same subsequence to the norm of the limiting Gaussian vector, strictly positive almost surely. Thus \(Nh\to\infty\) in probability while \(G/N\to g>0\).

**Coverage.** Away from a tie between boundary maximizers, \(u_\theta^TSu_\theta/G\) has, up to an additive constant, the cosine-squared profile with its unique constrained maximizer approaching the positive-arc maximizer of the pilot. Except when this maximizer equals the fixed interior truth, one side of the truth has an evidence maximum separated from the global one by a strictly positive fraction of \(G\). Full prior support supplies positive mass on an interior interval arbitrarily close to the global maximizer, including a maximizer at a boundary. The smaller one-sided raw probability at the truth is then exponentially small in \(G\). Formally one may exclude arbitrarily small neighborhoods of the pilot angles giving a tie or that interior truth. On the remaining compact sets the separation is uniform, and finitely many such prior intervals give a common positive lower mass. The excluded whole-line angles have arbitrarily small limiting probability. For the constrained pilot only interior truth requires exclusion, since endpoint pilots themselves have unique constrained maximizers.

The raw central ranks are fixed. For the powered report the same exponential bound holds with \(hG\to\infty\). For the normal-score warp, if the lower nominal rank is \(b\in(0,1/2)\), its raw rank is \(\Phi(\Phi^{-1}(b)/\sqrt h)\), whose negative logarithm is \(O(1+1/h)=o_p(G)\) by the Gaussian tail bounds. Thus even its remote ranks exclude that exponentially small side, giving the claimed zero coverage. Every subsequence with bounded \(\Delta/\sqrt N\) has a further convergent \(d\), which proves the result throughout the stated window.

On the identifying arrays of §2 these two specific reports *do* shrink around the truth. T42 gives \(\lambda_1(B)-\lambda_2(B)\sim_p\Delta\), while \(\lambda_1(B)=O_p(N+\Delta)\). Hence \(hG\to\infty\) using \(G\ge c'N\) for growing \(N\), and \(h\to1\) for bounded \(N\). The same gap bounds apply outside fixed neighborhoods of the truth. This says nothing about nominal coverage on identifying arrays, which requires control of errors on the shrinking interval's own scale.

An exactly observed all-zero series contributes zero weighted scatter and returns its starting law. Gaussian readings with zero *true* signal are random measured pairs. It is their random fitted direction, made sharply preferred by the second fit, that causes the false precision above.

## 4. The same conclusions for other monotone residual weights

**Corollary.** Both §2 and §3 hold for every fixed, bounded, strictly positive, Lipschitz \(w\) that is even, nonincreasing in \(|z|\), and nonconstant under a standard Gaussian. The function and its tuning constants are held fixed as \(N\) changes. Under these hypotheses, changing the taper alone cannot fix the weak-signal coverage of the three reports.

To see the strict gap, take independent \(Z,Z'\sim N(0,1)\). Then
\[
 2g=E[(w(Z)-w(Z'))(Z'^2-Z^2)]>0 .
\]
The product is nonnegative by monotonicity and strictly positive with positive probability by nonconstancy. The uniform law of large numbers and signal-perturbation argument in §3 use boundedness and the Lipschitz property. One can divide \(w\) by its finite supremum in the scatter bounds and restore this constant at the end. All three reports retain the same limiting coverage. Concentration here is on the closed angular arc. An endpoint selected in that limit corresponds to a mass approaching zero or infinity.

For the identifying proof the low rows and the bound on the bad high rows go through with these same two constants. For the good high rows, replace the special cap-weight estimate by a consequence of the *first* fit. The true direction is a feasible competitor under either pilot domain. Hence
\[
 \sum_i(v_p^TY_i)^2\le\sum_i(v_0^TY_i)^2=\sum_i\zeta_i^2 .
\]
For \(N\to\infty\), at most \((1+o_p(1))N/4\) rows can have \(|v_p^TY_i|>2\). If at least half of the rows have low amplitude then \(D_{\rm low}\ge gN/2+o_p(N)\) and all good high contributions are nonnegative. Otherwise at least \((1/4-o_p(1))N\) good high rows have \(|v_p^TY_i|\le2\). Each of these has weight at least \(w(2)>0\) and length eventually above one, so their combined gap is at least a fixed positive multiple of \(N\) with probability tending to one. This restores the needed lower bound in both cases. The angular bound \(|C_i|=o(1)D_i\) uses only the direction of good rows. For bounded \(N\), tightness of the largest fitted residual and positivity of \(w\) on compact sets give \(1/\min_iw_i=O_p(1)\) as before.

Constant weights sit exactly at \(g=0\); the argument producing an order-\(N\) null gap no longer applies. A weight tending to a constant with sample size, or a separate rule for retaining and reporting uncertainty, is outside this fixed-weight corollary.

**Finite checks.** [Tests for this exploration](../tests/test_e12_per_reading.py) exercise the population moments, the two pilot domains, the connection with O14's existing study, identifying allocations, and the weighted and unweighted null-gap scales. The asymptotic claims rest on the proofs.

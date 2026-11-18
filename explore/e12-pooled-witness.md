# A family of pooled scalar laws on one instrument

## Claim and standard mathematics

Take independent Gaussian scalar force–acceleration readings from one instrument with positive channel standard deviations, independent channels, fixed instrument ratio \(s=\sigma_F/\sigma_a\), and fixed interior true mass. The signals can have any deterministic signs and sizes. This section constructs a family with exact measured-zero neutrality, bounded informative evidence on a specified class of incompatible shrinking-noise data, full-rate consistency, and asymptotically nominal central intervals on every identifying array. Two members also widen the interval on the 5 and 6 example. They give different finite answers. The other premises, the pure-noise reporting promise, joint composition and selection of a unique member remain open.

The underlying fit is Gaussian orthogonal regression. The interval argument uses the local normal form of a generalized posterior and independently establishes its match to the fitted angle's sampling law. The project already derived the scalar scatter and score limits and the reference spectral precision in [T42](../definition/theorems/T42.md) and [T43](../definition/theorems/T43.md). We use those results and prove that a family of continuous pooled discounts preserves their calibration. For the general distinction between concentration and calibrated generalized-posterior uncertainty, see [Miller, 2021](https://jmlr.org/papers/v22/20-469.html).

## A pooled law and its report

Write \(Y_i=(F_i/\sigma_F,a_i/\sigma_a)\), \(\theta=\arctan(m/s)\), \(u_\theta=(\sin\theta,\cos\theta)\) and \(v_\theta=(\cos\theta,-\sin\theta)\). For the scatter \(B=\sum Y_iY_i^\top\), let \(L\ge l\ge0\) be its eigenvalues, \(H=L-l\), and \(\hat\theta\) any principal line on the **whole projective circle**. The least residual used by the discount is instead
\[
 Q=\min_{0\le\theta\le\pi/2}v_\theta^\top B v_\theta .
\]
The closure includes both axes. With positive off-diagonal scatter \(Q=l\). With a nonpositive off-diagonal it is the smaller of the two diagonal entries. In particular \(0\le Q\le L\).

Choose constants \(b,\alpha>0\) to name a member of the family. They are choices in this construction. No measurement premise selects them. Define
\[
 u=e^{-bt/2},\qquad h_b(t)=(1-u)\,[1+(2+1/b)u],\qquad
 C=\sum_i h_b(\|Y_i\|^2).
\]
At \(Q>0\) set
\[
 d=\min\{1,\ (C/Q)(L/Q)^\alpha\},\qquad
 K=d\, H^2/L .
\]
Use \(d=1\) when \(Q=0\) and \(K=0\) when \(L=0\). For any fixed proper prior \(\pi\) on positive angles the law is
\[
 P(d\theta\mid Y)\ \propto\ \pi(d\theta)
           \exp[-K\sin^2(\theta-\hat\theta)/2].
\]
Report its median (use the lower median on nonunique-median laws) and its central equal-tailed interval, with \(m=s\tan\theta\). The coverage result below requires a positive continuous prior density in a neighborhood of the true angle. Reporting this particular law's central interval is part of the witness. It leaves open a different law with a separately calibrated report. The weak-signal factor \(H/L\) is T43's spectral control. The additional factor \(d\) limits confidence on incompatible shrinking-noise data.

**Exact finite properties.** The kernel is positive and bounded, so every finite input gives a proper law. On measured all-zero data it returns the prior. An observed zero row changes neither \(C\) nor \(B\). Even an appended row tending to zero changes the law continuously in total variation: away from \(Q=0\) this follows by continuity, at a nonzero scatter with \(Q=0\) the discount equals one in a neighborhood, and at \(H=0\) the kernel tends uniformly to one. Eigenvector ambiguity at a tie is harmless. The summaries add across blocks: retaining \(B\) and \(C\) allows one to recompute exactly the result of pooling the raw rows for a fixed member. This is same-object raw-data pooling. It provides no physical composition law for different masses or independent routes.

Permutation of the rows, simultaneous sign reversal of a measurement pair, and consistent changes of physical units preserve the angle law. Channel interchange reflects the angle if the prior is transformed with it. For a continuous strictly positive prior the median and central interval transform accordingly. No assumption of a population of pushes enters the computation. Whether L19's evidence requirements permit this use of the row lengths remains open.

## Fixed readings as both noise scales shrink

Multiply each fixed whitened row by \(\sqrt{t}\) with \(t\to\infty\), holding \(s\) fixed. Subscript 0 refers to the unscaled scatter. On any fixed shape with \(Q_0>0\), \(C\to n_+\), the number of nonzero rows, since \(h_b(0)=0\) and \(h_b(x)\to1\) as \(x\to\infty\). The expression inside the minimum eventually falls below one. Therefore
\[
 K\longrightarrow {n_+\over Q_0}
        \left({L_0\over Q_0}\right)^\alpha {H_0^2\over L_0}.
\]
For \(H_0>0\) this is finite and strictly positive. The limiting kernel has bounded and nonzero evidence contrast across the positive mass domain. A full-support prior thus retains spread and receives directional information. Convergence of the kernel is uniform and the laws converge in total variation. For an interior rank-one compatible shape \(Q_0=0\) and \(K=tL_0\), so a prior with positive continuous density at that angle concentrates on the compatible mass.

**Two limits of the base construction matter for literal T30.** On a positive scalar Gram \(B_0=cI\) every candidate has the same positive leftover. Every global scalar power of Gaussian quadratic evidence is constant there, including this construction. For a rank-one axis Gram, every interior mass has positive leftover, yet \(Q_0=0\) on the closed arc. With a full-support prior the law concentrates toward the corresponding endpoint as \(t\) grows. Thus bounded *and* informative evidence for all fixed readings with positive leftover at every interior mass is stronger than the proved shape statement. In contrast, a negatively sloped non-axis rank-one Gram does have \(Q_0>0\) and is covered by the theorem.

## Full-rate point consistency and interval coverage

Consider any deterministic triangular array \(Y_i=r_i u_{\theta_0}+Z_i\) with independent \(Z_i\sim N(0,I_2)\). Write
\[
 \Delta=\sum_i r_i^2,\qquad V=N+\Delta,\qquad I=\Delta^2/V.
\]
The only growth hypothesis is \(I\to\infty\). This includes arbitrarily sparse signal and fixed-\(N\) amplification. For any fixed \(0<p<1\), the central probability-\(p\) interval of the witness has sampling coverage tending to \(p\); the point converges in probability to the true mass.

**1. The discount tends to one.** For \(Y\sim N(\mu,I_2)\), completing the square gives
\[
 E\exp(-b\|Y\|^2/2)=(1+b)^{-1}
                \exp[-b\|\mu\|^2/(2(1+b))].
\]
Expanding \(h_b\) in \(u,u^2\) therefore gives
\[
 E h_b(\|Y\|^2)=1+b^{-1}
 \{\exp[-b\|\mu\|^2/(2(1+b))]-
   \exp[-b\|\mu\|^2/(1+2b)]\}\ \ge 1.
\]
Equality holds at zero signal; the expected count per row also tends to one as its signal grows. As \(h_b\) is bounded, independence gives \(C\ge N-O_p(\sqrt N)\) uniformly in the signals. At the true interior line the total transverse squared residual is \(\chi_N^2\), so \(Q\le N+O_p(\sqrt N)\). Since \(L/Q\ge1\), \(d\to1\) in probability whenever \(N\to\infty\).

On a bounded-\(N\) subsequence with \(I\to\infty\), one of the signal norms diverges (take the largest in each row of the array). Hence \(L\to\infty\) and \(C\ge1-o_p(1)\), whereas \(Q=O_p(1)\). Since \(\alpha>0\), the cap yields \(d\to1\) here too. Any sequence of \(N\)'s has a further bounded or diverging subsequence, which proves the claimed convergence for all identifying arrays. The multiplier \((L/Q)^\alpha\) is one way of meeting the bounded-\(N\) case even when the realized transverse residual is large.

**2. The centre and precision have matching limits.** \(E B=NI_2+\Delta u_0u_0^\top\) and \(E\|B-E B\|_F^2=6V\). Thus \(H/\Delta\to1,\ L/V\to1\) in probability. In the true parallel/perpendicular basis let \(X_i=r_i+\xi_i,\ W_i=\zeta_i\), with the two standard normal noises independent. The off-diagonal scatter \(S=\sum_i X_iW_i\), conditional on \(X\), is Gaussian with variance \(\sum X_i^2\). This variance has mean \(V\) and variance \(2N+4\Delta\), so \(S/\sqrt V\Rightarrow N(0,1)\). The difference of the parallel and perpendicular diagonal entries is \(\Delta+O_p(\sqrt V)\). The two-dimensional eigenvector formula \(2(\hat\theta-\theta_0)=\operatorname{atan2}(2S,B_{\parallel\parallel}-B_{\perp\perp})\), with the representative close to truth, now gives
\[
 \sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1),
 \qquad K/I=d(H/\Delta)^2(V/L)\longrightarrow1
\]
in probability for the second limit. No leverage condition enters.

**3. The quantiles have the required sampling coverage.** With probability tending to one, the fitted angle lies in a closed neighborhood inside the positive arc and \(K\to\infty\) in probability. At \(\theta=\hat\theta+x/\sqrt K\) the kernel tends to \(e^{-x^2/2}\) uniformly on bounded \(x\), and the prior density tends to its strictly positive value at \(\theta_0\). On a small fixed neighborhood the density is bounded and \(\sin^2(\theta-\hat\theta)\ge c(\theta-\hat\theta)^2\) for some \(c>0\). Outside that neighborhood the contribution of the proper prior is exponentially small compared with the normalizing scale \(1/\sqrt K\). These bounds give the standard normal local limit and, for each fixed rank \(\gamma\in(0,1)\), angle quantile
\[
 \hat\theta+\{\Phi^{-1}(\gamma)+o_p(1)\}/\sqrt K .
\]
Part 2 now gives nominal repeated-sampling coverage at the two central ranks and consistency at rank one half. The monotone physical mass map preserves coverage. No finite-sample, uniform near-null, or independent-route coverage guarantee follows.

## Distinct finite laws and the ruling's two scenes

For \(F=(5.5,5.6)\) or \(F=(5,6)\), equal accelerations of one, \(\sigma_F=0.1,\ \sigma_a=0.002\), use \(b=1\). In this table the prior is uniform in noise angle, hence half-Cauchy with \(c=s=50\). This explicitly differs from the centre \(c=5\) in O14's scratch calculation. The discounts depend only on the readings.

| \(\alpha\) | close \(d\) | far \(d\) | close 95% mass half-width | far 95% mass half-width |
|---|---:|---:|---:|---:|
| 1 | 1 | 1 | 0.13944 | 0.13944 |
| 1/8 | 1 | 0.12841 | 0.13944 | 0.38914 |
| 1/16 | 1 | 0.07210 | 0.13944 | 0.51934 |

The exponent one fails the concrete discount direction. The smaller exponents both widen the far pair, yet they disagree about the amount while meeting every clause proved above. More generally \(L/Q\to\infty\) on these same fixed force readings if only \(\sigma_a\to0\) while \(\sigma_F\) stays fixed. Here \(Q\) tends to the nonzero force-only residual sum in noise units and \(d\to1\) for every positive \(\alpha\). Thus passing the finite scenes does not supply a rule depending only on gap/\(\sigma_F\) uniformly over acceleration precision. This one-channel limit is distinct from common shrinkage of both standard deviations in the deterministic theorem.

The count also matters at fixed \(\alpha\). With rows \((3,3)/\sqrt2,(2.5,-2.5)/\sqrt2\), \(L=9,\ Q=6.25,\ H=2.75\), \(\alpha=1/8\), the precisions are \(0.29605\) for \(b=1\) and \(0.28186\) for \(b=2\). Exact null mean, neutrality, continuous count, fixed-shape limiting tilt and identifying-array calibration all survive the change.

## A coverage obstruction for persistent gap-only widening

A stronger interpretation of the two-reading example conflicts with exact regular coverage for an explicitly delimited class of reports. Take two independent rows with the same deterministic signal \(A u_0\), \(A\to\infty\), and Gaussian errors. Set \(W_i=v_0^\top Z_i\), \(Z=(W_1+W_2)/\sqrt2\) and \(G=(W_1-W_2)/\sqrt2\). The orthogonal rotation gives **independent** standard normals \(Z,G\). Part 2 of the proof, or direct expansion of the two-dimensional eigenvector equation, gives \(\sqrt2 A(\hat\theta-\theta_0)=Z+o_p(1)\). Also \(v_{\hat\theta}^\top Y_i=W_i-A(\hat\theta-\theta_0)+o_p(1)\), so the fitted transverse residuals tend in probability to \(G/\sqrt2,-G/\sqrt2\) in this coupling. Hence the between-reading gap remains random at the same noise-unit scale even at arbitrarily high signal.

Suppose now that the reported central level \(p\) interval has endpoints
\[
 \hat\theta\ \pm\ {z\,g(G)\over\sqrt2 A}+o_p(1/A),
 \qquad z=\Phi^{-1}((1+p)/2).
\]
Assume a finite measurable multiplier with \(g\ge1\), and \(g>1\) on a set of positive standard-normal probability. This includes a persistent widening depending on the magnitude of the fitted gap when the multiplier is continuous almost everywhere and the displayed endpoint expansion holds. Independence and continuity of \(Z\) give limiting coverage
\[
 P(|Z|\le z\,g(G))=E[\,2\Phi(z\,g(G))-1\,]>p .
\]
For example, widening by 50% only when \(|G|>2\) adds \(2\Phi(-2)[2\Phi(1.5z)-1-p]\) to the nominal probability. Even the exact Gaussian location experiment with two readings has this mean–difference identity. The slope experiment approaches it on this identifiable slice.

Thus for this efficient centre and symmetric, never-narrower local intervals, persistent gap-only inflation produces **overcoverage**. This applies to the stronger target of asymptotic equality to the nominal level along every identifying array. A conservative coverage requirement, a changed centre, a different endpoint shape, or widening that vanishes in this experiment needs a different argument. The positive witness above uses the last possibility. This does not settle property N's intended promise.


## Why an exactly unbiased continuous per-row count is unavailable

Suppose bounded continuous \(h:\mathbb R^2\to\mathbb R\) obeys \(h(0)=0\), and its expected contribution is required to be exactly one for every amplitude and interior mass under the declared Gaussian instrument. The possible mean vectors contain the open positive quadrant. Therefore
\[
 \int(h(y)-1)e^{-\|y\|^2/2}e^{\mu^\top y}\,dy=0
\]
on an open set of real \(\mu\). Boundedness and the Gaussian tail give an entire transform in both coordinates. Apply the one-variable identity theorem successively on an open rectangle in the quadrant to obtain zero for all complex \(\mu\). Taking \(\mu\) purely imaginary gives a zero Fourier transform of the integrable signed density. Injectivity gives \(h=1\) almost everywhere, and continuity gives \(h(0)=1\), a contradiction. This is Gaussian location-family completeness for a count that uses one observed row and the declarations, independently of an unknown mass or a series fit. Counting nonzero rows obeys exact measured-zero neutrality and has exact mean one per Gaussian observation, but jumps at a measured zero. Counts depending jointly on the full sample are outside this argument.

**Finite checks.** [Witness implementation](e12_scalar_witness.py) and [tests](../tests/test_e12_pooled_witness.py) cover Gaussian moments by quadrature, zero padding, the constraints on the positive arc, the incompatible ray and its boundary exceptions, the finite scenes, channel interchange, and fixed-count amplification. The limiting claims above use analytic arguments.

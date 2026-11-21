# Which uncertainty laws preserve addition exactly?

The two labelled positive masses are \(X,Y\). On an exact shared-acceleration reading, contraction with the reciprocal acceleration gives each mass. The common-acceleration sum then reads \(X+Y\). The common-force composition reads \(H=(X^{-1}+Y^{-1})^{-1}\). These are degree-one functions of the same pair. This suggests asking when the entire law of a composite can be obtained by applying the operation to two scale parameters. The answer holds for any positive probability law, including a half-Cauchy.

## A converse for a common shape

Fix a law \(S>0\) and scales \(a,b>0\). Allow arbitrary dependence between \(X\overset d=aS\) and \(Y\overset d=bS\). Each of the following three conditions implies the other two.

1. \(X+Y\overset d=(a+b)S\).
2. \(X/a=Y/b\) almost surely.
3. \(H\overset d=(a^{-1}+b^{-1})^{-1}S\).

**Proof.** Write \(U=X/a,\ V=Y/b,\ w=a/(a+b)\). Both \(U,V\) have law \(S\). For any fixed \(t>0\), strict convexity of \(f(z)=e^{-tz}\) gives
\[
 E f(wU+(1-w)V)\ \le\ wE f(U)+(1-w)E f(V)=E f(S).
\]
Equality holds precisely on \(U=V\) almost surely. Condition 1 implies equality. Conversely condition 2 gives the two asserted identities pathwise with one draw \(U\). To go from 3 back to 2, apply the same proof to \(1/X,1/Y\), which have scales \(1/a,1/b\) and common shape \(1/S\). Taking reciprocals in 3 gives condition 1 for that pair. Every expectation used in the argument is bounded. In particular **one Laplace value at one positive \(t\)**, together with the two specified marginals, suffices for either converse.

There is a quantitative form. Put \(Z=(X+Y)/(a+b)\) and \(c_w=\min(w,1-w)\). Then
\[
 E e^{-tS}-E e^{-tZ}
 \ \ge\ c_w\,E\big[(e^{-tU/2}-e^{-tV/2})^2\big].
\]
Indeed for \(w\le1/2\), set \(A=e^{-tU},B=e^{-tV}\). Weighted arithmetic–geometric mean gives \(A^wB^{1-w}\le 2w\sqrt{AB}+(1-2w)B\). Subtract from \(wA+(1-w)B\). Exchange the variables for \(w\ge1/2\). At equal scales the bound is an identity. A single transform can be insensitive to differences far out in the tail.

For the project's positive half-Cauchy family, the declared centre is the scale (and median). Common-shape additivity holds exactly when the two scale parameters multiply the same unit half-Cauchy draw. Independently drawing the two scales gives this property only for a point mass. Along a common ray, every positively homogeneous reported point composes for both operations. This includes the median and each quantile at the same rank.

## All quantiles, with different marginal shapes

Let \(Q_X(p)=\inf\{x:F_X(x)\ge p\}\) for \(0<p<1\), and similarly for \(Y\). With any two fixed positive marginals the following are equivalent: (i) every quantile adds, (ii) the pair has the law of \((Q_X(V),Q_Y(V))\) for one uniform rank \(V\), and (iii) every quantile composes under reduced mass. Condition (ii) is called comonotonicity. A curved increasing relation can satisfy it; proportional marginal shapes give precisely the rays above.

For the forward proof from (ii), both composite functions of the rank are nondecreasing and left-continuous, hence give the composite quantiles, including at atoms. Conversely let \(A=e^{-tX},B=e^{-tY}\) for any fixed \(t>0\). With the marginals fixed,
\[
 E[AB]=\int_0^1\!\int_0^1\Pr(A>s,B>u)\,ds\,du
 \ \le\int_0^1\!\int_0^1\min\{\Pr(A>s),\Pr(B>u)\}\,ds\,du.
\]
The right side is attained by shared ranks, because the two upper-tail events are nested for every pair of thresholds. Condition (i) gives exactly the same sum law as in that coupling, so the two sides are equal. Their nonnegative integrand difference vanishes almost everywhere, and right continuity of upper-tail probabilities in each threshold gives equality everywhere by approaching thresholds from above. This is the joint upper-tail law of a shared-rank pair. Reversing both coordinates back yields (ii). For (iii), equality of the full reduced-mass law with the shared-rank benchmark implies equality in law of \(1/X+1/Y\) with the corresponding benchmark for the reciprocals. Both reciprocals reverse order together. Apply the same argument to them. This is the positive-mass version of the established characterization by the distribution of the sum; the bounded transform accommodates the project's infinite-mean marginals. [Sum-distribution characterization](https://www.sciencedirect.com/science/article/abs/pii/S0167668710000739)

A single point is weaker. Put equal probability on \((1,1),(1,4),(4,1)\). The marginal medians are 1 and 1 and the total median is 5. Conversely a two-point anti-ordered pair with unequal weights can have all its medians compose simply because one atom carries a majority. An individual match of medians does not diagnose a common scale. [Point readout](../problems/O2.md)

## What the superposition is

Write \(T=X+Y,\ r=X/T\). Every positive joint law has a unique allocation law \(\eta(dr)\) and, for \(\eta\)-almost every \(r\), a conditional law \(\mu_r(dT)\). For any bounded test function \(f\),
\[
 \int f(x,y)\,\mu(dx,dy)=\int\!\int f(rT,(1-r)T)\,\mu_r(dT)\,\eta(dr).
\]
This is disintegration on ordinary Borel spaces. The radial law can vary with the ray. Both mass operations act on every component by scaling: \(X=rT,\ Y=(1-r)T,\ H=r(1-r)T\). The mass law reported by each of these observables is a linear mixture of the corresponding pushforwards. The components have a precise invariance under common dilation; a fixed probability distribution on a ray need not itself be invariant under dilation. A nonlinear point such as the median is taken *after* mixing and loses the rankwise addition guarantee.

An exact measurement of \(r\) selects a ray and leaves its conditional radial law. An exact measurement of \(T\) selects the line of constant total, generally meeting many rays. For a concrete pair of preparations take equal weights on \((T,r)=(4,1/4),(8,3/4)\), or exchange their allocations. Total 4 in the first preparation gives \(X=1\), while in the second it gives \(X=3\). Their unconditioned distributions of \(T\) and \(r\) separately are identical. Which component a measurement selects depends on the instrument and on what was known jointly beforehand. [Measurement theorem](e12-mass-measurement.md) · [Which composite laws identify the preparation](e12-mass-projections.md) · [Exact checks](../tests/test_e12_mass_rays.py)

**What changes after a reading?** Reweighting by any strictly positive, finite likelihood with finite normalizer makes the posterior and prior mutually absolutely continuous. Being concentrated on a ray is therefore unchanged. So is comonotonicity: it is equivalently the absence, almost surely, of a discordant pair \((X-X')(Y-Y')<0\) for two independent copies, and the corresponding product measures are equivalent. For this equivalence, nested marginal threshold events give shared ranks. If instead two threshold events are not nested up to null sets, each opposite difference event has positive probability, producing two discordant copies. Thus a noisy positive-likelihood update can move arbitrarily close to one ray, while exact selection from a prior crossing multiple rays needs a sharp observation or a limit.

A sharp total \(T=t\) gives a further point-level identity even when it leaves several allocations unresolved. On this conditional law \(Y=t-X\), so the set of medians of \(Y\) is the reflection of the set of medians of \(X\) about \(t/2\). **If the two marginal medians are unique, their sum is exactly the measured total.** For a total known only to lie in \([t-\epsilon,t+\epsilon]\) almost surely, monotonicity of quantiles applied to \(t-\epsilon-X\le Y\le t+\epsilon-X\) yields \(|\operatorname{med}X+\operatorname{med}Y-t|\le\epsilon\) when the medians are unique. An exact reduced mass gives the dual identity by applying reflection to \(1/X+1/Y=1/H\) and inverting the unique medians. Thus a genuine sharp reading can restore median composition without selecting a unique ray. Requiring every rank to compose at a fixed total forces both masses to be constant, since otherwise the constraint \(Y=t-X\) creates discordance.

# Which Gaussian evidence factor loses least on matching preparations?

Fix a mass candidate and whiten the declared force and acceleration channels. Its possible exact readings form a linear graph \(G\subset\mathbb R^p\) of dimension \(k\). Write \(P\) for the orthogonal projection onto \(G\), \(R=I-P\), and \(K(y)=\|Ry\|^2\). For the two-force allocation test \(Ry\) can be expressed through the standardized force contrast. For a one-body mass it is the covector residual after fitting the candidate drive. The proof uses this same decomposition in any dimension. [Graph construction](../definition/theorems/T36.md) · [Force contrast](e12-testing-a-mass-ray.md)

## A common denominator for specified reference readings

A specified Gaussian reference preparation with exact signal \(b\) predicts density \(p_b(y)=(2\pi)^{-p/2}\exp[-\|y-b\|^2/2]\). Its mean is stipulated before observing \(y\). Seek a positive finite candidate factor \(D_G(y)\) independent of which reference is stipulated. Ask that the comparison \(p_b(Y)/D_G(Y)\) have expectation at most one at **each** candidate drive \(a\in G\). Such a nonnegative ratio is called an e-value against \(G\). Independent comparisons multiply, with the same expected-value bound conditional on earlier readings if the next reference is chosen using only past information. No distribution over actual drives enters. The admissible references must include every \(b\in G\). The result below works even if every \(b\in\mathbb R^p\) is allowed. [Grünwald, de Heide and Koolen on e-values and logarithmic growth](https://doi.org/10.1093/jrsssb/qkae011)

There is a natural place to compare losses. When the reference and the actual candidate preparation match, its density really is \(p_a\). How large can the expected **log** comparison to the common denominator be in that case? This is an additional performance criterion. It compares rules on a correctly specified Gaussian preparation, and the log adds under independent repetitions.

**Theorem.** Among all measurable common denominators satisfying these expectation bounds, the unique maximizer of \(E_a\log[p_a(Y)/D_G(Y)]\), for *every* \(a\in G\) simultaneously, is
\[
 D_*(y)=2^{-k/2}(2\pi)^{-p/2}\exp[-K(y)/2].
\]
Maximizing at just one matching drive, even zero, forces this same choice up to a Lebesgue-null set. The maximum is \(k(\log 2-1)/2\). The resulting mass-dependent factor is the familiar Gaussian profile. The numerical coefficient here calibrates the comparisons with actual reference densities.

**Validity.** Resolve \(Y=(U,V)\) into \(G,G^\perp\) in orthonormal coordinates. At \(a\), these are independent \(N(a,I_k)\) and \(N(0,I_{p-k})\). For \(b=(h,d)\) in the same coordinates,
\[
 \frac{p_b(Y)}{D_*(Y)}
 =2^{k/2}\exp[-\|U-h\|^2/2]\exp[V^\top d-\|d\|^2/2],\qquad
 E_a[p_b(Y)/D_*(Y)]=\exp[-\|a-h\|^2/4].
\]
This follows by two Gaussian integrals. The mean reaches one at \(a=Pb\), the candidate drive closest to the stipulated reference. At \(b=a\) the log comparison has the claimed expectation since \(E_a\|U-a\|^2=k\).

**Optimality.** For any competitor put \(s=D_*/D_G\). Using validity only when reference and actual candidate both equal \(h\in G\) yields
\[
 E_{U\sim N(h,I_k/2),\,V\sim N(0,I_{p-k})}s(U,V)\ \le 1 .
\]
Indeed \(p_h/D_*\) times the \(p_h\) density is exactly this narrower normal product. Average these inequalities over \(h\sim N(a,I_k/2)\). Convolution adds the two tangent variances, giving \(E_a s(Y)\le1\). This integration combines known inequalities for all fixed \(h\) and imposes no random drive on the experiment. Jensen gives \(E_a\log s\le0\). Equality requires \(s=1\) almost surely under the full-support \(p_a\), proving uniqueness. A quantitative form follows from \(E_a[s-1-\log s]\le-E_a\log s\). A log-score shortfall at most \(\epsilon\) therefore implies
\[
 \Pr_a\{|\log(D_G/D_*)|\ge t\}
 \le \epsilon/(t-1+e^{-t}),\qquad t>0 .
\]
This bounds changes on matching Gaussian data. Calibration under that sampling law alone does not bound a rule uniformly on rare deterministic conflicts.

**A characterization using only comparisons.** The proof establishes the stronger inequality \(E_a[D_*/D_G]\le1\) at each drive for every admissible denominator \(D_G\). Thus \(D_*\) can itself be put in the reference position against **any** admissible denominator, even though it integrates to infinity along the candidate graph when \(k>0\). It is the unique admissible factor with this property. If a second \(D\) has it, testing the two against each other gives \(E_a[D/D_*]\le1\) and \(E_a[D_*/D]\le1\). The nonnegative function \(x+x^{-1}-2\) then has zero mean, forcing \(D=D_*\) almost everywhere. Requiring a factor that can serve as this universal reference is an alternative operational premise to expected-log maximization. On an independent block a factor passing all fixed-mass comparisons still requires a separate interpretation before it can be normalized across masses.

The feasible denominators themselves admit an exact transform test. Let \(s=D_*/D_G\) and \(M(x,d)=E[s(U,V)]\) under independent \(N(x,I_k/2)\) and \(N(d,I_j)\), \(j=p-k\). Completing the same squares gives \(E_a[p_{(h,d)}/D_G]=\exp[-\|a-h\|^2/4]M((a+h)/2,d)\). For every transverse mean \(d\) achieved by a reference, let \(H_d\) be the set of achievable tangent means at that \(d\). Uniform validity is equivalent to
\[
          M(x,d)\le \exp[\operatorname{dist}(x,H_d)^2]\qquad\text{at every }x\in G.
\]
Indeed vary \(a=2x-h\), then take the infimum over \(h\in H_d\). For all translations as allowed references this becomes \(M\le1\) everywhere. For the physical reference family \(H_0=G\), the slice used in the optimality proof. This tests the calibration of a proposed change of working factor and makes the allowed preparations explicit.

If the denominator may instead depend on the reference, \(p_{Pb}\) yields the ratio \(\exp[V^\top d-\|d\|^2/2]\). It has mean one under **every** \(a\in G\) and expected log \(\|d\|^2/2\) under \(p_b\). No valid ratio can exceed that expectation because Jensen applied to its quotient by \(p_b/p_{Pb}\) gives the same bound under \(p_b\). At each specified alternative the common denominator incurs the fixed log penalty \(k(1-\log2)/2\) relative to this reference-dependent optimum. Using one denominator permits the reference numerator to cancel in comparisons between mass candidates. It is an operational assumption with a measurable cost.

## What mere validity still allows

Write \(j=p-k\). For each \(0\le\lambda\le1\) the common denominator
\[
 D_\lambda(y)= (2\pi)^{-p/2}2^{-k/2}(2-\lambda)^{-j/2}
                    \exp[-\lambda K(y)/2]
\]
is valid against every Gaussian reference of the declared noise covariance. The same factorization and completion of the square give
\[
 E_a[p_b(Y)/D_\lambda(Y)]
 =\exp[-\|a-h\|^2/4-(1-\lambda)\|d\|^2/(2(2-\lambda))].
\]
Every member has unit mean when the specified reference equals the actual candidate preparation. Yet \(\lambda=0\) leaves a mass prior untouched and \(\lambda=1\) gives the full profile tilt. The matching expected log shortfall relative to \(D_*\) is \(j[1-\lambda-\log(2-\lambda)]/2\), strictly positive when \(j>0\) and \(\lambda<1\). Thus maximizing the expected log in the theorem is stronger than demanding sharp arithmetic-mean calibration, even at every matching reference. These e-comparisons use declared simple reference Gaussians. Representing the inferred law as a Bayesian posterior under a proper distribution of true pushes is a separate requirement. [Gaussian-mixture representation test](e12-conditioning-nuisance.md)

## The consequence on conflicting readings

Independent pushes on the same candidate are a single product Gaussian with graph the direct sum of the individual graphs. Each coordinate of the drive list can be any deterministic value, so the theorem applies to the **entire block** without first assuming that final mass laws multiply. The maximizing divisor equals the product of the individual \(D_*\). The rank and dimension normalization is independent of mass for fixed declarations. If the candidate update is obtained by multiplying the declared prior by its common divisor and normalizing, the theorem selects \(\pi(dm)\exp[-\sum_i K_i(m)/2]\). Returning a law by this normalization is an explicit additional premise. Candidate-wise e-value bounds themselves only concern evidence against fixed candidates.

This choice fails the E12 spread test on two incompatible positive scalar readings. With equal independent errors \(\sigma\), take \((F,A)=(2,1),(1,1)\). Their exact force-to-acceleration ratios disagree. In unit-noise coordinates the scatter is \(\sigma^{-2}\left(\begin{smallmatrix}5&3\\3&2\end{smallmatrix}\right)\). It has an interior principal mass angle \(\hat\theta\) and eigenvalue gap \(3\sqrt5/\sigma^2\). The normalized tilt is proportional to \(\exp[-3\sqrt5\sin^2(\theta-\hat\theta)/(2\sigma^2)]\). Under any positive continuous angular prior it concentrates at \(\hat\theta\) as \(\sigma\downarrow0\): outside a fixed neighborhood it is exponentially smaller than on a smaller neighborhood of positive prior mass. Thus exact maximization of matching expected log comparison, a single denominator across reference preparations, and interpreting that denominator as the final block tilt cannot together supply the required spread. The first two premises do give an independent derivation of the *working* Gaussian factor and an exact calibration baseline for modifying it. [Earlier composition obstruction](e12-product-obstruction.md) · [Numerical and analytic checks](../tests/test_e12_residual_bets.py)


## Arbitrarily small log loss suffices to retain spread

There is no positive lower bound on the additional matching log loss required for the fixed-conflict spread condition. Keep a block with residual dimension \(j=p-k>0\) independent of the candidate. Write \(K_m\) for its sum of squared whitened residuals and \(Q=\inf_{m>0}K_m\). For any budget \(\epsilon>0\), choose \(A_j\ge j\) large enough that
\[
 \tfrac12 E[X{\bf1}_{\{X>A_j\}}]\le\epsilon,\qquad X\sim\chi^2_j.
\]
The left side is \(j\,\Pr(\chi^2_{j+2}>A_j)/2\) by integrating the gamma density. Set \(d=\min(1,A_j/Q)\), with \(d=1\) at \(Q=0\), and replace the exponent \(-K_m/2\) of the optimal block denominator by \(-dK_m/2\). The constant prefactor is the same. This new denominator is at least \(D_*\) everywhere, so all reference tests retain their expectation bound. At every matching true mass and drive, \(Q\le K_m\sim\chi^2_j\). The shortfall from the optimum in expected log comparison is
\[
 E_m[(1-d)K_m/2]\le\tfrac12 E_m[K_m{\bf1}_{\{K_m>A_j\}}]\le\epsilon .
\]
The same budget works at every block length by using the corresponding dimension. If padding with exactly zero observed rows must be neutral for the normalized mass law, omit these rows from the dimension used to select \(A_j\), returning the prior if no rows remain. A nonsingular Gaussian reading equals zero with probability zero under every preparation. The expectation bounds therefore persist.

For a fixed conflicting shape on the scalar positive arc with common noise shrinking to zero, write \(K_m=K_m^0/\sigma^2\) and assume \(Q^0=\inf_m K_m^0>0\), including the endpoint limits. Then eventually \(dK_m=A_j K_m^0/Q^0\) exactly. The continuous bounded costs on that compactified arc give a proper positive tilt with positive interval width; it is informative whenever the shape's residual depends on mass. In the two-reading example above the candidate at the best fit and any distinct interior candidate have limiting odds strictly between zero and infinity for a two-point prior with positive weights. At the exact optimum their odds diverge. Different finite \(A_j\) above the required threshold give different spreads on the same conflict at the same log-loss budget.

For the common isotropic scalar experiment the relaxed rule still concentrates on each identifying array of [T42](../definition/theorems/T42.md). Indeed \(d^{-1}=\max(1,Q/A_j)=O_p(1)\) by \(Q\le K_{\theta_0}\), \(A_j\ge j\) and Markov's inequality. The spectral calculation gives a consistent principal mass angle and an eigenvalue gap tending to infinity in probability; multiplying the gap by \(d\) retains divergence. The normalized angular tilt concentrates at the same limit. This does not establish calibration of its quoted intervals. Exact expected-log optimality can recover the working Gaussian profile. Imposing spread excludes attainment of that optimum under the final-tilt interpretation, and every strictly positive tolerance permits it to be approached. This leaves a real finite-data choice after the covector and exact noise distribution have been fixed.

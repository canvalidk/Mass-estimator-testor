# A partial proof for mass from scalar readings

**Result.** On one scalar Gaussian instrument, a law of positive mass, its median and a companion interval can jointly satisfy measured-zero neutrality, consistency at the sharp identifying threshold, nonvanishing uncertainty on every fixed incompatible series as the errors shrink, and a finite-sample coverage guarantee. Two explicit constructions satisfy all four and give different laws for the same inputs. This proves existence and failure of uniqueness for these stated requirements. The table at the end identifies the extra requirements on which the project's uniqueness question turns. [Project target](../GOAL.md) · [Proof and outstanding premises](../definition/proof.md)

## Start with the reciprocal pairing

Take one physical spatial vector space with metric \(g\). A nonzero acceleration \(a\) gives a normalized covector \(a^\vee=g(a,\cdot)/g(a,a)\), so \(a^\vee(a)=1\). The normalization alone leaves many possible covectors in multiple dimensions. Requiring the construction to respect all spatial orthogonal transformations selects this one. Indeed the subgroup fixing \(a\) contains the reflection reversing \(a^\perp\), so the selected covector must vanish there; its value on \(a\) fixes the remainder. In 3D the half-turn about \(a\) already suffices. Spatial dimension and the physical symmetries assumed matter. [Original reciprocal-vector argument](../docs/mass_as_force_paired_with_reciprocal_acceleration_2026-10-10.md) · [29 October ruling](../rulings.md#r-2026-10-29-3--the-reading-space-form-is-kept-central)

Write the linear inertia relation as \(F^i=M^i{}_j a^j\). Under a change of spatial frame the components of \(M\) transform by conjugation. Physical invariance of the *same inertia* under all spatial rotations in 3D forces \(M=m\,\mathrm{Id}\). In 2D, allowing all linear maps, a reflection as well as rotations supplies the corresponding conclusion. Positivity of \(m\) is a separate physical premise. For an exact nonzero acceleration the scalar then follows by contraction:
\[
                 a^\vee_i F^i=m\,a^\vee_i a^i=m .
\]
In contrast the uncontracted tensor \(F^i a^\vee_j\) on such a pair is \(m\) times the projector onto the acceleration line. A single push measures the action on that line; physical isotropy supplies the extension to other directions. [T71](../definition/theorems/T71.md)

For measured, nonzero \(F,a\), keep both \(a^\vee(F)\) and \(F^\vee(a)\). Their product is \(\cos^2\phi\) for their spatial angle \(\phi\); their signs distinguish alignment from anti-alignment. The contraction gives the longitudinal coefficient even off the common line. What the scalar inertia hypothesis leaves unexplained is a vector:
\[
 F-m a=(a^\vee(F)-m)a+F_\perp,\qquad F_\perp=F-a^\vee(F)a .
\]
For a full spatial reading with independent, spatially isotropic channel errors, the two components are orthogonal and the declared variance of each component of \(F-ma\) is \(\sigma_F^2+m^2\sigma_a^2\). Thus the candidate residual is derived directly from the pairing:
\[
 K(m)=\frac{\|a\|_g^2\,[m-a^\vee(F)]^2+\|F_\perp\|_g^2}
             {\sigma_F^2+m^2\sigma_a^2}.
\]
For example, with unit channel errors, \(F=a=(2,2)\) gives both pairings 1, whereas \(F=(4,1),\,a=(1,4)\) gives both \(8/17\). In the latter case one is a mass and the other an inverse mass, so they disagree as reciprocal answers. Direct contraction gives \(K_{\rm second}(m)-K_{\rm first}(m)=9\) for every \(m\). Normalizing the two Gaussian evidence curves against the same prior discards precisely this constant. The reciprocal pairings already distinguish these readings before fitting any distribution. [Same example as two scalar readings and as one joint vector reading](e12-pooled-information.md)

In general the physical spatial metric fixes the reciprocal and the instrument covariance fixes the size of the observed residual. With full subspace, set \(D_m=[I,-mI]\), \(y=(F,a)\). Then \(D_m y\) is the residual vector, \(D_m\Sigma D_m^\top\) its covariance, and their quadratic contraction is T36's \(K(m)\). This form also works at a measured \(a=0\), where \(a^\vee\) itself is undefined. The Gaussian law below is one way of assigning relative weight to candidate masses using this geometry. The geometric derivation specifies the residual and exposes what normalization can erase. The derivation of a unique law from the permitted contractions, the instrument and the composition premises is still open. The scalar construction below measures how much freedom survives the particular requirements we tested. [T36](../definition/theorems/T36.md) · [T32](../definition/theorems/T32.md)

## What the reciprocal construction proves next

**A series has a reciprocal at each candidate mass.** Move the candidate by a small amount. In each reading the fitted acceleration determines a tangent \(t_i\) to the moving graph in noise units. Apply its reciprocal covector to the unexplained residual \(r_i\). The result \(b_i=t_i^\vee(r_i)\) is the local mass change that best absorbs the discrepancy. Requiring all readings to move by the same amount gives a unique correction \(\bar b=\sum h_i b_i/\sum h_i\), \(h_i=\|t_i\|^2\), when the total tangent is nonzero. The residual splits exactly into within-reading misalignment, disagreement between these signed corrections, and the common displacement. This works for full correlated covariances and declared subspaces. At a fixed true candidate the three squared lengths have independent chi-square laws, conditionally on the fitted graph projections. In the common-shift tangent experiment, translation covariance and exact one-sided conditional coverage at every rank uniquely give the normal law \(N(\bar b,1/H)\) for the displacement, where \(H=\sum h_i\). In the corresponding separate-shift experiment the disagreement sum \(J=\sum h_i(b_i-\bar b)^2\) gives the unique test maximizing worst-direction detection power at a stated false alarm rate. [Proof and block-combination formula](e12-tangent-pairings.md) · [Uniqueness of the local diagnostic](e12-tangent-diagnostic.md)

**What that forces for the mass law.** Using the unmodified noise-metric pairing of the common correction with a displacement as the differential of log evidence integrates uniquely to the ordinary Gaussian law times the chosen prior. But that differential discards conflict common to every candidate. One might try to repair it using a universal decreasing response to the absolute candidate cost. Every such response either collapses on a fixed conflicting pair along some shrinking-noise sequence or fails to concentrate on repeated measurements at a fixed nonzero push. This leaves a concrete obligation for a full proof: derive how the other invariant information, such as the orthogonal disagreement, enters the relative evidence. [Law-level derivation and impossibility theorem](e12-reciprocal-evidence.md)

Three independent consequences test the physical scope of that claim.

* On one isotropic instrument the squared norm of the exterior product of stacked force and acceleration separates into within-reading wedge terms and between-reading tensor terms. For locally compatible readings the between term is a squared difference of masses times the two acceleration energies. No alignment of the spatial frames of different readings is needed. [Exterior-product identity](e12-tangent-pairings.md)
* Choosing a common reciprocal direction independently and isotropically forces a Cauchy law for the projected ratio. Its median on positive mass is the ratio of the vector lengths. Applied to stacked noisy readings it can converge to \(\sqrt{5/2}\) when the true mass is 2. Allowing observation-dependent directions instead can represent any equivariant law. [Direction theorem and limits](e12-random-reciprocal.md)
* For a separable instrument, correlation changes the angular aperture of the positive scalar graphs to \(\arccos(-\rho)\). The old worst-case mean residual \(k-\tfrac12\) survives down to correlation \(-1/2\), sharp in dimension one. Over all admissible covariances the sharp mean is \(k\). Further, two distinct observations on some instruments have identical *absolute* noise-residual curves at every mass yet different physical reciprocal pairings. [Aperture](e12-positive-aperture.md) · [Sharp scalar threshold](e12-aperture-mean.md) · [Correlated-instrument information theorem](e12-covector-information.md)

## Experiment and quantifiers

A reading has two scalar components with independent Gaussian errors of positive standard deviations \(\sigma_F,\sigma_a\). The ratio \(s=\sigma_F/\sigma_a\) is common to the series. Set \(Y_i=(F_i/\sigma_F,a_i/\sigma_a)\), \(\theta=\arctan(m/s)\), \(u_\theta=(\sin\theta,\cos\theta)\) and \(v_\theta=(\cos\theta,-\sin\theta)\). The positive mass domain is the **open** arc \(0<\theta<\pi/2\). Use a fixed proper prior with continuous positive density on this arc and transform it with the units and channel swap. A reciprocal-symmetric prior gives its declared centre as its median.

* In a **sampling experiment**, \(Y_i=r_i u_{\theta_0}+Z_i\), where \(\theta_0\) is fixed in the interior, \(r_i\) are arbitrary deterministic signed amplitudes, and the independent \(Z_i\) are standard bivariate Gaussians. Write \(\Delta=\sum r_i^2,\ V=N+\Delta,\ I=\Delta^2/V\). Identification below means \(I\to\infty\), including fixed \(N\) with growing signal and arrays where one row carries all the signal. Sharpness refers to the uniform Gaussian energy class of T53. At **random true null**, every \(r_i=0\) and the observed rows still contain noise. [T42](../definition/theorems/T42.md) · [T53](../definition/theorems/T53.md)
* In a **fixed-data test**, take a finite list \(x_i\) in noise coordinates at a reference scale, then put \(Y_i=\sqrt t\,x_i,\ t\to\infty\). The observed force and acceleration stay put and both declared errors shrink at the same rate. *Incompatible* means no finite positive \(m\) fits every row. An axis row alone is incompatible in this sense. The all-measured-zero list fits every \(m\) and contributes neutral evidence. The exact-input branch with declared zero error has its own definition. [T30's whole-series condition](../definition/theorems/T30.md) · [T1](../definition/theorems/T01.md)

Orthogonal regression gives the pilot direction. The sampling proof below uses the local normal approximation to a law defined by a data-dependent loss and separately computes the pilot's sampling distribution. The second calculation is what justifies a coverage claim. [Miller on generalized posteriors](https://jmlr.org/papers/v22/20-469.html) · [Yu, Wang and Samworth on eigenvector perturbation](https://arxiv.org/abs/1405.0680)

## 1. A three-reading decision on the per-reading rule

Set the accelerations to \(1,1,1\), the forces to \(0.2,1,1.4\), and both errors to \(\sigma\). In reference units the scatter is \(\left(\begin{smallmatrix}3&2.6\\2.6&3\end{smallmatrix}\right)\) with eigenvalues \(5.6,0.4\). Thus the unique best positive line is \(m=1\) and the total leftover there is \(0.4/\sigma^2\). The three individual leftovers are \(0.32/\sigma^2,0,0.08/\sigma^2\). [O14](../problems/O14.md) · [T36](../definition/theorems/T36.md)

O14's shared-fit per-reading rule freezes \(w_i=\min(1,1/K_i(1))\) with \(w_i=1\) at zero. Here \(K_i(m)=(F_i-m)^2/[\sigma^2(1+m^2)]\). The outer two terms \(w_i K_i(m)\) are uniformly bounded over positive mass for small \(\sigma\). The middle one remains
\[
 \frac{(1-m)^2}{\sigma^2(1+m^2)}.
\]
Outside any fixed angular neighborhood of \(m=1\), this term is at least \(c/\sigma^2\). Inside an angular interval of radius \(\sigma\) about the fit, all three terms are bounded and the prior gives mass of order \(\sigma\). Their weighted Gaussian law therefore concentrates at 1. Every fixed central interval collapses despite the incompatible readings. One vanishing residual is enough to give this counterexample.

For comparison, O14's *simple pooled* rule assigns all rows \(w=\min(1,3\sigma^2/0.4)\). Once uncapped, every mass contrast is independent of \(\sigma\). With a half-Cauchy starting law of centre 1 the original scratch calculation gives the following central 95% intervals. [Original calculation](../tools/leftover_at_series_best_mass.py) · [Reproduction checks](../tests/test_e12_scalar_proof.py)

| Common error | Per-reading weight | Simple pooled weight |
| --- | --- | --- |
| 0.1 | [0.803, 1.348] | [0.497, 2.011] |
| 0.01 | [0.973, 1.029] | [0.497, 2.011] |

This comparison alone selects neither rule. Counting rows in the simple pooled numerator makes an appended exactly zero row change the law once the discount is active. In a different experiment, the per-reading rule **is** point-consistent on every identifying array above, with either a positive-arc or whole-line pilot. Yet in the whole weak-signal window \(\Delta=O(\sqrt N)\), that rule fits noise around its pilot and produces an eigenvalue gap of order \(N\). Coverage of its raw central interval tends to zero at each fixed interior mass. At the true null with this taper the gap divided by \(N\) tends to \(0.33326\ldots\). Two proposed unweighted-Gram adjustments also fail there. [Per-reading proof and scope](e12-per-reading.md)

## 2. A joint scalar existence theorem

**Theorem.** For every \(b,\alpha>0\), the construction below supplies a raw law, a point and a report with the following properties.

1. The raw law is proper and continuous in the rows in total variation, neutral under appending a zero, and equals the prior on all-zero data. Row order, individual row sign reversals, consistent changes of units and channel interchange have their expected actions.
2. Every nonempty fixed incompatible shape gives a limiting raw law with a bounded, strictly positive, nonconstant density relative to the prior. Its central intervals retain positive width. Every fixed nonzero compatible interior shape gives concentration at its common mass.
3. On every identifying array the fitted angle has error asymptotically \(N(0,1/I)\) in the standardized sense. The centered local raw law tends to the matching normal in total variation in probability. The median is consistent and the raw central \(p\)-interval has sampling coverage tending to \(p\), for each fixed \(0<p<1\).
4. A companion interval covers the truth with probability at least \(p\) in every finite sampling experiment and tends to nominal coverage on every identifying array. On a random true-null block it covers the whole positive domain with probability in \([p^2,p]\) for every \(N\ge1\), tending to \(p\) as \(N\) grows.

**Construction.** Write \(B=\sum_iY_iY_i^\top\), \(L\) for its largest eigenvalue, \(H\) for its eigenvalue gap, \(\hat\theta\) for its principal line on the whole circle modulo \(\pi\), and \(Q=\min_{[0,\pi/2]}v_\theta^\top Bv_\theta\). For \(h_b(x)=(1-z)[1+(2+1/b)z]\), \(z=\exp(-bx/2)\), put
\[
 C=\sum_i h_b(\|Y_i\|^2),\quad
 d=\min\{1,(C/Q)(L/Q)^\alpha\},\quad K=dH^2/L.
\]
Use \(d=1\) at \(Q=0\) and zero log evidence at \(L=0\). Define
\[
 \begin{aligned}
 K'&=\min\{K,1+\min(B_{11},B_{22})^2\},\\
 U(\theta)&=-\frac{\sum_i\log\!\left(1+(v_\theta^\top Y_i)^2/(1+\|Y_i\|^2)\right)}
                       {1+\sum_i\|Y_i\|^2/(1+\|Y_i\|^2)},\\
 \ell(\theta)&=-\tfrac12 K'\sin^2(\theta-\hat\theta)
       +e^{-H^2/L}\frac{Q}{1+Q}U(\theta).
 \end{aligned}
\]
Return the normalized law \(\pi(d\theta)e^{\ell(\theta)}\), its median in mass, and the interval described below. The cap handles an axis reading. The last term can distinguish finite observed configurations whose scatter has equal eigenvalues. The choices here are explicit mathematical witnesses. Each ingredient can be recovered from the individual readings’ relative Gaussian mass-evidence curves. L19 may still restrict how those curves are combined. [Information retained by each curve](e12-pooled-information.md) [Pooled construction](e12-pooled-witness.md) · [Executable log evidence](e12_scalar_completion.py)

**Finite data.** Since \(\log(1+x)\le x\) for \(x\ge0\), \(-1\le U\le0\). Every finite tilt is bounded. A zero row contributes exactly zero to each sum. The weight is continuous at a nonzero \(Q=0\) scatter because \(C>0\) there and \(d=1\) nearby. At an eigenvalue tie \(K'\le H\) removes the ambiguous fitted angle. At \(B=0\), \(K'\le L\to0\) and \(Q\le L\). These bounds give uniform convergence in angle at the exceptional points and hence continuity of the normalized law. Channel interchange reflects \(\theta\) into \(\pi/2-\theta\) and interchanges the diagonals of \(B\). The terms are symmetric under that swap and row order, even in each row, and expressed in noise units.

**Fixed-data proof.** Under the scaling above, \(B=tB_0\), \(Q=tQ_0\), \(L=tL_0\), \(H=tH_0\). Write \(n_*\) for the number of nonzero rows. Then \(C\to n_*\). The following cases exhaust the nonzero fixed shapes.

* If \(H_0>0\) and \(Q_0>0\), both diagonal entries are positive and the cap is eventually inactive. The second term tends uniformly to zero. The limiting precision is
  \[
    \frac{n_* H_0^2}{Q_0L_0}\left(\frac{L_0}{Q_0}\right)^\alpha>0 .
  \]
  The limiting sine-squared tilt is bounded and nonconstant on the positive arc, also when the fitted line has negative slope.
* If \(Q_0=0\), all nonzero rows lie along a single line in the closed positive arc, because each squared residual must vanish. On an axis the smaller diagonal is zero, so \(K'=\min(tL_0,1)\). This gives finite nonzero limiting evidence about the direction. In the interior both diagonals grow linearly, \(K'=tL_0\) eventually, and the compatible mass receives all limiting probability in every neighborhood.
* If \(H_0=0\), then \(B_0=L_0\,\mathrm{Id}\) with \(L_0>0\), \(K'=0\), and the multiplier of \(U\) tends to one. Write each nonzero line's angle as \(\psi_i\) modulo \(\pi\). Uniformly on the closed positive arc the limit of \(U\) is
  \[
    -\frac{1}{1+n_*}\log\prod_{i:x_i\ne0}\bigl(1+\sin^2(\theta-\psi_i)\bigr).
  \]
  The product's coefficient of \(e^{2in_*\theta}\) is \((-1/4)^{n_*}e^{-2i\sum_i\psi_i}\ne0\). A trigonometric polynomial with that coefficient cannot be constant on an open arc. Thus even a balanced scatter gives nonconstant limiting evidence in this construction.

Uniform log-evidence convergence in the incompatible cases gives total variation convergence. The limiting density is bounded above and below relative to a positive full-support prior. This proves both the spread and the information claims for a finite fixed list. No uniform lower bound on evidence contrast across *different* shapes is asserted.

**Identifying-array proof.** The pooled calculation gives \(L/V\to1,\ H/\Delta\to1,\ d\to1\) in probability and \(\sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1)\), using the local representative of the projective angle. The score at truth is conditionally normal given the parallel projections, so the calculation permits even one dominant reading. The diagonal scatter entries are \(N+\Delta u_{0j}^2+O_p(\sqrt V)\). For fixed interior \(\theta_0\), each is at least a fixed positive fraction of \(V\) with probability tending to one. Since \(K\le L=O_p(V)\), the cap is inactive with probability tending to one. The remaining log-evidence perturbation is uniformly bounded by \(\exp(-H^2/L)=o_p(1)\). Normalizing two kernels with a uniform log difference \(o_p(1)\) changes their laws by \(o_p(1)\) in total variation, including after the local rescaling. The point and interval assertions follow from the pooled local normal limit. [Detailed estimates and quantile proof](e12-pooled-witness.md)

**Report proof.** Set \(z_p=\Phi^{-1}((1+p)/2)\). Take the interval hull of the raw central \(p\)-interval and all masses satisfying \((v_\theta^\top B u_\theta)^2\le z_p^2u_\theta^\top B u_\theta\). The latter score test accepts the true mass with probability \(p\) on every nonempty Gaussian experiment, conditionally on its parallel readings. Hull containment gives the finite guarantee. On identifying arrays the two sets have the same first-order endpoints, \(\hat\theta+(\mp z_p+o_p(1))/\sqrt I\). At the random null the two boundary tests determine whether the hull is the whole domain. Each passes with marginal probability \(p\), and their joint probability has the stated bounds and limit. For an external independent point combined by the specified inverse-squared log-width rule, a null block changes the point with probability at most \(1-p\) at every sample size. This bound allows positive influence at arbitrarily large \(N\). Every connected interval containing the entire score set strictly overcovers each interior truth in finite nonempty Gaussian experiments. [T45](../definition/theorems/T45.md) · [Coverage, null influence and connectedness proofs](e12-report.md)

## 3. Which part is forced?

Within the class of quadratic angle laws centered on the same Gaussian pilot, suppose precision divided by \(I\) has a deterministic limit \(\rho\in(0,\infty)\) along an identifying array. The local central interval has sampling coverage tending to \(2\Phi(z_p/\sqrt{\rho})-1\). Requiring asymptotically nominal coverage at any one interior level forces \(\rho=1\). This pins down the first-order precision in that class. It allows distinct finite answers even with the starting law and all physical inputs held fixed.

With the project's half-Cauchy prior of centre 5, accelerations \(1,1\), \(\sigma_F=0.1,\sigma_a=0.002\), set \(b=1\). The following are **raw central 95% intervals** for the completed laws, before taking the score hull. Both laws widen for the far pair.

| Exponent | Forces 5.5 and 5.6 | Forces 5 and 6 |
| --- | --- | --- |
| \(1/8\) | [5.410, 5.689] | [5.104, 5.882] |
| \(1/16\) | [5.410, 5.689] | [4.968, 6.007] |

There is also an exact proof of nonuniqueness. On the incompatible series in §1, \(L_0/Q_0=14\) and \(H_0>0\), so the limiting precision in §2 is a fixed positive multiple of \(14^\alpha\). Distinct positive exponents give distinct limiting laws with the same prior. If two rules agreed on every finite series along this path, their limits would agree. Thus the theorem admits an **uncountable family of distinct laws** even without the numerical comparison. Two members also give the requested direction of width change on the project's finite example. Neither exponent has been established as an instrument input. [Computational checks](../tests/test_e12_scalar_proof.py)

Even the choice of point can be tested on a single reachable law. The geometric mean of its quartiles and its median both respect units, channel swap and a reciprocal prior centre on the measured zero. They agree to first order on every identifying array. For one compatible measurement they give 0.747047 and 0.763997 respectively. Connected reports satisfying the same coverage bounds can also differ. [Alternative reports, readouts and their proofs](e12-report.md)

**The data needed for a discount.** The two series \(((2,2),(2,2))\) and \(((4,1),(1,4))\) have exactly the same pooled relative Gaussian evidence about mass at every common error scale. Their scatter matrices differ by \(9\,\mathrm{Id}\), which disappears on normalization. The first series agrees on mass 1. The second reads 4 and \(1/4\). No subsequent use of the normalized aggregate alone, even a nonlinear one, can distinguish them. A Gaussian sampling argument strengthens this to incompatibility of law concentration on all identifying arrays with positive limiting angular spread on every fixed conflict. The individual scalar curves distinguish the series. For a single joint 2D vector reading even its own relative curve can hide a conflict. Keeping that curve and its minimum residual restores every feature invariant under a common orthogonal change of spatial frame, for the full-subspace isotropic instrument. [Minimal counterexample and sampling proof](e12-pooled-information.md)

| Project question | What this branch establishes | What remains |
| --- | --- | --- |
| O14, whole-series leftover | The per-reading rule fails the three-reading test. The completed law passes the stated fixed-data test, including an axis and equal eigenvalues. | Which reading information L19 licenses, and which physical meaning fixes a single discount. An insistence that every Gaussian contrast can only shrink prevents informative evidence at a balanced scatter. |
| L11, G22/G25, identifying data | Both constructions of §2 recover the true scalar mass and give locally calibrated raw intervals under the stated information condition. | Covariance shapes, vector subspaces, multiple instruments and the full-type quantifiers. Literal unrestricted L11 already has an impossibility theorem. [T41](../definition/theorems/T41.md) |
| G28, uncertainty and G6, the point | The companion report has the finite guarantee above. A proper law, median and report can be supplied together. | What a measurement must promise on pure noise, whether the report must be an interval of the raw law, and what selects the median under ratio-scale symmetries. |
| T12 and G27, composition | Multiplication of **final** independent factors with bounded log contrast on every fixed conflicting list prevents even singleton law concentration. The result also holds for Gaussian observations. | Raw multiplication followed by recomputation of a block-level correction remains viable. The joint law of different bodies or routes and its composition protocol remain open. [Proof and precise hypotheses](e12-product-obstruction.md) |

A precise next proof obligation is to characterize which statistics of a conflicting series may change its evidence under L19 and the intended composition rule. The two exponents above and the equal-eigenvalue example give concrete cases on which a proposed characterization must give a verdict. A full uniqueness proof also needs the physical reporting and joint-composition premises.

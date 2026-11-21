# A finite discount curve from the conditional likelihood ratio

Use \(n\) independent scalar force–acceleration readings on one instrument with a common known positive definite covariance. In whitened coordinates each row of \(Y\) has mean \(r_i u_{\theta_0}\) and identity noise, where all signed \(r_i\) are arbitrary. The positive mass hypotheses form an open arc of the projective line. Use length-angle \(\theta\) in the whitened plane. Put \(G=Y^\top Y\), \(L=\lambda_{\max}G\), \(H=L-\lambda_{\min}G\), and \(s=u_\theta^\top G u_\theta\). The drop in twice the profile log evidence relative to the best line on the **whole projective circle** is
\[
 h=L-s=H\sin^2(\theta-\hat\theta).
\]
The ordinary score test can accept both eigenlines because the derivative vanishes at both. The profile drop distinguishes the best line from the orthogonal worst line. This is the Gaussian conditional likelihood-ratio construction of Moreira specialized to the reciprocal mass experiment. The issue here is what relative mass weight its finite-sample ranks could justify. [Moreira, 2003](https://doi.org/10.1111/1468-0262.00438) · [T45 score test](../definition/theorems/T45.md)

## Exact calibration at each candidate

Hold fixed the vector \(U=Yu_\theta\), and let \(V\) be the transverse readings. At a true candidate \(V\sim N(0,I_n)\) independently of \(U\). For \(s=\|U\|^2>0\), let \(X=(U^\top V)^2/s\) and \(J=\|V\|^2-X\). They are independent with laws \(\chi_1^2,\chi_{n-1}^2\), taking \(J=0\) for \(n=1\). A random drop is at most \(t>0\) iff the two-dimensional Gram matrix in this basis is bounded above by \((s+t)I\). Its determinant inequality gives
\[
                 X/t + J/(s+t)\le 1 .
\]
Since \(s+h=L\) for the observed drop, the exact candidate tail probability has one common curve for the entire observed block:
\[
 p_L(h)=E_J\,2\bar\Phi\!\left(\sqrt{h(1-J/L)_+}\right).
\]
Use \(p_L(0)=1\) and return 1 on the zero scatter. At any fixed true candidate \(p_L(h)\) is uniform conditionally on \(U\) for almost every \(U\). The expectation uses a known reference chi-square, independent of the unknown drives. On a realized nonzero scatter \(p_L(h)\) is strictly decreasing in \(h\). In particular its high-rank regions follow the principal line alone on the full circle, with the physical mass arc subsequently restricting which directions can be reported. [Calculator](e12_conditional_likelihood_ratio.py) · [Checks](../tests/test_e12_conditional_likelihood_ratio.py)

## What weight would reproduce all these ranks on the Gaussian scale?

Make one explicit interpretation choice. Assign to a candidate with observed two-sided rank \(p\) the same **relative density to the peak** as a scalar unit Gaussian at that rank. This fixes
\[
 z_L(h)=\Phi^{-1}(1-p_L(h)/2),\qquad
 K_L(h)=e^{-z_L(h)^2/2},\qquad
 P(dm\mid Y)=\frac{K_L(h(m))\pi(dm)}{\int K_L(h(m'))\pi(dm')}.
\]
This is a proper positive law for every starting probability \(\pi\). Here \(p\) is a null tail probability; it acquires the displayed density interpretation only under the extra matching condition. The Gaussian rank is two-sided and the factor is relative to \(h=0\), even when that optimum is outside the positive arc. The earlier T47 used Gaussian *quantiles of a supplied law* for width calibration. This construction uses conditional *sampling ranks* of the data at each candidate. [T47](../definition/theorems/T47.md)

**The discount varies with distance from the best line.** Write \(D_L(h)=z_L(h)^2/h\) for \(h>0\). It lies in \((0,1]\) and has limit \([E\sqrt{(1-J/L)_+}]^2\) at zero. For \(n=1\) it is identically 1. For every finite \(L>0,\ n\ge2\), \(D_L(h)\) strictly decreases with \(h\). Hence one scalar multiplying the Gaussian log factor cannot reproduce all conditional ranks on a nonconstant cost curve. To prove this let \(W=\sqrt{(1-J/L)_+}\) and put \(x=\sqrt h\). Then \(\Phi(z_L(x^2))=E\Phi(xW)\). For \(b>1\) the function \(v\mapsto\Phi(b\Phi^{-1}(v))\) is strictly concave on \([1/2,1)\), since its derivative is \(b\exp[-(b^2-1)(\Phi^{-1}(v))^2/2]\). Jensen gives \(z_L(b^2x^2)<b z_L(x^2)\) when \(W\) is nonconstant. The limit at zero follows by a first-order expansion of \(\Phi\) under the bounded expectation.

## What this law passes and where it fails

**Identifying arrays.** At a fixed interior true angle set \(\Delta=\sum r_i^2\), \(I=\Delta^2/(n+\Delta)\). For any deterministic sequence with \(I\to\infty\), the law centered at the best angle and rescaled by \(\sqrt I\) converges in total variation in probability to \(N(0,1)\), for a bounded angular prior density continuous and positive at truth. Its central intervals have asymptotically nominal coverage. This includes the weak-signal identifying regime and different drive sizes on every row.

Indeed in the true axes the two Gram diagonals are \(n+\Delta+O_p(\sqrt{n+\Delta})\) and \(n+O_p(\sqrt n)\), and the cross term is \(O_p(\sqrt{n+\Delta})\). Thus \(L=n+\Delta+O_p(\sqrt{n+\Delta})\), \(H=\Delta+O_p(\sqrt{n+\Delta})\). For \(\theta=\hat\theta+x/\sqrt I\) with bounded \(x\), \(h(1-J/L)_+\to x^2\) in probability over the data and independent reference \(J\). Boundedness and continuity give \(K_L(h)\to e^{-x^2/2}\) in probability, uniformly on compact \(x\) by monotonicity on either side of the mode. For tail control take the data event \(L\ge n+3\Delta/4,\ L\le 2(n+\Delta),\ H\ge\Delta/2\), of probability tending to one. The chi-square exponential moment bound gives \(\Pr(J>n+\Delta/2)\le e^{-cI}\) for an absolute \(c>0\) (apply the Chernoff bound at order \(\Delta/(n+\Delta)\)). On its complement \(h(1-J/L)\ge c'I\sin^2(\theta-\hat\theta)\). Gaussian tail bounds give \(p_L(h)\le C e^{-c''I\sin^2(\theta-\hat\theta)}\). Also \(e^{-z^2/2}\le C\sqrt{2\bar\Phi(z)}\) for \(z\ge0\) by the Gaussian lower tail bound. The fixed physical arc has length strictly less than \(\pi\), so these estimates give an integrable Gaussian envelope on the rescaled arc. Normalization proves the local limit. The usual principal-line fluctuation \(\sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1)\) gives coverage. [T42 identifying limit](../definition/theorems/T42.md)

**Fixed conflict.** For the incompatible unit-noise rows \((4,1),(1,4)\) scaled by \(\sqrt t\), \(L=25t,\ H=16t\). The best line stays at mass 1, and for every competing fixed angular displacement the law concentrates away from it as \(t\to\infty\). More precisely, outside any fixed neighbourhood of the best line the factor decays exponentially, by splitting \(J\) at \(L/2\). Within angular distance \(t^{-1/2}\), \(z_L^2\le h\) keeps the factor bounded below. Thus the law puts probability tending to one at the compromise mass. It fails the project's requirement that arbitrarily precise contradictory readings retain uncertainty.

**Exactly zero readings.** Appending a zero leaves \(G,L,H\) fixed and increases the degrees of freedom of \(J\). Coupling the new \(J\) to the old one plus an independent \(\chi^2_1\) strictly increases \(p_L(h)\) for \(h>0\). For a positive definite scatter with a principal line inside the physical arc and unequal eigenvalues, the peak factor is still 1, so the normalized law changes for any full-support prior. Omitting exact zero rows recovers measured-zero neutrality without affecting the sampling theorem almost surely. The resulting rule is discontinuous when a nonzero row approaches zero with the other rows fixed. Any continuous law agreeing almost everywhere with the displayed rank-matching rule at each length inherits its failure of measured-zero neutrality, by full support of the Gaussian observation law. This pinpoints the premise conflict for this particular operational definition of relative weight. Neither the finite discount nor its failure applies to all ways of using the covector. [Other covector rank construction](e12-conditional-ranks.md) · [What two-body consistency requires](e12-joint-score-geometry.md)

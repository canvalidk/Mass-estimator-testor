# Four numbers suffice to join the readings, while the law stays underdetermined

Fix one instrument with independent Gaussian errors, isotropic spatial channels, a fixed positive error ratio \(s=\sigma_F/\sigma_a\), full spatial subspace and common dimension \(k\ge1\). The observed \(k\times2\) matrix \(Y_i\) has columns \(F_i/\sigma_F,a_i/\sigma_a\). Each declared reading may have its own true acceleration vector of any size or sign. For an interior candidate put \(\theta=\arctan(m/s)\), \(u=(\sin\theta,\cos\theta)\), \(v=(\cos\theta,-\sin\theta)\). Fix the same prior with positive continuous angular density in every comparison.

**Theorem.** Four additive real numbers suffice for an associative block state that retains the full Gram, has exact measured-zero neutrality, produces informative finite-spread limits on every fixed incompatible shape, and gives consistent, locally calibrated raw laws on every identifying Gaussian array. The same state and join work for a continuum of different readouts. Those readouts agree on **every singleton** and have the **same limiting law for each fixed conflict**. Four is minimal among continuous Euclidean states retaining the Gram on the scalar instrument, by the [lower bound](e12-joining-gram.md). The dimension \(k\) is a known instrument declaration.

## The fourth number uses each reading's own mass residual

The Gram \(B=\sum_iY_i^\top Y_i\) retains the two squared lengths and their inner product. For \(f_i=F_i/\sigma_F,\ b_i=a_i/\sigma_a\), the T36 residual measured in noise units is \(K_i(\theta)=\|f_i\cos\theta-b_i\sin\theta\|^2\). Channel interchange fixes the candidate \(\theta=\pi/4\), or \(m=s\). Evaluating each reading at this common mass gives the additive statistic
\[
 R_i=K_i(\pi/4)=\tfrac12\|f_i-b_i\|^2,\quad c_k=k/(1-2^{-k/2}),\quad
 h_k(x)=c_k(1-e^{-x/2}),\qquad C=\sum_i h_k(R_i).
\]
The factor \(e^{-R_i/2}\) is the maximized Gaussian likelihood at \(m=s\) divided by the unrestricted maximum for that observation. A reading exactly compatible with \(s\) contributes zero for every observed magnitude. The expected contribution is \(k\) at true mass \(s\) for *any* fixed signal. The sum retains how the residual at that one mass is divided among the readings. The response \(h_k\) is a choice of rule. Evaluating it requires only the declared instrument and each reading's own mass residual. The rule estimates no distribution of push sizes and transfers no signal estimate between readings. This places the fourth scalar within the type of width information that L19 permits. L10 supplies no signal-size law, and the rulings do not select this response, the pilot, the precision cap or the tie term below. [T36](../definition/theorems/T36.md) · [L10 and L19](../definition/premises.md) · [L19's width ruling](../rulings.md#r-2026-09-30-5--ruling-1-of-1157-l19s-scope-the-birge-move-is-allowed)

A join of blocks is simply \((B_A+B_B,C_A+C_B)\). Under an independent orthogonal change of the spatial basis of any reading, \(Y_i\mapsto O_iY_i\), each term is unchanged. A consistent change of units leaves standardized inputs unchanged. Interchanging the channels swaps the diagonal entries, reflects the angle, and preserves \(R_i\). Taking two spatial components of one observed vector and declaring them to be separate pushes can change \(C\) even when it leaves \(B\) unchanged. In dimension \(k\ge2\) the relative mass curve can lose the additive floor needed to evaluate \(R_i\). For example the vector reading \(f=b=(2,2)\) and the reading \(f=(4,1),\,b=(1,4)\) have identical relative curves and \(R_i=0,9\), respectively. The absolute T36 residual from the declared data retains that distinction. [What the individual mass curves recover](e12-pooled-information.md)

## One family with fixed singleton behavior

Let \(L\ge l\) be the eigenvalues of \(B\), let \(H=L-l,\ T=\operatorname{tr}B,\ Q=\min_{\theta\in[0,\pi/2]}v^\top Bv\), and let \(\hat\theta\) be its principal line on the full projective circle. For \(B\ne0\) define
\[
 d=\min\{1,(C/Q)(L/Q)^{1/8}\},\qquad
 K=\min\{dH^2/L,\ 1+\min(B_{11},B_{22})^2\},
\]
with \(d=1\) when \(Q=0\). Use the normalized evidence relative to the fixed prior with log factor
\[
 \ell_\tau(\theta)=
 -\frac{K}{2(1+\tau a)}\sin^2(\theta-\hat\theta)
 -e^{-H^2/L}\frac{Q}{1+Q}\frac{(v^\top n_B)^2}{1+T},
 \qquad n_B=(\sqrt{B_{11}},\sqrt{B_{22}}),
\]
where \(\tau\ge0\) is fixed for a given readout and
\[
 a=\frac{4\det B}{(1+T)^3}\frac{C-h_k(R)}{1+C},
 \qquad R=\sum_i R_i=(T-2B_{12})/2.
\]
Return the prior for \(B=0\). The squared exterior area detects failure of the two stacked channel vectors to be collinear. The factor \(C-h_k(R)\) is positive exactly when at least two readings have positive residuals at \(m=s\), since \(h_k(x)+h_k(y)-h_k(x+y)=c_k(1-e^{-x/2})(1-e^{-y/2})\). Both a positive exterior area and a split residual are necessary for the parameter to have any effect. In particular \(a=0\) for every singleton, even when its force and acceleration vectors are nonparallel. Also \(0\le a\le1/(1+T)\), by \(4\det B\le T^2\).

The last angular term uses the two observed channel lengths, with a gate that vanishes at large spectral precision. It is bounded in magnitude by one. The choice to favor their ratio when the Gram is balanced is explicit. For \(B=tI\) the Gram itself has no preferred **line in channel space**. Arbitrary rotations mixing force and acceleration are not part of the physical spatial frame group used here. If one imposes that extra symmetry with a rotation-invariant domain and prior, it would prohibit choosing a preferred line at such a scatter. The positive mass arc, channel roles and spatial rotations are distinct assumptions. [Reciprocal division and physical isotropy](e12-scalar-proof.md#start-with-the-reciprocal-pairing)

At an eigenvalue tie the first term has coefficient zero. As \(B\to0\), \(K\le L\) and \(Q\le L\) force both terms to zero uniformly. At \(Q=0,L>0,C>0\) one has \(d=1\) in a neighborhood. At \(Q=C=0,L>0\), all \(R_i\) vanish. Nearby, \(C\ge h_k(R)\ge R\ge Q\), since \(h_k'(0)=c_k/2>1\) and \(Q\) minimizes the residual over an arc containing \(\pi/4\). Again \(d=1\) in a neighborhood. Elsewhere all terms are continuous and any principal-vector sign cancels. This proves continuity of the normalized raw laws in total variation, propriety, measured-zero neutrality and covariance under the transformations above.

## Fixed conflict and repeated measurements

Scale a fixed nonempty list by \(\sqrt t\). Then \(B=tB_0\), \(C\to c_k n_*\) for \(n_*\) readings with \(K_i(\pi/4)>0\) in the original list, and \(a\to0\) whenever \(B_0\ne0\). If \(Q_0>0\), at least one such reading exists since \(R_0\ge Q_0\). If \(Q_0,H_0>0\), the limiting precision is
\[
 c_k n_*\,\frac{H_0^2}{Q_0L_0}(L_0/Q_0)^{1/8}.
\]
It is finite and positive, the channel cap is eventually inactive and the last angular term vanishes. If \(Q_0=0\), the channels are collinear along a single line in the closed positive arc. An axis gives a limiting precision of one with zero gate. For an interior compatible line the precision diverges and the law concentrates there. If \(H_0=0\), \(B_0=L_0 I\) and the entire log factor tends to \(-\sin^2(\theta-\pi/4)\). These exhaust all nonzero shapes, including negatively sloped lines. On every incompatible shape the log factor converges uniformly to a bounded nonconstant function. Its law retains positive angular variance and differs from the prior. All \(\tau\) give precisely the same limits.

For \(N\) independent vector readings put \(p=kN\), flatten their spatial coordinates, write the \(p\) whitened mean pairs as \(r_j u_{\theta_0}\), and put \(\Delta=\sum_j r_j^2,\ I=\Delta^2/(p+\Delta)\). Assume a fixed interior mass and \(I\to\infty\). At the reference candidate the whitened transverse vector \((f_i-b_i)/\sqrt2\) is \(Z_i+\nu_i\), where \(Z_i\sim N(0,I_k)\) and \(\nu_i\) is its deterministic mean. Gaussian integration gives
\[
 E h_k(R_i)=c_k[1-2^{-k/2}e^{-\|\nu_i\|^2/4}]\ge k.
\]
Because \(k\) is fixed, independence and boundedness give \(C\ge p-O_p(\sqrt p)\) if \(N\to\infty\). At truth, the transverse residual is \(\chi^2_p\), so \(Q\le p+O_p(\sqrt p)\) and \(d\to1\). On bounded-\(N\) subsequences with identifying signal, \(L\to\infty\) and \(Q=O_p(1)\). Even at \(\theta_0=\pi/4\), where increasing the signal leaves \(R_i\) central, \(1/C\le1/h_k(R_1)\) is tight. Indeed the translated Gaussian density is bounded by \((2\pi)^{-k/2}\), so the probability that \(\|Z_1+\nu_1\|\le\epsilon\) is at most this bound times the volume of the radius-\(\epsilon\) ball, uniformly in the signal. Thus \(Q^{9/8}/(C L^{1/8})\to0\) in probability and \(d\to1\) also in this case.

The flattened Gaussian calculation gives \(H/\Delta\to1,\ L/(p+\Delta)\to1\), and \(\sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1)\). Both channel diagonals are of order \(p+\Delta\) at the fixed interior angle. The cap is inactive with probability tending to one, \(a=o_p(1)\), and the precision \(K/(1+\tau a)\) divided by \(I\) tends to one. The tie term vanishes uniformly. Expanding the centered sine-squared kernel on the \(1/\sqrt I\) scale gives a standard normal local law in total variation in probability. Thus the mass median is consistent and central intervals have asymptotically nominal coverage. The existing conditional-score hull also supplies at least nominal finite-sample coverage. [Spectral and quantile estimates](e12-pooled-witness.md) · [Score-hull proof](e12-report.md)

## A finite experiment separates the laws

Use scalar readings \((F,a)=(2,1),(1,2)\) in unit noise. Then \(B=\left(\begin{smallmatrix}5&4\\4&5\end{smallmatrix}\right)\), \(L=9,\ H=8,\ Q=1\). Both reference residuals are \(1/2\), so \(C=2c_1(1-e^{-1/4})>1\), \(K=64/9\) and \(a>0\). The angle at mass 1 is the principal line. At mass \(1/2\) the squared angular sine from that line is \(1/10\). The log evidence ratio favoring 1 over \(1/2\), for \(\tau_1\) minus the same log ratio for \(\tau_2\), is exactly
\[
 \frac{K}{20}\left[(1+\tau_1a)^{-1}-(1+\tau_2a)^{-1}\right].
\]
It is nonzero for every distinct pair of parameters. Normalization and the common tilt cancel. With the half-Cauchy mass prior of centre 1, the angle prior is uniform and channel interchange gives the same median 1. The equal-tail 95% mass interval is \([0.113638,8.799860]\) for \(\tau=0\) and \([0.102600,9.746570]\) for \(\tau=64\). On the project's O14 example, both accelerations are one, \(\sigma_F=0.1,\sigma_a=0.002\), and the prior centre is 5. With \(\tau=0\), forces 5.5 and 5.6 give \([5.409620,5.688507]\), while forces 5 and 6 give \([5.287793,5.708992]\). These are raw equal-tail intervals before the score hull. [O14](../problems/O14.md) · [Implementation](e12_compact_join.py) · [Finite checks](../tests/test_e12_compact_join.py)

The independent-block join, every singleton answer, the regular local precision and each incompatible limit all agree across the family. A finite joint readout remains free under these hypotheses. Agreement across the family on singletons has a different content from agreement with the flat-displacement Gaussian baseline in L10. If L10 requires the *final* singleton factor to equal that baseline up to mass-independent normalization, this family fails that requirement. For instance, on the singleton \((2,1)\) the cap gives angular precision 2, whereas the flat Gaussian factor has precision 5. Here \(Q=0\) and the tie term vanishes. The fourth scalar uses the reading's own allowed residual and introduces no signal-size distribution. Deriving the entire readout under L10 and L19, heterogeneous instrument covariances, physical multi-body composition and the project's reporting premises remains open.

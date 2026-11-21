# Joining the readings at every candidate mass

Consider the fixed isotropic Gaussian instrument of the [four-number construction](e12-joining-laws.md). It has independent channels, a full spatial subspace of fixed dimension \(k\ge1\), and positive noise ratio \(s=\sigma_F/\sigma_a\). The true acceleration of each reading is an arbitrary fixed vector. Fix a proper prior with positive continuous density on the open positive angle arc. This construction uses each reading's residual at *every* shared candidate mass. It proves the same finite ambiguity in the mass law after removing the choice of a reference mass for the fourth statistic.

## An additive curve of measured disagreement

The bounded transform below is the Welsch–Leclerc robust loss applied to the transverse residual norm. Minimizing its sum is equivalent to maximizing the sum of Gaussian residual kernels. We use its minimum to set the uncertainty of the mass law. Its high-precision limit and the two distinct finite readouts are derived below. [Barron, robust loss family, eq. (7)](https://openaccess.thecvf.com/content_CVPR_2019/papers/Barron_A_General_and_Adaptive_Robust_Loss_Function_CVPR_2019_paper.pdf)

In noise units \(Y_i=[f_i,b_i]=[F_i/\sigma_F,a_i/\sigma_a]\). At \(\theta=\arctan(m/s)\), let \(v_\theta=(\cos\theta,-\sin\theta)^\top\), \(B_i=Y_i^\top Y_i\), and
\[
 K_i(\theta)=v_\theta^\top B_i v_\theta=\|f_i\cos\theta-b_i\sin\theta\|^2,\qquad
 h(x)=c_k(1-e^{-x/2}),\qquad c_k=\frac{k}{1-2^{-k/2}} .
\]
These are the reading's own T36 residual and a bounded function of that residual. A block \(D\) stores the Gram \(B=\sum_i B_i\) and the continuous curve \(S_D(\theta)=\sum_i h(K_i(\theta))\) on the *closed* arc \([0,\pi/2]\). Both add when blocks are joined. Form
\[
 C=\min_\theta S_D(\theta),\qquad Q=\min_\theta v_\theta^\top Bv_\theta .
\]
At any candidate the first curve sums separately bounded residual contributions, so that even a very large residual contributes at most \(c_k\) per reading. Its minimum is zero exactly when the observed list admits a common fit on the closed arc. Indeed \(h(x)+h(y)-h(x+y)=c_k(1-e^{-x/2})(1-e^{-y/2})\ge0\). Thus \(S_D(\theta)\ge h(\sum_i K_i(\theta))\) and \(C\ge h(Q)\). If \(Q=0\), at its minimizing angle every individual residual vanishes and \(C=0\). For a singleton \(C=h(Q)\), including a noncollinear vector reading. Minimization occurs *after* joining the curves. The previous construction retains four additive numbers by sampling each reading at a fixed candidate. This construction supplies an exact additive function state and makes no claim of a minimal dimension for that state. [T36](../definition/theorems/T36.md) · [Information retained by the Gram](e12-joining-gram.md)

## The mass law and its fixed-data limits

Let \(L\) be the largest eigenvalue of \(B\), \(H\) its eigenvalue gap, \(T=\operatorname{tr}B\), and \(\hat\theta\) its principal line on the whole projective circle. Put \(n=(\sqrt{B_{11}},\sqrt{B_{22}})^\top\). For \(L>0\) define
\[
 d=\min\{1,(C/Q)(L/Q)^{1/8}\},\quad
 P=\min\{dH^2/L,\,1+\min(B_{11},B_{22})^2\},\quad
 a=\frac{4\det B}{(1+T)^3}\frac{C-h(Q)}{1+C}.
\]
Use \(d=1\) at \(Q=0\). For each fixed \(\tau\ge0\) multiply the prior by the exponential of
\[
 -\frac{P}{2(1+\tau a)}\sin^2(\theta-\hat\theta)
 -e^{-H^2/L}\frac{Q}{1+Q}\frac{(v_\theta^\top n)^2}{1+T}
\]
and normalize. At \(L=0\) return the prior. This is neutral under appending a measured zero and invariant under independent changes of orthonormal spatial frame and consistent changes of units. Interchanging channels reflects the evidence in \(\pi/4\). Continuity of each minimum follows from uniform continuity on the compact arc under changes of a finite list. At \(Q=0\), \(h'(0)=c_k/2>1\) gives \(C\ge h(Q)\ge Q\) nearby, hence \(d=1\) there. At a spectral tie \(P=0\); at zero both angular terms go uniformly to zero. The raw laws are proper and continuous in total variation. The nonnegative coefficient \(a\le1/(1+T)\) vanishes on every singleton, so all \(\tau\) give exactly the same singleton laws.

Scale a fixed shape by \(Y_i\mapsto\sqrt t\,Y_i\). The Gram and every untransformed cost scale by \(t\). Write \(B_0,Q_0,L_0,H_0\) for their original values. The bounded residual minimum has the limit
\[
 C(t)\longrightarrow c_k\,j_*,\qquad
 j_*=\min_\theta \#\{i:K_i(\theta)>0\}.
\]
At an angle attaining this count each positive term saturates, giving the upper bound. For the lower bound, take an arbitrary sequence \(t\to\infty\) and pass to a converging subsequence of minimizing angles. Every residual positive at the limit angle is bounded away from zero along that subsequence and contributes \(c_k\) in the limit. There are at least \(j_*\) of them. In particular, on a conflict with \(Q_0>0\), \(j_*\ge1\). It counts the fewest entire readings one must remove to obtain an exact common fit on the closed arc.

* For \(Q_0,H_0>0\), the limiting precision is \(c_k j_* H_0^2(L_0/Q_0)^{1/8}/(Q_0L_0)>0\). Both channel diagonals are positive, so the cap is eventually inactive. The second angular term vanishes. The limiting sine-squared log factor is bounded and nonconstant even for a negatively sloped principal line.
* For \(Q_0=0,L_0>0\) all readings fit along a single channel line in the closed positive arc. At an interior fit \(P=tL_0\) eventually and the law concentrates there. At a boundary axis one channel vanishes on the whole list, so \(P=\min(tL_0,1)\). This gives finite informative spread when the exact fit occurs only at zero or infinite mass. Treating this case as an interior compatible mass would give the wrong limit.
* For \(H_0=0,L_0>0\), \(B_0=L_0I\), \(P=0\), and the limiting log factor is \(-\sin^2(\theta-\pi/4)\).

On every fixed incompatible shape the log factor converges uniformly, and the limiting law differs from the prior and has positive spread. Every \(\tau\) gives the same limit, since \(a\to0\) on each nonzero shape. The choice of angular term on a balanced Gram uses the physical channel axes. Physical rotations act on the left of \(Y_i\) and leave this Gram unchanged. [Background on the axis and tie cases](e12-scalar-proof.md)

## Uniform control on identifying arrays

For \(N\) independent vector pushes write \(p=kN\). The \(p\) individual spatial-coordinate pairs have independent standard Gaussian errors and mean pairs \(r_j(\sin\theta_0,\cos\theta_0)\) at a fixed interior mass. Put \(\Delta=\sum_jr_j^2\), \(I=\Delta^2/(p+\Delta)\), and assume only \(I\to\infty\). For every \(\theta\) the transverse vector is standard Gaussian plus some deterministic mean \(\nu_i(\theta)\). Gaussian integration yields
\[
 E h(K_i(\theta))=c_k[1-2^{-k/2}e^{-\|\nu_i(\theta)\|^2/4}]\ge k.
\]
Taking the minimum over a continuum of angles requires controlling sampling error uniformly. The following bound is independent of all these means:
\[
 E\sup_\theta\left|\sum_i\{h(K_i(\theta))-Eh(K_i(\theta))\}\right|
 \ \le\ 2c_k\sqrt{2N\log(8N+2)} .
\]
To prove it, introduce an independent copy for each reading with its own corresponding mean. Jensen's inequality bounds the left side by the expected supremum of the sample-minus-copy sum. Each difference is symmetric under swapping that pair, so independent unbiased random signs \(\epsilon_i\) can multiply the differences. The triangle inequality gives at most \(2E\sup_\theta|\sum_i\epsilon_i h(K_i(\theta))|\).

Condition on the observed sample and write each \(h(K_i)\) as the integral of \(\mathbf1_{\{h(K_i)>z\}}\) for \(0<z<c_k\). For each fixed \(z\) this is a threshold on \(K_i\), a constant plus a sinusoid in \(2\theta\). Each nonconstant term can cross that threshold at most twice on the closed arc. Constant terms give a fixed indicator. Between and at the at most \(2N\) crossing points there are at most \(4N+1\) distinct indicator vectors, including ties. For each vector \(w\), \(E_\epsilon e^{\lambda\sum_i\epsilon_i w_i}\le e^{N\lambda^2/2}\) since \(\cosh t\le e^{t^2/2}\). Bounding the exponential of the maximum by the sum over the vectors and their negatives, then optimizing \(\lambda\), bounds its expected absolute maximum by \(\sqrt{2N\log(8N+2)}\). Integrating over \(z\) proves the displayed bound. This argument needs independent readings and allows each true push to be different.

For \(N\to\infty\), the bound and Markov's inequality give \(C\ge p-o_p(N)\). The true-mass residual is \(\chi^2_p\), so \(Q\le p+O_p(\sqrt p)\). Since \(L\ge Q\), \(d\to1\) in probability. If \(N\) stays bounded, \(Q=O_p(1)\) and \(L\to\infty\) in probability. On \(Q>0\), \(Q/C\le Q/h(Q)\), uniformly bounded on every bounded interval including the right limit at zero. Thus \((Q/C)(Q/L)^{1/8}\to0\) in probability and again \(d\to1\). These two subsequence arguments cover arbitrary identifying arrays.

The [spectral proof](e12-pooled-witness.md) on the \(p\) independent coordinate pairs gives \(L/(p+\Delta)\to1\), \(H/\Delta\to1\) in probability and \(\sqrt I(\hat\theta-\theta_0)\Rightarrow N(0,1)\). Both channel diagonals grow at order \(p+\Delta\) at an interior truth, so the cap is inactive with probability tending to one. Also \(a=o_p(1)\), and \(P/[(1+\tau a)I]\to1\). The final bounded angular term vanishes uniformly. The same local sine-squared expansion therefore gives a standard normal limit in total variation in probability for the law of \(\sqrt I(\theta-\hat\theta)\). The median is consistent and raw central intervals have asymptotically nominal coverage. The [conditional-score hull](e12-report.md) retains its finite-sample lower coverage guarantee.

## The disagreement statistic for general T36 instruments

Independent readings can have different positive-definite noise covariances and declared spatial subspaces of dimensions \(1\le k_i\le d\), with \(d\) bounded. Use a single **physical** mass \(m\) across instruments. Whiten each observation and project it onto the span \(W_i\) of the two whitened channel copies of its declared subspace. Call that projection \(z_i\). In homogeneous mass coordinates \(m=x/y\), let \(A=xE_i+yF_i\), where \(E_i,F_i\) are bases for the force and acceleration copies formed from the same subspace basis. The two copies are linearly independent jointly since whitening is invertible. Thus \(A\) has rank \(k_i\) at every projective mass, including both physical endpoints. With \(G=A^\top A,\ u=A^\top z_i\), the T36 cost is
\[
 K_i(m)=\|z_i\|^2-u^\top G^{-1}u.
\]
The projection onto \(W_i\) discards a mass-independent squared length outside that declared span. Inside \(W_i\) the residual has dimension \(k_i\). At each candidate its noise is standard normal in that residual space, possibly with a deterministic mean. Therefore \(E h_{k_i}(K_i(m))\ge k_i\) for all candidates and all true signals, with \(h_{k_i}(x)=c_{k_i}(1-e^{-x/2})\) as above. At the true mass the sum of the untransformed costs has law \(\chi^2_p\) for \(p=\sum_i k_i\).

The same uniform-bound argument applies. At a fixed level \(t\), the sign of \(K_i(m)-t\) is the sign of the homogeneous polynomial
\[
  (\|z_i\|^2-t)\det G-u^\top\operatorname{adj}(G)u .
\]
Its degree is at most \(2k_i\), and \(\det G>0\) on the projective line. Unless identically zero it has at most \(2k_i\) distinct real projective zeros, including a possible root at infinity. At a fixed layer \(0<z<c_{\max}=\max_i c_{k_i}\), the threshold for reading \(i\) is \(h_{k_i}^{-1}(z)\) if \(z<c_{k_i}\), otherwise its indicator is always zero. An identically zero threshold polynomial also gives a fixed indicator. Across the closed physical arc there are at most \(2p\) threshold-equality points and at most \(4p+1\) indicator vectors including all ties. The ghost-sample and random-sign calculation above therefore yields
\[
 E\sup_m\left|\sum_i\{h_{k_i}(K_i(m))-E h_{k_i}(K_i(m))\}\right|
 \le 2c_{\max}\sqrt{2N\log(8p+2)} .
\]
Since \(d\) is bounded, the minimum of the transformed sum is at least \(p-o_p(N)\) as \(N\to\infty\). The minimum \(Q\) of the untransformed sum is at most its realized value at truth, with the \(\chi^2_p\) benchmark above.

Each block can retain the two additive functions \(\sum_i K_i(m)\) and \(\sum_i h_{k_i}(K_i(m))\) on the common physical arc. The latter has zero minimum exactly at an exact common fit on that arc. Scaling a fixed list of \(z_i\) by \(\sqrt t\) with the graph families held fixed gives the limiting transformed minimum
\[
       \min_{m\in[0,\infty]}\ \sum_{i:K_i(m)>0}c_{k_i}.
\]
Continuity on the compact projective arc gives the same minimizing-subsequence proof as above, now with the contribution of each saturated reading weighted by its own \(c_{k_i}\). These conclusions extend the disagreement statistic. The spectral pilot, cap, tie term and calibrated mass law above were proved on the fixed isotropic instrument. Calibration of an entire readout for general T36 instruments remains open. [T36 geometry](../definition/theorems/T36.md) · [Finite checks](../tests/test_e12_residual_join.py)

## Two positive readings distinguish the laws

The numerical implementation brackets the minimum on the full closed arc. On each interval it adds the individual minimum costs after transforming them. It also bounds the sum below by its linear interpolation minus an upper curvature bound times the squared interval width divided by eight. Splitting the interval with the lowest bound controls the objective gap to the specified tolerance, subject to floating-point accuracy. This matters on conflicting high-precision readings with narrow minima at several masses.

For \((f,b)=(2,1),(1,2)\), \(B_{11}=B_{22}=5,\ B_{12}=4,\ L=9,\ H=8,\ Q=1\). Here the minimizing angle of the transformed residuals can be found exactly. The two costs are \(5/2-2\sin(2\theta)\pm(3/2)\cos(2\theta)\). Minimizing their sum after applying \(h\) is equivalent to maximizing
\[
 e^q\cosh(3\sqrt{1-q^2}/4),\qquad q=\sin(2\theta)\in[0,1].
\]
The logarithmic derivative for \(q<1\) is at least \(1-9q/16>0\), using \(\tanh x\le x\). The optimum is \(\theta=\pi/4\), where both residuals are \(1/2\). Thus \(C=2h(1/2)>h(1)>1\) for \(k=1\), \(a>0\), and \(P=64/9\). The difference, for two readouts, of their log evidence ratios favoring \(m=1\) over \(m=1/2\) is exactly
\[
 \frac{P}{20}\{(1+\tau_1a)^{-1}-(1+\tau_2a)^{-1}\}
\]
in unit channel noise. The ratios differ whenever \(\tau_1\ne\tau_2\). With the half-Cauchy mass prior of centre 1 both medians are 1 by channel symmetry. The raw equal-tail 95% mass intervals are \([0.113638,8.799860]\) at \(\tau=0\) and \([0.102600,9.746570]\) at \(\tau=64\), identical here to the fixed-candidate construction because their minimizing and reference angles coincide.

An asymmetric pair changes the reported point as well. With readings \((3,1),(1,2)\) in unit noise and the same prior, numerical integration gives a median of \(1.528611\) and raw equal-tail 95% interval \([0.307038,15.258568]\) at \(\tau=0\). At \(\tau=64\) the median is \(1.510203\) and the interval is \([0.270048,16.247667]\). On O14's two cases, take both accelerations as one, \(\sigma_F=0.1,\sigma_a=0.002\), the half-Cauchy prior of centre 5, and \(\tau=0\). Forces 5.5 and 5.6 give the raw interval \([5.40962,5.68851]\), while forces 5 and 6 give \([5.19840,5.79408]\). [O14](../problems/O14.md) · [Implementation](e12_residual_join.py) · [Finite checks](../tests/test_e12_residual_join.py)

The only common argument supplied to the individual residuals is a candidate mass. The width statistic estimates no reading's signal size, direction or spread for use on another. It fits the type of residual information permitted by L19, including the minimal number of incompatible readings in the fixed-data limit. The choice to transform and minimize those residuals, the spectral pilot, the cap and the angular term on balanced data still require physical justification. Agreement within the family on every singleton also leaves the L10 flat-displacement baseline unresolved. On a scalar unit-noise singleton \((2,1)\) the cap gives precision 2, whereas the flat Gaussian factor has precision 5. Requiring that factor as the **final normalized singleton law** excludes the displayed family. An absolute residual on a spatial reading also retains its mass-independent floor, which normalized relative mass curves can discard. Thus the stated ambiguity concerns the hypotheses proved here. A unique law under the full project premises, including physical composition and reporting, remains open. [L10 and L19](../definition/premises.md) · [Vector reading example](e12-pooled-information.md)

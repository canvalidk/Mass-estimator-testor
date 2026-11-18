# What is lost when the relative Gaussian curves are pooled?

**Proposition.** In the scalar common-instrument model, each reading's relative Gaussian evidence curve determines its whitened scatter exactly. Adding those curves as log factors loses one scalar from the pooled scatter. There are two strictly positive two-reading experiments with identical pooled curves at every common error scale. One has an exact common positive mass. The other has a fixed conflict. Consequently any function only of the pooled curve, the prior, the sample size and the declared instrument must give the same law in both experiments. The distinction survives arbitrary nonlinear post-processing.

## The missing scalar

Write a whitened reading \(y=(f,a)\), \(B=\sum yy^\top\). Up to a constant the pooled Gaussian log evidence in the mass angle is
\[
 -\tfrac12v_\theta^\top Bv_\theta
 =-\tfrac14\{(B_{11}-B_{22})\cos 2\theta-2B_{12}\sin 2\theta\}
 \quad\pmod{\text{constants}}.
\]
Since \(1,\cos 2\theta,\sin 2\theta\) are linearly independent on the open mass arc, two curves agree modulo a constant **if and only if** their scatter matrices differ by a scalar multiple of the identity. Equivalently the retained quantity is the complex number \(T=\sum_i(f_i+ia_i)^2\). Adding \({\rm tr}(B)\) to the record suffices to recover the scatter. Normalizing against any one known positive prior preserves exactly the same equivalence. [Gaussian evidence and the spectral calculation](../definition/theorems/T42.md)

For a *single* row its scatter has rank at most one. Its trace is \(\sqrt{(f^2-a^2)^2+4f^2a^2}=f^2+a^2\), already recoverable from the two coefficients. In particular the measured zero is recognizable. The individual relative curves hence determine every rank-one scatter and each row up to simultaneous sign reversal. They determine \(C,K',Q,U\) in the completed construction, including the energy normalizations. An outside population of signal amplitudes is unnecessary for this recovery. This establishes an exact statement about the information used. A restriction on **how** those curves may affect a different reading's contribution would be a further condition. [L19](../definition/premises.md) · [Scalar construction](e12-scalar-proof.md)

## Two readings suffice

Use a common unit reference error in both channels and the same fixed prior. Consider
\[
 D=((2,2),(2,2)),\qquad E=((4,1),(1,4)).
\]
The respective scatter matrices are
\[
 B_D=\begin{pmatrix}8&8\\8&8\end{pmatrix},\qquad
 B_E=\begin{pmatrix}17&8\\8&17\end{pmatrix}=B_D+9\,\mathrm{Id}.
\]
The duplicate readings agree exactly on mass 1. The other two read 4 and \(1/4\). Multiplication of both scatters by any \(t>0\) preserves equality of the relative curves. Thus even recording the *entire path* of pooled relative evidence over common noise scales cannot distinguish these experiments. Requiring concentration on \(D\) and persistent spread on \(E\) is enough to rule out every such curve-only update. The singleton recovery above proves minimality in reading count.

The three-reading O14 example has the same feature. Three copies of \(\sqrt{13/15}(1,1)\) give a scatter differing from that of \((0.2,1),(1,1),(1.4,1)\) by \(0.4\,\mathrm{Id}\), with the latter having the larger trace. Both examples keep the number of observations and all declarations fixed. [Finite checks](../tests/test_e12_scalar_proof.py)

## A single vector reading and its missing residual

For one full \(k\)-dimensional reading with the same independent isotropic errors on each spatial coordinate, set \(f=F/\sigma_F,\ a=A/\sigma_a\). Its candidate cost is
\[
 K(\theta)=\|f\cos\theta-a\sin\theta\|^2 .
\]
The same two coefficients describe its relative curve, now with the channel scatter \(B=\left(\begin{smallmatrix}\|f\|^2&f\cdot a\\ f\cdot a&\|a\|^2\end{smallmatrix}\right)\). At \(k\ge2\) the rank-one identity need not hold. Take the two displays in the preceding example **as single 2D readings**, with force and acceleration both \((2,2)\) in the first, force \((4,1)\) and acceleration \((1,4)\) in the second. Both declarations now specify one joint reading with the same full subspace and the same covariance. At each fixed reference mass the second cost exceeds the first by 9, identically in mass. Their minimum residuals are respectively 0 and 9. Thus even retaining each normalized reading curve loses the contradiction *within a single vector reading*. The argument uses its declared joint Gaussian cost; it does not treat its two coordinates as separately declared readings. [T32](../definition/theorems/T32.md) · [T36](../definition/theorems/T36.md) · [Vector-reading check](../tests/test_e12_scalar_proof.py)

One additional scalar, the minimum residual \(Q\) on the closed positive arc, recovers the entire curve \(K\): subtract the minimum of any chosen representative of the relative cost, then add \(Q\). This recovers the trace, hence the full channel scatter. Conversely that scatter already classifies a vector pair up to simultaneous orthogonal changes of its spatial frame. Indeed, equality of the two Gram matrices makes the map sending the two vectors to their counterparts a well-defined isometry on their spans, and it extends to an orthogonal map in the ambient dimension. This characterizes exactly the information in one such reading modulo the frame. It supplies no law of mass by itself. The project's T70 compares another loss involving pooling **full** scatters across differently structured experiments. Here dimensions, declarations and count are held fixed, and normalization alone erases the trace. [T70](../definition/theorems/T70.md)

## The obstruction also applies under Gaussian sampling

Exact agreement on deterministic data and consistency in probability on Gaussian observations are different premises. Here one can prove a sampling obstruction as well. Use \(N=2\), fixed prior, equal channel errors \(\sigma\), and let the update be any measurable function of the pooled relative curve and the declarations. Measure spread by the variance of the **angle** of its returned law. This variance is bounded because the arc is bounded. Suppose every fixed incompatible list retains positive \(\liminf\) angular variance as \(\sigma\downarrow0\). In particular, the bounded-contrast property in the main theorem implies this condition. Then the returned law cannot concentrate in probability for every fixed interior signal pair in the Gaussian experiment, even with two equal pushes of growing strength in noise units.

Identify a raw observation with a complex number \(z=F+ia\). At scale \(\sigma\), the observed curve is encoded by \(\sigma^{-2}\sum_jz_j^2\). Let \(g_\sigma(t)\) be the angular variance output at **raw** statistic \(t=\sum_j z_j^2\), absorbing the known scale into the function. Take \(w\) uniform in a compact rectangle \(R\) in the open positive quadrant. Each \(w\) represents a possible true positive mass \(\operatorname{Re}(w)/\operatorname{Im}(w)\). For every such \(w\), there exists a *fixed* incompatible positive list with statistic \(w^2\). For example, for sufficiently small fixed \(\epsilon>0\) use
\[
 z_\pm=\frac{w}{\sqrt2}\left(\sqrt{1+\epsilon^2}\ \pm i\epsilon\right).
\]
Their squares sum to \(w^2\); their directions differ. Compactness inside the quadrant permits the same \(\epsilon\) for every \(w\). Thus \(\liminf_{\sigma\downarrow0}g_\sigma(w^2)>0\) pointwise.

Now take the two true means **both** \(w/\sqrt2\), with independent isotropic Gaussian noise of scale \(\sigma\). The orthogonal sum and difference in row index preserve the sum of squares, so the random raw statistic has the mixture law of
\[
 T_\sigma=(w+\sigma Z_1)^2+\sigma^2 Z_2^2 ,
\]
with independent standard planar Gaussian \(Z_1,Z_2\). This law converges in total variation to that of \(w^2\). Indeed the law of \(w+\sigma Z_1\) tends in total variation to the uniform law on \(R\) by continuity of translations in \(L^1\). Squaring is a measurable map and cannot increase that distance. Also \(w^2\) has a planar density since squaring is injective on \(R\) with nonzero Jacobian. Translating that density by the independent \(\sigma^2Z_2^2\) changes it by \(o(1)\) in \(L^1\), by translation continuity and bounded convergence. These two bounds give the assertion by the triangle inequality.

Take any sequence \(\sigma_j\downarrow0\). Boundedness of \(g_\sigma\), total variation convergence and Fatou's lemma now give
\[
 \liminf_j E\,g_{\sigma_j}(T_{\sigma_j})
 \ \ge\ E\,\liminf_j g_{\sigma_j}(w^2)>0 .
\]
If the returned law concentrated in probability at each fixed true \(w\)'s angle, its angular variance would tend to zero in probability and in mean, and bounded convergence over \(w\) would make the left side zero. This contradiction allows data-dependent functions and pointwise conflict bounds. It concerns law concentration. A spread-out law can still have an accurate median.

This proof identifies which compression fails on one instrument. The trace repairs the algebraic loss in the pooled scatter. It need not recover all individual curves or settle whether scatter alone should supply informative evidence at every balanced configuration. The [duplicate-reading obstruction](e12-product-obstruction.md) addresses a distinct requirement that final factors multiply even when arbitrary per-row information is available.

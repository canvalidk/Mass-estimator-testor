# What a reciprocal residual curve retains about a vector reading

Start with the candidate force residual \(r(t)=f-tx\). For every spatial covector \(b\), the observed disagreement in its direction is \(b^\top r(t)\). Its squared size in its own declared variance is \((b^\top r)^2/(b^\top V(t)b)\). Consequently the strongest directional disagreement is
\[
 K(t)=\sup_{b\ne0}\frac{(b^\top r(t))^2}{b^\top V(t)b}
     =r(t)^\top V(t)^{-1}r(t).
\]
Cauchy–Schwarz in the covariance metric proves this identity. At nonzero residual equality requires \(b\) proportional to \(V(t)^{-1}r(t)\). Thus one candidate mass selects a most sensitive testing covector from the reading and its covariance. It varies with the candidate. This is the dual version of the projection in [T36](../definition/theorems/T36.md). Here we ask the inverse question. Given these extrema at all candidate masses, can one reconstruct the original physical pairings \(a^\vee(F),F^\vee(a)\)?

## Full spatial subspace and declared covariance

Work with one full \(k\)-dimensional reading. Choose fixed positive force and acceleration units \(u_F,u_a\), let \(f=F/u_F,\ x=a/u_a,\ t=m u_a/u_F\), and transform the joint covariance by the same change of units. Use the same spatial basis in both channels. Write its blocks as \(\Sigma=\left(\begin{smallmatrix}A&B\\B^\top&D\end{smallmatrix}\right)>0\), giving the matrix polynomial
\[
 V(t)=A-t(B+B^\top)+t^2 D .
\]
Gram pairings of observed physical vectors use their physical spatial metric \(g\), which in orthonormal axes is the identity. Effective noise coordinates can carry a different spatial metric. Multiplying a dimensionless forward pairing by \(u_F/u_a\) restores mass units, and reversing the ratio restores inverse mass for the reverse pairing.

The **relative curve** is \(K(t)\) up to one additive constant on the entire open positive mass interval. Equivalently it is the normalized raw Gaussian mass law with the same known positive starting law. The **absolute curve** is \(K(t)\) itself. The relative curve plus its infimum on positive mass specifies the absolute curve. The following result concerns exact functions of mass. Inverting residues from noisy evaluations on a grid can be ill-conditioned. No noise-free estimate of the underlying true vectors is proposed.

## Reconstruction for a simple quadratic spectrum

**Theorem.** Suppose the determinant \(h(z)=\det V(z)\) has \(2k\) simple complex roots. For fixed declared \(\Sigma\), the relative curve determines a real reading up to \(k\) independently choosable signs on explicitly defined pairs of real coordinates. If every pair is nonzero, exactly \(2^k\) readings give that curve. The absolute curve can remove some or all of these ambiguities apart from global sign.

The spectral objects are those of the standard quadratic eigenvalue problem. The proof uses its resolvent and an elementary adjugate argument. The observations below identify the equivalence classes of *readings* for this experiment. [Tisseur and Meerbergen, The quadratic eigenvalue problem](https://eprints.maths.manchester.ac.uk/466/)

1. \(V(t)\) is positive definite for every real \(t\), by its covariance representation. The leading block \(D\) is positive definite, so \(h\) has degree \(2k\). Its roots come in nonreal conjugate pairs. Choose the \(k\) upper-half-plane roots \(z_j\) and a nonzero complex covector \(w_j\) with \(V(z_j)w_j=0\), using coordinates in the dual spatial basis. At a simple root the nullity is one, since a larger nullity would make every cofactor zero and force \(h'(z_j)=0\). The symmetric adjugate is then \(\gamma_j w_jw_j^\top\) with \(\gamma_j\ne0\). Transpose here is bilinear, without complex conjugation.
2. The rational continuation of \(K\) has numerator \(r(z)^\top{\rm adj}(V(z))r(z)\) and denominator \(h(z)\). Its residue at \(z_j\) is \(\gamma_j [w_j^\top(f-z_jx)]^2/h'(z_j)\), including a zero residue if the pole cancels. Equality of real relative curves on the positive interval gives equality of all residues by polynomial identity. Conversely, equal residues at all roots imply that the difference is a polynomial. It is constant, since each curve tends to \(x^\top D^{-1}x\) at infinity.
3. Set \(c_j=w_j^\top(f-z_jx)\). The map from real readings to these \(k\) complex contractions is a real-linear isomorphism. In its kernel the vector rational function \(V(z)^{-1}(f-zx)\) has removable poles at each upper and lower root by the rank-one adjugate formula and conjugation. It is then polynomial and tends to zero at infinity as \(O(1/z)\), so it vanishes identically. Hence \(f=x=0\). Real dimensions on the two sides are both \(2k\). The residue determines \(c_j^2\), and each independent sign for a nonzero \(c_j\) is achievable via the inverse linear map. This proves the theorem. Independent complex sign choices refer to simultaneous signs on each pair \((\operatorname{Re} c_j,\operatorname{Im} c_j)\). The actual \(f,x\) are real.

Simple roots occur on an open dense subset of positive-definite joint covariances. A repeated root requires the discriminant polynomial of \(h\) to vanish. The polynomial is nonzero for \(D=I,\ B=0\), \(A\) diagonal positive with distinct entries. The simple-root condition therefore holds generically in the unrestricted type. Repeated roots are physically important, including the spatially isotropic instrument.

**Exactly what adding one absolute residual recovers.** Write \(\xi=(\operatorname{Re}c_1,\operatorname{Im}c_1,\ldots,\operatorname{Re}c_k,\operatorname{Im}c_k)\). By the isomorphism, \(K(\infty)=x^\top D^{-1}x=\xi^\top G\xi\) for a known real symmetric \(G\). Join indices \(i\ne j\) when the \(2\times2\) block \(G_{ij}\ne0\). Changing the normalization of any \(w_j\) applies an invertible coordinate change within its own block and preserves these connections.

A sign matrix \(S=\operatorname{diag}(s_1I_2,\ldots,s_k I_2)\) preserves this quadratic form for all data iff \(s_i s_j=1\) at each connection. The allowed signs are constant on each connected component. For any other sign choice the exceptional readings satisfying \(\xi^\top(S^\top GS-G)\xi=0\) lie in a proper algebraic hypersurface. There are finitely many sign choices. Outside their union and the zero-mode loci, the absolute curve leaves **exactly** \(2^b\) readings when there are \(b\) connected components. In particular, for a given instrument with a connected graph, the whole reading and its physical reciprocal pairings are generically recoverable up to simultaneous sign from the absolute curve. Both pairings are unchanged by that global sign. This conclusion holds almost surely under any nonsingular Gaussian sampling law. It need not hold on every fixed conflicting series.

## Two fixed instruments

**Cross-channel correlation can remove sign ambiguity when absolute residual is retained.** Take the full two-dimensional instrument
\[
 A=\begin{pmatrix}2&.3\\.3&3.7\end{pmatrix},\quad
 B=\begin{pmatrix}.4&-.15\\.8&.1\end{pmatrix},\quad
 D=\begin{pmatrix}1&.1\\.1&.8\end{pmatrix}.
\]
The smallest eigenvalue of the joint covariance is about 0.629. The upper roots are near \(0.474586+1.428980i\) and \(-0.025219+2.020063i\). For \(f=(1,-2),x=(.7,1.2)\), reversing just the first contraction reconstructs approximately \(f'=(-1.049856,-5.130173)\), \(x'=(-2.654676,.429211)\). Their absolute residual curves differ by the constant \(5.551836\). Thus this instance has four readings with the same relative curve, and their two different absolute curves each leave just the global sign pair. The theorem proves the constancy. The finite check reconstructs the alternatives using a generalized linear eigenvalue problem and verifies the offset at multiple positive masses against direct physical residuals.

**Even the absolute curve can lose the physical reciprocal.** Set \(g=I\) and \(S=\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)\). Use \(D=SS^\top,\ A=S\,\operatorname{diag}(1,4)S^\top,\ B=0\). A common spatial transform by \(S^{-1}\) makes the noise matrices diagonal. In those coordinates reverse the second entry of both vectors. In the *same original physical axes*, the two readings are
\[
 (f,x)=((3,2),(2,1)),\qquad (f',x')=((-1,-2),(0,-1)).
\]
Both absolute residuals equal \((1-t)^2/(1+t^2)+(2-t)^2/(4+t^2)\) at every mass. Yet the ordered forward and reverse physical pairings are \((8/5,8/13)\) and \((2,2/5)\). The sign flip in the diagonal noise coordinates is an isometry for both covariance blocks. Conjugated back to physical space it is a shear followed by a reflection and does not preserve \(g\). This is an equality at fixed covariance and fixed physical metric, with two different observed readings. It exhibits exactly the ambiguity that is lost if the noise metric silently replaces the spatial metric.

There is a concrete symmetry choice here. The active map \(J=S\operatorname{diag}(1,-1)S^{-1}\), applied to both channels while the declarations stay fixed, sends each latent graph point \((t h,h)\) to another point on the same graph. It preserves the channel covariances. The Gaussian sampling family at each fixed mass is therefore invariant under this map with a relabeling of the unknown push. Requiring an update to respect all such symmetries of the statistical experiment would force identical output on these two readings. Euclidean physical frame covariance alone does not impose that requirement, since \(J^\top gJ\ne g\). A passive coordinate change would transform \(g\) as well. Neither physical frame covariance nor the existence of the reciprocal alone requires a different answer on the two readings. The example locates a choice of invariance premise with observable consequences.

## Repeated isotropic groups

For uncorrelated channels, transform the same spatial coordinates of \(f,x\) to make \(D=I\) and \(A\) diagonal by congruence. Group equal diagonal entries \(s_j^2\) into eigenspaces and use the Euclidean norm in these noise coordinates. Polynomial division gives
\[
 K(t)=\sum_j\|x_j\|^2+
 \sum_j\frac{\|f_j\|^2-s_j^2\|x_j\|^2-2t(f_j\cdot x_j)}{s_j^2+t^2}.
\]
Distinct denominators are coprime. The relative curve therefore specifies each difference \(d_j=\|f_j\|^2-s_j^2\|x_j\|^2\) and contraction \(e_j=f_j\cdot x_j\). In a one-dimensional group \(\sqrt{d_j^2+4s_j^2e_j^2}=\|f_j\|^2+s_j^2\|x_j\|^2\), recovering both squared norms. In a group of dimension at least two the Gram of \((f_j,s_jx_j)\) can be incremented by \(cI_2\) for any \(c>0\) while preserving \(d_j,e_j\); the enlarged Gram is realizable by two vectors in that group. Its continuous ambiguity is removed by one group trace. For a single spatially isotropic group this is the E12 loss of absolute residual under normalization. The [two-reading comparison](e12-pooled-information.md) treats that case.

**Consequence for the proof.** The covariance metric and the physical metric have distinct measurable consequences. If the intended law uses the physical reciprocal pairings, even retaining every absolute Gaussian residual can be insufficient on the full reading type. Conversely, for a simple spectrum and a connected graph, these curves generically suffice to recover all simultaneous-sign-invariant uses of the reading. The existence of directional testing covectors alone supplies no probability measure on possible masses. A selection principle must specify which contractions may enter the update, and with which physical and instrument tensors. [Numerical checks](../tests/test_e12_covector_information.py) · [Main scalar proof](e12-scalar-proof.md)

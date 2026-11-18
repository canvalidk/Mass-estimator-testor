# Dividing by a shared reciprocal direction

**Question.** For a single spatial pair \((F,a)\), can symmetry select a law for \(m=b(F)/b(a)\)? The same covector \(b\) appears in both contractions, so changing its sign or length leaves the ratio unchanged. The reciprocal of \(a\) along this direction is \(b/b(a)\). The project's October argument already observed that varying the direction through a noncollinear pair produces every mass. Here are the resulting law under independent isotropic directions, two uniqueness statements with different assumptions, and an identifying experiment that this law fails. [Original argument, §8.3](../docs/mass_as_force_paired_with_reciprocal_acceleration_2026-10-10.md)

## 1. The law from a direction selected independently of the pair

Let \(g\) be the Euclidean spatial metric in dimension \(d\ge2\). Use a fixed physical unit system. Take the probability law of the unoriented covector direction \([b]\) to be independent of \(F,a\) and invariant under \(O(d)\). This measure is unique: averaging any continuous function over rotations of one line gives the same value for every starting line by transitivity, and integrating this identity fixes its integral under any invariant measure. Equivalently take \(b=g(Z,\cdot)\) for a standard isotropic Gaussian \(Z\); only its line enters. This auxiliary \(Z\) samples directions, with no assumption on the distribution of the true push. The standard probability input is the ratio of two centered jointly normal variables. [Marsaglia (1965)](https://doi.org/10.1080/01621459.1965.10480783)

First take noncollinear \(F,a\), write
\[
 A=g(a,a),\quad B=g(a,F),\quad C=g(F,F),\quad D=AC-B^2>0 .
\]
In an orthonormal basis starting along \(a\), then along \(F-(B/A)a\), the terms in the quotient are driven by two independent standard normals \(Z_1,Z_2\). Thus
\[
 M=\frac{b(F)}{b(a)}=\frac BA+\frac{\sqrt D}{A}\frac{Z_2}{Z_1},
 \qquad f_M(m)=\frac{\sqrt D}{\pi(A m^2-2Bm+C)} .
\]
The slope is finite almost surely. Its location \(B/A=a^\vee(F)\) is the forward reciprocal pairing. Its scale is the transverse length divided by \(\|a\|\). The denominator is exactly \(\|F-ma\|^2\). There is no fitted exponent. The Cauchy law describes random physical directions *conditional on the observed vectors*. Turning it into the law returned by a noisy mass measurement is a further choice.

**Positivity and point.** Write \(r=\sqrt{C/A}\) and \(\rho=B/\sqrt{AC}\). The probability of a positive quotient is \(p_+=\tfrac12+\arcsin(\rho)/\pi\). Condition on \(M>0\). The dimensionless \(X=M/r\) has density
\[
 f_{+}(x)=\frac{\sqrt{1-\rho^2}}{\arccos(-\rho)}
           \frac{1}{1-2\rho x+x^2},\qquad x>0 .
\]
Since \(x^{-2}f_+(1/x)=f_+(x)\), its unique median is 1. The positive mass median is exactly \(\|F\|/\|a\|\), including for an obtuse noncollinear pair. Channel interchange sends this law to \(1/M\) since that equality holds for each usable direction.

**Boundary cases.** At \(F=m_0a\), \(a\ne0\), the unconditioned quotient is exactly \(m_0\) almost surely, including in dimension one. Positive conditioning gives the point mass for \(m_0>0\) and has zero conditioning probability for \(m_0\le0\). The interior formula does give a unique *continuous limit* at a nonzero antiparallel pair: as \(\rho\downarrow-1\) its positive density for \(X\) tends to \(1/(1+x)^2\). This follows by the limit of the normalizing constant and dominated convergence for \(\rho\le0\). It is an extension of the formula to the boundary. If \(a=0\) and \(F\ne0\), the projective quotient is infinity almost surely. If \(F=0\) and \(a\ne0\), it is zero. At \(F=a=0\) the construction supplies no direction quotient. Scaling a given noncollinear pair toward zero also leaves its ratio law fixed. Distinct such rays have different limiting laws, so no data-independent assignment at zero makes this update continuous there. A physical instrument can return a prior on measured zero; that requires additional structure.

## 2. Which symmetry would force this law?

Express the two channels in fixed dimensionless reference units, and write
\[
 T=[F\ a],\qquad G=T^\top gT=\begin{pmatrix}C&B\\B&A\end{pmatrix}.
\]
Consider a rule taking only a full-rank pair to a probability on the **whole real projective line**. Spatial orthogonal invariance means the output factors through \(G\), since two pairs with the same Gram are related by an orthogonal transformation. Suppose in addition that applying *any* \(H\in GL(2,\mathbb R)\) to the channels changes the output by its projective action. Explicitly \(T'=T H^\top,\ G'=HGH^\top\), and \([b(F'):b(a')]=H[b(F):b(a)]\). Then the preceding unconditioned law is the **unique** such rule.

Indeed, at \(G=I_2\) the orthogonal channel group preserves the input Gram and acts transitively on projective slopes, so the output there must be the uniform angle law. Every positive-definite \(G\) is \(H I_2 H^\top\); two choices of \(H\) differ on the right by an orthogonal factor. Equivariance fixes \(P_G=H_*P_I\). Taking the direction of \(H Z\), \(Z\sim N(0,I_2)\), realizes that transport, whose slope density is exactly the formula in §1. For a nonzero rank-one Gram a continuous extension from full rank gives its supported line (couple the Gaussian representatives through \(G^{1/2}Z\)). No rule on the zero Gram alone can obey this full group covariance: orthogonal invariance would require uniform projective direction there, and diagonal channel stretching changes that distribution.

This is a conditional uniqueness result for **mixing force and acceleration channels**. Rotating physical space acts on the other side of \(T\). A passive coordinate change transforming the covariance, prior and candidate family as well is also a different premise. General channel mixing changes which cone represents positive mass and can introduce a pole at a physical mass. Independent changes of units and the channel swap by themselves do not imply this larger group. Conditioning the result on positivity follows only after designating the positive physical cone.

**Conversely, equivariance of a data-dependent direction need not constrain the resulting mass law at all.** For a prescribed \(m>0\) define
\[
 b_m(v)=(m,1)G^{-1}T^\top g v .
\]
It satisfies \(b_m(F)=m,\ b_m(a)=1\). It is the unique minimum-dual-norm covector obeying those two equations: its Riesz representative lies in the span of the two observed vectors, and the orthogonal complement adds norm without changing either constraint. It transforms with the spatial frame. Interchanging the channels and replacing \(m\) by \(1/m\) gives the same direction, with covector scaled by \(1/m\). Positive changes of units have the corresponding action: for \(F'=uF,\ a'=v a\), use \(m'=u m/v\) and \(b'_{m'}=b_m/v\). Consequently sampling \(m\) from *any* specified equivariant law and taking \([b_m]\) realizes that very law as a reciprocal-direction quotient. This proves an exact distinction between exogenous isotropic sampling of \([b]\) and spatial covariance alone. The converse is representational, so it supplies no method for choosing the target law.

**When the channel-noise tensor is retained.** Suppose the instrument carries a positive-definite \(2\times2\) channel covariance \(S\) with isotropic spatial blocks \(S_{uv}\,\mathrm{Id}\). A change of dimensionless channels \(H\) now transports *both* \(G\) and \(S\) by congruence. For a law on the whole real projective line depending only on these two tensors, full channel covariance has a different implication. It permits an explicit classification. Whiten \(S\) and take orthogonal eigenaxes of \(S^{-1/2}GS^{-1/2}\) with ordered eigenvalues \(\lambda_1\ge\lambda_2\). With unequal eigenvalues its stabilizer consists exactly of the sign flips of the two eigenaxes. Thus for each ordered eigenvalue pair any law for the acute angle from the first eigenline, lifted by equal probabilities for the two signs of that angle, gives an equivariant rule. At angle zero or \(\pi/2\) the two projective lines coincide. Every such rule arises in this way: diagonalization reduces each orbit to that ordered representative, and the sign flips are the only ambiguity in transporting its law back. A measurable assignment in the eigenvalues gives a measurable rule. Consequently covariance of the declarations does leave freedom even on the full slope line.

At \(G=tS\), including \(t=0\), the whitened stabilizer grows to all of \(O(2)\). The whole-line law there must equal the law of the line of a \(N(0,S)\) channel vector, for every \(t\). Equivalently, an intrinsic **relative evidence** computed only from these tensors, defined on the whole projective line, channel-equivariant up to an overall constant and carried separately from the prior, cannot distinguish slopes at this balanced configuration. For a positive bounded continuous factor this follows by normalizing against the unique stabilizer-invariant probability and applying transitivity. A prior or a positive cone taken as a further structural input can reduce the stabilizer. These are constraints under the stated channel symmetry. Physical spatial rotations by themselves do not force flat relative evidence there.

## 3. Testing this direction law as a mass measurement

On the branch's two single vector readings, \(F=a=(2,2)\) gives mass 1 for every usable direction. For \(F=(4,1),\,a=(1,4)\), the full direction law is Cauchy with location \(8/17\) and scale \(15/17\); its positive median is 1 and its width is positive. Thus the covector calculation preserves the distinction that disappears when the Gaussian evidence is normalized. An orthogonal pair with equal lengths in noise units gives the half-Cauchy of scale one after conditioning. With that common starting law its output on this conflicting shape coincides with the starting law. [Two-reading obstruction](e12-pooled-information.md)

There is also a sharp failure on identifying data. For the following experiment apply the same construction to two vectors obtained by stacking scalar readings across the series. This is a proposed extension; the rotation group in row index is larger than the group's physical action on an individual spatial reading. Let the \(i\)th force be \(2+\epsilon_i\) and acceleration \(1+\delta_i\) with all errors independent standard Gaussian. The true mass is 2 and the signal per reading is fixed and nonzero. The strong law gives \((A,B,C)/N\to(2,2,5)\) almost surely. Thus the fitted direction law tends to Cauchy with location 1 and scale \(\sqrt{3/2}\) before positive conditioning. Its positive median tends to \(\sqrt{5/2}\), and its positive central intervals retain nonzero width. The scalar energy grows as \(\Delta=5N\), so \(I=\Delta^2/(N+\Delta)=25N/6\to\infty\). This is an identifying Gaussian sequence on which even the direction-law median is inconsistent. The failure needs neither assumed random pushes nor a simulation. It names the omitted quantity: the raw pairings incorporate error energy at the same order as physical signal energy on this sequence.

## 4. Changing the directional metric cannot give the Gaussian local shape

A choice of positive-definite spatial metric, including a choice computed from the observations, still produces a conditional Cauchy location and scale when the auxiliary direction is uniform in that metric. Consider any sequence of these positive-conditioned Cauchy laws that concentrates at a fixed interior mass \(m_0\). Then their medians \(r\to m_0\). In the dimensionless expression in §1 the shape must have \(\rho\to1\), since any subsequence with \(\rho\to\rho_*<1\) has a nondegenerate limiting density, including the stated limit at \(-1\). Thus \(\mu=r\rho\to m_0\), \(\gamma=r\sqrt{1-\rho^2}\to0\), and the discarded Cauchy probability to the left of zero tends to zero. These statements also hold in probability for random laws concentrating in probability, by taking almost surely converging subsubsequences.

For any fixed interior rank \(u\), the centered quantile in units of \(\gamma\) consequently tends to \(\tan(\pi(u-1/2))\). If \(W_p\) is the central probability-\(p\) width, this forces
\[
 \frac{W_{0.95}}{W_{0.50}}\ \longrightarrow\ \frac{\tan(0.475\pi)}{\tan(0.25\pi)}
 =12.706\ldots .
\]
Changing to the smooth mass-angle coordinate at an interior truth leaves the limiting ratio unchanged. A locally Gaussian law instead has ratio \(\Phi^{-1}(0.975)/\Phi^{-1}(0.75)=2.906\ldots\). This rules out the *full* locally Gaussian shape of the regular witness for any repair that only changes the metric within this direction construction. It does not prohibit a coverage calibration at one level, or frequentist coverage obtained from additional randomization.

The geometric law therefore supplies a conditional answer to what the direction sweep means. Converting its spread over possible projections into uncertainty about the *one* physical mass requires a rule that uses the declared instrument noise. [Finite and rotation checks](../tests/test_e12_random_reciprocal.py)

# What an updated mass law must remember

Fix one scalar instrument, its noise ratio and a positive continuous angular prior. Write \(y_i=(F_i/\sigma_F,a_i/\sigma_a)\). Two blocks are *equivalent for future readings* for a specified output if appending any common finite list gives identical outputs, including appending the empty list. Any state from which the output can be updated exactly must distinguish inequivalent blocks. A state that identifies precisely the equivalent blocks captures the minimal information in this sense. Its numerical encoding need not be unique.

## Tracking the reciprocal disagreement tensor

For \(B=\sum_i y_i y_i^\top\), the candidate-independent squared exterior norm on one instrument is \(\det B\). After observing row \(x\),
\[
 \det(B+xx^\top)-\det B=x^\top\operatorname{adj}(B)x .
\]
The two-by-two determinant identity holds also for a singular or zero block. Thus two blocks are equivalent for the future exterior norm iff their **whole** scatter matrices agree. Indeed, the empty extension determines \(\det B\). For \(x=(1,t)\) with \(t=1,2,3\) the increments give three values of the quadratic polynomial \(B_{22}-2tB_{12}+t^2 B_{11}\), fixing all three entries. Conversely every extension adds its scatter to \(B\). These positive scalar probes keep the measurement declarations fixed.

This is a predictive reason to retain more than the current norm. The same three numbers can be organized as the ordinary relative Gaussian mass curve and the current squared exterior norm. The curve gives \(B_{11}-B_{22}\) and \(2B_{12}\). Positivity gives \(\operatorname{tr} B=\sqrt{(B_{11}-B_{22})^2+4B_{12}^2+4\det B}\). Together these recover the scatter. The retained state suffices for every function of the pooled scatter, including Gaussian updates and recomputed tensor norms. [Exterior-product derivation](e12-tangent-pairings.md) · [Information in the relative curve](e12-pooled-information.md)

If one tracks both within-reading and between-reading **exterior** squared norms, the corresponding minimal predictive state is \((B,L)\), \(L=\sum_i\det B_i\), with \(B_i\) the per-reading channel Gram. The present within norm equals \(L\), the between norm is \(\det B-L\), and future blocks add their scatters and their individual determinant sums. Knowledge of both present and future norms distinguishes \(B\) by the same three probes, which can be declared vector readings with collinear components. For two full vector readings the extra number can vary independently: for any positive-definite \(B\), individual Grams \(sB,(1-s)B\) give \(L=(s^2+(1-s)^2)\det B\). Varying \(B\) and \(0<s<1/2\) yields an open four-dimensional state family. For scalar rows \(L=0\).

## The reciprocal diagnostics across all candidate masses

At a **fixed** candidate, the reciprocal tangent decomposition already has a finite additive state. There is a different answer if the update must retain its within, between and common residual energies at **every** positive candidate mass. Take the full vector reading with independent spatially isotropic channels, common instrument ratio \(s\), and its two-channel Gram in noise units
\[
 B_i=\begin{pmatrix}A_i&C_i\\C_i&D_i\end{pmatrix}>0 .
\]
This full-rank class requires spatial dimension at least two. It permits arbitrary relative spatial frames and all energies \(A_i,D_i>0\) and correlations with \(C_i^2<A_i D_i\). Write \(t=m/s=\tan\theta,\ u=(\sin\theta,\cos\theta),\ v=(\cos\theta,-\sin\theta)\). The candidate-parallel energy is \(h_i=u^\top B_i u\); the score is \(s_i=v^\top B_i u\). The part of the residual orthogonal to the reading's own tangent is
\[
 W_i(\theta)=v^\top B_i v-\frac{s_i^2}{h_i}=\frac{\det B_i}{h_i}.
\]
These are the quantities in the [reciprocal tangent proof](e12-tangent-pairings.md). Here all \(h_i\) are positive. Define \(x_i=-C_i/A_i,\ y_i=\sqrt{\det B_i}/A_i>0\). The whole curve has the rational form
\[
 \frac{W(t)}{1+t^2}=\sum_i
 \frac{A_i y_i^2}{(t-x_i)^2+y_i^2}.
\]
Thus exact knowledge on an open positive mass interval determines every distinct shape \((x_i,y_i)\) and the total force energy \(\sum_{j:(x_j,y_j)=(x_i,y_i)} A_j\) at each shape. Indeed the rational function has a pole at \(x_i+i y_i\) with residue \(y_i\sum A_j/(2i)\); a common pole has contributions of the same sign in its imaginary part. Conversely the shape and its total energy give the group's scatter by multiplying that energy by
\[
 \begin{pmatrix}1&-x_i\\-x_i&x_i^2+y_i^2\end{pmatrix}.
\]
These statements concern **exact** functions. Pole recovery from a noisy or sparsely sampled curve is a separate numerical problem.

This yields a necessary and sufficient **predictive state for the three energy curves** on these full-rank blocks: group proportional Gram matrices and retain their sums. Necessity follows already for the empty extension from \(W\). Sufficiency follows because those groups determine \(W\), the total scatter \(B\), the common energy \(T=(v^\top B u)^2/(u^\top B u)\), and \(J=v^\top Bv-W-T\). Future full-rank readings add to the groups or create new groups. This comparison is of the energies. Quoting conditional chi-square ranks across varying experiment declarations may additionally require the reading counts and dimensions. The candidate-independent exterior state \((B,L)\) above does not generally preserve the candidate-dependent split.

For example, take two readings with Grams \(\operatorname{diag}(3,1),\operatorname{diag}(1,3)\). Compare the two Grams with diagonal entries \((2,2)\) and off-diagonal entries respectively \(+1,-1\). Both blocks have \(B=4I\), \(L=6\), the same reading dimensions, and the same instrument. At \(t=1\) the first gives \((W,J,T)=(3,1,0)\) and the second gives \((4,0,0)\). The pooled Gaussian residual and both candidate-independent exterior energies agree exactly. [Numerical check](../tests/test_e12_joining_state.py)

**Including all singular readings.** For a rank-one \(B_i\), its contribution to \(W\) vanishes wherever \(h_i>0\). If \(h_i=0\) at an interior candidate \(t_i>0\), the local tangent vanishes, so that reading contributes its entire residual \(\operatorname{tr}B_i\) to \(W(t_i)\). There is at most one such positive candidate per nonzero rank-one reading. The zero matrix contributes zero everywhere. The exact predictive state for the three energy curves on **all** common-ratio positive-semidefinite blocks therefore consists of the grouped sums of full-rank Grams, the pooled rank-one scatter, and the sum of traces at each positive singular candidate. For necessity, off those finitely many points the rational within curve identifies all full-rank groups. At each exceptional point the excess over its continuous limit identifies the corresponding trace sum. The total residual \(W+J+T=v^\top Bv\) on an open arc identifies the whole pooled scatter since its numerator is a quadratic in \(t\). Subtracting the full-rank groups gives the pooled rank-one scatter. For sufficiency all three parts update by addition. They give the within curve including its exceptional points and the total scatter, which supplies \(T\) with value zero when its denominator vanishes. This covers parallel, antiparallel, zero-channel and zero observations for positive finite candidate masses. In particular positive scalar rows alone need just their pooled scatter for these three energy curves. [Local zero-tangent convention](e12_tangent_pairings.py)

There is a precise cost in state size. For \(n\) readings with distinct \(x_i\), ordered by \(x_i\), the \(3n\) parameters \((A_i,x_i,y_i)\), with \(A_i,y_i>0\) and \(x_i\) unrestricted in sign, vary on an open set. These Grams are continuously realized using two spatial coordinates per reading, for example \(f_i=(\sqrt{A_i},0),\,a_i=\sqrt{A_i}(-x_i,y_i)\). On this family the within curve determines every parameter. Thus any continuous real state permitting the full geometric split to be evaluated across masses has dimension at least \(3n\), by invariance of domain as below. This lower bound follows from the reciprocal diagnostic itself, without any choice of mass law. It leaves open a mass law using a coarser part of the tensor information. [Finite checks](../tests/test_e12_joining_state.py)

## The current completed law hides information needed at the next reading

Here the output is specifically the [completed E12 law](e12-scalar-proof.md) for any fixed positive \(\alpha,b\), its median, and its score-hull report at any fixed confidence level. Fix \(u=(1,1)/\sqrt 2\). Two blocks with exactly two nonzero, strictly positive observations are
\[
 D=(u,2u),\qquad E=(\sqrt{5/2}\,u,\sqrt{5/2}\,u).
\]
Both give \(B=5uu^\top\), positive-arc minimum residual \(Q=0\), and the same undamped spectral precision after the fixed cap. The entire angular law agrees because the nonquadratic term is multiplied by \(Q/(1+Q)\). The median and every fixed-level score hull agree. In particular these blocks also give the same current determinant, absolute Gaussian cost and candidate-wise total tangent correction. Nevertheless appending a *single strictly positive row* on any other line (e.g. \((1,2)\)) gives different completed mass laws. This involves the same instrument, same prior and same count on each side. It concerns the returned law itself after appending, strengthening the [earlier hull-only example](e12-report.md). The distinction follows algebraically below.

**Theorem for the E12 completion.** Two finite scalar blocks are equivalent for its future completed mass laws iff their unordered multisets of nonzero rank-one matrices \(y_i y_i^\top\) coincide, with multiplicities. Independent row signs and measured zero padding are immaterial. The same characterization holds if the outputs also include the median and the score hull. This is a theorem about this witness on exact data and exactly specified functions. It places no general requirement on a different rule for returning mass.

**Proof.** The forward computation depends only on that multiset, also after a common extension. For the converse, append two nonparallel strictly positive rows on lines disjoint from all observed block lines in either experiment. They guarantee full-rank scatter and a strictly positive gate for both unions. One added row suffices if both old blocks have some nonzero row and the marker direction is disjoint from all old lines. Treat the marker rows and their lengths as known.

Write \(y_i=\sqrt{e_i}(\sin\psi_i,\cos\psi_i)\) up to sign, \(c_i=e_i/(1+e_i)\in(0,1)\), and use all nonzero rows of the union for the index \(i\). Up to an angle-independent constant the E12 log factor on this union has the form
\[
 \ell(\theta)= -\tfrac12 K'\sin^2(\theta-\hat\theta)
                -\lambda\sum_i\log[1+c_i\sin^2(\theta-\psi_i)],
 \qquad \lambda=
 \frac{e^{-H^2/L}Q}{(1+Q)(1+\sum_i c_i)}>0 .
\]
Here \(L,H\) are the top eigenvalue and its gap. The nonquadratic term carries information discarded when the gate vanishes. Two equal normalized laws with the same positive prior give identical log factors up to a constant on the open positive arc, by continuity. Put \(z=e^{2i\theta},\,q_i=e^{2i\psi_i}\). Each argument of the logarithm becomes the rational function
\[
 g_i(z)=1+c_i/2-\tfrac14 c_i(z/q_i+q_i/z).
\]
Its only root strictly inside the unit disk is \(z_i=q_i r_i\), where
\[
 r_i=\frac{\sqrt{1+c_i}-1}{\sqrt{1+c_i}+1}.
\]
The root is simple, nonzero and has phase \(q_i\). The other root \(q_i/r_i\) lies outside the disk. The map \(e_i\mapsto c_i\mapsto r_i\) is strictly increasing. A pole location \(z_i\) therefore identifies both the row energy and the real unoriented row line. The derivative \(g_i'/g_i\) has residue one at this root. The derivative in \(z\) of the quadratic spectral term has possible poles only at zero or infinity. Hence \(d\ell/dz\) is rational, with residue \(-\lambda\) times the row multiplicity at each distinct \(z_i\). Equality along the arc implies equality of the rational functions by the polynomial identity theorem. The known single occurrence of a marker yields \(-\lambda\) from its residue even when the two datasets have different total energies. Dividing the remaining residues by this one recovers all multiplicities. Removing the markers gives identical original multisets. This also proves the strictly positive example above. Mathematically positive gates may be extremely small, and the argument concerns identifiability, without a claim of stable numerical inversion. [Executable evidence formula](e12_scalar_completion.py)

**Continuous state lower bounds.** Any fixed-length real state that continuously encodes each \(n\)-row block and allows every future completed law to be evaluated needs at least \(2n\) coordinates in the unrestricted scalar experiment. In fact, restrict to positive rows with strictly increasing first components; this is a nonempty open subset of \(\mathbb R^{2n}\). There, the unordered outer-product multiset determines the ordered block. An injective continuous map from this open set into \(\mathbb R^d\) with \(d<2n\) contradicts invariance of domain after embedding \(\mathbb R^d\) as a coordinate hyperplane in \(\mathbb R^{2n}\). Even supplying the **present law, median and every fixed-level report in addition** still requires at least \(n-1\) continuous real coordinates of supplementary information on aligned blocks with fixed total energy. Their strictly ordered positive energies with fixed sum form an open simplex slice of that dimension, all giving identical current outputs. A future marker identifies the full allocation. No dimensional lower bound is asserted for discontinuous encodings. [Hatcher, *Algebraic Topology*, Theorem 2B.3](https://pi.math.cornell.edu/~hatcher/AT/AT-doublepage.pdf)

## One conditional restriction beyond this choice of law

Suppose a finite-dimensional real summary is additive over rows, acts by a continuous **linear representation** under rotations of the whitened force–acceleration channel plane, and the resulting data factor on the full projective mass line is equivariant under the same rotations. A positive continuous factor can be normalized against the uniform projective measure to remove multiplicative constants. Some balanced incompatible lists then give a constant factor. To see this, a compact rotation group has an invariant Euclidean inner product (average any positive inner product). The complexified commuting rotations diagonalize jointly, giving a finite set of integer Fourier frequencies; equivalently each character of the circle is \(e^{ik\phi}\). Choose \(N\ge3\) larger than all absolute nonzero frequencies. The sum of the summaries of \(N\) equal-energy rows at angles \(2\pi j/N\) has only its rotation-fixed part, by summing roots of unity. The scatter is isotropic and the rows admit no common mass. The factor must now be invariant under every rotation and hence constant. This argument imposes a rotation action **mixing the two channels** and does not follow merely from rotating physical space. Keeping the physical cone, prior or instrument axes as additional directional information can change the argument. [Channel and spatial symmetries](e12-random-reciprocal.md)

The three-coordinate result for the exterior norm and the unbounded state requirement for the E12 witness distinguish two concrete uses of the same data. The latter reveals an implementation cost of choosing a law with information in arbitrarily high angular harmonics. Equating that cost with a physical requirement, or selecting a different finite-state mass law by fiat, would go beyond either theorem. [Finite checks](../tests/test_e12_joining_state.py)

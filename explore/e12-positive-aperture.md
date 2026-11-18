# How instrument correlation changes the cost of a positive scalar inertia

**Result for G5.** Fix the full spatial dimension \(k\ge1\). With independent channels T38 gives \(k-\tfrac12\) as the supremum of the *mean* least residual compatible with a positive mass. On a separable instrument with nonnegative force–acceleration noise correlation, the same number is sharp, and a half-and-half chi-square mixture gives the sharp **entire upper tail**. Negative correlation can already invalidate that mean for one physical dimension. Over all positive-definite declared covariances the sharp worst-case mean is \(k\) and the sharp tail is \(\chi^2_k\). These are noise benchmarks for a scalar graph fitted to the observations. A probability law on mass needs an additional argument. [T38](../definition/theorems/T38.md) · [G5](../definition/proof.md)

## The covariance determines the aperture

Use the exact reciprocal construction \(F=m a\) with \(m>0\). Begin with channel covariance \(C=\left(\begin{smallmatrix}\sigma_F^2&\rho\sigma_F\sigma_a\\\rho\sigma_F\sigma_a&\sigma_a^2\end{smallmatrix}\right)\), \(|\rho|<1\), and common positive-definite spatial factor \(S\), with total covariance \(C\otimes S\). Spatial whitening acts identically on both vectors and preserves \(F=m a\). Channel whitening maps \((m,1)\) to \(C^{-1/2}(m,1)^\top\). The angle between the two *oriented* endpoint rays corresponding to zero and infinite mass has
\[
 \cos\beta=\frac{(C^{-1})_{12}}{\sqrt{(C^{-1})_{11}(C^{-1})_{22}}}=-\rho,
 \qquad 0<\beta<\pi .
\]
The positive coefficients select the arc of length \(\beta\) between these rays, viewed on the projective circle of circumference \(\pi\). Taking the shorter projective distance between endpoints would incorrectly replace an obtuse aperture by its complement. A spatial rotation has no effect on the aperture. The dual channel metric is what measures it.

In these noise coordinates the \(2\times k\) observation is \(Y=u_{\theta_0}r^\top+Z\), with \(\theta_0\) in the aperture, a deterministic \(r\in\mathbb R^k\), and independent standard Gaussian entries of \(Z\). Everything here also describes a series of independent scalar or isotropic spatial readings with the same whitened channel arc and parametrization, by adjoining their columns; \(k\) then counts all independent spatial components. Write \(B=YY^\top\), \(l\le L\) for its eigenvalues, and \(\psi\) for its principal line. If \(A\) is the closed mass arc and \(d_\pi\) its projective angular distance, the infimum of the squared residual over positive mass is
\[
 Q_\beta=l+(L-l)\sin^2 d_\pi(\psi,A).
\]
Using the closure leaves the infimum unchanged. Arbitrary SPD full covariances are covered by the universal bound at the end, while the aperture formula uses separability. [Reading-space residual](../definition/theorems/T36.md)

The mean threshold actually extends past the right angle. A reflection argument gives \(k-\tfrac12\) for all \(\rho\ge-1/2\). For scalar readings this mean fails for every \(\rho<-1/2\), while the upper-tail threshold is \(\rho=0\). [Proof of the mean threshold](e12-aperture-mean.md)

## A sharp comparison when the aperture is at least a right angle

For a right-angle aperture, the following tail version uses T38's eigenline proof. At fixed \(\|r\|\), the eigenvalues' distribution is independent of the orientation of the true line. Conditional on them, the principal line has a density even and decreasing in its projective distance from the true line. Indeed, starting with central Gaussian noise and averaging its mean-shift likelihood ratio over right orthogonal transformations multiplies the uniform angle density by a factor proportional to
\[
 E\exp\{\|r\|\,U_1\sqrt{\,l+(L-l)\cos^2(\psi-\theta_0)}\},
\]
where \(U_1\) is a coordinate of a uniform unit vector in \(\mathbb R^k\). This is increasing in the square root by sign symmetry (also for \(k=1\)). For any threshold \(c\), the rejection set \(Q_\beta>c\) is either empty, the full circle, or the central part of the complementary arc. The overlap of a symmetric decreasing angle density with this latter interval is largest when the true line is closest to its centre. This follows by integrating the overlaps with concentric superlevel arcs of the density. Among orientations in \(A\), the two endpoints are therefore least favourable, for **any** aperture and at each threshold. The same reduction gives the largest mean at an endpoint for fixed strength. [T38's mean argument](../problems/O13.md) · [Noncentral Wishart reduction](../definition/theorems/T69.md)

For a right angle and a boundary truth, write the measured column-vectors along the true channel and its normal as \(X=r+\xi\) and \(Y_\perp=\zeta\), with independent standard Gaussian errors. Keeping the first observed vector fixed, fit the other along its nonnegative ray. This is a subset of the allowed graph fits. With \(e=X/\|X\|\), the squared residual of this construction is
\[
 D=\|\zeta-e(e^\top\zeta)_+\|^2
   =\|(I-ee^\top)\zeta\|^2+(e^\top\zeta)_-^2 .
\]
Here \(v_+=\max(v,0)\) and \(v_-=\max(-v,0)\). The true minimizing distance is at most \(D\). Conditional on \(X\), the two projected Gaussian components are independent, and the sign of \(e^\top\zeta\) is fair and independent of its magnitude. Hence \(D\) has the mixture law
\[
 B_k=\tfrac12\,\mathcal L(\chi^2_{k-1})+\tfrac12\,\mathcal L(\chi^2_k),
 \qquad \chi^2_0=0 .
\]
For \(\beta\ge\pi/2\), at each true line there is a right-angle subarc of \(A\) containing that line. Fitting on the larger arc can only lower the residual. This proves \(P(Q_\beta>c)\le P(B_k>c)\), writing \(B_k\) also for a draw from the mixture.

For every fixed \(0<\beta<\pi\), the strong-signal limit at an endpoint is this same mixture. In a frame with the endpoint at angle zero take \(X=n e_1+\xi,\;Y_\perp=\zeta\). The principal line has signed angle \(\zeta_1/n+o(1/n)\) almost surely in this coupling. The least unrestricted residual tends to \(\|\zeta\|^2-\zeta_1^2\) by the two-by-two eigenvalue formula. If \(\zeta_1>0\) the principal line eventually lies in \(A\); if \(\zeta_1<0\) the nearest admissible line is the zero-angle endpoint, giving residual \(\|\zeta\|^2\). Taking interior true angles \(o(1/n)\) preserves the limit and keeps all true masses positive and finite. A candidate at truth gives a \(\chi^2_k\) upper bound in law, hence uniform integrability. For \(c>0\) this proves sharpness of the survival bound when \(\beta\ge\pi/2\), and taking expectations gives the sharp mean \(k-\tfrac12\). At \(k=1,c=0\) the sign of the fitted line gives the mixture's atom directly. This is a finite-experiment envelope with a limiting least-favourable sequence. [Boundary likelihood-ratio precedent](https://doi.org/10.1214/aoms/1177728725)

## An acute aperture can cost more than half an extra dimension

At zero true signal the principal line is uniform on the whole projective circle, independent of its eigenvalues. The complement has length \(\pi-\beta\). Integration from its ends gives
\[
 E_0 Q_\beta
 =k-\frac{\beta+\sin\beta}{2\pi}\,E_0(L-l).
\]
Indeed \(E_0\operatorname{tr}B=2k\) and the angle average of \(\sin^2 d_\pi(\psi,A)\) is \((\pi-\beta-\sin\beta)/(2\pi)\). For \(k=1\) the expected gap is 2. With \(\rho=-\sqrt3/2\) the aperture is \(\pi/6\), and
\[
 E_0Q_\beta=1-\frac{\pi/6+1/2}{\pi}
            =\frac56-\frac1{2\pi}\simeq0.67418>\tfrac12 .
\]
The truth can have any positive mass on this zero-push experiment. By continuity of the Gaussian shift and uniform integrability, nonzero sufficiently weak pushes also exceed one half. This disproves a covariance-independent reading of T38's constant. Its original independent-channel theorem is unaffected.

There is a stronger tail counterexample for **every** \(\rho<0\) in dimension one. With standardized observations \((X,Y)\), the positive graphs fill the same-sign quadrants and their axes. On opposite signs the minimum squared Mahalanobis distance is \(\min(X^2,Y^2)\). Each axis belongs to the closure and its squared distance is the square of its standardized normal coordinate. Any path into the same-sign region crosses an axis. At the boundary mean \((R,0)\), for \(u>0\),
\[
 P(Q_\beta>u^2)-\Phi(-u)
 =P(X<-u,Y>u)-P(X\le u,Y<-u).
\]
The first probability is asymptotic to \(\Phi(-R-u)\). Indeed \(Y\mid X=x\) has mean \(\rho(x-R)\) and variance \(1-\rho^2\), so \(P(Y>u\mid X=x)\to1\) uniformly over \(x<-u\). The subtracted probability is bounded above by
\[
 \Phi(u-R)\,\Phi\!\left(\frac{-u-\rho(u-R)}{\sqrt{1-\rho^2}}\right).
\]
For fixed negative \(\rho\), the second Gaussian tail here decays as \(\exp[-c_\rho R^2+O(R)]\), which dominates the \(\exp[O(R)]\) ratio of the first factor to \(\Phi(-R-u)\). Thus the difference is positive for sufficiently large finite \(R\). By continuity it stays positive for a nearby strictly positive true mean \((R,b)\). Since \(P(B_1>u^2)=\Phi(-u)\), the right-angle tail envelope fails for every acute aperture even though the strong-signal boundary **limit** is the same. This tail argument does not determine the sharp mean for each negative correlation separately.

## What can be forced for an arbitrary instrument?

For any SPD instrument and any fixed true positive mass, the profiled residual \(Q_+\) is no larger than the squared Gaussian residual at truth. On T36's declared span the latter has the \(\chi^2_k\) distribution. This follows by whitening and projecting onto the \(k\)-dimensional complement of the true graph. The same argument on a direct sum uses the total rank. Conversely, take separable instruments with \(\rho\downarrow-1\) and zero true push. In whitened coordinates couple all measurements by one standard Gaussian matrix and rotate the aperture to have a fixed endpoint. As \(\beta\downarrow0\) the minimum tends almost surely to the residual at that endpoint, which has law \(\chi^2_k\). The endpoint residual bounds the coupled minima, so the means approach \(k\). Thus the upper bound is sharp for means and for every positive tail threshold when the supremum includes the declared covariance. No singular covariance is used in any finite experiment.

For a *fixed* acute separable instrument we have also reduced the sharp mean and each sharp tail to a supremum over one scalar, the signal strength at one endpoint. An exact representation uses mutually independent \(A\sim\chi'^2_k(\|r\|^2)\), \(T\sim N(0,1)\), \(W\sim\chi^2_{k-1}\). Evaluate the aperture minimum with scatter \(\left(\begin{smallmatrix} A&\sqrt{A}T\\\sqrt{A}T&T^2+W\end{smallmatrix}\right)\). This follows by decomposing the normal channel into its projection onto \(X\) and its orthogonal complement. The fixed-instrument optimization for arbitrary nonseparable covariances remains open. Choosing a calibrated tail probability, an expected cost, or an evidence weight are distinct operations. A dimension-based noise benchmark also changes under padding the experiment with an extra declared dimension; automatically using it to weight the returned law would owe the project's exact-zero test. [T36](../definition/theorems/T36.md) · [Checks](../tests/test_e12_aperture.py)

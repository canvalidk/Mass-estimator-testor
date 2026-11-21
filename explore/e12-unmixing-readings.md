# Can the original instrument recover a superposition from its readings?

A single force–acceleration reading with an arbitrary push has no declared sampling law on mass alone. There is nevertheless an identification result on the **same sensor channels**. Suppose an entire repeatable preparation is specified by its predictive observation law. Its mass and drive may have any dependence; the calculation requires no model for that dependence. The new input is the full predictive law, in contrast to one finite record of unrelated unknown pushes. [Single-record obstruction](e12-conditioning-nuisance.md)

## The reciprocal identifies everything except a zero push

Let \(m=(m_1,\ldots,m_k)\) be a vector of positive labelled masses and \(a\in\mathbb R^d\) their shared true spatial acceleration. An exact prepared signal is
\[
             z=(m_1a,\ldots,m_ka,a).
\]
Let \(\Lambda\) denote the unknown distribution of \(z\) in the preparation and let \(P\) be the observable law after adding independent centered Gaussian error of fixed, known positive definite covariance \(\Sigma\) to these channels. The signal can include an atom at zero. Each individual Gaussian error vector may have correlations across all force and acceleration channels.

**Theorem.** The full law \(P\) determines the joint mass-and-drive law on \(a\ne0\) and the probability \(q\) of zero drive. If no other mass information is provided, all compatible mass laws are exactly
\[
                   \mu_{\rm active}+q\eta ,
\]
where \(\mu_{\rm active}\) is an identified subprobability of total weight \(1-q\) and \(\eta\) ranges over probability laws on positive labelled mass vectors.

**Proof.** The multivariate characteristic functions obey
\[
        \widehat\Lambda(t)=\widehat P(t)\exp(t^\top\Sigma t/2).
\]
The Gaussian transform is everywhere nonzero. Characteristic functions identify arbitrary probability laws, so the given \(P\) fixes \(\Lambda\), without moment assumptions. Away from \(a=0\), pair every force vector with the **same reciprocal covector**
\[
 a^\vee(f)=\langle a,f\rangle/\|a\|^2,\qquad
 M(z)=\big(a^\vee(f_1),\ldots,a^\vee(f_k)\big).
\]
This inverts the exact-signal map, independently of the spatial coordinate frame. Push the restriction of \(\Lambda\) forward by \((M,a)\) to recover the joint on all active pushes. Every zero drive maps to the single signal \(z=0\), giving \(q=\Lambda(\{0\})\). Attaching any mass law \(\eta\) to that atom leaves the same exact-signal and observation laws. This proves sharpness. If the acquisition sums the labelled forces before recording them, the active acceleration and total mass are recovered. Their allocation between bodies is then unobserved. [Gaussian deconvolution](https://arxiv.org/abs/1907.11024)

**The law after seeing a noisy reading is also identified up to precisely this ambiguity.** For a particular observation \(y\) let \(\phi_\Sigma\) be the Gaussian noise density and \(p(y)=\int\phi_\Sigma(y-z)\Lambda(dz)>0\). Define the subprobability measure
\[
 A_y(B)=\frac{\int_{z\ne0}{\bf1}_{\{M(z)\in B\}}\phi_\Sigma(y-z)\,\Lambda(dz)}
                   {p(y)},\qquad
 q_y=\frac{q\phi_\Sigma(y)}{p(y)}.
\]
In the canonical Bayes version given by the Gaussian density, exactly the posterior laws \(A_y+q_y\eta\) are compatible with \(P\). One and the same \(\eta\) applies to every possible outcome. This follows by conditioning each component against its Gaussian response; the undriven component has the same density for every mass. If pushes are almost surely nonzero, \(P\) therefore determines the complete mass posterior at each reading, even though the drive itself was unrecorded. If the marginal law \(\pi\) of masses in this **same preparation** is separately known, it fixes even the missing mass assignment: compatibility is precisely \(\pi\ge\mu_{\rm active}\) as measures, and for \(q>0\), \(q\eta=\pi-\mu_{\rm active}\). When \(q=0\) the compatibility condition requires equality. A declared starting uncertainty law need not describe the mass frequencies of a repeatable preparation. Equating them is an experimental premise. [Measurement framework](e12-mass-measurement.md) · [Exact finite-mixture checks](../tests/test_e12_unmixing_readings.py)

**A formula using the declared preparation prior directly.** If the supplied \(\pi\) is compatible, write the same answer as
\[
 p(y)\Pi_y(dm)=\phi_\Sigma(y)\,\pi(dm)+
 \int_{z\ne0}[\phi_\Sigma(y-z)-\phi_\Sigma(y)]\,\delta_{M(z)}(dm)\,\Lambda(dz).
\]
This follows by adding and subtracting the zero-signal likelihood on the active component. Tested against any bounded continuous mass function, the integral's integrand extends continuously to zero by setting it to zero, because its likelihood difference vanishes there. Hence if compatible preparation marginals converge weakly to \(\pi\) and their predictive observation laws converge weakly to \(P\), the corresponding posteriors converge weakly at each fixed \(y\), **even when \(q>0\)**. Indeed deconvolution gives weak convergence of the exact-signal laws as shown below, and the displayed integrands are bounded and continuous at every physical limiting signal. This bypasses the possibly discontinuous operation of deciding how much observed probability sits at exactly zero.

There is a second direct expression. If the predictive density is available near \(y\), its nearby values determine the conditional moment-generating function of the *exact signal* \(Z\):
\[
 E[e^{t^\top Z}\mid Y=y]
 =e^{t^\top y+t^\top\Sigma t/2}\frac{p(y+\Sigma t)}{p(y)} .
\]
Completing the square inside the Gaussian mixture integral proves the identity and finiteness for every real \(t\). A neighbourhood of 0 determines this conditional law by uniqueness of moment-generating functions. Apply the same reciprocal on its active part. Thus deconvolution is one proof of identification; the local predictive density also encodes the full conditional signal law. The formula uses population density values, which a single reading cannot supply.

## Infinite repeats and small drives

Under independent repeats from one such preparation with \(q=0\), there exists a **pointwise weakly consistent** estimator of the mass law and, for any fixed fresh reading \(y\), of its conditional law. For an existence proof take a countable family of finite rational mixtures of signals at positive rational masses and rational nonzero acceleration vectors. It is weakly dense among the allowed exact-signal laws, by approximation of the latent law and continuity of the signal map. At sample size \(n\), choose the first candidate whose Gaussian convolution has distance from the empirical observation law within \(1/n\) of the countable infimum, in any metric for weak convergence. This choice is measurable, for example with the bounded-Lipschitz metric. Empirical weak convergence and the triangle inequality give convergence of the selected convolved laws to \(P\) almost surely. Dividing characteristic functions and using the continuity theorem gives weak convergence of their signal laws to \(\Lambda\). The reciprocal inverse is continuous on the signal space away from \(a=0\), which has full limiting probability. Pushing forward gives the asserted convergence of mass laws. Multiplying first by the bounded continuous function \(\phi_\Sigma(y-z)\) and dividing by its positive limiting mean gives the corresponding conditional result. This construction establishes existence without claiming a practical rate or an estimator for a finite list of arbitrary preparations.

The weak-drive limit forbids **uniform** recovery. Choose two distinct fixed mass vectors and the same nonzero acceleration \(\epsilon a_0\). Their mass laws are point masses at different locations, while their observational Gaussians have total variation \(2\Phi(d_\epsilon/2)-1\), where \(d_\epsilon\) is their whitened mean separation, proportional to \(\epsilon\). This expression follows by thresholding the scalar whitened projection at the midpoint of the means. With \(n\) independent repeats the distance is \(2\Phi(\sqrt n\,d_\epsilon/2)-1\). Drives with \(\epsilon=o(n^{-1/2})\) keep these experiments indistinguishable even as \(n\) grows. Thus the positive result states what a fixed repeatable superposition encodes. The original unrestricted task gives a single series, a declared starting law and error covariances. It promises no population from which to learn the predictive observation law. The repeatability assumption is additional to L10 and L19. [Reading premises](../definition/premises.md) · [Information limits on arbitrary reading lists](../definition/theorems/T41.md)

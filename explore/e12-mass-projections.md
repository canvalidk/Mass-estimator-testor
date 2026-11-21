# Which measurements determine a joint law of two positive masses?

Write \(X,Y>0\) for the two labelled masses with an arbitrary joint law. In each experiment a declared scalar observable \(g\) gives the law of \(g(X,Y)\). Across experiments these may be supplied as **separate marginal laws**. Two measurements made on the *same realization*, with their pairing retained, provide the joint law of two observables and can carry much more information. The distinction matters even with exact instruments. Here the masses and standards are expressed numerically in fixed units.

## Paired measurements select a point; separate laws can leave an ambiguity

A total \(T=X+Y\) and fraction \(r=X/T\) taken on the same pair give \((X,Y)=(rT,(1-r)T)\). A total alone fixes a line of possibilities. A fraction alone fixes a ray. Across a population or an uncertainty law, the two *separate* distributions leave their dependence open. For example compare equal mixtures of \((T,r)=(4,1/4),(8,3/4)\) and \((4,3/4),(8,1/4)\). Both have the same total and fraction distributions. At total 4 one assigns force-bearing mass \(X=1\), the other \(X=3\).

There is a similar pair of strictly physical, labelled composite measurements. Attach a known standard \(b\) to the first body and form the reciprocal of its reduced mass with the second body:
\[
 W_b=\frac1{X+b}+\frac1Y .
\]
For \(0\le b<c\), two exact readings on the same realization give \(W_b-W_c=(c-b)/[(X+b)(X+c)]>0\). Thus
\[
 X=-\tfrac12(b+c)+\sqrt{(b-c)^2/4+(c-b)/(W_b-W_c)},\qquad
 Y=[W_b-(X+b)^{-1}]^{-1}.
\]
A reading of \(W_b\) and an independently listed distribution of \(W_c\) do not give this difference on individual pairs. These constructions use known mass additions and the common-force composition of the corresponding positive masses, under the declared protocols. They concern observable laws that are available to the observer; the existing premises do not automatically deliver all of them. [Known-standard composition](e12-known-mass.md) · [Earlier joint-law result](e12-joining-bodies.md)

## A finite list of separate smooth measurements always has a blind direction

**Theorem.** Choose any finite list of smooth real observables \(g_1,\ldots,g_J\) on the open positive mass quadrant. There exist two distinct smooth, strictly positive, full-support joint densities that have exactly the same pushforward probability for each of the \(J\) observables. The two may both be exchangeable, including when the selected instruments distinguish body labels. Smoothness and the continuum of admissible mass pairs matter; a known finite support or parametric model changes this claim.

**Proof.** In a small open ball inside the quadrant look for a differential operator \(P=\sum_{|\alpha|\le m}c_\alpha(x,y)\partial^\alpha\). Require \(P(h\circ g_j)=0\) for every smooth univariate \(h\). The chain rule expands this into at most \(m+1\) homogeneous linear conditions per observable on the coefficients \(c_\alpha\), one for each \(h^{(k)}\circ g_j,\ 0\le k\le m\). There are \((m+1)(m+2)/2\) coefficient functions. Choose \(m=2J\), so there are more unknown coefficients than equations. The coefficient matrix is smooth. Near a point of maximal rank in the ball it has constant rank. Solving an invertible rank-size minor after setting one free coefficient to one gives a smooth nonzero null vector on that patch. The zero-rank case permits any free coefficient. Multiplication by a compactly supported cutoff leaves all the conditions intact.

This gives a nonzero local operator \(P\). Some smooth compactly supported \(\psi\) has nonzero formal adjoint image \(q=P^*\psi\). Indeed one can prescribe the derivatives of a test function \(v\) at a point where a coefficient is nonzero to make \(Pv\ne0\); then choose \(\psi\) supported where \(Pv\) has constant nonzero sign. Integration by parts gives \(\int h(g_j)q=0\) for all \(j,h\). In particular \(\int q=0\), and smooth scalar test functions determine the pushforward measures. Adding and subtracting a small multiple of this smooth compactly supported signed density from a strictly positive smooth base yields the claimed laws. For exchangeability choose the ball disjoint from its coordinate-swapped image, include every \(g_j\) and its swapped version before constructing \(q\), and add the swapped copy of \(q\). Its disjoint supports prevent cancellation. Use an exchangeable base.

Consequently even access to the full marginal output distributions for the parts, total, reduced mass, ratio and any *finite* set of additive or reciprocal standards cannot identify every positive joint law. The same ambiguity persists if a noisy scalar instrument at setting \(j\) has a known response kernel depending on the pair only through \(g_j\). The identical pushforwards enter the same kernel integral. A joint force–acceleration record whose distribution depends on other components can convey further information. For an exact finite example with equal part, total and reduced laws, see the [exchangeable rational construction](e12-joining-bodies.md#what-if-both-physical-composition-laws-were-known-globally).

## One adjustable labelled standard suffices

**Theorem.** Suppose the observer has the separate laws of \(W_b\) for every \(b\) in a set with an accumulation point strictly inside \((0,\infty)\), or just for \(b=1,2,\ldots\) in a chosen mass unit. These laws uniquely determine the joint law of \((X,Y)\). Neither part law nor any finite moment of the masses is required. Knowing the reduced-mass laws themselves gives the same information by reciprocation.

**Proof.** Fix \(t>0\). The observations supply
\[
 H_t(b)=E\exp[-t/Y-t/(X+b)].
\]
For complex \(b\) in the open right half-plane this is bounded and holomorphic: the integrand is analytic and has magnitude at most one, so integration on compact triangles and Morera's theorem apply. In the accumulation-point case, equality for two candidate laws extends to all real \(b>0\). In both cases the limit \(b\to\infty\) along the available settings gives \(H_t(\infty)=E e^{-t/Y}\).

Put \(k_t(u)=\sum_{n\ge1}(-1)^{n+1}t^n u^{n-1}/[n!(n-1)!]\). This analytic kernel has \(k_t(0)=t\). Expanding the exponential and integrating powers of \(u\) gives, for every real \(z>0\),
\[
 1-e^{-t/z}=\int_0^\infty e^{-zu}k_t(u)\,du .
\]
The integral of the sum of absolute terms is \(e^{t/z}-1\). Consequently
\[
 H_t(\infty)-H_t(b)
   =\int_0^\infty e^{-bu} k_t(u)\,L(u,t)\,du,\qquad
 L(u,t)=E[e^{-uX-t/Y}].
\]
For two candidate laws with the same data, take the difference of this identity at \(b=1,2,\ldots\). The finite signed measure with density \(e^{-u}k_t(u)[L_1(u,t)-L_2(u,t)]\), pushed through \(u\mapsto e^{-u}\), has every nonnegative integer moment zero. Polynomial approximation on \([0,1]\) makes this signed measure zero. The functions are continuous and the zeros of \(k_t\) are isolated, hence \(L_1(u,t)=L_2(u,t)\) for every \(u>0\). This works for each \(t>0\). The resulting bivariate Laplace transforms determine the law of \((X,1/Y)\): the mixed integer values give all joint polynomial moments of \((e^{-X},e^{-1/Y})\) in \([0,1]^2\), including axis moments by bounded convergence. Inverting both coordinates recovers \((X,Y)\). The proof establishes uniqueness with ideal predictive laws. It supplies no stable reconstruction from finitely many noisy histograms. [Lin and Dou, Laplace–Stieltjes uniqueness from countable values](https://doi.org/10.1080/02331888.2021.1893728)

There is a simple check on the labelling requirement. On three equally likely pairs \((1,2),(2,3),(3,1)\), compare their transposes. All individual marginal, total and ordinary reduced-mass distributions agree. Adding one unit to the *first* mass before reciprocal composition gives \(\Pr(W_1=1)=1/3\) for the first law and zero for its transpose. A standard added only to the total, or indiscriminately averaged over the body labels, cannot reveal the same orientation. [Exact and numerical checks](../tests/test_e12_mass_projections.py)

The settings can also be randomized within one declared protocol. Choose \(B\) independently of the underlying masses from a known continuous distribution with positive density on a bounded interval inside \((0,\infty)\), and retain the **setting together with the reading**. Equality of joint predictive laws for \((B,W_B)\) entails equality of \(H_t(b)\) for almost every \(b\) in that interval for each fixed \(t\), so the same analytic argument applies. This is a two-coordinate observed record. One realized pair of numbers does not determine an arbitrary latent probability law.

## What precision would identifying all the settings require?

The one-standard sweep has no inverse that is uniformly continuous in total variation on smooth full-support joint densities, even when its input discrepancy is the supremum of total variation over **every** \(b\ge0\). To see this, work in \(x=X,\ z=1/Y\). Choose a nonnegative smooth nonzero bump \(\phi\) supported in \((1,2)^2\), and take the signed density
\[
 q_N(x,z)=\partial_z\{\phi(x,z)\sin[N(x+z)]/N\}.
\]
It has integral zero, is uniformly bounded, and its \(L^1\) norm stays bounded below as the integer \(N\) grows. Indeed the leading term is \(\phi\cos[N(x+z)]\); its absolute integral tends to \((2/\pi)\int\phi\) by periodic averaging, and the derivative-of-bump term has norm \(O(1/N)\). Add and subtract the same small constant multiple of \(q_N\) from a smooth strictly positive density. This leaves a fixed positive lower bound on the total variation between the two joint laws. Returning to the \(Y\) coordinate preserves that distance.

The density difference at \(W_b=w\) is proportional to \(\int q_N(x,w-1/(x+b))\,dx\) with unit Jacobian. The phase inside the integral is \(N[x+w-1/(x+b)]\). Its derivative divided by \(N\) is \(1+(x+b)^{-2}\ge1\) and all higher derivatives are uniformly bounded on the support for \(b\ge0\). The same bounds hold there for all fixed-order derivatives of the amplitudes. Integration by parts in \(x\), with no endpoint terms, gives \(O_k(N^{-k})\) for every fixed integer \(k\), uniformly in \(b,w\). These densities vanish outside \(1\le w\le3\). Thus the supremum of the output total variations tends to zero while the joint total variation stays positive. This example concerns uniform robustness over a large class; it is compatible with exact identification. Further independent instrument noise can only reduce the total variation of each output.

## Comparison with linear measurements

For weighted lumps \(X+sY\), positive slopes \(s\) accumulating anywhere inside \((0,\infty)\) also determine the joint law. For each \(u>0\), \(E e^{-uX-usY}\) is bounded and holomorphic for \(\operatorname{Re}s>0\). Agreement on the given slopes gives equality on the entire half-line by the identity theorem. The bivariate Laplace transform then identifies the joint law by the same compact moment argument. Even slopes \(s=n,\ n\ge1\), suffice by themselves: for fixed \(u>0\), these are the positive integer moments of the finite measure of \(e^{-uY}\) weighted by \(e^{-uX}\). That measure has no atom at zero, so moments of strictly positive order determine it (multiply its difference by the coordinate before using polynomial approximation). This strengthens the earlier all-integer-weights argument. Here the coefficients act on the **same two masses**. Physically replacing an integer coefficient by fresh independently uncertain copies gives a different transform. Finite smooth observables still have the ambiguity proved above. For compactly supported joint laws, every infinite set of distinct planar slopes suffices by polynomial interpolation of bounded moments. No unbounded-support claim about every possible sparse sequence is required here. [Cuesta-Albertos, Fraiman and Ransford, *A sharp form of the Cramér–Wold theorem*](https://personales.unican.es/cuestaj/CramerWold.pdf)

For a calibrated Gaussian instrument observing \(g(X,Y)\) plus independent noise of known variance, its predictive marginal identifies that same scalar pushforward by dividing characteristic functions by the nonvanishing Gaussian transform. This observation also applies to the reciprocal-standard experiment if that is what the instrument measures. Exact composite-law predictions here are premises to be supplied or derived. Merely forming them by pushing forward an already chosen joint law makes each choice consistent with itself and selects no unique update. Repeated observations of one fixed unknown pair and draws from a population with a specified joint law have different sampling interpretations.

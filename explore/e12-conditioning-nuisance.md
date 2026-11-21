# When noisy force–acceleration readings can act by ordinary conditioning

For an exact set observation, the square root of a density is multiplied by the indicator of the observed set and normalized in \(L^2\). A noisy instrument with declared sampling density \(k(y\mid m)\) instead multiplies it by \(\sqrt{k(y\mid m)}\) and normalizes. This follows by disintegrating the joint measure \(\pi(dm)k(y\mid m)dy\). For conditionally independent outcomes the multipliers compose by multiplication. On hypotheses with positive prior odds, the posterior odds equal the prior odds times the ratio of those sampling densities. One may test this implication on a two-point prior, or by continuity on the interior of a continuous full-support mass prior. Agreement of event updates alone supplies no probabilities for the instrument outcomes. Here we ask when the project's measured pair really supplies such a kernel on mass.

## A raw outcome with an unrestricted push

On a common unit-noise scalar instrument write \(Y\sim N(r u_\theta,I_2)\), where \(\theta=\arctan m,\ u_\theta=(\sin\theta,\cos\theta)\), and \(r\) is an arbitrary fixed signed drive amplitude in graph arclength units. Consider any candidate-independent data reduction or randomized instrument \(B(dz\mid y)\). For it to supply an exact **mass-only sampling kernel** on this whole experiment, the distribution of \(Z\) must be the same at every allowed drive for a fixed mass. If zero is allowed, all mass lines contain that same zero mean. So this kernel must be identical at every mass. If zero is excluded but can be approached, the Gaussian laws converge in total variation to the common zero-mean law, and data reduction can only decrease total variation distance. The same argument works for vector observations with any fixed nonsingular covariance and for independent finite lists permitting jointly vanishing pushes.

The scalar conclusion still holds if the allowed drives at each mass contain an interval bounded away from zero, or any set with a finite accumulation point. For each output event let \(f(y)=B(\text{event}\mid y)\in[0,1]\). Its expectation against \(N(r u_\theta,I)\) is real analytic in \(r\): relative to the zero-mean Gaussian density the integrand has multiplier \(\exp(r u_\theta^\top y-r^2/2)\), uniformly dominated on compact complex \(r\)-sets by an integrable function. Constancy on the allowed set forces constancy for all \(r\), including zero. The pivot \(K_\theta(Y)=\|(I-u_\theta u_\theta^\top)Y\|^2\) has a nuisance-free noise law at *each true candidate*. Its definition uses the candidate, whereas the event \(Z\) above must be the same measured outcome for all candidates. Inverting a pivot and specifying a joint outcome kernel are different mathematical operations. [T36](../definition/theorems/T36.md) · [T45](../definition/theorems/T45.md)

A controlled drive \(r\ne0\) supplies a sampling kernel directly. If only its noise-unit magnitude \(r_0>0\) is known and either sign is possible, \(YY^\top\) has the same law at both signs and expectation \(I+r_0^2u_\theta u_\theta^\top\), which distinguishes different mass lines. This illustrates exactly where the accumulation assumption enters. In an apparatus that controls *physical acceleration*, it is the known graph mean \((ma,a)\) that enters the Gaussian likelihood. A declared proper law of the unknown drive also defines a marginal kernel by integration. L10 and L19 do not supply that law in the unrestricted reading type. [Premises](../definition/premises.md)

## A four-observation test of Bayesian representability

First consider a weaker way to ask for an ordinary posterior. Draw the true mass \(\Theta\) from the fixed starting law \(\pi\) and hold the push \(r\) fixed at an unknown value. Suppose the returned law \(P_Y\) is required to be **reliable conditional on what it reports** under every such fixed drive. Precisely, for each mass set \(B\),
\[
 E_{\pi,r}[{\bf1}_{\{\Theta\in B\}}\mid P_Y]=P_Y(B)\quad\text{almost surely}.
\]
A true posterior has this property in its own joint sampling model by iterated expectation. At \(r=0\), \(Y\) and \(\Theta\) are independent, so the left side equals \(\pi(B)\) almost surely. Applying this to a countable determining family of intervals proves \(P_Y=\pi\) almost surely under the full-support zero-mean Gaussian. Every nonzero-drive observation law is absolutely continuous with respect to that Gaussian, so the same readout must equal \(\pi\) almost surely at *every* mass and drive. If drives approaching zero are included, apply the reliability identity multiplied by bounded functions of the report and pass to the limit in total variation. Thus a nontrivial mass law cannot be a reliable ordinary posterior in every one of these nuisance worlds. This impossibility uses calibration of the law conditional on its report; finite coverage of a quoted interval and local asymptotic coverage ask for different things. If there is a prior distribution on drives, posterior reliability is averaged over that distribution, so the zero-drive step need not hold separately. [T27 and its change of scope under L10](../definition/theorems/T27.md)

Allow even an arbitrary **mass-specific** proper probability law \(\nu_\theta\) for the latent drive, declared before the observed pair. Write \(\phi_2\) for the standard normal density. The raw data density becomes
\[
 p_\theta(y)=\phi_2(y)\,M_\theta(u_\theta^\top y),\qquad
 M_\theta(s)=\int e^{rs-r^2/2}\,\nu_\theta(dr).
\]
This follows by expanding the squared Gaussian distance. The transform is finite and analytic for every real \(s\), including for heavy-tailed drives. For two distinct interior masses let \(p=u_1^\top y,\ q=u_2^\top y\); these are independent linear coordinates on the **data plane** (their noises need not be statistically independent). Their log likelihood ratio is \(R(p,q)=\log M_1(p)-\log M_2(q)\). Therefore
\[
 R(p_1,q_1)-R(p_1,q_0)-R(p_0,q_1)+R(p_0,q_0)=0 .
\]
The test concerns odds, so multiplying all candidates' factors by an arbitrary positive data-dependent function or using different fixed prior odds leaves it unchanged. This is necessary for ordinary posterior conditioning under *any* such proper Gaussian push mixture. It does not require a population shared between readings.

The [earlier E12 laws](e12-joining-residuals.md) agree on each strictly positive singleton. Their log odds between candidate masses 1 and 2, relative to starting odds, are
\[
 R=\frac{d}{2}\left[\frac{(F+A)^2}{2}-\frac{(2F+A)^2}{5}\right],\qquad
 d=\min\{1,(1+\min(F^2,A^2)^2)/(F^2+A^2)\}.
\]
Use the unnormalized projection coordinates \(P=F+A,\ Q=2F+A\), which just rescale the two axes above. Four physical inputs give the following exact comparison.

| \(P\) | \(Q\) | \((F,A)\) | E12 log odds |
| ---: | ---: | --- | ---: |
| 4 | 6 | \((2,2)\) | \(2/5\) |
| 4 | 7 | \((3,1)\) | \(-9/50\) |
| 5 | 6 | \((1,4)\) | \(53/170\) |
| 5 | 7 | \((2,3)\) | \(27/20\) |

The mixed difference is \(2751/1700\). Hence neither E12 witness can equal ordinary conditioning of the full observed pair under the specified Gaussian experiment with any proper candidate-specific drive distributions. The discrepancy uses positive measured forces, positive measured accelerations and positive candidate masses. E12 was constructed as a generalized evidence rule. This calculation locates the extra measurement premise that would exclude it. [Checks against the existing readout](../tests/test_e12_conditioning_nuisance.py)

## Exactly which common functions of residual cost are representable?

For a single physical reading with \(k\ge1\) whitened spatial components, suppose one uses the same positive continuous evidence factor \(\exp g(K_\theta(Y))\) at two distinct positive masses, with \(K_\theta=\|Yv_\theta\|^2,\ v_\theta=(\cos\theta,-\sin\theta)\) and \(Y\) the \(k\times2\) matrix of observed channel pairs. Conditional on an arbitrary signed drive vector \(r\in\mathbb R^k\), the mean is \(r u_\theta^\top\) and all \(2k\) noise components are independent unit normals. The factor is assumed to reproduce the raw Gaussian-mixture likelihood **ratios**, on all data for the two candidates, for some proper pre-data laws \(\nu_\theta(dr)\). A mass-independent data factor is allowed. Then necessarily
\[
 g(x)=g(0)-\lambda x/2,\qquad 0\le\lambda<1,\qquad
 \nu_\theta=N\left(0,\frac{\lambda}{1-\lambda}I_k\right)
\]
in the noise-unit coordinates of the candidate graph at both masses. The zero value denotes a point mass at zero drive. Conversely these choices reproduce the claimed odds. At \(\lambda=1\), the flat profile corresponds to improper graph arclength as in T48. It has no proper push law in this representation. The result concerns factors that are a *common function of cost alone*, not all possible updates or all return laws on mass.

**Proof of the shape restriction.** First restrict \(Y\) to one active spatial coordinate. The multivariate Gaussian factorization above still makes the log odds separable in \(p=u_1^\top y,\ q=u_2^\top y\). Write \(c=\cos(\theta_2-\theta_1)>0,\ s=\sin(\theta_2-\theta_1)\ne0\), and \(G(t)=g(t^2)\). Expressing the two transverse projections using \(p,q\) gives log odds
\[
 G((q-cp)/s)-G((cq-p)/s).
\]
The two arguments are independent linear coordinates. Holding the second at zero shows \(G-G(0)\) is real analytic because the Gaussian-mixture odds were analytic. Its mixed derivative in \(p,q\) is \((-c/s^2)[G''((q-cp)/s)-G''((cq-p)/s)]\). Separability forces this to vanish. Thus \(G''\) is constant, and \(G\) is even, giving the asserted linear form for \(g\) with a real coefficient.

**Proof of the drive-law restriction.** On the entire \(k\)-component experiment put \(p=Yu_1,\ q=Yu_2\). The same expansion gives \(p_\theta(Y)/\phi_{2k}(Y)=M_\theta(Yu_\theta)\) with \(M_\theta(x)=\int\exp(r^\top x-\|r\|^2/2)\nu_\theta(dr)\). All vector pairs \(p,q\) are attainable. Since \(K_\theta(Y)=\|Y\|_F^2-\|Yu_\theta\|^2\), the exact odds identity forces \(M_1(p)e^{-\lambda\|p\|^2/2}=M_2(q)e^{-\lambda\|q\|^2/2}\) for all \(p,q\). Both sides are the same positive constant. Hence the probability measure proportional to \(e^{-\|r\|^2/2}\nu_\theta(dr)\) has moment-generating function \(e^{\lambda\|x\|^2/2}\). Its covariance forces \(\lambda\ge0\), and transform uniqueness identifies \(N(0,\lambda I_k)\). Undoing the tilt gives a finite measure precisely when \(\lambda<1\), and normalization gives the stated Gaussian. Direct integration proves the converse. [Gaussian MGF and uniqueness, Guntuboyina §3.2](https://www.stat.berkeley.edu/~aditya/resources/FullLectureNotes201AFall2019.pdf)

Choosing \(\lambda\) in this class specifies a signal variance on the candidate graph. L10 withholds that physical input. The transform and converse sharpen T48's distinction between normalized residual noise and a full sampling kernel. They neither invalidate its improper calculation nor derive a compulsory drive distribution from Newton's law. Exact event conditioning on a state with a declared measurement law is available; passing from the enlarged \((m,r)\) state to mass alone is a substantive extra step.

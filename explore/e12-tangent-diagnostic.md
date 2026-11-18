# A uniquely selected test of disagreement in the Gaussian tangent experiment

The [reciprocal tangent theorem](e12-tangent-pairings.md) gives signed local mass corrections and their noise scales at a supplied candidate mass. Here is one question for which these data force an answer. Once a false alarm rate is stated, which test has the greatest guaranteed power against disagreement of a specified size when its direction among readings is unknown?

## A unique probability law for the local displacement

The conditional common-tangent experiment below also has a positive law-level result. Under a common shift \(d\), the pooled reciprocal correction \(\bar b\) has distribution \(N(d,1/H)\), \(H=\sum h_i>0\). All orthogonal residual components are independent of it. Hold the fitted graph projections and those components fixed. Require the reported law for a displacement on the **real tangent line** to shift by \(c\) when the observed common correction shifts by \(c\). Require its upper endpoint at each rank \(0<u<1\) to cover the true \(d\) with conditional probability \(u\), for every \(d\). Then there is exactly one such law:
\[
                      \mathcal Q=N(\bar b,1/H).
\]
Indeed equivariance gives its rank-\(u\) quantile as \(\bar b+c_u\) for fixed ancillary data. Conditional coverage is \(\Phi(\sqrt H\,c_u)\), so \(c_u=\Phi^{-1}(u)/\sqrt H\). These quantiles determine the law and the displayed normal attains precisely those coverages. The same argument works for one active reading. Independent blocks combine by summing their tangent precisions and precision-weighted corrections. This is conditional calibration in a linear location experiment; it does not assert a normal posterior on positive mass for finite moves across the curved graphs.

## The uniquely selected disagreement test

Condition on the fitted graph projections and keep the \(n\ge2\) nonzero tangents. In the conditional Gaussian *tangent shift experiment*, the whitened transverse residual of reading \(i\) has mean \(d_i t_i\) and identity covariance. With \(h_i=\|t_i\|^2\), the observed reciprocal correction \(b_i=\langle t_i,r_i\rangle/h_i\) gives independent
\[
 y_i=\sqrt{h_i}b_i\sim N(\sqrt{h_i}d_i,1).
\]
Equality of the local shifts means \(d_i=d\) for every reading and an arbitrary common \(d\). Put \(w_i=\sqrt{h_i}\), \(P=I-ww^\top/(w^\top w)\), and \(J=\|Py\|^2=\sum h_i(b_i-\bar b)^2\). At the true candidate in the original measurement model \(d=0\) and this distribution holds exactly. Nonzero \(d_i\) describe shifts **along the candidate tangents**. A finite move to a new physical mass generally follows a curved graph and has different conditional means.

For each fixed size \(0<\alpha<1\), the test \(J>\chi^2_{n-1,1-\alpha}\) has size \(\alpha\) for every common \(d\). At each unknown common component and every fixed contrast magnitude \(\delta=\|P E y\|>0\), it uniquely maximizes the minimum power over all directions of that contrast, among all measurable tests controlling size at every common \(d\). Uniqueness is up to null sets. It has the same power at every such direction,
\[
  P\bigl(\chi'^2_{n-1}(\delta^2)>\chi^2_{n-1,1-\alpha}\bigr).
\]
No invariance restriction on the competing tests is needed. The premise here is that power in the *worst* direction at a fixed metric size is what the diagnostic must optimize. Declaring a specific alternative direction would change the problem. [Neyman and Pearson, 1933](https://doi.org/10.1098/rsta.1933.0009)

**Proof.** Fix the common component \(d w\). Under both null and all alternatives \(d w+\delta v\), \(v\) a unit vector in \(w^\perp\), the component along \(w\) has the same law, independently of the contrast. Any within-reading transverse residual left over from this tangent experiment has the same independent law as well. Average the alternatives over the unit sphere of \(w^\perp\), writing \(q=n-1\). Their likelihood ratio to the null is
\[
 e^{-\delta^2/2} E\,e^{\delta\,\|Py\| V_1}
\]
for a uniform unit vector \(V\) in \(\mathbb R^q\). Sign symmetry replaces the exponential by a hyperbolic cosine, making the ratio strictly increasing in \(\|Py\|\) for positive norm. This also holds for \(q=1\), where the sphere consists of two points. Neyman–Pearson gives the unique most powerful level-\(\alpha\) test of this mixture, the stated threshold. Its power depends only on \(\delta\). Any other test's minimum power is no greater than its mixture power, which cannot exceed the threshold test's power. Equality would require equality in Neyman–Pearson, forcing the same rejection region almost everywhere. Positive Gaussian densities make the exceptional sets common. The uniform distribution over alternatives is a proof of the bound, not a population law for the pushes.

For the example \((4,1),(1,4)\) at mass 1 and unit channel errors, \(J=9\) with one degree of freedom. The exact upper-tail probability at this candidate is \(2\Phi(-3)\simeq0.00270\). The common signed mass correction is zero. The same numbers in one declared planar reading give no between-reading contrast, and a within-reading residual of 9 with one degree of freedom. These statements refer to a **specified** candidate mass. Scanning many masses and quoting the most striking value would require another calibration.

This selects a way to *test* local disagreement if the stated decision objective is required. A tail probability here is not a probability that the candidate mass is true. Using it as a multiplier in a mass law would introduce a new assumption. [Finite checks](../tests/test_e12_tangent_diagnostic.py)

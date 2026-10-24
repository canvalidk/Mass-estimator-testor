# Weighting, repeated measurements, and mass uncertainty

**Third working contribution — 15 September 2026.** This follows the
[Gaussian-conditioning investigation](gaussian_conditioning.md). It examines
what could justify a measure, whether a proposed repeated-trial extension
learns the right mass, and what the uncertainty intervals actually achieve.
Newton II is assumed throughout.

## Findings

- The standardized distance-tube measure is exactly the **joint Jeffreys
  measure** for the Gaussian compatible-pair model. It has a precise
  interpretation in terms of local distinguishability of sampling laws.
- Pooling separately fitted first moments can converge to the **wrong mass**
  under repeated measurements of one unchanged latent pair. The earlier
  readout audit's general endorsement of that procedure needs qualification.
- Uniformly calibrated, always-finite mass intervals are impossible over
  unrestricted positive masses when arbitrarily weak excitation is allowed.
- A sampling study shows large differences between the three measures' 95%
  interval coverage. A narrow conditional interval can have very poor coverage
  at particular fixed true masses.
- We now have an exact local information benchmark. At strong compatible
  signals with isotropic errors, the measured magnitude ratio already has
  the efficient first-order influence. Using direction does not by itself
  establish an information advantage in that regime.

These findings concern different levels of the procedure. None contradicts
the composition theorem for the readout **given one probability law**.

## 1. A statistical interpretation of the tube measure

For known, parameter-independent covariance \(\Sigma\), take the sampling
model

$$
Y\mid\vartheta\sim N(\Phi(\vartheta),\Sigma),\qquad
\Phi(f,\alpha,\mathbf u)=(f\mathbf u,\alpha\mathbf u).
$$

In local coordinates for \(\vartheta=(f,\alpha,\mathbf u)\), let
\(J=D\Phi\). The score and Fisher information are

$$
\nabla_\vartheta\log p(Y\mid\vartheta)=J^T\Sigma^{-1}(Y-\Phi),
\qquad I(\vartheta)=J^T\Sigma^{-1}J.
$$

The Jeffreys rule assigns reference volume \(\sqrt{\det I}\,d\vartheta\).
This definition and the relevance of experimental design are discussed in
[Sun and Berger, §2.1](https://www2.stat.duke.edu/~berger/papers/OB-sequential.pdf).
The calculation above makes it exactly the compatible surface volume in the
\(\Sigma^{-1}\) metric. Therefore, for independent isotropic channels,

$$
d\Pi_J\propto
\left(\frac{f^2}{\sigma_F^2}+\frac{\alpha^2}{\sigma_a^2}\right)^{(d-1)/2}
df\,d\alpha\,d\Omega.
$$

The omitted factor \((\sigma_F\sigma_a)^{-1}\) is constant in the parameters.
This is the second contribution's tube measure. Full covariance requires
the full determinant, rather than separate scalar standardization.

Its interpretation is concrete: nearby Gaussian sampling laws have

$$
D_{\rm KL}(p_\vartheta\|p_{\vartheta+\delta})
=\tfrac12\delta^T I(\vartheta)\delta+o(\|\delta\|^2).
$$

The reference measure assigns equal weight to equal local information-volume
cells. This is a warrant for a particular convention. It does not establish
that the joint prior is optimal for the single quantity mass, or that its
credible intervals have a specified sampling coverage.

### The choices remain substantive in mass coordinates

Write \(m=f/\alpha\), \(\mathbf v=\alpha\mathbf u\), with
\(\alpha=\|\mathbf v\|\). Relative to \(dm\,d^d\mathbf v\), the three
reference densities are, up to parameter-independent constants,

| Measure | Density in \((m,\mathbf v)\) |
|---|---|
| Flat magnitudes | \(\alpha^{2-d}\) |
| Joint Jeffreys / covariance-metric tube | \(\alpha(m^2/\sigma_F^2+1/\sigma_a^2)^{(d-1)/2}\) |
| Angular tolerance | \(m^{d-1}\alpha^d\) |

The Jacobian is \(df\,d\alpha\,d\Omega=\alpha^{2-d}dm\,d^d\mathbf v\).
These expressions show how each construction weights mass and excitation
together. The measures are improper before observing data, so integrating out
their infinite nuisance volume does not produce a proper marginal mass prior.
The posteriors studied here are proper.

Calling a measure symmetric or coordinate-invariant cannot select it alone:
any measure transforms consistently when its Jacobian is retained. The
Jeffreys rule adds a particular information metric. Conversely, separately
declaring a density flat again after changing coordinates changes the measure.

## 2. What consistent updating requires

For repeated observations of the **same latent pair**, any one fixed reference
measure \(\Pi\) can be updated coherently:

$$
dP(\vartheta\mid Y_1,\ldots,Y_N)
\propto \left[\prod_i L_i(\vartheta)\right]d\Pi(\vartheta).
$$

Multiplication is associative and order-independent. This requirement alone
does not choose among flat, angular, or tube measures.

For repeated identical-covariance measurements, the full vector sample mean
with covariance \(\Sigma/N\) retains all the mean-parameter likelihood
information. Rebuilding the tube reference at \(\Sigma/N\) merely multiplies
it by a constant, so that special experiment does not distinguish the choices.

If the force-to-acceleration precision ratio changes, rebuilding the tube
reference can change its shape. For example, in \(d=3\), a first experiment
with inverse channel variances (precisions) \((1/\sigma_F^2,1/\sigma_a^2)=(1,1)\)
has weight \(f^2+\alpha^2\). After a second
experiment with precisions \((9,1)\), the joint-design weight is
\(10f^2+2\alpha^2\). Those are different priors.

For two observed-zero results with this combined likelihood, retaining the
first reference gives \(7/(11\sqrt5)=0.284590\), while choosing the
combined-design reference gives \(1/\sqrt5=0.447214\), in the chosen mass
unit. The difference is prior redesign. Ordinary updating with either prior
chosen once remains coherent. A planned-design reference can be chosen at
the outset and retained during updating.

Repeated observations of one pair and trials with different excitations but
one mass are different models. In the latter, the latent state is
\((m,\mathbf v_1,\ldots,\mathbf v_N)\); a prior or excitation model must be
specified there. Single-trial composition does not supply it.

## 3. A concrete failure of separately fitted pooling

The managed readout audit, §6, recommends

$$
\widehat m_{\rm pool}=\frac{\sum_i\overline f_i}{\sum_i\overline\alpha_i}
$$

from separately fitted single-trial distributions, describing its acceleration
weighting as right. The following counterexample shows that this is not a
generally consistent way to learn a common mass.

Use a known physical direction, independent Gaussian channels, flat positive
magnitude measures, and repeatedly measure one true pair
\(f=2,\alpha=1\), with \(\sigma_F=\sigma_a=1\). The true mass is 2.
For a scalar reading \(y\), the fitted positive latent mean is

$$
h_\sigma(y)=y+\sigma\frac{\phi(y/\sigma)}{\Phi(y/\sigma)},
$$

where \(\phi,\Phi\) are the standard normal density and CDF. With \(Z\sim N(0,1)\),

$$
E[h_1(2+Z)]=2.142871302725,\qquad
E[h_1(1+Z)]=1.414715858842.
$$

The law of large numbers therefore gives

$$
\boxed{\widehat m_{\rm pool}\longrightarrow1.514700842104\ne2.}
$$

This persists with unlimited independent measurements. Each separate fit
has a positivity correction; averaging those fitted means retains unequal
corrections in the two channels.

Combining the likelihoods first gives

$$
\widehat m_{\rm joint}
=\frac{h_{1/\sqrt N}(\overline F)}
       {h_{1/\sqrt N}(\overline a)}\longrightarrow2.
$$

The first construction fits a separate pair for each reading and pools its
moments. The second learns the common pair from all the data. The composition
proof addresses propagation within one joint law, and does not equate these
two inference operations. The example challenges the general pooling
recommendation, not the validity of the single-trial readout.

More generally in this fixture, put \(s=\sigma_F/\sigma_a\),
\(p=f/\sigma_F\), \(q=\alpha/\sigma_a\), and
\(H(t)=E[h_1(t+Z)]\). The separate-fit limit is \(sH(p)/H(q)\).
For positive signals it lies strictly between \(m\) and \(s\) when
\(m\ne s\). To see this, write \(H(t)=t+b(t)\), where
\(b(t)=E[\phi(t+Z)/\Phi(t+Z)]\) is positive and decreasing.
Then \(H(t)/t\) decreases, while \(H(t)\) increases because a truncated
normal's posterior mean increases with its reading. This bounds the ratio
between 1 and \(p/q\), in the appropriate order.

## 4. A limit on any uncertainty guarantee

Let an interval procedure return \(C(Y)=[L(Y),U(Y)]\) for positive mass.
At true zero excitation, \(\mathbf a^*=\mathbf0\) and
\(\mathbf F^*=m\mathbf a^*=\mathbf0\). The sampling law \(P_0\) is the
same for every \(m>0\).

If \(U<\infty\) almost surely, its coverage obeys

$$
P_0\{m\in C(Y)\}\le P_0\{U(Y)\ge m\}\longrightarrow0
\quad\text{as }m\to\infty.
$$

A strictly positive lower endpoint analogously fails as \(m\downarrow0\).
For uniform 95% coverage at null over all positive masses, an interval would
need \(P_0(U=\infty)\ge0.95\) and \(P_0(L=0)\ge0.95\). It must therefore
contain the whole positive mass domain on at least 90% of null datasets.

The obstruction survives exclusion of the exact boundary. In fixed units,
take \(m_n\to\infty\) and \(\|\mathbf a_n\|=m_n^{-2}\). Both true vectors
tend to zero. For fixed nonsingular Gaussian covariance, these sampling laws
converge to \(P_0\) in total variation, so probabilities of all interval events
differ from their null probabilities by a quantity tending to zero. Finite
intervals consequently cannot maintain a positive uniform coverage guarantee
over all positive masses and arbitrarily small nonzero excitations.

This does not assert that every fixed mass undercovers, and it does not apply
unchanged to a bounded mass domain or a declared lower bound on excitation.
Nor does it prohibit a conditional 95% credible interval: that is a probability
statement under the chosen posterior, rather than a uniform sampling guarantee.
It tells us that tuning cannot make these two promises universally coincide.

## 5. A first sampling calibration

The experiment uses true compatible pairs with an unknown common direction
in three dimensions, independent unit-variance Gaussian channels, and
\(s=\sigma_F/\sigma_a=1\). Let \(r=m/s\) and
\(t=\|\mathbf a^*\|/\sigma_a\). The true standardized means are
\((rt,0,0)\) and \((t,0,0)\). The inference receives the noisy vectors,
their covariance, and the selected measure; it does not receive the true mass
or direction.

For each of 12 scenarios, 2,048 independent datasets were generated. The same
noise draws were reused across scenarios and measures to make comparisons
paired. Each construction produces an equal-tail interval containing 95% of
its conditional mass distribution. Coverage below is the fraction of datasets
on which that interval contains the **fixed true mass**.

| True \(m/s\) | True acceleration SNR \(t\) | Flat coverage | Tube coverage | Angular coverage |
|---:|---:|---:|---:|---:|
| 1 | 0 — null | 99.95% | 99.12% | 99.90% |
| 1 | 4 | 95.12% | 94.38% | 96.19% |
| 4 | 0.25 | 96.83% | 93.12% | 51.86% |
| 4 | 1 | 97.02% | 96.29% | 79.98% |
| 16 | 0 — null | 73.19% | 65.43% | 0 / 2,048 |
| 16 | 0.25 | 93.31% | 92.33% | 7.96% |
| 16 | 1 | 97.31% | 97.27% | 74.56% |
| 16 | 4 | 95.07% | 95.07% | 93.55% |

The numerical script prints all 12 scenarios; its optional JSON also includes
Wilson intervals for each coverage proportion, median point error, and RMS
log point error. With 2,048 trials, the largest binomial standard error is
1.10 percentage points. The large differences above substantially exceed
that sampling uncertainty. Zero observed successes do not prove a zero
population coverage probability.

There are corresponding point effects. At \(m/s=4,t=1\), the median returned
point divided by truth is 0.709 for flat, 0.742 for tube, and 0.499 for angular.
At \(m/s=16,t=4\), those ratios are 0.995, 0.995, and 0.892. The geometric
construction affects practical behavior well beyond the exact zero-reading
examples.

The angular construction's narrow conditional intervals have poor coverage
on several of these fixed-parameter experiments. Flat and tube also fail a
uniform 95% interpretation, as the null theorem requires. This is a targeted
calibration study, not a ranking over all experiments, priors, or loss functions.

## 6. A precise benchmark for the information question

Write the true acceleration vector as \(\mathbf v\). In the model
\(Y\sim N((m\mathbf v,\mathbf v),\Sigma)\), eliminate the unknown
\(\mathbf v\) from the Fisher information by taking its Schur complement.
The resulting information about mass is

$$
I_{m\cdot v}=\mathbf v^T V(m)^{-1}\mathbf v,
\quad
V(m)=\Sigma_{FF}+m^2\Sigma_{aa}
-m(\Sigma_{Fa}+\Sigma_{aF}).
$$

One derivation uses the contrast \(K=[I,-mI]\). It annihilates the nuisance
mean derivative \([mI;I]\), and satisfies \(K\Sigma K^T=V(m)\) and
\(K[\mathbf v;0]=\mathbf v\). The projection of the mass score off the
nuisance score gives the displayed result.

For independent isotropic channels,

$$
I_{m\cdot v}=\frac{\|\mathbf v\|^2}{\sigma_F^2+m^2\sigma_a^2},
\qquad
I_{\log m}^{-1}=\frac1{\mathrm{SNR}_F^2}+\frac1{\mathrm{SNR}_a^2}.
$$

Here \(\mathrm{SNR}_F=m\|\mathbf v\|/\sigma_F\) and
\(\mathrm{SNR}_a=\|\mathbf v\|/\sigma_a\).

This quantifies the weaker-channel bottleneck. It is an exact local information
calculation, not a claim that the finite-sample ratio distribution has finite
variance or that any tested estimator attains this information bound.

At an identified nonzero state \(\mathbf v=\alpha\mathbf u\), the measured
magnitude ratio has first-order error

$$
\frac{\|\widetilde{\mathbf F}\|}{\|\widetilde{\mathbf a}\|}-m
=\frac{\mathbf u\cdot\epsilon_F-m\mathbf u\cdot\epsilon_a}{\alpha}
+O(\|\epsilon\|^2).
$$

Under isotropic noise, this is the efficient first-order influence. The
perpendicular components still affect finite-data direction inference and
global likelihood shape, but they do not add a first-order mass score after
the nuisance direction is accounted for. An advantage over the magnitude
ratio therefore needs a specified finite-signal, covariance, uncertainty,
or reuse comparison.

For independent isotropic channels and invariant inference with no preferred
direction, the three observed
quantities \(\|\widetilde F\|^2,\|\widetilde a\|^2,
\widetilde F\cdot\widetilde a\) retain the full geometry up to a common
orthogonal transformation. A norm-only procedure discards the third. This
is an invariance statement, not sufficiency for estimating the oriented
latent vector and not a numerical measure of mass information lost.

## 7. Consequences and next questions

The investigation has progressed from a generic free weighting to concrete
choices with concrete consequences:

1. **Single-trial measure:** joint Jeffreys has a warrant; mass-focused
   inference may justify a different nuisance treatment. Its name alone
   does not choose it.
2. **Repeated learning:** use a joint model for the quantities actually shared
   between measurements. Separately fitted moment pooling has no general
   consistency guarantee and fails the explicit example above.
3. **Uncertainty promise:** distinguish conditional mass probability from
   repeated-sampling coverage, and state the operating range for any
   calibration claim. The weak-excitation obstruction is structural.
4. **Information comparison:** assess performance relative to the available
   information under a stated model and regime. Counting vector components
   or comparing interval widths does not answer that question.

These findings constrain further tuning. A useful next target is the
mass-focused joint model for repeated trials, with explicit excitation
assumptions and an uncertainty contract that acknowledges weak identification.

## Additional correction to the earlier readout audit

Its §3 says the two directional regression readouts bracket the ratio of
means. This is not general. Equiprobable pairs \((f,\alpha)=(1,1),(9,2)\) give

$$
\frac{E[f]}{E[\alpha]}=\frac{10}{3}
<\frac{E[f\alpha]}{E[\alpha^2]}=\frac{19}{5}
<\frac{E[f^2]}{E[f\alpha]}=\frac{82}{19}.
$$

Only the ordering between the two regression readouts is guaranteed by
Cauchy–Schwarz when the second moments exist. This qualifies the general
bracketing statement; it does not dispute its displayed numerical examples
or the later rigorous composition theorem.

## Reproduction, verification, and sources

[weighting_checks.py](weighting_checks.py) reproduces the pooling counterexample,
Fisher identities, and sampling study. It uses Python and NumPy. The default
run checks the algebraic/numerical fixtures and all 12 sampling scenarios;
`--checks-only` omits the Monte Carlo study. Files are written only if an
`--output` path is supplied.

~~~powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/weighting_checks.py'
~~~

Checks include:

- Gauss–Hermite orders 64, 128, 256 for the pooling limit.
- Fisher-volume and efficient-information checks for 60 random cases in
  dimensions 1–3, including comparison with the existing local-uncertainty code.
- Analytic null CDFs, reciprocal symmetry, and the prior contribution's
  independently verified mass integrator.
- Quadrature refinement from 96 radial nodes and 64 nodes per angular half
  to 128 and 96. In each scenario, a fixed 32 observations and the 32 nearest
  a coverage boundary were checked; overlaps were counted once. Largest
  relative point change was \(5.94\times10^{-8}\), largest absolute CDF change
  \(3.37\times10^{-10}\); coverage classifications agreed on that subset.
- An independent Cartesian-wedge CDF calculation on nine actual sampled
  observations agreed within \(7.23\times10^{-11}\), including weak and
  strong signals. These are checked-case errors, not a certified bound over
  every simulated observation.

The probability obstruction and inconsistency result rest on proofs; Monte
Carlo illustrates particular consequences. Fixed seed: 20260915. Machine
results, source hashes, independent reviews, and development probes remain
under `.tools/mass_equation_20260915/weighting/`.

Managed main was resolved once to
[8782a293297330f7a354d87ebe62792018fa8866](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
the same snapshot as the previous contributions. The catalog was consulted.
This pass read the [readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md)
in full, particularly §§3, 5–6, and the
[working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md)
on repeated trials. The [main investigation](README.md) supplies the prior
source map and governing user steering. The local scope was the active
estimator, practical formulas, earlier investigation scripts, and the new
checks. Selection targeted weighting, pooling, common mass, calibration, and
Fisher information. The external source above supports the Jeffreys definition;
this was not an estimator-priority search. No managed source or baseline
estimator was changed, and no rejection gate was introduced.

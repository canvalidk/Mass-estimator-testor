# Zero readings can inform shared candidates without changing the mass marginal

**Record — 28 September 2026.** Written by Codex at the user's request.

**Finding.** An all-zero experiment can leave the mass law equal to its starting law while updating shared candidate size. That update can affect mass once a nonzero reading is combined with it. An unchanged mass marginal is therefore not sufficient evidence that a reading supplied no information.

A second distinction is essential: omitting a structurally absent coordinate, conditioning an uncertain coordinate on a precise zero, and rebuilding a lower-dimensional default model are different operations. Their equality needs an additional compatibility premise.

## 1. What the rotation story does and does not establish

Rotating two aligned observed vectors onto one axis makes their observed transverse components zero. It does not establish that the corresponding aligned candidate vectors have no transverse components.

In the common-noise version of ours27, the standardized conditional candidate law is $N(A_\theta,I)$. If $A_\theta$ has zero transverse coordinates, those candidate coordinates still have $N(0,1)$ distributions. Directional uncertainty is present. Flat integration makes the transverse contribution to the **mass evidence** constant; it does not remove the candidate uncertainty.

Thus the precise issue is whether the chosen candidate measure should allow that uncertainty to feed back into mass through shared size or direction. Rotation symmetry alone does not prohibit such feedback.

If a frame is selected using the observations, its resulting zero components are not additional independent observations. A pointwise rotation of the original data, covariance and candidate model is valid; retaining the original experiment is essential when interpreting the result.

## 2. A shared-size example

Take one common force–acceleration noise pair, independent standardized component errors, and $M>0$ stacked candidate components. Put

$$
x=F/\sigma_F,\qquad y=a/\sigma_a,\qquad
s=\sigma_F/\sigma_a,\qquad m=s\tan\theta,
\qquad 0<\theta<\pi/2.
$$

Define the projected and residual observations

$$
A_\theta=x\sin\theta+y\cos\theta,\qquad
B_\theta=x\cos\theta-y\sin\theta,
$$

$$
H_\theta=\|A_\theta\|^2,\qquad
Q_\theta=\|B_\theta\|^2,\qquad
T=H_\theta+Q_\theta=\|x\|^2+\|y\|^2.
$$

Consider the illustrative candidate family

$$
v\mid t\sim N(0,(t-1)I_M),\qquad t\geq1.
$$

The observed parallel coordinates then have variance $t$, and the residual coordinates have variance one. Apart from data-only constants,

$$
L(m,t)=t^{-M/2}
\exp\left(-\frac{Q_\theta}{2}-\frac{H_\theta}{2t}\right).
$$

For fixed $m$, the Fisher information of this variance experiment is $I_{tt}=M/(2t^2)$. Its one-parameter Fisher length measure is proportional to $dt/t$. Using that measure gives

$$
p(m,t\mid D)\propto
\pi(m)t^{-M/2-1}
\exp\left(-\frac{Q_\theta}{2}-\frac{H_\theta}{2t}\right).
$$

Here $\pi$ is a proper mass starting density. The posterior is proper for finite data: the data factor is bounded and the tail $t^{-M/2-1}$ is integrable for $M>0$.

This is a conditional reference construction inside a declared Gaussian family. It does not derive that family or uniquely select a mass definition. Also, $t\geq1$ is not the full multiplicative scale group, so an unqualified Haar-measure claim would be inappropriate. The construction was already explored in the earlier learned-size investigation; its existence is not a new discovery.

Set $\beta=1-1/t$. Then

$$
p(\beta\mid m,D)\propto
(1-\beta)^{M/2-1}\exp(\beta H_\theta/2),
\qquad 0<\beta<1,
$$

and

$$
v\mid\beta,m,D\sim N(\beta A_\theta,\beta I_M).
$$

Equivalently, the marginal mass evidence is

$$
e^{-T/2}\int_0^1
(1-\beta)^{M/2-1}e^{\beta H_\theta/2}\,d\beta.
$$

## 3. The zero reading updates size

For an all-zero dataset, $H_\theta=Q_\theta=0$ for every mass. Consequently,

$$
p(m\mid D=0)=\pi(m),
\qquad
\beta\mid D=0\sim\operatorname{Beta}(1,M/2),
\qquad
E[\beta\mid D=0]=\frac{2}{M+2}.
$$

There is no preferred mass angle, but there is information about the shared size state.

More generally, appending $k$ independent, noisy zero components multiplies the joint density by

$$
t^{-k/2}=(1-\beta)^{k/2}.
$$

With active nonzero data, $m$ and $\beta$ are generally dependent through $H_\theta$. Changing the size distribution can therefore change the mass marginal. Integrating out size and retaining only $p(m)=\pi(m)$ would lose information needed for this later update.

Multiplying likelihoods and then integrating the **same joint state** is order independent. This is the practical consequence of the distinction established in the [axis-splitting record](ours27_axis_splitting_does_not_force_independent_candidates_2026-09-28.md).

The result does not say that every zero reading must alter mass. In this example an isolated all-zero dataset explicitly does not. Nor does it establish a universally neutral zero-reading rule for heterogeneous instruments; the displayed symmetry belongs to the declared common-noise experiment.

## 4. Omission and precise conditioning have different limits

Let one uncertain transverse candidate coordinate and its observation be

$$
b=\delta\sqrt{t-1}\,\xi,\qquad
Y=b+\epsilon\eta,
$$

where $\xi,\eta$ are independent standard Gaussians. At the observed value $Y=0$, remove the likelihood constant independent of $t$. The remaining factor is exactly

$$
\left[1+\left(\frac{\delta}{\epsilon}\right)^2(t-1)\right]^{-1/2}.
$$

This gives three distinct comparisons:

| Comparison | Effect on the shared size law |
|---|---|
| $\delta/\epsilon\to0$ | The factor tends to one: the candidate coordinate becomes structurally negligible relative to measurement resolution. |
| $\delta=\epsilon$, with $\epsilon\to0$ | The factor remains $t^{-1/2}$ for every positive $\epsilon$. |
| Fixed $\delta>0$, with $\epsilon\to0$ | The observation increasingly conditions an unknown physical coordinate; it is not an omitted input. |

Weak collapse of a physical candidate coordinate to zero is therefore insufficient for statistical irrelevance when the measurement error shrinks at the same rate. The standardized experiment still distinguishes values of $t$.

At exact structural absence, an **independent** added noise coordinate contributes a parameter-free factor. If its noise is correlated with retained measurement errors, however, observing it can inform those errors. Structural neutrality then requires either independence or the appropriate conditional treatment of the retained covariance.

A lower-dimensional conditional restriction of the original joint measure also need not equal a fresh lower-dimensional default. Resetting the default forgets how the zero was established.

## 5. The previous continuity failure remains real

The distinctions above do not make the scalar reference law pass the stronger test “arbitrarily precise measured zero equals a fresh lower-dimensional default.”

For example, use a uniform starting law in $\theta$, unit active-axis errors, and aligned active observations $F=3$, $a=1$. Add two zero components whose force and acceleration errors both have standard deviation $\epsilon$. In this noise-relative candidate model the mass law stays three-dimensional for every $\epsilon>0$, although the physical candidate lengths change.

Use the existing ratio-of-mean-sizes readout. With $D_\epsilon=\operatorname{diag}(1,\epsilon,\epsilon)$ in three dimensions, it is

$$
\widehat m
=\frac{E[\sin\theta\,\|D_\epsilon v\|]}
       {E[\cos\theta\,\|D_\epsilon v\|]}.
$$

The expectations use the same joint candidate law. The one-dimensional case has only the active entry of $D_\epsilon$.

| Experiment | Mass readout |
|---|---:|
| One active dimension | 2.2064082220 |
| Three dimensions, $\epsilon=1$ | 2.0184357974 |
| Three-dimensional model, limit $\epsilon\to0$ | 2.0696736242 |

The limiting result differs from a fresh one-dimensional fit. The channel-noise ratios, and hence the uniform-angle starting law, are unchanged in this comparison. The discrepancy cannot be dismissed as a change in $\pi$.

These values were checked by deterministic integration; increasing the mass-law quadrature order from 200 to 400 changed its density by less than $4\times10^{-14}$. The exact likelihood factor above establishes the persistence of the shared-state information independently of the numerical example.

## 6. What this licenses

The [existing candidate-size and null-readings record](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours27/ours27_candidate_size_and_null_readings_2026-09-28.md) already distinguishes known-empty from measured-empty axes within ours27. The new result here is the shared-state update and the explicit separation of structural, conditioning and model-reset limits.

Measured-zero neutrality is a property of a specified estimator; it has not been derived from rotation or the ratified shared-quantity argument. Permitting shared-size feedback would therefore require an explicit change to the stronger neutrality requirement, if that requirement is retained as a defining premise.

The scalar example demonstrates the distinction. It is not a recommended final law: the [repair-limits record](ours27_limits_of_scalar_and_tensor_size_repairs_2026-09-28.md) supplies a separate regular-regime counterexample to its general width calibration.

## 7. Provenance and status

Managed context was read at **canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8**, selected through its catalog. Alongside the linked null-readings record, relevant sources are the [full-covariance formulation](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md), the [learned-push-scale proposal](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_push_scale_2026-09-27.md), and its [status as a hint](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_omega_is_a_hint_2026-09-27.md).

Local verification is in .tools/ours27_direction_2026_09_28/shared_size_probe.py and measure/null_conditioning_audit.py within the same scratch directory. The derivations above are self-contained.

This local record has not been submitted to VD-docs. No new mass equation, starting-law uniqueness claim or persistence argument is adopted here.

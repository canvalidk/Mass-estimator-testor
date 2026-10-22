# A Median-Based Estimator for Inertial Mass

## 1. Purpose

Newton's second law relates the true net force, true inertial acceleration, and inertial mass of one material object:

$$
\mathbf F^*=m\mathbf a^*,
\qquad
m>0.
$$

In an experiment, however, we do not normally receive the exact values $\mathbf F^*$ and $\mathbf a^*$. We receive measurements

$$
\widetilde{\mathbf F},
\qquad
\widetilde{\mathbf a},
$$

with uncertainties

$$
\Sigma_F,
\qquad
\Sigma_a.
$$

The earlier errors-in-variables estimator asks:

> Which positive mass permits the **single best corrected force–acceleration pair**?

That is a useful question, but it discards most of the uncertainty distribution. A different approach is:

> For each possible positive mass, how much total support does the complete measurement model assign to it?

Once that distribution over positive mass has been constructed, the reported mass can be taken to be its **median**.

The estimator is therefore not

$$
\widehat m=\arg\min Q(m),
$$

but

$$
\boxed{
P(m\leq \widehat m_{50}\mid \text{measurements})=\frac12.
}
$$

The result $\widehat m_{50}$ is the posterior median of inertial mass.

---

## 2. Why this is different from choosing the best corrected pair

Suppose the measured vectors are exactly zero:

$$
\widetilde{\mathbf F}=\mathbf0,
\qquad
\widetilde{\mathbf a}=\mathbf0.
$$

A best-fit-pair procedure observes that, for every positive mass, it can choose

$$
\mathbf a^*=\mathbf0,
\qquad
\mathbf F^*=\mathbf0.
$$

Then

$$
\mathbf F^*=m\mathbf a^*
$$

is satisfied exactly for every $m>0$.

Consequently the minimum residual is

$$
Q_{\min}(m)=0
$$

for every mass.

From this, however, one may conclude only:

> Every positive mass possesses a perfectly fitting latent point.

It does **not** follow that every mass has the same total evidential support.

A candidate mass may permit a large region of plausible latent force–acceleration states, while another candidate mass may permit only a very narrow region around the exact zero point.

The median estimator preserves this information.

---

# 3. Measurement model

Assume the instruments obey

$$
\widetilde{\mathbf F}
=
\mathbf F^*+\boldsymbol\epsilon_F,
$$

$$
\widetilde{\mathbf a}
=
\mathbf a^*+\boldsymbol\epsilon_a,
$$

with

$$
\mathbf F^*=m\mathbf a^*.
$$

For Gaussian measurement errors,

$$
\boldsymbol\epsilon_F\sim
\mathcal N(\mathbf0,\Sigma_F),
$$

$$
\boldsymbol\epsilon_a\sim
\mathcal N(\mathbf0,\Sigma_a).
$$

For a proposed mass $m$ and proposed underlying acceleration $\mathbf a^*$, the likelihood of obtaining the measurements that actually occurred is therefore

$$
L(m,\mathbf a^*)
\propto
\exp\left[
-\frac12
\left(
(\widetilde{\mathbf F}-m\mathbf a^*)^T
\Sigma_F^{-1}
(\widetilde{\mathbf F}-m\mathbf a^*)
+
(\widetilde{\mathbf a}-\mathbf a^*)^T
\Sigma_a^{-1}
(\widetilde{\mathbf a}-\mathbf a^*)
\right)
\right].
$$

This is essentially the same measurement model used by the earlier errors-in-variables estimator.

The difference comes in what we do with it.

---

# 4. The old operation: maximize over the hidden state

The old estimator effectively performs

$$
L_{\mathrm{profile}}(m)
=
\max_{\mathbf a^*}
L(m,\mathbf a^*).
$$

It asks:

> For this mass, what is the single most favourable underlying acceleration?

Then it compares masses according to those best cases.

This operation is called **profiling** over the latent variable.

It throws away everything except the highest point of the likelihood for each $m$.

Schematically:

```text
all possible latent accelerations
          |
          v
take only the best one
          |
          v
compare masses
```

This is why the exact zero-zero observation becomes problematic.

Every mass has the same perfect point:

$$
\mathbf a^*=0.
$$

The profile likelihood therefore sees no difference between them.

---

# 5. The new operation: integrate over the hidden state

The median-based approach instead keeps the entire region.

We introduce an **experimental or latent-state model**

$$
\pi(\mathbf a^*),
$$

describing which underlying accelerations were plausible in the experiment before observing the measurement.

Then calculate

$$
L_{\mathrm{marg}}(m)
=
\int
L(m,\mathbf a^*)
\pi(\mathbf a^*)\,d\mathbf a^*.
$$

This is the marginal likelihood for mass.

Its interpretation is:

> Add together the support from every underlying acceleration that could plausibly have produced the observed measurements under this candidate mass.

Schematically:

```text
candidate mass m
     |
     v
all possible latent accelerations a*
     |
     +-- very plausible
     +-- somewhat plausible
     +-- barely plausible
     +-- implausible
     |
     v
add the support from all of them
     |
     v
total support for m
```

Two masses can therefore possess the same maximum likelihood while having different marginal likelihoods.

One may have a broad region of compatible latent states.

Another may have only one narrow spike.

The profile estimator treats those cases as equal if their peaks have the same height.

The marginal estimator does not.

---

# 6. Turn the support into a distribution over positive mass

We must also state what masses were possible before the measurement. Let

$$
\pi(m)
$$

be a prior or experimentally licensed distribution over

$$
m>0.
$$

Then

$$
\boxed{
p(m\mid
\widetilde{\mathbf F},
\widetilde{\mathbf a})
\propto
\pi(m)
L_{\mathrm{marg}}(m).
}
$$

Equivalently,

$$
\boxed{
p(m\mid
\widetilde{\mathbf F},
\widetilde{\mathbf a})
\propto
\pi(m)
\int
L(m,\mathbf a^*)
\pi(\mathbf a^*)\,d\mathbf a^*,
\qquad m>0.
}
$$

This produces an entire probability distribution over positive inertial mass rather than one point.

Negative masses never appear because the domain is

$$
m>0.
$$

---

# 7. The median estimator

Once the mass distribution has been obtained, define

$$
\widehat m_{50}
$$

by

$$
\boxed{
\int_0^{\widehat m_{50}}
p(m\mid
\widetilde{\mathbf F},
\widetilde{\mathbf a})\,dm
=
\frac12.
}
$$

Half of the posterior probability lies below $\widehat m_{50}$, and half lies above it.

Thus

$$
\boxed{
\widehat m_{\mathrm{median}}
=
\operatorname{Median}
\left[
p(m\mid\text{data})
\right].
}
$$

The estimator can also report an interval, for example

$$
[m_{0.16},m_{0.84}]
$$

or

$$
[m_{0.025},m_{0.975}],
$$

rather than pretending that the point median contains all available information.

---

# 8. Why use the median?

The median is not being selected because it is universally the one true statistic.

It corresponds to a particular estimation rule.

The posterior mean minimizes expected squared error:

$$
E[(\widehat m-m)^2].
$$

The posterior median minimizes expected absolute error:

$$
E[|\widehat m-m|].
$$

For distributions with a large high-mass tail, these can behave very differently.

Imagine schematically that most of the distribution sits around $5\,\mathrm{kg}$, but a tiny probability extends to enormously large masses.

Then one might have something qualitatively like

```text
mode       ≈ 5 kg
median     ≈ 5 kg
mean       >> 5 kg
```

because the mean is pulled by the distant tail.

The median scarcely moves.

This is particularly attractive in the present problem because ratios and near-zero denominators naturally generate long mass tails.

---

# 9. How different measurements produce different median estimates

The important feature of this estimator is that it need not behave identically in every low-signal case.

The **experimental context matters**.

---

## 9.1 Strong aligned force and acceleration

Suppose the experiment measures approximately

$$
\widetilde F=20\ \mathrm N,
\qquad
\widetilde a=4\ \mathrm{m\,s^{-2}},
$$

with uncertainties small compared with those values.

Then almost all plausible corrected states lie near

$$
F^*\approx20\ \mathrm N,
\qquad
a^*\approx4\ \mathrm{m\,s^{-2}}.
$$

Newton II therefore requires

$$
m\approx\frac{20}{4}=5\ \mathrm{kg}.
$$

The posterior becomes narrow around $5\,\mathrm{kg}$.

Consequently,

$$
\text{mode}
\approx
\text{median}
\approx
\text{mean}
\approx
5\ \mathrm{kg}.
$$

In the strong-signal limit, the choice among these summaries becomes largely irrelevant.

The exact-law result reappears.

---

# 10. Nonzero force but unresolved acceleration

Suppose instead that a controlled experiment establishes approximately

$$
F^*\approx50\ \mathrm N,
$$

while the acceleration instrument returns

$$
\widetilde a=0
$$

with an absolute uncertainty

$$
\sigma_a=0.1\ \mathrm{m\,s^{-2}}.
$$

Newton II predicts

$$
a^*=\frac{50}{m}.
$$

A small mass predicts an easily detectable acceleration.

For example,

$$
m=10\ \mathrm{kg}
\quad\Rightarrow\quad
a^*=5\ \mathrm{m\,s^{-2}},
$$

which is wildly inconsistent with a reading near zero at this precision.

A much larger mass predicts a much smaller acceleration:

$$
m=500\ \mathrm{kg}
\quad\Rightarrow\quad
a^*=0.1\ \mathrm{m\,s^{-2}}.
$$

That is entirely plausible.

The posterior distribution therefore shifts toward **large masses**.

Instead of a symmetric point determination, the result may look like

```text
very small mass        strongly disfavoured
moderate mass          somewhat disfavoured
large mass             plausible
very large mass        increasingly hard to distinguish
```

The most useful output may consequently be a lower bound such as

$$
m>m_{\min}
$$

at a declared probability level.

The median may also exist, but its numerical value can depend strongly on the permitted upper mass range or prior because the experiment has mainly established that the object is **heavy**, rather than exactly how heavy.

---

# 11. Nonzero acceleration but unresolved force

Reverse the experiment.

Suppose an acceleration is independently known to be appreciable,

$$
a^*\approx5\ \mathrm{m\,s^{-2}},
$$

but the force channel reports approximately

$$
\widetilde F=0.
$$

Newton II predicts

$$
F^*=ma^*.
$$

Large masses predict large forces and therefore become inconsistent with the force non-detection.

Small masses predict forces close to zero.

The posterior therefore moves toward **small masses**.

So:

```text
known force + unresolved acceleration
        -> pushes mass upward

known acceleration + unresolved force
        -> pushes mass downward
```

The same numerical-looking “zero” measurement therefore has different meaning depending on how the experiment was generated.

---

# 12. What happens at measured $0/0$?

Now consider

$$
\widetilde{\mathbf F}=\mathbf0,
\qquad
\widetilde{\mathbf a}=\mathbf0.
$$

The profile estimator says every positive mass has a perfect latent solution:

$$
\mathbf F^*=\mathbf0,
\qquad
\mathbf a^*=\mathbf0.
$$

The median estimator asks a stronger question:

> Apart from that exact point, what other latent states could plausibly have produced the observation for each candidate mass?

The answer depends on the experimental model.

If the experiment made appreciable acceleration states plausible, then large $m$ implies that even small latent accelerations would generate appreciable force:

$$
F^*=ma^*.
$$

A zero force reading may therefore penalize large masses.

If the experiment instead made appreciable force states plausible, then small $m$ predicts appreciable accelerations:

$$
a^*=\frac{F^*}{m}.
$$

A zero acceleration reading may penalize small masses.

If the experiment genuinely expected

$$
F^*=a^*=0
$$

almost certainly, then zero-zero may carry almost no mass information.

Therefore:

$$
\boxed{
\text{a }0/0\text{ reading does not have one universal interpretation.}
}
$$

Its information about mass is conditional on what latent states the experiment could reasonably have generated.

---

# 13. A simple mathematical illustration of the $0/0$ effect

Take a one-dimensional experiment and suppose the latent acceleration before measurement has a Gaussian experimental distribution

$$
a^*\sim\mathcal N(0,\tau_a^2).
$$

The measurement model is

$$
\widetilde a=a^*+\epsilon_a,
\qquad
\widetilde F=ma^*+\epsilon_F,
$$

where

$$
\epsilon_a\sim\mathcal N(0,\sigma_a^2),
\qquad
\epsilon_F\sim\mathcal N(0,\sigma_F^2).
$$

For an observed zero-zero pair,

$$
\widetilde F=0,
\qquad
\widetilde a=0,
$$

integrating over $a^*$ gives a mass-dependent likelihood proportional to

$$
\boxed{
L(0,0\mid m)
\propto
\frac{1}{
\sqrt{
m^2\tau_a^2\sigma_a^2+
\sigma_F^2(\tau_a^2+\sigma_a^2)
}
}
}.
$$

This is **not constant in $m$**.

So although every mass possesses the exact latent solution $a^*=0$, the total probability associated with the observed zero-zero result varies with mass.

In this particular experimental model, larger masses are increasingly penalized.

That happens because the experiment allowed nonzero accelerations before measurement. A large mass would turn those into large forces, making simultaneous zero force and zero acceleration less probable.

A differently designed experiment can produce the opposite dependence.

That is why the experimental excitation model belongs in the estimator.

---

# 14. The special role of $\sigma_F/\sigma_a$

The ratio

$$
\boxed{
m_\sigma=\frac{\sigma_F}{\sigma_a}
}
$$

has units of mass and marks the scale at which the two measurement noise floors correspond under Newton II:

$$
\sigma_F=m_\sigma\sigma_a.
$$

It is therefore a natural **apparatus mass scale**.

There is also a separate result concerning the positive stationary-root statistic previously studied for repeated symmetric noise-only trials.

Under the isotropic symmetric noise model, swapping normalized force and acceleration noise transforms the dimensionless root as

$$
r\rightarrow\frac1r.
$$

This reciprocal symmetry places the median at

$$
r=1,
$$

hence

$$
\boxed{
\operatorname{median}(M_{\mathrm{root}})
=
\frac{\sigma_F}{\sigma_a}.
}
$$

That result is real, but it must not be confused with the posterior-median estimator defined in this note.

They are different constructions:

```text
old null simulation:
    calculate one positive stationary root per noise realization
    -> take median across repeated realizations

posterior-median estimator:
    observe the actual data
    -> construct p(m | data)
    -> take the median of that mass distribution
```

The first has the $\sigma_F/\sigma_a$ median by reciprocal symmetry.

The second has whatever median follows from the complete measurement, excitation, and prior model.

If one wishes the posterior estimator also to reproduce

$$
\sigma_F/\sigma_a
$$

as its calibrated null median, that becomes an additional design requirement on the model rather than an automatic theorem.

---

# 15. Anti-aligned measurements

Suppose

$$
\widetilde{\mathbf F}\cdot
\widetilde{\mathbf a}<0.
$$

The raw measurements point in directions incompatible with an exact positive-scalar Newton-II pair.

The old stationary-point method can produce a positive algebraic root even though that root may correspond to a maximum rather than a minimum of the residual objective.

The posterior method does not need to manufacture a negative mass or select that root.

Instead it asks:

> For each $m>0$, how plausible is it that the underlying true vectors were codirectional but the measurement errors produced the observed anti-aligned pair?

If the uncertainties are tiny, this event is extremely improbable for every mass.

The posterior therefore says chiefly:

```text
the data are badly compatible with the Newton-II measurement model
```

rather than pretending to have a precise mass.

If the uncertainties are large, substantial corrections may be plausible and a broad positive-mass distribution can remain.

Thus an anti-aligned trial can yield:

```text
broad mass estimate
weak information
large model-discrepancy warning
```

rather than a negative mass.

---

# 16. Several trials

The method becomes much more useful when several independent trials share one mass.

For trials $i=1,\ldots,N$,

$$
\mathbf F_i^*=m\mathbf a_i^*.
$$

The joint posterior becomes

$$
p(m\mid D)
\propto
\pi(m)
\prod_i
\left[
\int
L_i(m,\mathbf a_i^*)
\pi_i(\mathbf a_i^*)
\,d\mathbf a_i^*
\right].
$$

Then define

$$
\boxed{
\int_0^{\widehat m_{50}}
p(m\mid D)\,dm
=
\frac12.
}
$$

One weak or nearly perpendicular trial no longer controls the result.

A collection of strong trials can concentrate the mass distribution tightly.

As the evidence improves,

$$
p(m\mid D)
\longrightarrow
\text{a narrow distribution around }m_0.
$$

Then

$$
\operatorname{mean}(m)
\approx
\operatorname{median}(m)
\approx
\operatorname{mode}(m)
\approx
m_0.
$$

The disagreement among estimators is therefore greatest precisely where the experiment carries the least information.

---

# 17. What additional inputs this estimator requires

The increased information does not come for free.

The operation requires more than

$$
\widetilde{\mathbf F},
\qquad
\widetilde{\mathbf a},
\qquad
\Sigma_F,
\qquad
\Sigma_a.
$$

It also needs to know something about the experiment that produced the latent states.

At minimum:

```text
measured force
measured acceleration

force uncertainty
acceleration uncertainty

Newton-II constraint
    F* = m a*

latent-state / excitation model
    what true accelerations or forces were plausible?

positive-mass domain

mass prior or experimentally licensed mass range
```

This is not arbitrary statistical decoration.

Without the excitation model, the same measured zero-zero pair cannot tell us whether:

- a nonzero force was attempted and acceleration failed to appear;
- a nonzero acceleration was expected and force failed to appear;
- neither quantity was expected to differ from zero.

Those are physically different experiments.

---

# 18. The role of the prior on mass

A posterior median requires a probability distribution over $m$, so a prior or experimentally permitted mass range is unavoidable.

This needs care.

For example,

```text
uniform prior in m
```

and

```text
uniform prior in log(m)
```

are not the same assumption.

Low-information measurements can therefore produce prior-sensitive medians.

That is not necessarily a defect. It is an honest statement that the experiment itself did not determine the mass sharply.

A strong experiment overwhelms reasonable prior differences.

A weak experiment does not.

Thus the estimator should report both:

```text
posterior median
posterior interval
```

and, where relevant,

```text
sensitivity to the chosen prior/excitation model
```

rather than disguising a weak determination as a precise object property.

---

# 19. Possible result type

The operation could therefore be represented as

$$
\operatorname{MassMedianEstimate}:
$$

```text
MeasuredForceData
× MeasuredAccelerationData
× ForceUncertaintyModel
× AccelerationUncertaintyModel
× ExperimentalExcitationModel
× PositiveMassPrior
    ->
MassPosteriorResult
```

where the result may contain

```text
median mass
credible interval
posterior distribution
residual/model-compatibility statistic
measurement policy
provenance
```

For example:

```text
MassPosteriorResult:
    median = 5.03 kg
    68% interval = [4.91, 5.16] kg
    model compatibility = good
```

or:

```text
MassPosteriorResult:
    median = prior-sensitive
    95% lower bound = 240 kg
    upper tail = poorly constrained
    model compatibility = good
```

or:

```text
MassPosteriorResult:
    mass distribution = broad
    model compatibility = poor
    cause = strongly anti-aligned measurements
```

The result does not need to pretend that every experiment produces one equally meaningful point number.

---

# 20. Comparison with the previous estimator

| Feature | Best-fit / profile estimator | Posterior-median estimator |
|---|---|---|
| Uses measurement uncertainties | Yes | Yes |
| Enforces Newton II on latent vectors | Yes | Yes |
| Uses one best latent pair | **Yes** | No |
| Uses all plausible latent states | No | **Yes** |
| Requires excitation/latent-state model | Not necessarily | **Yes** |
| Requires mass prior/range | No | **Yes** |
| Naturally restricted to $m>0$ | Can be | **Yes** |
| Handles long tails robustly | Poorly if summarized by mean/root | **Median is robust** |
| Can return intervals | Yes, with extra analysis | **Naturally** |
| Literal $0/0$ automatically flat? | Yes under profile residual | **Not necessarily** |
| Anti-aligned data | Can produce misleading stationary root | **Produces broad/low-compatibility mass distribution** |
| Strong-data limit | Exact Newton-II ratio | Exact Newton-II ratio |

---

# 21. The conceptual difference in one picture

The old estimator asks:

```text
For each m:

    find the single best hidden pair
        (F*, a*)

    and ask:
        how close can I get to the measurements?
```

The median estimator asks:

```text
For each m:

    consider every plausible hidden pair
        (F*, a*)

    weight each by:
        measurement likelihood
        experimental plausibility

    add all of that support

Then:

    compare all positive masses
    construct the mass distribution
    take its 50% point
```

The difference is therefore:

$$
\boxed{
\text{best compatible point}
\quad\text{versus}\quad
\text{total compatible probability}.
}
$$

---

# 22. Main conclusion

The median-based estimator does not define mass by dividing two measured vectors, nor by choosing the single corrected pair with the smallest residual.

It constructs a distribution over positive inertial mass from:

$$
\text{Newton II}
+
\text{measured force}
+
\text{measured acceleration}
+
\text{uncertainties}
+
\text{experimental excitation model}
+
\text{mass prior}.
$$

It then reports the halfway point:

$$
\boxed{
P(m\leq\widehat m_{50}\mid D)=\frac12.
}
$$

The method therefore gives different answers in different informational regimes:

```text
strong aligned signal
    -> narrow mass distribution near F/a

known force + unresolved acceleration
    -> evidence toward large mass / lower bound

known acceleration + unresolved force
    -> evidence toward small mass / upper region suppressed

zero-zero measurement
    -> answer depends on what the experiment could have produced

strongly anti-aligned measurements
    -> broad estimate and/or model incompatibility

many good trials
    -> sharply concentrated common mass
```

The essential change is that uncertainty is no longer treated merely as a penalty attached to one optimal correction.

It becomes part of the information from which the mass distribution itself is constructed.
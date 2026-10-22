# The Disagreement Rule for Inertial-Mass Estimation

## Why $m=\mathrm{undefined}$ is the correct end of the $0/0$ trace, but not of the $F/0$ trace

**Status:** working argument note  
**Project:** VD Newton / inertial-mass estimation  
**Purpose:** state the rule that separates **disagreement with Newton II** from **agreement without mass identification**, and derive the required behavior of the inertial-mass estimator from that distinction.

---

## 1. Thesis

Newton II requires, for one material object in one relevant inertial-frame-qualified account,

$$
\mathbf F=m\mathbf a,
\qquad m>0.
$$

When force and acceleration are used to determine mass, two questions must be asked in order:

1. **Agreement:** Is there at least one positive scalar $m$ for which the supplied force and acceleration satisfy Newton II?
2. **Identification:** If so, do those values select one unique positive scalar?

These questions generate two fundamentally different failures:

```text
no positive mass fits
    -> disagreement with Newton II

positive masses fit, but no unique one is selected
    -> agreement without identification
```

The estimator is required to handle the first problem on the empirical branch. It is not required to erase the second.

Therefore:

```text
measured disagreement with Newton II
    -> invoke the estimator

exact zero-zero agreement without identification
    -> allow m = undefined to end the trace
```

The estimator is not a general machine for preventing every undefined result. Its core Newtonian role is to prevent ordinary measurement disagreement from being mistaken for the nonexistence or undefinedness of inertial mass.

---

## 2. The compatibility set

For a supplied force–acceleration pair, define the set of positive masses compatible with Newton II:

$$
\mathcal C(\mathbf F,\mathbf a)
=
\{m>0:\mathbf F=m\mathbf a\}.
$$

For vector force and acceleration, exactly three structures can occur.

### 2.1 One compatible mass

If $\mathbf F$ and $\mathbf a$ are nonzero and codirectional, then

$$
\mathcal C(\mathbf F,\mathbf a)=\{m_0\}
$$

for one positive scalar

$$
m_0=\frac{|\mathbf F|}{|\mathbf a|}.
$$

The pair agrees with Newton II and identifies one mass.

### 2.2 Every positive mass is compatible

If

$$
\mathbf F=\mathbf0,
\qquad
\mathbf a=\mathbf0,
$$

then

$$
\mathbf0=m\mathbf0
$$

for every $m>0$. Therefore

$$
\mathcal C(\mathbf0,\mathbf0)=(0,\infty).
$$

The pair agrees with Newton II, but it does not identify one mass.

### 2.3 No positive mass is compatible

If exactly one vector is zero, or if the two nonzero vectors are not codirectional, then

$$
\mathcal C(\mathbf F,\mathbf a)=\varnothing.
$$

The pair disagrees with Newton II.

This trichotomy is the foundation of the estimator rule.

---

## 3. The rule: repair disagreement, preserve non-identification

The required behavior can be stated compactly:

> **When the supplied empirical values disagree with Newton II, the mass-determination operation must estimate an underlying compatible relation rather than allow the raw disagreement to become the mass result. When the supplied exact values already agree with Newton II but do not identify one mass, the operation must preserve that non-identification rather than manufacture a mass estimate merely to avoid an undefined output.**

This gives two different meanings to an otherwise similar-looking arithmetic failure.

```text
F/0 is undefined because the pair does not satisfy the law.
0/0 is undefined even though the pair does satisfy the law.
```

The first undefined result appears **before compatibility has been secured**.

The second appears **after compatibility has been secured but uniqueness has failed**.

Only the first calls for Newton-II reconciliation on the empirical branch.

---

## 4. Why $F/0$ must not terminate as $m=\mathrm{undefined}$

Use $F/0$ as shorthand for

$$
\mathbf F\neq\mathbf0,
\qquad
\mathbf a=\mathbf0.
$$

The exact law requires

$$
\mathbf F=m\mathbf0=\mathbf0.
$$

But the supplied force is nonzero. Thus

$$
\nexists m>0:\mathbf F=m\mathbf a.
$$

The pair is not merely unhelpful. It disagrees with Newton II.

If one naively continues with

$$
m=\frac{\mathbf F}{\mathbf0},
$$

then vector division fails and the trace may appear to produce

$$
m=\mathrm{undefined}.
$$

But on an empirical trial, that is the wrong end of the trace.

Why? Because it would convert this chain

```text
independently measured force and acceleration disagree with the exact law
```

into this conclusion

```text
the material object's inertial mass is undefined
```

without any intermediate argument.

That is precisely the inference ordinary Newtonian practice refuses to make. A force channel returning a nonzero value while an acceleration channel returns zero is ordinarily treated as a measurement-resolution, calibration, model, address, interval, or uncertainty problem—not as an immediate demonstration that the object has no inertial mass.

The empirical trace must therefore continue:

```text
DEMAND: inertial-mass(A)

INPUTS:
    measured net-force F~ != 0
    measured inertial-acceleration a~ = 0

CHECK:
    no positive scalar exactly relates the reported central values

RESULT OF RAW CHECK:
    measured disagreement with Newton II

DO NOT RETURN:
    m[A] = undefined

INSTEAD:
    invoke declared joint estimation policy
    estimate latent F* and a*
    require F* = m a*, m > 0
    retain corrections, uncertainty, and goodness of fit

RETURN:
    estimated mass result
    or model-incompatibility result if credible reconciliation fails
```

The estimator does not merely replace an inconvenient division. It prevents the measurement disagreement from being mis-typed as the mass value.

### Exact-input rider

If $\mathbf F\neq0$ and $\mathbf a=0$ are explicitly stipulated as exact Newtonian quantities, then no estimator should rescue them. The correct result is an exact contradiction:

```text
NoPositiveMassFits
```

The estimator is licensed because empirical outputs are fallible estimates of latent Newtonian quantities. It is not licensed to rewrite explicitly exact premises.

---

## 5. Why $0/0$ should terminate as $m=\mathrm{undefined}$

Now take

$$
\mathbf F=\mathbf0,
\qquad
\mathbf a=\mathbf0.
$$

The law gives

$$
\mathbf0=m\mathbf0
$$

for every positive mass.

There is no disagreement to repair.

The compatibility check does not merely pass for one value. It passes for all values in the positive-mass domain:

$$
\mathcal C(\mathbf0,\mathbf0)=(0,\infty).
$$

The mass-determination route then attempts to identify which member of this set is the object's mass. The scalar readout reaches

$$
m=\frac00.
$$

This is undefined, and here the undefined result is correct.

It expresses exactly what the law says about the trial:

```text
Newton II is satisfied.
No unique positive mass is selected.
Therefore this route returns no mass value.
```

The trace may therefore end:

```text
DEMAND: inertial-mass(A)

INPUTS:
    exact net-force = 0
    exact inertial-acceleration = 0

CHECK AGREEMENT:
    0 = m0 for every m > 0
    PASS

CHECK IDENTIFICATION:
    is exactly one m selected?
    NO

READOUT:
    0/0 = undefined

RETURN:
    m[A] from this trial = undefined

HALT REASON:
    law-compatible non-identification
```

Nothing has gone wrong with Newton II. The trace has reached an information boundary.

In this scoped sense,

$$
\boxed{m=\mathrm{undefined}}
$$

is the right answer: the demanded mass value is undefined **from this force–acceleration trial**.

It should not be read as the stronger claim that the material object cannot possess inertial mass in any account. It means that this route does not return one.

---

## 6. The same arithmetic word hides two different trace causes

Both $F/0$ and $0/0$ may be labelled “undefined” by elementary arithmetic. But the VD trace must preserve why the undefined result occurred.

| Case | Newton-II agreement? | Positive masses fitting | Why the quotient fails | Correct trace behavior |
|---|---:|---:|---|---|
| $\mathbf F\neq0,\ \mathbf a=0$ | No | None | The supplied values contradict the law | Empirical: estimate; exact: contradiction |
| $\mathbf F=0,\ \mathbf a\neq0$ | No | None | The supplied values contradict the law | Empirical: estimate; exact: contradiction |
| Nonzero non-codirectional pair | No | None | No positive scalar relates the vectors | Empirical: estimate; exact: contradiction |
| $\mathbf F=0,\ \mathbf a=0$ | Yes | Every $m>0$ | The law does not select one scalar | Return $m=\mathrm{undefined}$ from the trial |
| Nonzero codirectional pair | Yes | Exactly one | No failure | Return the unique mass |

The relevant distinction is therefore not:

```text
defined quotient
versus
undefined quotient
```

It is:

```text
law disagreement
versus
law agreement without identification
```

Arithmetic comes after that classification.

---

## 7. What the rule requires of the estimator

The disagreement rule constrains the estimator's behavior.

### 7.1 The estimator must be triggered by empirical law disagreement

If the independently measured central vectors do not admit a positive scalar relation, the estimator must not allow the raw codirectionality or zero-denominator error to stand as the mass result.

It must instead infer a latent pair

$$
\mathbf F^*,
\qquad
\mathbf a^*,
$$

and a positive scalar $m$ satisfying

$$
\mathbf F^*=m\mathbf a^*.
$$

### 7.2 The estimator must preserve the original disagreement

Because the latent pair satisfies Newton II by construction, the estimator must retain the difference between

$$
(\widetilde{\mathbf F},\widetilde{\mathbf a})
$$

and

$$
(\mathbf F^*,\mathbf a^*).
$$

Otherwise it would install the law and then report the installed agreement as evidence.

A suitable result therefore retains:

```text
estimated mass
latent compatible force and acceleration
measurement corrections
uncertainty
residual or goodness of fit
estimation policy
provenance
```

### 7.3 The estimator must not manufacture a mass merely to eliminate $0/0$

When an exact zero-zero pair already satisfies Newton II for every positive mass, there is no compatibility failure for the estimator to repair.

A procedure whose only justification is reconciliation should therefore preserve the non-identifying result:

```text
m from this trial = undefined
```

Forcing a positive point estimate here would answer a new question by adding new assumptions. It would not be required by Newton II itself.

### 7.4 The estimator must preserve the exact-data boundary

The estimator may reconcile empirical results because those results represent uncertain access to latent quantities.

It must not reconcile explicitly exact incompatible inputs. Those produce contradiction, not estimation.

### 7.5 The estimator must be allowed to report failed reconciliation

The estimator must not make Newton II unfalsifiable. If the measured disagreement is too large to be credible under the declared uncertainty and measurement model, the proper output is not a confident mass. It is a poor-fit or model-incompatibility result.

Thus:

```text
ordinary disagreement within measurement standing
    -> estimated mass

persistent or excessive disagreement
    -> Newton-II/model compatibility warning
```

---

## 8. The resulting mass-determination operation

A compact representation is:

$$
\operatorname{MassFromForceAcceleration}(\mathbf F,\mathbf a)
=
\begin{cases}
\dfrac{|\mathbf F|}{|\mathbf a|},
& \mathcal C(\mathbf F,\mathbf a)=\{m_0\},\\[10pt]
\mathrm{undefined},
& \mathcal C(\mathbf F,\mathbf a)=(0,\infty),\\[6pt]
\operatorname{Estimator}(\widetilde{\mathbf F},\widetilde{\mathbf a},\Pi),
& \mathcal C(\widetilde{\mathbf F},\widetilde{\mathbf a})=\varnothing
  \text{ and the inputs are empirical},\\[6pt]
\mathrm{contradiction},
& \mathcal C(\mathbf F,\mathbf a)=\varnothing
  \text{ and the inputs are exact}.
\end{cases}
$$

Here $\Pi$ stands for the declared uncertainty model, reconciliation policy, and other licensed estimator inputs.

In result-type form:

```text
MassDeterminationResult :=
    ExactMass(m)
  | EstimatedMass(m, latent_pair, uncertainty, residual, policy, provenance)
  | UndefinedMassFromTrial
  | ExactNewtonIIContradiction
  | PoorModelCompatibility
```

The crucial result is that `UndefinedMassFromTrial` and `ExactNewtonIIContradiction` are not the same output.

---

## 9. Why this supports $m=\operatorname{Estimator}(F,a)$ on the empirical disagreement branch

Suppose the empirical mass operation were simply

$$
m=\frac{\widetilde{\mathbf F}}{\widetilde{\mathbf a}}.
$$

Most independent vector measurements would fail exact codirectionality. The operation would therefore return a codirectionality or division error on almost every trial.

If that error were allowed to carry directly into the mass headword, almost every experiment would conclude that inertial mass is absent or undefined, and almost every experiment would disagree with Newton II.

But that is not how Newtonian experimental reasoning treats ordinary measurement disagreement. It treats the reported vectors as uncertain access to underlying quantities expected to satisfy

$$
\mathbf F^*=m\mathbf a^*.
$$

Therefore, on the empirical disagreement branch,

$$
\boxed{
\widehat m
=
\operatorname{Estimator}
(\widetilde{\mathbf F},\widetilde{\mathbf a},\Pi)
}
$$

is not optional numerical decoration. It is the operation required to prevent measurement disagreement from being mistaken for the mass result.

The zero-zero case shows the limit of that conclusion. The estimator is constitutive where empirical values disagree with the law. It is not thereby licensed to replace every law-compatible non-identifying trace with a point estimate.

---

## 10. Boundary: a richer statistical question is possible

This note derives the behavior of the **Newton-II reconciliation estimator**: the estimator required because independently measured force and acceleration ordinarily fail to form an exact law-compatible pair.

A richer probabilistic analysis may ask a further question:

> Given the complete uncertainty, excitation, and prior model, how much evidential support does a measured zero-zero outcome assign to different masses?

Such a model may extract information from a measured zero-zero event. That is a separate, stronger inferential proposal. It is not needed to repair disagreement with Newton II, because zero-zero contains no such disagreement.

Therefore the core argument should keep the two claims separate:

```text
Core reconciliation claim:
    disagreement requires estimation

Stronger probabilistic claim:
    even a law-compatible null observation may carry model-dependent information
```

The first claim determines why $F/0$ cannot end at `m = undefined` on the empirical branch.

The second is optional downstream machinery and does not change why exact $0/0$ may correctly terminate there.

---

## 11. Main conclusion

The estimator's behavior is governed by whether the force–acceleration inputs **disagree with Newton II** or merely **fail to identify a unique mass while agreeing with it**.

For $\mathbf F\neq0$ and $\mathbf a=0$:

$$
\nexists m>0:\mathbf F=m\mathbf a.
$$

The raw values disagree with Newton II. On an empirical trace, `m = undefined` is therefore premature and misleading. The trace must invoke estimation, preserve the disagreement as a residual, and return either a supported estimated mass or a model-compatibility failure.

For $\mathbf F=0$ and $\mathbf a=0$:

$$
\forall m>0:\mathbf0=m\mathbf0.
$$

The values agree with Newton II completely, but select no unique mass. The inversion reaches

$$
\frac00=\mathrm{undefined},
$$

and that undefined result carries through correctly:

$$
\boxed{
m\text{ from this trial}=\mathrm{undefined}.
}
$$

Thus the correct governing rule is:

> **Do not repair every undefined quotient. Repair empirical disagreement with the law. Preserve undefinedness when the law is satisfied but the demanded value is not identified.**

---

## Source basis

This note develops a new synthesis from the current project arguments in:

- `why_inertial_mass_must_be_estimated.md`
- `why_reconciliation_cannot_be_upstream_of_E10-informing(1).md`
- `why_the_mass_entry_must_be_long_joint_measurement_d3(1).md`
- `median_based_estimator_inertial_mass-informing.md`

The disagreement/non-identification rule is the new organising contribution of this note.

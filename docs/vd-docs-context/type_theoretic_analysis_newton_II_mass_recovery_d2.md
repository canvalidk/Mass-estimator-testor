# Type-Theoretic Analysis of Newton-II Mass Recovery

## Purpose

Newton II relates three addressed quantities for one material object:

$$
\mathbf F_{\mathrm{net}} = m\mathbf a,
\qquad
m>0.
$$

The three directions of use are not type-theoretically identical.

Two directions are total scalar–vector operations. The mass direction takes two vectors, but whether a mass exists depends on a relation between their values that is not guaranteed by their ordinary types.

---

## 1. Net-force determination

Given a positive inertial mass and an inertial-acceleration vector,

$$
\operatorname{NetForce}:
\mathrm{PositiveMass}
\times
\mathrm{AccelerationVector}
\to
\mathrm{ForceVector},
$$

$$
(m,\mathbf a)
\mapsto
m\mathbf a.
$$

Every correctly typed input pair produces a force vector. This is a total operation.

---

## 2. Inertial-acceleration determination

Given a net-force vector and a positive inertial mass,

$$
\operatorname{Acceleration}:
\mathrm{ForceVector}
\times
\mathrm{PositiveMass}
\to
\mathrm{AccelerationVector},
$$

$$
(\mathbf F_{\mathrm{net}},m)
\mapsto
\frac{1}{m}\mathbf F_{\mathrm{net}}.
$$

Every correctly typed input pair produces an acceleration vector. This is also a total operation.

---

## 3. The naïve mass signature is too broad

A naïve signature would be:

$$
\operatorname{Mass}:
\mathrm{ForceVector}
\times
\mathrm{AccelerationVector}
\to
\mathrm{PositiveMass}.
$$

This cannot be a total function.

The input types establish that the two supplied quantities are vectors of the correct physical kinds. They do not establish that the particular vector values satisfy

$$
\mathbf F_{\mathrm{net}}=m\mathbf a
$$

for some positive scalar $m$.

The missing requirement is relational:

$$
\mathbf a\neq\mathbf 0
\quad\text{and}\quad
\exists m>0:
\mathbf F_{\mathrm{net}}=m\mathbf a.
$$

For positive inertial mass, the nonzero force and acceleration vectors must therefore be codirectional.

---

## 4. Mass recovery as a checked function

The exact mass operation is better represented as:

$$
\operatorname{MassExact}:
\mathrm{ForceVector}
\times
\mathrm{AccelerationVector}
\to
\mathrm{MassRecoveryResult}.
$$

The operation first checks the relation between the vector values.

Possible results include:

```text
ExactMass(m)
ZeroZeroUnderdetermined
ZeroNonzeroInconsistency
NoncodirectionalInconsistency
```

The procedure is:

```text
receive force and acceleration vectors
    -> check zero/nonzero status
    -> check whether the nonzero vectors are codirectional
    -> if compatible, recover the positive scalar
    -> otherwise return the relevant non-mass result
```

For a compatible nonzero pair,

$$
m
=
\frac{F_{\parallel}}{a_{\parallel}}
=
\frac{|\mathbf F_{\mathrm{net}}|}{|\mathbf a|}.
$$

The magnitude ratio is the final scalar readout after compatibility has been established.

---

## 5. Refined-type formulation

The compatibility check may instead be represented by a refined input type:

$$
\mathrm{CodirectionalNonzeroPair}
=
\left\{
(\mathbf F,\mathbf a)
\;\middle|\;
\mathbf a\neq\mathbf 0
\land
\exists m>0:\mathbf F=m\mathbf a
\right\}.
$$

The trace then separates into two functions:

$$
\operatorname{CheckMassPair}:
\mathrm{ForceVector}
\times
\mathrm{AccelerationVector}
\to
\operatorname{Result}
\left(
\mathrm{CodirectionalNonzeroPair},
\mathrm{MassPairFailure}
\right),
$$

followed by

$$
\operatorname{ReadMass}:
\mathrm{CodirectionalNonzeroPair}
\to
\mathrm{PositiveMass}.
$$

In this formulation, the mass readout is total on the refined domain. The additional checks are performed while establishing that the supplied pair belongs to that domain.

---

## 6. Partial-function formulation

The same structure may be written compactly as a partial function:

$$
\operatorname{MassExact}:
\mathrm{ForceVector}
\times
\mathrm{AccelerationVector}
\rightharpoonup
\mathrm{PositiveMass}.
$$

This states that the function is defined only for compatible nonzero vector pairs.

For the lesson and the VD trace, the result-valued formulation is more informative because it preserves the reason the operation did not return a mass.

---

## 7. Measured-vector extension

For measured vectors, exact codirectionality will generally not hold.

The two measured vectors alone do not determine how the discrepancy should be corrected. A mass-estimation operation therefore requires an additional estimation policy:

$$
\operatorname{MassEstimate}:
\mathrm{MeasuredForceVector}
\times
\mathrm{MeasuredAccelerationVector}
\times
\mathrm{EstimationPolicy}
\to
\mathrm{MassEstimateResult}.
$$

Newton II supplies the exact relation that the underlying vectors must satisfy:

$$
\mathbf F_{\mathrm{net}}^{*}
=
m\mathbf a^{*}.
$$

The estimation policy supplies the rule for selecting an aligned underlying pair from the imperfect measurements.

This measured-vector extension is separate from the exact checked mass-recovery function.

---

## 8. Main conclusion

```text
PositiveMass × AccelerationVector
    -> ForceVector

ForceVector × PositiveMass
    -> AccelerationVector

ForceVector × AccelerationVector
    -> checked MassRecoveryResult
```

The first two operations are total because their input types are sufficient to define the output.

The mass operation is checked because the ordinary vector types do not encode the value-level relation required for a positive inertial mass to exist.

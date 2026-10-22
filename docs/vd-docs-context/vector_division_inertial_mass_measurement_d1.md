# Welp… Vector Division Is Not Defined

## How do we actually infer inertial mass from $\mathbf F=m\mathbf a$?

At school, Newton’s second law looks almost too easy:

$$
\mathbf F=m\mathbf a.
$$

Rearrange it:

$$
m=\frac{F}{a}.
$$

Then force and acceleration are introduced as vectors, and a small disaster appears:

> **Welp… vector division is not defined.**

If the force and acceleration vectors are exactly parallel and point the same way, the problem can be reduced to one dimension. Choose the common direction and divide the corresponding scalar components:

$$
m=\frac{F_\parallel}{a_\parallel}
 =\frac{|\mathbf F|}{|\mathbf a|}.
$$

But real measurements will almost never produce exactly parallel vectors. Suppose an experiment reports

$$
\mathbf F=(10,0)\ \mathrm N,
\qquad
\mathbf a=(4,3)\ \mathrm{m\,s^{-2}}.
$$

What mass has been measured?

Several plausible-looking calculations give different answers:

$$
\frac{|\mathbf F|}{|\mathbf a|}=2.0\ \mathrm{kg},
$$

$$
\frac{\mathbf a\cdot\mathbf F}{\mathbf a\cdot\mathbf a}
=1.6\ \mathrm{kg},
$$

$$
\frac{\mathbf F\cdot\mathbf F}{\mathbf F\cdot\mathbf a}
=2.5\ \mathrm{kg}.
$$

A componentwise average in the displayed coordinate system gives yet another value:

$$
\frac12\left(\frac{F_x}{a_x}+\frac{F_y}{a_y}\right)
=\frac12\left(\frac{10}{4}+\frac03\right)
=1.25\ \mathrm{kg}.
$$

These are not different notations for the same operation. They are different estimators, based on different assumptions about what went wrong when the vectors failed to align.

The question is therefore not merely algebraic. It is a measurement problem.

---

## 1. The exact law and the measurement problem are different

Newton’s second law says that the **true** net force and **true** acceleration of a Newtonian material object satisfy

$$
\mathbf F^*=m\mathbf a^*,
$$

where $m$ is one positive scalar.

The stars matter. Instruments do not return the inaccessible exact values $\mathbf F^*$ and $\mathbf a^*$. They return measured values such as

$$
\widetilde{\mathbf F},
\qquad
\widetilde{\mathbf a},
$$

which are affected by calibration errors, sensor noise, trajectory reconstruction, numerical differentiation, frame errors, timing errors, omitted force contributions, and other imperfections.

If two **exactly specified** vectors are unaligned, then there is no scalar $m$ satisfying

$$
\mathbf F=m\mathbf a.
$$

That is not negotiable. Multiplication by a scalar can change a vector’s length and, if the scalar is negative, reverse its direction. It cannot rotate the vector through an arbitrary angle.

But measured vectors are not required to lie exactly on the model. A fitted parameter can exist even when no individual measurement satisfies the model perfectly. Saying that mass cannot be estimated because the raw vectors are not exactly parallel would be like saying that the slope of a fitted straight line cannot exist because the experimental points are not perfectly collinear.

The real task is:

> Find the positive scalar $m$ for which plausible corrections to the measured force and acceleration make the data agree with $\mathbf F^*=m\mathbf a^*$ as well as possible.

The phrase **“as well as possible”** must then be defined.

---

## 2. Why the obvious answers are not universally correct

### 2.1 Divide the magnitudes

The most tempting answer is

$$
\widehat m_{\mathrm{mag}}
=\frac{|\widetilde{\mathbf F}|}{|\widetilde{\mathbf a}|}.
$$

This is rotationally invariant and always positive when both vectors are nonzero. It is also the correct ratio of the two measured magnitudes.

Its problem is that it discards all directional information.

If

$$
\widetilde{\mathbf F}=(10,0),
\qquad
\widetilde{\mathbf a}=(0,5),
$$

then the force and acceleration are perpendicular, yet the formula returns

$$
\widehat m_{\mathrm{mag}}=2\ \mathrm{kg}.
$$

If they are antiparallel, it returns the same positive value. Thus the formula estimates a **magnitude-response ratio**, but it is not automatically the best estimate of the scalar coefficient in the full vector law.

It becomes appropriate only when direction has deliberately been excluded from the measurand, or under a special measurement model that makes this ratio the preferred estimator.

### 2.2 Divide the components and average

Another idea is

$$
m_x=\frac{F_x}{a_x},
\qquad
m_y=\frac{F_y}{a_y},
\qquad
m_z=\frac{F_z}{a_z},
$$

followed by an average.

For one vector pair, this is not acceptable as a fundamental procedure.

First, it becomes unstable whenever an acceleration component is close to zero. Second, it depends on the arbitrary choice of coordinate axes.

For the example

$$
\mathbf F=(10,0),
\qquad
\mathbf a=(4,3),
$$

the average of the two component ratios is $1.25\ \mathrm{kg}$. Rotate the axes by $45^\circ$, without changing the physical vectors, and the component ratios become approximately

$$
1.43\ \mathrm{kg}
\qquad\text{and}\qquad
10.0\ \mathrm{kg},
$$

whose average is about

$$
5.71\ \mathrm{kg}.
$$

A physical mass estimate cannot change merely because the analyst rotated the graph paper.

Taking the Euclidean norm of the component ratios is no repair. Even in a perfect two-dimensional trial, where both component ratios equal the true mass $m$, it returns

$$
\sqrt{m^2+m^2}=\sqrt2\,m.
$$

In three dimensions it would return $\sqrt3\,m$.

There is an important exception: **three genuinely separate experiments performed along three physically selected orthogonal directions** may each produce a legitimate one-dimensional mass estimate. Those estimates can then be combined statistically. A simple arithmetic mean is justified only when they are independent and equally precise; otherwise a weighted mean or, preferably, one joint fit should be used. That is a multi-trial estimation procedure, not componentwise division of one vector observation.

### 2.3 Project force onto acceleration

One can choose the scalar that minimizes

$$
\left\|\widetilde{\mathbf F}-m\widetilde{\mathbf a}\right\|^2.
$$

The solution is

$$
\widehat m_{F\mid a}
=
\frac{\widetilde{\mathbf a}\cdot\widetilde{\mathbf F}}
{\widetilde{\mathbf a}\cdot\widetilde{\mathbf a}}.
$$

This is a perfectly respectable least-squares estimator **if acceleration is treated as effectively exact and all relevant error is assigned to force**. It asks:

> Holding the measured acceleration fixed, which scalar makes $m\mathbf a$ come closest to the measured force?

It is not universally correct because measured acceleration is generally uncertain too.

### 2.4 Project acceleration onto force, then invert

One can instead fit

$$
\mathbf a\approx\beta\mathbf F,
\qquad
\beta=\frac1m.
$$

Minimizing acceleration residuals gives

$$
\widehat\beta
=
\frac{\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}}
{\widetilde{\mathbf F}\cdot\widetilde{\mathbf F}},
$$

and therefore

$$
\widehat m_{a\mid F}
=
\frac{\widetilde{\mathbf F}\cdot\widetilde{\mathbf F}}
{\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}}.
$$

This is appropriate **if force is treated as effectively exact and all relevant error is assigned to acceleration**.

It differs from the previous projection because ordinary least squares is asymmetric: regressing $F$ on $a$ is not the same problem as regressing $a$ on $F$ and taking the reciprocal.

The fact that force may have been deliberately applied first does not make it error-free. Causal order does not determine regression direction. Instrument uncertainty does.

### 2.5 Take the method with the smallest ordinary residual

This still does not solve the problem. A force residual is measured in newtons, while an acceleration residual is measured in metres per second squared. They cannot be added or compared as bare Euclidean lengths without introducing a scale that says how much one newton of correction should cost relative to one metre per second squared of correction.

The measurement uncertainties provide that scale.

---

## 3. There is no unique estimate “before uncertainty”

This is the central result.

Suppose an estimator is built only from the two measured vectors, is unchanged by a rotation of the coordinate axes, and has the dimensions of mass. Then it must have the general form

$$
\widehat m
=
\frac{|\widetilde{\mathbf F}|}{|\widetilde{\mathbf a}|}
\,g(\cos\theta),
$$

where $\theta$ is the angle between the vectors and $g$ is some dimensionless function.

Exact Newton II requires only

$$
g(1)=1,
$$

because all sensible methods must return $|\mathbf F|/|\mathbf a|$ when the vectors point in exactly the same direction.

But the law does not specify $g(c)$ for $c\neq1$.

The candidate formulas correspond to different choices:

$$
g(c)=1
\quad\text{(magnitude ratio)},
$$

$$
g(c)=c
\quad\text{(force residuals only)},
$$

$$
g(c)=\frac1c
\quad\text{(acceleration residuals only)}.
$$

Infinitely many other choices are possible.

Therefore the measured vectors and the exact law do not determine a unique mass estimate once the measurements fall away from the exact-law surface.

Choosing how to extend the definition away from perfect alignment is already choosing which deviations are considered more plausible. In other words:

> **The point estimator and the uncertainty model cannot be separated.**

Uncertainty is not decorative error bars added after the “real” mass has already been calculated. For inconsistent measured inputs, it helps define what the best estimate is.

---

## 4. The standard answer: errors-in-variables fitting

The standard statistical framework is called **errors-in-variables regression**. Closely related names include **Deming regression**, **total least squares**, and **orthogonal-distance regression**.

The basic idea is simple:

1. Both force and acceleration measurements may be wrong.
2. Introduce unknown corrected values representing the quantities the instruments were trying to measure.
3. Require the corrected values to satisfy Newton’s law exactly.
4. Choose the mass that requires the smallest total correction, where every correction is measured relative to its own uncertainty.

For trial $i$, let

$$
\widetilde{\mathbf F}_i,
\qquad
\widetilde{\mathbf a}_i
$$

be the measured vectors. Let

$$
\mathbf a_i^*
$$

be the unknown corrected acceleration, so the corrected force is constrained to be

$$
\mathbf F_i^*=m\mathbf a_i^*.
$$

Let

$$
\Sigma_{F,i},
\qquad
\Sigma_{a,i}
$$

be the covariance matrices describing the measurement uncertainty of force and acceleration. Under independent Gaussian measurement errors, estimate the mass by solving

$$
\boxed{
\widehat m
=
\underset{m>0,\{\mathbf a_i^*\}}{\operatorname{argmin}}
\sum_i
\left[
(\widetilde{\mathbf F}_i-m\mathbf a_i^*)^{\!T}
\Sigma_{F,i}^{-1}
(\widetilde{\mathbf F}_i-m\mathbf a_i^*)
+
(\widetilde{\mathbf a}_i-\mathbf a_i^*)^{\!T}
\Sigma_{a,i}^{-1}
(\widetilde{\mathbf a}_i-\mathbf a_i^*)
\right].
}
$$

This formula looks formidable, but its meaning is straightforward:

- the first term charges for correcting the force measurement;
- the second term charges for correcting the acceleration measurement;
- a correction along a precise direction is expensive;
- a correction along a noisy direction is cheap;
- the selected mass is the positive scalar for which the total standardized correction is smallest.

The covariance matrices solve the units problem. The force and acceleration residuals are not compared as raw numbers; each is divided by the uncertainty scale appropriate to that measurement and direction. A large angular mismatch is therefore highly costly when both directions are known precisely, but may be quite ordinary when the directional uncertainties are large.

If force and acceleration errors are correlated, they are represented together by a full joint covariance matrix rather than by two independent matrices.

Under normally distributed errors, this uncertainty-weighted fit is a maximum-likelihood estimate. “Best” then has a precise meaning:

> It is the value of $m$ under which the measurements that actually occurred are most probable, given the calibrated error behaviour of the measurement procedures.

The same fit supplies more than a point estimate. It can also provide:

- an uncertainty interval for $m$;
- standardized residuals showing where the data disagree with the law;
- a test of whether one common scalar mass adequately describes all trials;
- evidence of offsets, omitted force contributions, frame errors, timing errors, or direction-dependent response.

For non-Gaussian measurement errors or strong outliers, the likelihood or loss function should be changed accordingly. The principle remains the same: the estimator must follow the measurement model.

---

## 5. The earlier formulas are limiting cases, not rival laws

The projection formulas reappear naturally from the general fit.

### Acceleration uncertainty is negligible

If $\Sigma_a$ is tiny compared with $\Sigma_F$, then acceleration is effectively fixed. The fit reduces to

$$
\widehat m
=
\frac{\mathbf a\cdot\mathbf F}{\mathbf a\cdot\mathbf a}.
$$

### Force uncertainty is negligible

If $\Sigma_F$ is tiny compared with $\Sigma_a$, then force is effectively fixed. The fit reduces to

$$
\widehat m
=
\frac{\mathbf F\cdot\mathbf F}{\mathbf F\cdot\mathbf a}.
$$

### Both are uncertain

The answer depends on their relative uncertainties and correlations. Neither projection is privileged in advance.

### Only magnitudes are being fitted

If the declared measurand is the ratio of response magnitudes, then

$$
\widehat m=\frac{|\mathbf F|}{|\mathbf a|}
$$

is exactly the quantity being estimated. But that is a different measurement statement from fitting the coefficient in the complete vector equation.

### Several orthogonal trials

If three independent directional trials return $m_x,m_y,m_z$ with equal uncertainty, their arithmetic mean is a sensible estimate of a common mass. If their uncertainties differ, use an inverse-variance weighted mean. Better still, fit the original force and acceleration measurements from all trials simultaneously.

Their spread is not merely annoying noise. It tests the claim that inertial response is adequately represented by one direction-independent scalar.

---

## 6. What if the vectors point in opposite directions?

Suppose

$$
\widetilde{\mathbf F}=(5,0),
\qquad
\widetilde{\mathbf a}=(-5,0).
$$

The magnitude ratio gives $1\ \mathrm{kg}$. The signed projection estimators give a negative result. Which conclusion is correct depends on the uncertainties.

- If both vectors are very precise, no plausible small correction makes them point the same way. The positive-scalar Newtonian model is strongly inconsistent with the trial.
- If both vectors have uncertainties comparable to their measured magnitudes, corrected vectors that point the same way may remain plausible. The mass estimate will then be broad and weak.

The correct response is not to take an absolute value automatically, nor to declare immediately that the object has no mass. It is to fit the positive-mass model under the stated uncertainties and report both the resulting mass interval and the remaining incompatibility.

A failed fit says that the measurements and model are not jointly adequate. Possible causes include a reversed sensor axis, a wrong coordinate transformation, an omitted force, a non-inertial frame, mismatched times, the wrong tracked object, underestimated uncertainty, or genuinely non-scalar response.

---

## 7. Why multiple trials matter

A single unaligned force–acceleration pair contains too little information to diagnose much. A serious mass experiment should normally use several trials:

$$
(\widetilde{\mathbf F}_1,\widetilde{\mathbf a}_1),
\ldots,
(\widetilde{\mathbf F}_n,\widetilde{\mathbf a}_n).
$$

The same positive scalar $m$ is fitted across all of them.

Good experimental design varies the magnitude of the applied force and, where relevant, its direction. This does three jobs:

1. estimates the common mass more precisely;
2. distinguishes random measurement scatter from systematic error;
3. tests whether one isotropic scalar response actually survives changes of direction.

In a real analysis it can also be useful to test

$$
\mathbf F=\mathbf b+m\mathbf a
$$

before forcing $\mathbf b=\mathbf0$. A nonzero fitted offset may reveal sensor zero error or an omitted approximately constant force. The physical Newtonian model still predicts zero intercept once net force and acceleration have been correctly constructed; allowing an intercept temporarily is a diagnostic, not a redefinition of the law.

---

## 8. What should students be taught?

Students should not be taught that vectors can simply be divided.

They should be taught two connected ideas at different levels.

### Idealized mechanics problems

When the problem supplies exact, aligned force and acceleration vectors, choose their common axis and use

$$
m=\frac{F_\parallel}{a_\parallel}.
$$

The equality of the other components is part of the law’s consistency, not a separate vector-division operation.

### Actual measurements

When force and acceleration are measured, mass is not obtained by a unique algebraic quotient from one imperfect pair. It is estimated as a common positive slope across the data, with uncertainty in both force and acceleration included in the fit.

A useful student-facing statement is:

> **Inertial mass is the common positive scalar that makes the measured force and acceleration trials agree with $\mathbf F=m\mathbf a$ after the smallest adjustments permitted by their measurement uncertainties.**

Or, more briefly:

> **Mass is the uncertainty-weighted common slope of force against acceleration.**

The full practical algorithm is:

1. identify one object, frame, and time qualification for every force–acceleration pair;
2. measure several force and acceleration trials;
3. determine the uncertainties and correlations of both measurements;
4. fit one positive scalar using errors-in-variables or orthogonal-distance regression;
5. report the mass estimate and its uncertainty;
6. inspect the residuals to decide whether one scalar mass is an adequate model.

This is not needless complication added to Newton’s law. It is what the simple equation becomes when it is used as an experimental measurement rule rather than as an exact textbook constraint.

---

## 9. The concise conclusion

For exact aligned vectors,

$$
m=\frac{|\mathbf F|}{|\mathbf a|}
$$

is harmless because every sensible estimator agrees.

For measured unaligned vectors, there is no unique mass value before a measurement-error model is specified. Magnitude division, the two projection formulas, and directional averaging answer different questions or assume different locations for the error.

The standard general answer is:

$$
\boxed{
\text{fit one positive }m\text{ to }\mathbf F^*=m\mathbf a^*
\text{ by minimizing uncertainty-weighted corrections to both measured vectors.}
}
$$

That is an errors-in-variables or orthogonal-distance-regression problem.

---

## Sources and further reading

1. **Joint Committee for Guides in Metrology, _Evaluation of measurement data — Guide to the expression of uncertainty in measurement_ (JCGM 100:2008).** The international metrology guide for constructing measurement models and propagating uncertainties and covariances.  
   <https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf>

2. **Paul T. Boggs and Janet E. Rogers, “Orthogonal Distance Regression.”** A NIST treatment of orthogonal-distance regression as maximum-likelihood estimation for measurement-error models with normally distributed errors.  
   <https://docs.scipy.org/doc/external/odr_ams.pdf>

3. **Thomas W. Bartel, Sara Stoudt, and Antonio M. Possolo, “Force Calibrations using Errors-in-Variables Regression and Monte Carlo Uncertainty Evaluations” (NIST, 2016).** A direct force-metrology application showing why uncertainty in both applied force and measured response belongs in the regression.  
   <https://www.nist.gov/publications/force-calibrations-using-errors-variables-regression-and-monte-carlo-uncertainty>

4. **Derek York, Norman M. Evensen, Margarita López Martínez, and Jonás De Basabe Delgado, “Unified equations for the slope, intercept, and standard errors of the best straight line,” _American Journal of Physics_ 72, 367–375 (2004).** A widely used treatment of straight-line fitting when both coordinates have varying and possibly correlated errors.  
   <https://doi.org/10.1119/1.1632486>

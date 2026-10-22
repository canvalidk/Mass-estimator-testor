# ours26 with a full noise covariance

**2026-09-27.** Record, filed at the user's request ("Go ahead", on Claude's
proposal to record the unrecorded 27 September results by thread), written by
Claude.

**Where the results come from.** The equations in §§2–5 come from a completion
message VD-Newton wrote on 27 September (not filed). Claude re-derived and
checked each one (§7). The point in §2 that the combination rule forces the
determinant to cancel is Claude's.

**Relations.** Extends `ours26_per_axis_noise_2026-09-26.md` (its §8 lists
correlated errors as owed). Binds nothing.

---

## 1. Setting

Stack a reading's measured components as $\mathbf y=(\tilde{\mathbf F},\tilde{\mathbf a})$,
with a known joint Gaussian covariance

$$
\Sigma=\begin{pmatrix}\Sigma_{FF}&\Sigma_{Fa}\\\Sigma_{aF}&\Sigma_{aa}\end{pmatrix}>0 .
$$

This allows different SDs per component, correlation between axes, and
correlation between the force and acceleration errors. ours26 is the case of
diagonal $\Sigma$.

## 2. The law

Write $B_m=\begin{pmatrix}mI\\I\end{pmatrix}$, so the true pair is $B_m\mathbf a^\ast$. For each
candidate mass define

$$
\mathbf r_m=\tilde{\mathbf F}-m\tilde{\mathbf a},\qquad
V_m=\Sigma_{FF}-m(\Sigma_{Fa}+\Sigma_{aF})+m^2\Sigma_{aa},\qquad
Q(m)=\mathbf r_m^\top V_m^{-1}\mathbf r_m .
$$

**Integrating the true acceleration out.** With a push weight flat in shape,
$w(m)\,d\mathbf a^\ast$, a reading's likelihood is

$$
w(m)\,\big|B_m^\top\Sigma^{-1}B_m\big|^{-1/2}\,e^{-Q(m)/2}.
$$

- **The determinant depends on $m$.** In a random 3D case with correlated
  channels it is 2.13, 0.580 and 0.0462 at $m$ = 0.3, 1.7 and 5.
- **The combination rule forces it to cancel.** Nothing attached to a reading
  may depend on $m$ (lemma T7), so $w(m)\propto|B_m^\top\Sigma^{-1}B_m|^{1/2}$.
  With diagonal $\Sigma$ this weight is exactly flat in each axis's noise-unit
  push, which is ours26's weighting. So the extension is not a free convention:
  given a flat shape, its $m$-dependence is forced.
- **Result:**

$$
\boxed{\;p(m)\propto\pi(m)\,e^{-Q(m)/2}.\;}
$$

**Not re-derived:** the flat *shape* itself when axes are correlated. Lemma T6
got it from axis splitting on independent axes; with correlated axes a reading's
axes are not independent readings, so that argument has to be restated.

## 3. The true acceleration at a fixed mass

Given $m$ and the readings, the true acceleration is Gaussian,
$\mathbf a^\ast\,\vert\,m\sim N(\mathbf b_m,K_m)$, with

$$
\mathbf b_m=\tilde{\mathbf a}-(\Sigma_{aF}-m\Sigma_{aa})V_m^{-1}\mathbf r_m,\qquad
K_m=\Sigma_{aa}-(\Sigma_{aF}-m\Sigma_{aa})V_m^{-1}(\Sigma_{Fa}-m\Sigma_{aa}).
$$

These are the flat-weight posterior mean $(B_m^\top\Sigma^{-1}B_m)^{-1}B_m^\top\Sigma^{-1}\mathbf y$
and covariance $(B_m^\top\Sigma^{-1}B_m)^{-1}$, written by conditioning. With
one noise pair, $\mathbf b_m=(m\tilde{\mathbf F}+s^2\tilde{\mathbf a})/(s^2+m^2)$,
ours26's value.

## 4. The readout

Within a candidate $\mathbf F^\ast=m\,\mathbf a^\ast$, so only the acceleration
weight is needed:

$$
\boxed{\;\widehat m=\frac{\int_0^\infty m\,A(m)\,\pi(m)\,e^{-Q(m)/2}\,dm}{\int_0^\infty A(m)\,\pi(m)\,e^{-Q(m)/2}\,dm},\qquad
A(m)=E\,\lVert N(\mathbf b_m,K_m)\rVert .\;}
$$

**Computing $A$ without sampling.** From
$\sqrt X=\frac1{2\sqrt\pi}\int_0^\infty(1-e^{-tX})\,t^{-3/2}\,dt$ and, for
$\mathbf z\sim N(\mathbf b,K)$,

$$
E\,e^{-t\lVert\mathbf z\rVert^2}=\det(I+2tK)^{-1/2}\exp\!\big(-t\,\mathbf b^\top(I+2tK)^{-1}\mathbf b\big),
$$

one integral in $t$ gives $A$. This is ours26's Laplace-transform method with a
full covariance.

**Any proper $\pi$ gives a mass:** $A(m)=O(1/m)$ as $m\to\infty$ and $A$ stays
bounded as $m\to0$ (`ours26_starting_law_2026-09-27.md` §2).

## 5. Combination

For readings with independent errors that share one covariance $\Sigma$, keep

$$
S=\sum_i\mathbf y_i\mathbf y_i^\top,\qquad n=\text{number of readings}.
$$

- The law: $\sum_iQ_i(m)=\operatorname{tr}\big(V_m^{-1}L_mSL_m^\top\big)$,
  $L_m=(I,\,-mI)$.
- The readout: each $\mathbf b_{m,i}=C_m\mathbf y_i$ with
  $C_m=(B_m^\top\Sigma^{-1}B_m)^{-1}B_m^\top\Sigma^{-1}$, so
  $\sum_i\mathbf b_{m,i}\mathbf b_{m,i}^\top=C_mSC_m^\top$. With the stacked
  covariance $K_m\otimes I_n$, that is all the stacked length's law needs.
- **Combining is adding** $(S,n)$. Groups with different covariances carry
  separate states, as ours26's noise pairs do. With one isotropic noise pair the
  state reduces to ours25b's $(P,Q,D,M)$.
- **Correlation between readings** needs their joint covariance; batches then
  cannot simply be added.

## 6. What is open

- The flat shape under correlated axes (§2).
- Force–acceleration correlated errors have never been run in the tester; this
  record is algebra and a spot check, not an implementation.
- ours26 §8's sketch (a correlation that is the same in every reading, with the
  channels sharing principal axes, reduces to ours26 by a rotation) is a special
  case of this form.
- Owed: a tester implementation, with worlds that have correlated errors.

## 7. Checks

Run by Claude in its sandbox on 27 September (not kept), for a random 3D case
with a full, force–acceleration-correlated $\Sigma$:

- the integrated exponent equals $Q(m)$ to $5\times10^{-15}$, and the
  determinant factor is as quoted in §2;
- $\mathbf b_m$ and $K_m$ match the flat-weight posterior to $10^{-15}$;
- $A$ by the $t$-integral, 3.53443, against 3.53483 from $2\times10^6$ samples;
- $mA(m)$ tends to a constant as $m\to\infty$ (§4);
- diagonal $\Sigma$ gives ours26's Beale sum and conditional mean exactly;
- with five readings, $\sum_iQ_i$ and $\sum_i\mathbf b_i\mathbf b_i^\top$ rebuilt
  from $S$ to $2\times10^{-15}$.

# Exact mass curves, combination, and the estimator's role

**Working contribution 11 — 17 September 2026.** User-requested consolidation
of the discussion about repeated readings, additivity, the exact \(J(h)\)
reduction, and combination without replacing uncertainty curves by error bars.
This is a local working account. It does not adopt a new fusion policy and has
not been submitted to VD-docs. Newton II remains an assumed law throughout.

The central conclusion is that **exact combination is possible once a joint
measurement model is specified**. What generally fails is representing the
answer as another member of the original single-reading family, or retaining
enough information in just the mass point and one spread parameter. Those are
different claims. The proofs below establish the narrower failures, not an
impossibility of exact inference or of every finite representation.

The exact \(J(h)\) formula was already derived in
[contribution 7, §3](operational_warrant.md#3-an-exact-reduction-that-reveals-how-the-mass-laws-differ).
It is reproduced and developed here so that the combination account is
self-contained. The distinction between a probability law, its scalar
readout, and its sampling calibration runs through all the results.

## 1. What the estimator returns, and what additivity means

Write a compatible latent pair as

$$
\mathbf F^*=f\mathbf u,\qquad \mathbf a^*=\alpha\mathbf u,
\qquad f,\alpha>0,\quad \|\mathbf u\|=1.
$$

The observations, their uncertainty model, and a chosen reference measure
produce a joint law \(P\) on these quantities. The original point readout is

$$
\widehat m=R(P)=\frac{E_P[f]}{E_P[\alpha]}.
\tag{1}
$$

There is also an induced mass variable \(M=f/\alpha\), with distribution
\(p(m)\). Equation (1) is generally **not** \(E_P[M]\). In the original
finite-data Gaussian construction, the latter diverges even though (1) is
finite. Quantiles and the spread of log mass remain available.

The composition argument concerns operations on quantities in an already
specified joint law:

$$
R(f_1+f_2,\alpha)=R(f_1,\alpha)+R(f_2,\alpha),
\tag{2}
$$

$$
\frac1{R(f,\alpha_1+\alpha_2)}
=\frac1{R(f,\alpha_1)}+\frac1{R(f,\alpha_2)}.
\tag{3}
$$

Equation (2) uses a common uncertain acceleration and additive force
magnitudes in the stipulated codirectional composition. Equation (3) uses a
common force magnitude and additive relative-acceleration magnitudes. Its
mechanical example is the reduced mass of two bodies subject to equal and
opposite forces; it does not say that their total material mass is their
reduced mass.

With the stated domain, law-dependence, calibration, and regularity premises,
these requirements select the ratio of means. The
[local uniqueness audit](README.md#2-a-stronger-uniqueness-theorem-with-substantial-premises)
also shows how some regularity follows on the full positive integrable
domain. The managed proofs are listed in §11.

Crucially, (2)–(3) do not construct \(P\), decide which latent variables are
shared between experiments, select a prior or reference measure, or establish
coverage of uncertainty intervals. They are not a rule for combining
separately fitted experimental results.

For example, partitioning one law into cases with probabilities \(w_k\)
gives, with \(F_k=E[f\mid k]\) and \(a_k=E[\alpha\mid k]\),

$$
\widehat m=\frac{\sum_k w_k F_k}{\sum_k w_k a_k}
=\frac{\sum_k w_k a_k(F_k/a_k)}{\sum_k w_k a_k}.
\tag{4}
$$

This is exact regrouping within one law. It does not justify applying the
same weights to independently constructed posteriors.

## 2. Repeating one pair and observing different pairs

If all measurements concern the same stable true vector pair, independent
Gaussian readings \(y_i=(\mathbf F_i,\mathbf a_i)\), with known covariance
matrices \(\Sigma_i\), have the exact sufficient reduction

$$
\Sigma_*^{-1}=\sum_i\Sigma_i^{-1},\qquad
y_*=\Sigma_*\sum_i\Sigma_i^{-1}y_i.
\tag{5}
$$

One can then apply the estimator to this refined pair using the chosen
reference. Equal covariances give the sample mean and covariance
\(\Sigma/N\). This is the simple scenario discussed in the conversation.
It still requires excitation: a truly zero force and acceleration do not
identify mass merely because they are measured more precisely.

For different true force–acceleration pairs, the shared object is instead
the mass. A typical model has a different unknown excitation vector
\(\mathbf v_i\) for each trial:

$$
\mathbf F_i^{\rm obs}=m\mathbf v_i+\epsilon_{F,i},\qquad
\mathbf a_i^{\rm obs}=\mathbf v_i+\epsilon_{a,i}.
\tag{6}
$$

Averaging these into one pair can cancel excitation or lose information.
The correct joint treatment depends on the nuisance vectors, dependence
between readings, and reference assumptions.

Nor is summing separately fitted means automatically equivalent to refining
the pair first. The explicit scalar positive-normal example in
[contribution 3](weighting_and_calibration.md) has true force 2,
acceleration 1, and unit channel errors. Separately fitting each noisy
reading and pooling the fitted moments tends to approximately
\(2.1428713/1.4147159=1.5147008\), rather than the true mass 2.
Combining the repeated-pair likelihoods first recovers the stable pair.

Inverse-variance weighting needs a separate qualification. For independent
Gaussian mass observations with known variances,

$$
\widehat m=\frac{\sum_i m_i/\tau_i^2}{\sum_i1/\tau_i^2},
\qquad \tau_*^2=\frac1{\sum_i1/\tau_i^2}
\tag{7}
$$

is exact for the Gaussian likelihood, and for an untruncated flat-prior
Gaussian posterior. It becomes an approximation when the actual non-Gaussian
mass curves are replaced by Gaussians. A positive-mass restriction or a
nonflat prior changes the posterior; correlation changes the weighting.
The original single-reading mass law has no finite ordinary mass variance
to insert into (7).

## 3. The exact J(h) result

### Assumptions and radial reduction

This closed form applies to **three-dimensional, independent isotropic
Gaussian channels**, uniform common direction, and flat positive magnitudes
\(df\,d\alpha\). The general correlated-covariance estimator is broader;
this particular reduction should not be applied to it unchanged.

Define

$$
U=\frac{\mathbf F^{\rm obs}}{\sigma_F},\qquad
V=\frac{\mathbf a^{\rm obs}}{\sigma_a},\qquad
s=\frac{\sigma_F}{\sigma_a},
$$

and introduce

$$
f=\sigma_F\rho\sin\theta,\qquad
\alpha=\sigma_a\rho\cos\theta,\qquad
\rho>0,\quad 0<\theta<\frac\pi2.
$$

Then \(m=s\tan\theta\), and the observed vectors enter through

$$
h(\theta)^2=\|\sin\theta\,U+\cos\theta\,V\|^2.
\tag{8}
$$

After discarding an observation-only factor, the likelihood is
\(\exp[-\rho^2/2+\rho\mathbf u\cdot(\sin\theta U+\cos\theta V)]\).
Its uniform direction average is

$$
e^{-\rho^2/2}\frac{\sinh(\rho h)}{\rho h}.
$$

The flat-magnitude Jacobian supplies one factor \(\rho\). Define

$$
J_k(h)=\int_0^\infty \rho^k e^{-\rho^2/2}
                 \frac{\sinh(\rho h)}{\rho h}\,d\rho.
\tag{9}
$$

The function called \(J(h)\) in the discussion is \(J_1(h)\):

$$
\boxed{J(h)=J_1(h)=\sqrt{\frac\pi2}\,
       \frac{e^{h^2/2}\operatorname{erf}(h/\sqrt2)}h},
\qquad J_1(0)=1.
\tag{10}
$$

For \(h>0\), cancel the radial factor in (9), express the hyperbolic sine
as two exponentials, and complete the squares:

$$
J_1(h)=\frac{e^{h^2/2}}h\int_0^h e^{-t^2/2}\,dt.
$$

This also gives a nonsingular representation and its small-\(h\) series:

$$
J_1(h)=\int_0^1 e^{h^2(1-u^2)/2}\,du
=1+\frac{h^2}3+\frac{h^4}{15}+O(h^6).
\tag{11}
$$

Related radial moments are exact too:

$$
J_2(h)=\sqrt{\frac\pi2}e^{h^2/2},\qquad
J_3(h)=(1+h^2)J_1(h)+1,\qquad
J_4(h)=(h^2+3)J_2(h).
\tag{12}
$$

For example, \(G(h)=hJ_1(h)\) satisfies \(G'=hG+1\), and
\(J_3=G''/h\); the values at zero follow by continuity. These identities
remove radial quadrature. Remaining one-dimensional integrals need not have
elementary antiderivatives to be exact mathematical definitions.

### The full mass curve, point, and conditional acceleration

Only four scalar data summaries are required for this isotropic mass curve:

$$
A=\|U\|^2,\quad B=\|V\|^2,\quad C=U\cdot V,\quad s,
\qquad A,B\geq0,\quad C^2\leq AB.
$$

With \(r=m/s\),

$$
h(m)^2=\frac{Ar^2+2Cr+B}{1+r^2},\qquad
Z_\theta=\int_0^{\pi/2}J_1(h(\theta))\,d\theta.
$$

The normalized density with respect to **\(dm\)** is

$$
\boxed{p(m)=\frac1{Z_\theta}\frac{s}{s^2+m^2}J_1(h(m))},
\qquad m>0.
\tag{13}
$$

The original point and the conditional acceleration mean are

$$
\widehat m=s\,
\frac{\int_0^{\pi/2}\sin\theta J_2(h(\theta))\,d\theta}
     {\int_0^{\pi/2}\cos\theta J_2(h(\theta))\,d\theta},
\tag{14}
$$

$$
\boxed{a(m):=E[\alpha\mid M=m]
=\sigma_a\cos\theta\frac{J_2(h(\theta))}{J_1(h(\theta))}},
\qquad \theta=\arctan(m/s).
\tag{15}
$$

The extra dimensional scale \(\sigma_a\) in (15) will matter when preserving
the ratio-of-means readout across trials.

Four parameters specify (13); they need not be recoverable uniquely from
that curve. If \(A=B=t,C=0\), then \(h^2=t\) is constant. Every such
dataset gives the same half-Cauchy mass curve and point \(s\), even though
its original data summaries differ.

The same formulas also preserve an earlier useful result: for equal
standardized lengths, \(h^2=t+C\operatorname{sech}\ell\), where
\(\ell=\log(m/s)\), so the mass law is symmetric in log mass and its point
is \(s\). Under the alternative tube reference the angular weight is
\(J_3\) rather than \(J_1\). Positive alignment makes that law narrower
in absolute log mass; negative alignment makes it wider; perpendicular
observations give identical curves. The proof and the distinct reference
assumptions remain in [contribution 7](operational_warrant.md).

## 4. What exact combination of mass curves looks like

### One explicit joint model

For trial \(i\), let the original law be

$$
dP_i=Z_i^{-1}L_i(f_i,\alpha_i,\mathbf u_i)
                         \,df_i\,d\alpha_i\,d\Omega_i.
$$

Transforming \(f_i=m_i\alpha_i\) gives
\(df_i\,d\alpha_i=\alpha_i\,dm_i\,d\alpha_i\). Consequently,

$$
p_i(m)=Z_i^{-1}\int \alpha_i
 L_i(m\alpha_i,\alpha_i,\mathbf u_i)\,d\alpha_i\,d\Omega_i.
\tag{16}
$$

Now declare independent measurement errors, a common mass, separate latent
pairs, and the joint reference

$$
d\Pi_N=dm\prod_{i=1}^N\alpha_i\,d\alpha_i\,d\Omega_i.
\tag{17}
$$

Integrating the nuisance quantities under (17) gives exactly

$$
\boxed{p_{1:N}(m)=
 \frac{\prod_{i=1}^Np_i(m)}{\int_0^\infty\prod_{i=1}^Np_i(t)\,dt}}.
\tag{18}
$$

The single-trial normalizers cancel. This operation is associative and
independent of processing order. The whole combined curve can be carried
forward and multiplied by the next curve. Retaining the factor list is an
exact finite representation for every finite number of trials, with storage
growing with that number.

Equation (17) is an explicit extension of the single-trial construction; it
is **not uniquely forced** by it or by readout additivity. One way to obtain
it is to start from independent single-trial laws and condition on vanishing
absolute mass differences, using a specified constant-width tolerance in
mass. Conditioning on equality of continuous quantities requires such a
prescription.

More generally, if individually fitted posteriors have proper priors
\(p_i(m)\propto\pi_i(m)\mathcal L_i(m)\), then independent evidence and
one desired batch prior give

$$
p_{1:N}(m)\propto\pi_{1:N}(m)
                        \prod_i\frac{p_i(m)}{\pi_i(m)}.
\tag{19}
$$

Multiplying posteriors without this correction can count the prior more than
once. The original flat-magnitude reference is improper, with infinite
nuisance marginal \(\int_0^\infty\alpha\,d\alpha\); it does not supply
a normalized mass prior by that integration. Equation (17), rather than an
unstated proper-prior interpretation, supplies the justification for (18).
The standard likelihood/posterior distinction is also set out in the
[PDG statistics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-statistics.pdf).

Shared calibration uncertainty or correlated errors break the factorization
used here. Individual marginal curves alone do not specify their dependence;
the joint uncertainty information must accompany them.

### A coordinate change versus a different combination rule

Let \(\ell=\log(m/m_0)\), and let \(g_i(\ell)=m p_i(m)\) be the density
with respect to \(d\ell\). Multiplying those densities and then converting
back to mass gives

$$
p_{1:N}^{\rm log}(m)\propto m^{N-1}\prod_i p_i(m).
\tag{20}
$$

This is a different model: equal-width matching in log mass replaces
equal-width matching in mass. It corresponds to multiplying the reference
(17) by \(m^{N-1}\).

Correctly expressing the original rule (18) in log coordinates instead gives
\(g_{1:N}(\ell)\propto\prod_i g_i(\ell)/m^{N-1}\). An ordinary change
of variables does not change the inference. Selecting a new flat matching
measure does. Reciprocal symmetry can motivate comparing rules, but by
itself does not prove a rule uniquely warranted or calibrated.

### A transparent zero-observation example

For exactly zero observed vectors and common \(s\), a single standardized
mass \(r=m/s\) has density \(q_1(r)=2/[\pi(1+r^2)]\). Two copies give

$$
q_2(r)=\frac4{\pi(1+r^2)^2}\quad\text{under (18)},\qquad
q_2^{\rm log}(r)=\frac{2r}{(1+r^2)^2}\quad\text{under (20)}.
\tag{21}
$$

Both are proper, exact, and different. For \(N\) copies the first rule gives

$$
q_N(r)=\frac{2\Gamma(N)}{\sqrt\pi\,\Gamma(N-1/2)}
                         (1+r^2)^{-N}.
\tag{22}
$$

The sum-magnitude readout defined in §7 equals \(s/2\) for the first density
in (21) and \(s\) for the second. The ordinary mass mean of the first is
\(2s/\pi\), illustrating yet another distinct readout.

These are statements about exactly zero *observations*. They are not a
sampling experiment with true zero excitation and random noisy readings.
Neither exact normalization nor a narrowing curve proves that informative
data were obtained. The separate calibration and null-excitation problems
are examined in [contribution 4](joint_inference_and_uncertainty.md).

## 5. A proof of nonclosure, with its scope made explicit

There is a simple tail obstruction to compressing (18) into another original
single-reading law. It applies beyond isotropic errors.

For finite observations and a positive-definite Gaussian covariance, let
\(H(f,\alpha)\) be the likelihood integrated over common direction under
the original flat-magnitude reference. It is positive, and a Gaussian
envelope controls its tails. The induced density is

$$
p(m)=\frac1Z\int_0^\infty\alpha H(m\alpha,\alpha)\,d\alpha
=\frac1{Zm^2}\int_0^\infty fH(f,f/m)\,df.
\tag{23}
$$

Dominated convergence therefore gives

$$
p(m)\sim\frac{C_\infty}{m^2}\quad(m\to\infty),
\qquad p(m)\to C_0\quad(m\downarrow0),
\qquad C_\infty,C_0>0.
\tag{24}
$$

The product of \(N\) such curves has tail \(m^{-2N}\). For \(N\geq2\)
it cannot equal any finite-data member of the original family, whose tail is
always \(m^{-2}\). The log-matching rule instead has tail
\(m^{-(N+1)}\) and behavior \(m^{N-1}\) at zero, so it also leaves the
family. A finite change of the original parameters cannot remove this
obstruction.

This proves nonclosure **for these combination rules and this family**.
It does not prove that no finite-dimensional joint inference exists.
For example, the different joint Jeffreys model in
[contribution 4](joint_inference_and_uncertainty.md) has a compact exact
reduction under a common instrument scale ratio: the dimension \(K=Nd\),
the ratio \(s\), and the three total standardized inner products
\(A_{\rm total},B_{\rm total},C_{\rm total}\). That is a different joint
law, with its own documented calibration limitations.

Nor can every such joint model be reconstructed from baseline mass curves
alone: the equal-norm perpendicular examples in §3 already give identical
curves for different original inner products.

## 6. A mass point and log spread are insufficient: an exact counterexample

Counting parameters suggests information loss, but does not prove it. The
following counterexample lies within the original isotropic Gaussian family.
Use reference units with \(s=1\), and set \(A=B=t,C=c\). Define the entire
function

$$
\mathcal H(z)=\int_0^1e^{z(1-u^2)/2}\,du;
\qquad \mathcal H(h^2)=J_1(h).
$$

The log-mass density is

$$
g_{t,c}(\ell)\propto\operatorname{sech}\ell\,
                \mathcal H(t+c\operatorname{sech}\ell).
\tag{25}
$$

It is symmetric, so its median mass is 1, its mean log mass is zero, and
the original ratio-of-means point is exactly 1. Compare

$$
(t,c)=(4,2),\qquad (t,c)=(8,c_*),
\tag{26}
$$

where \(c_*\) is the unique value in \((0,2)\) that makes the two log
variances equal. This definition, rather than a rounded decimal, specifies
the exact counterexample.

**Existence and uniqueness.** Write \(b(z)=(\log\mathcal H)'(z)\).
The integral representation makes \(b\) the mean of
\((1-u^2)/2\) under exponential tilting, so \(b>0\) and \(b'>0\).
For \(y=\operatorname{sech}\ell\), increasing \(c\geq0\) makes the
density ratio strictly increase with \(y\). Increasing \(t\) at fixed
\(c>0\) does too. Thus either change concentrates \(|\ell|\) and
strictly reduces its second moment. At \(c=0\) the density is proportional
to \(\operatorname{sech}\ell\), with variance \(\pi^2/4\), and

$$
\operatorname{Var}_{8,0}(\ell)
>\operatorname{Var}_{4,2}(\ell)
>\operatorname{Var}_{8,2}(\ell).
$$

Continuity and strict monotonicity give the claimed unique \(c_*\).

**The curves differ.** If the normalized densities were equal, then
\(\mathcal H(8+c_*y)=K\mathcal H(4+2y)\) on \(0<y\leq1\), for a
constant \(K\). Analytic continuation extends this identity to positive
\(y\). But
\(\mathcal H(z)\sim\sqrt{\pi/(2z)}e^{z/2}\) as \(z\to\infty\),
so equality would force \(c_*=2\), a contradiction.

These parameters are physically realizable observations: take
\(U=(\sqrt t,0,0)\) and
\(V=(c/\sqrt t,\sqrt{t-c^2/t},0)\), with equal component uncertainties.
Numerically,

$$
c_*\approx1.8625450752016093,\qquad
\operatorname{SD}(\log M)\approx1.308944317672899.
$$

The difference survives combination with the same third curve, whose
parameters are \((A,B,C,s)=(9,1,3,1)\), using (18):

| Quantity | First curve | Second curve |
|---|---:|---:|
| Original mass point | 1 exactly | 1 exactly |
| Original log standard deviation | 1.308944317673 | 1.308944317673 |
| Original probability \(M<1/2\) | 0.245522184873 | 0.245727551358 |
| Median after the same third reading | 1.143173028999 | 1.143364150981 |

This disproves sufficiency of the original mass point plus log standard
deviation for exact curve combination. It does not prove that every
conceivable uncertainty statistic, or an artificial real-number encoding
of a whole curve, is impossible. The useful conclusion concerns ordinary
measurement summaries and the specified uncertainty parameter.

## 7. An exact two-function output for a specified readout extension

The mass curve is enough to compute its own median, quantiles, and any
finite ordinary mass moments. To preserve a ratio-of-means construction,
one must also track how acceleration magnitude depends on candidate mass.
For each trial retain

$$
\bigl(p_i(m),a_i(m)\bigr),\qquad
a_i(m)=E[\alpha_i\mid M=m].
\tag{27}
$$

For the isotropic baseline, both functions are explicitly given by
(13) and (15). Under (17), nuisance laws are conditionally independent
given mass, and adding another trial does not change a trial's conditional
law at a fixed mass. Define

$$
A_{1:N}(m)=\sum_i a_i(m).
$$

If the declared extension is the ratio of summed expected magnitudes, then

$$
\boxed{\widehat m_{1:N}
=\frac{E[\sum_i f_i]}{E[\sum_i\alpha_i]}
=\frac{\int_0^\infty mA_{1:N}(m)p_{1:N}(m)\,dm}
       {\int_0^\infty A_{1:N}(m)p_{1:N}(m)\,dm}}.
\tag{28}
$$

This supplies an associative exact merger:

$$
\boxed{(p_A,A_A)\star(p_B,A_B)
=\bigl(\operatorname{normalize}(p_Ap_B),\ A_A+A_B\bigr)}.
\tag{29}
$$

Here \(A\) and \(B\) label disjoint batches. For the log-matching reference,
replace the first component by \(\operatorname{normalize}(m p_Ap_B)\).
The conditional acceleration functions still add because the reference
change depends only on mass.

Equations (28)–(29) are exact for this specified law and readout. Their
sum-magnitude choice is an additional across-trial convention, not physical
vector addition of experiments conducted at different times, and is not
forced by the original composition theorem. The stacked-norm construction
in contribution 4 is another possible extension. Choosing an ordinary
posterior mean instead is possible where it exists, but changes the
original readout.

### Why the mass curve alone cannot preserve this extension

Multiply both observed vectors and both channel uncertainties of one trial
by the same positive factor \(\lambda\), holding the units fixed.
Its \(A,B,C,s\), full mass curve, and original mass point are unchanged.
Its conditional acceleration function is multiplied by \(\lambda\).
The combined mass curve under (18) is therefore unchanged, while (28)
generally changes.

For a concrete original-family example, use these two trials, in fixed
consistent units:

| Trial | Observed force vector | Observed acceleration vector | \(\sigma_F\) | \(\sigma_a\) |
|---|---|---|---:|---:|
| 1 | \((2,0.2,0)\) | \((1,0,0)\) | 0.5 | 0.4 |
| 2 | \((5,0,0.1)\) | \((2.1,0.2,0)\) | 1 | 0.6 |

Their individual points are approximately 1.993026469631 and
2.373589604907. Equation (28) gives 2.008902161134. Scaling every input
of trial 1 by 10 leaves both individual mass curves and points unchanged,
but the combined readout becomes 1.997267503687.

This proves that mass curves and mass points do not suffice for (28).
The scaling example does not prove insufficiency if the original global
force and acceleration means are also supplied, since those means change
under scaling. The conditional function is a direct sufficient export;
recoverability from some alternative richer parameterization is a separate
question.

Outside the restricted Gaussian family, even global means can be
insufficient. Let mass be uniform on \(\{1,2,3\}\), and let conditional
acceleration be either \((1,1,1)\) or \((5/4,1/2,5/4)\), deterministically
at those masses. Both laws have \(E[\alpha]=1\), \(E[f]=2\), and the
same mass curve and point. Combine either with a third law uniform on
\(\{1,2\}\), with acceleration 1, by common-mass matching. The resulting
mass laws agree, but the summed-magnitude readouts are respectively
\(3/2\) and \(7/5\). This illustrates the information carried by the
conditional function without claiming an additional Gaussian-family proof.

## 8. Bessel functions and basis representations

The notation \(J_1(h)\) here denotes the integral (9). It is not the
standard ordinary Bessel function \(J_\nu\). The inner function
\(\sinh z/z\) *is* the modified spherical Bessel function \(i_0(z)\);
see [NIST DLMF, §10.49](https://dlmf.nist.gov/10.49).

Ordinary Bessel functions can form Fourier–Bessel orthogonal bases on a
specified interval with the appropriate weight and boundary conditions.
That concerns expansion of functions as sums. It does not imply that the
product of two basis functions is one basis function, or that arbitrary
products of our mass curves have only two parameters. See
[DLMF orthogonality](https://dlmf.nist.gov/10.22.E37) and
[Fourier–Bessel expansions](https://dlmf.nist.gov/10.23#iii).

There is even an exact product identity for the inner function:

$$
i_0(a)i_0(b)=\frac12\int_{-1}^1
i_0\!\left(\sqrt{a^2+b^2+2abt}\right)\,dt,
\qquad a,b\geq0.
\tag{30}
$$

For \(a,b>0\), substitute \(r=\sqrt{a^2+b^2+2abt}\); the right side is
\([\cosh(a+b)-\cosh(|a-b|)]/(2ab)=\sinh a\sinh b/(ab)\).
Zero arguments follow by continuity. The identity is a continuous mixture,
not collapse to one \(i_0\). The complete mass curves additionally contain
radial integration and the mass Jacobian.

A basis can nevertheless encode exact combination. For consistently
defined likelihood factors, map \(m>0\) to \(0<x<1\), for example by
\(x=(2/\pi)\arctan(m/m_0)\), and expand

$$
\log\mathcal L_i(m(x))=\sum_k c_{ik}\phi_k(x).
$$

Independence gives coefficient addition in the log likelihood:

$$
\log\mathcal L_{1:N}(m(x))
=\sum_k\left(\sum_i c_{ik}\right)\phi_k(x).
\tag{31}
$$

A full convergent expansion can be exact in its stated function-space
sense. Truncation is an approximation, and convergence of the reconstructed
normalized density and its tails needs checking. No special role for a
Bessel basis has been established here. Direct function evaluation,
one-dimensional quadrature, or a retained factor list is often simpler.

## 9. What this changes about the estimator's role

The investigation separates three operations:

1. Reconcile force and acceleration uncertainty with the assumed law and
   a declared reference measure.
2. Express the resulting information about mass as a full curve, with its
   dependence assumptions and provenance.
3. Select a scalar readout for the intended use.

The ratio-of-means theorem characterizes the third operation on a given
law. The \(J(h)\) result makes the first two operations substantially more
explicit in the isotropic case. Exact combination concerns a joint model
connecting several such outputs, not merely repeated application of the
scalar equation.

Thus the estimator can serve as a producer of reusable measurement
information. A point with an error bar is one presentation of that
information, and can be a useful approximation in a resolved regime; it is
not a generally sufficient accumulator. The two-function output in §7 is
a concrete exact accumulator for one chosen extension.

This does not establish that the current single-trial curve is a universal
likelihood, uniquely sufficient for every future model, or already
calibrated across repeated experiments. Raw observations and covariance
information remain useful when shared nuisance variables, calibration, or
the reference assumptions change. The choice of a justified joint model
and the assessment of its sampling behavior remain substantive work.

## 10. Implementation boundary and verification

No production estimator or fusion policy was changed for this account.
The current implementation in [estimator.py](../estimator.py) computes a
mass distribution on an angular quadrature grid and global expected
magnitudes. Its internal radial calculations contain the information needed
for conditional means, but its returned estimate does not explicitly export
the function (15). The practical readout in
[practical_formulas.py](../practical_formulas.py) is more compressed still.

An eventual implementation of (29) would need a common physical mass
coordinate, reference and dependence metadata, a curve representation, and
the conditional acceleration information for (28). Raw angular grid
indices from trials with different \(s_i\) cannot be multiplied as if they
represented the same mass. Equation (13) includes the required Jacobian.
Analytic equations may be exact while their finite grids and quadratures
remain numerical approximations.

The following verification material is retained in the workspace scratch
area rather than duplicated in the document collection:

- [Combination reproduction script](../../../.tools/mass_output_combination_audit_20260917/check_combination.py)
  and [saved checks](../../../.tools/mass_output_combination_audit_20260917/checks.json).
  Orders 256 and 512 reproduce the original estimator's three tested mass
  points to absolute differences below \(9.60\times10^{-14}\); normalizations
  are within about \(10^{-15}\) of 1. Staged and direct combinations differ
  in log density by at most \(1.43\times10^{-14}\). The zero-observation
  formulas and the scaling example in §7 are also reproduced.
- [Same-point, same-spread script](../../../.tools/mass_combine_audit_20260917/closure_same_spread.py),
  [initial output](../../../.tools/mass_combine_audit_20260917/closure_same_spread.json),
  and [refined output](../../../.tools/mass_combine_audit_20260917/closure_same_spread_refined.json).
  Composite Gauss–Legendre order 32 with log cutoff 48 was checked against
  order 64 with cutoff 52. The root, spreads, CDFs, and combined medians
  agree to the precision displayed in §6. The script also reports a density
  L1 difference; that nonsmooth absolute-value integral is less converged
  and is not used for the proof or the displayed numerical claims.

These are scoped numerical checks, not universal error certificates.
The tail and equal-spread arguments supply the mathematical proofs.
To rerun the second script, use its defaults and then
`--order 64 --cutoff 52`; the first script writes its own `checks.json`.

## 11. Source snapshot and coverage

Managed context was consulted through the GitHub connector at the single
VD-docs `main` commit
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866),
reconfirmed on 17 September 2026. `catalog.yaml` was the managed inventory.
The relevant canonical sources, all under
`Newton-analysis/_dscn_mass_estimator/`, were:

- [Full working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md).
- [Readout composition and uniqueness](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_readout_composition_uniqueness_2026-09-06.md).
- [Additional derivations](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md).
- [What forces the ratio-of-means readout](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md).
- [Law-assumed-correct steering](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md).

The active local investigation supplied the more recent analysis:
[the initial audit](README.md),
[Gaussian conditioning](gaussian_conditioning.md),
[weighting and calibration](weighting_and_calibration.md),
[joint inference](joint_inference_and_uncertainty.md),
[conditional uncertainty](conditional_uncertainty_and_design.md), and
[the exact radial reduction](operational_warrant.md), together with the two
implementation files and scratch reproductions linked above.

Selection followed the mass-estimator catalog entries and source links,
then the local investigation index and searches for the radial integrals,
composition, pooling, and uncertainty results. This is a consolidation of
the present discussion against those managed and active local sources,
not an exhaustive history of every mass-estimator document or a new
literature-priority claim. The new combination account and scratch checks
are local, unsubmitted work; the managed sources establish the baseline
they examine. No managed document was edited.

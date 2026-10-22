# The ours27 width defect and its connection to candidate geometry

**Record — 28 September 2026.** Written by Codex at the user's request.

**Finding.** Under full Gaussian measurement covariance, the flat ours27 evidence has an exact score-variance contribution that is absent from its expected curvature. In a regular repeated-experiment limit this can make its mass law too narrow. Separate exact identities show how changing candidate size and direction weighting changes both candidate moments and evidence curvature.

These results identify a mechanism and a benchmark. They do not select the unique mass equation, and they do not turn every weak-signal coverage failure into an algebraic error.

## 1. The experiment and the meaning of its candidates

For reading $i$, write

$$
y_i=\begin{pmatrix}F_i\\a_i\end{pmatrix},
\qquad
B_m=\begin{pmatrix}mI\\I\end{pmatrix},
\qquad
y_i\sim N(B_m a_i^*,\Sigma_i).
$$

The readings in this calculation are independent. Each known $\Sigma_i$ is positive definite and may include component and force–acceleration correlations. A supplied joint covariance can instead be handled by stacking the readings.

The pair $(ma_i^*,a_i^*)$ is an aligned candidate. Specifying $(m_0,a_i^*)$ in a repeated-experiment check means generating observations about that candidate to test the procedure. It does not posit a separate ontic mass that replaces the project's operational definition.

For trial mass $m$, define

$$
r_i=F_i-ma_i,
$$

$$
V_i=\Sigma_{FF,i}
-m(\Sigma_{Fa,i}+\Sigma_{aF,i})
+m^2\Sigma_{aa,i},
$$

$$
K_i=(B_m^\top\Sigma_i^{-1}B_m)^{-1},
\qquad
b_i=K_iB_m^\top\Sigma_i^{-1}y_i.
$$

The flat ours27 log evidence and conditional candidate law are

$$
\ell(m)=-\frac12\sum_i r_i^\top V_i^{-1}r_i,
\qquad
a_i^*\mid m,D\sim N(b_i,K_i),
$$

with mass law $p(m\mid D)\propto\pi(m)e^{\ell(m)}$. These are the quantities being audited; no determinant factor is added to this existing evidence.

## 2. Exact score and curvature identities

Evaluate at a specified candidate mass $m_0$. Let

$$
C_i=\Sigma_{aF,i}-m_0\Sigma_{aa,i}.
$$

Gaussian regression gives

$$
b_i=a_i-C_iV_i^{-1}r_i,
\qquad
b_i\sim N(a_i^*,K_i),
\qquad
r_i\sim N(0,V_i),
$$

and $b_i$ and $r_i$ are independent. Differentiating the displayed log evidence, including the dependence of $V_i$ on $m$, gives

$$
U(m_0)=\ell'(m_0)=\sum_i b_i^\top V_i^{-1}r_i.
$$

Define

$$
J=\sum_i a_i^{*\top}V_i^{-1}a_i^*,
\qquad
G=\sum_i\operatorname{tr}(K_iV_i^{-1}).
$$

Then the following are exact finite-experiment identities:

$$
E[U]=0,\qquad
\operatorname{Var}(U)=J+G,\qquad
-E[\ell''(m_0)]=J.
$$

The variance identity follows immediately by conditioning on the $b_i$:

$$
\operatorname{Var}(U\mid b_1,\ldots,b_n)
=\sum_i b_i^\top V_i^{-1}b_i.
$$

Taking its expectation gives $J+G$.

For the curvature identity, whiten a reading as $z_i=\Sigma_i^{-1/2}y_i$ and let $P_m$ be the orthogonal projector onto the columns of $\Sigma_i^{-1/2}B_m$. Apart from a data-only constant, its log evidence is $z_i^\top P_mz_i/2$. The trace of $P_m$ is constant. Differentiating $P_m^2=P_m$ twice gives

$$
P_mP_m''P_m=-2P_m(P_m')^2P_m.
$$

The candidate mean $\mu_i=\Sigma_i^{-1/2}B_{m_0}a_i^*$ lies in that subspace. Consequently the noise contributes zero expected curvature and

$$
-E[\ell_i'']
=\|P_{m_0}'\mu_i\|^2
=a_i^{*\top}V_i^{-1}a_i^*.
$$

The final equality is the squared orthogonal displacement of the candidate subspace under a change of $m$, expressed in the residual covariance.

## 3. The regular width consequence

Consider a sequence of experiments with a consistent interior profile maximum, a concentrating Hessian, a score central limit theorem, and no dominating reading. Suppose also that the fixed starting density is smooth and positive near the limiting mass, so its derivatives are negligible at the leading information scale.

The local expansion of the profile score and the local expansion of the mass law then give, respectively,

$$
V_{\mathrm{sampling}}
\sim\frac{J+G}{J^2},
\qquad
V_{\mathrm{law}}
\sim\frac1J.
$$

Hence

$$
\frac{V_{\mathrm{sampling}}}{V_{\mathrm{law}}}
\sim1+\frac GJ.
$$

The standard-deviation factor is $\sqrt{1+G/J}$, not $1+G/J$. The score and curvature identities are exact; these variance statements require the stated regular limit.

Here the sampling variance refers to the local profile/root center. The project's ratio-of-mean-sizes point has the same leading behavior only when its additional smooth weighting shifts that center by a negligible amount at the sampling scale. This record does not assert a universal finite-sample variance formula for that point.

In the common-noise angular coordinate, write $R^2$ for the squared standardized candidate length and $M$ for the number of candidate components. Then $J=R^2$ and $G=M$. Keeping $R^2/M=\lambda>0$ while $M$ grows gives a regular informative sequence but leaves the variance ratio at $1+1/\lambda$. Thus the defect is not restricted to a nonidentified zero-signal limit.

## 4. An exact conditional test benchmark

The same calculation gives

$$
Z(m)=
\frac{\sum_i b_i^\top V_i^{-1}r_i}
{\sqrt{\sum_i b_i^\top V_i^{-1}b_i}}.
$$

At the tested candidate mass, $Z$ is exactly standard normal conditional on the $b_i$, whenever its denominator is positive. Positive-definite measurement covariance makes the exceptional zero denominator a probability-zero event in the sampling experiment.

This is an exact test pivot for the declared Gaussian experiment. Inverting $|Z(m)|\leq z_{1-\alpha/2}$ supplies a confidence-set benchmark. It does not supply a probability law over mass: $Z(m)$ need not be monotone, may have multiple stationary branches, and its p-values cannot be normalized into posterior odds. For the literal all-zero observed input, the formula is undefined; it therefore cannot serve as the all-input mass procedure.

A numerical implementation was checked on a full correlated covariance with 250,000 draws at each of three candidate strengths. The nominal 95% acceptance fractions were 0.950724, 0.949104 and 0.949928. Rotation, unit-change, channel-swap and score-differentiation checks agreed within $10^{-8}$. These checks support the implementation; conditional Gaussianity supplies the proof.

## 5. Exact identities connecting size, direction and curvature

In the common-noise coordinates of the [zero-reading record](ours27_zero_readings_and_shared_candidate_information_2026-09-28.md), let

$$
A=x\sin\theta+y\cos\theta,\qquad
B=x\cos\theta-y\sin\theta,
\qquad H=\|A\|^2,\quad Q=\|B\|^2.
$$

Take a fixed isotropic candidate measure $\mu_M(dv)$ with finite Gaussian-damped moments through the orders used below. Define

$$
F_M(H)=
\int \exp\left(-\frac{\|v\|^2}{2}+A^\top v\right)\mu_M(dv),
\qquad
\kappa(H)=2\frac{d}{dH}\log F_M(H).
$$

Isotropy makes this integral a function of $H$. The evidence is $e^{-(H+Q)/2}F_M(H)$. Differentiating its logarithmic partition function with respect to $A$ gives

$$
E[v\mid A]=\kappa A,
$$

$$
\operatorname{Cov}(v\mid A)
=\kappa I+2\kappa'(H)AA^\top,
$$

$$
E[\|v\|^2\mid A]
=H\kappa^2+M\kappa+2H\kappa'(H).
$$

For $C=A^\top B$, the angular evidence derivatives are

$$
\ell'(\theta)=\kappa C,\qquad
\ell''(\theta)=\kappa(Q-H)+2\kappa'(H)C^2.
$$

At an interior maximum of $H$, let $h\geq q$ be the eigenvalues of the two-by-two Gram matrix of $x,y$. Then $H=h$, $Q=q$, $C=0$, and the evidence curvature is

$$
-\ell''=\kappa(h)(h-q).
$$

This concerns evidence curvature; a nonconstant starting density contributes its own derivatives.

For flat Cartesian measure, $\kappa=1$:

$$
v\mid A\sim N(A,I),\qquad E[\|v\|^2\mid A]=H+M.
$$

In a specified fixed-candidate experiment, $E[H]=R^2+M$, so this conditional candidate-size mean averages to $R^2+2M$. That is a diagnostic of the measure's size weighting. It does not establish that every posterior-variance term is erroneous or that $R$ is a hidden physical quantity the project must recover.

The distinction between flat radial measure $dr$ and flat Cartesian measure $dv=r^{M-1}dr\,d\omega$ is substantive. Merely changing integration order under the same measure cannot change the answer. Any proposed correction must identify an actual change of measure and justify it.

## 6. Fieller limits and avoidable under-width are different claims

Fieller-type inference can retain coverage by allowing unbounded or disconnected confidence sets. The obstruction concerns procedures that insist on bounded intervals throughout a parameter space with arbitrarily weak ratio identification. See the primary treatment by [von Luxburg and Franz](https://arxiv.org/pdf/0711.0198).

For the present experiment, a short argument isolates the obstruction. As both aligned candidate vectors shrink to zero, the data distribution approaches the same noise-only distribution $P_0$, whatever the candidate mass. Suppose an interval procedure has a finite upper endpoint $U(D)$ almost surely under $P_0$. Then

$$
P_0\{U(D)\geq m\}\longrightarrow0
\quad\text{as }m\longrightarrow\infty.
$$

For each such $m$, choose a sufficiently small positive candidate amplitude so that its Gaussian data distribution is arbitrarily close to $P_0$. Coverage at that candidate is bounded above by a quantity arbitrarily close to $P_0\{U(D)\geq m\}$. Uniform positive coverage over the whole weak-signal parameter space is therefore impossible for this always-finite-upper-endpoint procedure.

A proper mass probability law may still return finite interior quantiles at every input. Those quantiles do not thereby acquire uniform repeated-experiment coverage. Conversely, the impossibility result does not require failure at every low-signal cell.

The regular $J+G$ versus $J$ calculation concerns a distinct, avoidable calibration mismatch in informative sequences. Neither claim alone selects a new definition of mass.

## 7. Provenance and status

Managed context was read at **canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8**, selected through its catalog: the [full-covariance ours27 formulation](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md), [September 27 lemmas](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/mass_estimator_lemmas_2026-09-27.md), and [learned-scale hint](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_omega_is_a_hint_2026-09-27.md).

Local verification includes .tools/ours27_next_2026_09_28/full_covariance_conditional_score.py and .tools/ours27_direction_2026_09_28/geometry/size_direction_identities.py, with their saved results. The identities and proofs needed for this record are stated above.

This local record has not been submitted to VD-docs. Candidate-measure uniqueness, cross-trial persistence and the final equation remain open.

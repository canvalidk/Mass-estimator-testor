# The equation: ours30

Stated once. ours30 (R-2026-10-29-1) is ours29 (`../docs/noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md`
§3; lemma list d4: T8, T10, T30, T31) with one change: each reading's misfit weight reads its misfit at
its best **positive** mass, with numerator $k_i-\tfrac12$. Edit in place when the equation changes; git
holds the history. The equation note in VD-docs is cut from this file when Can asks for one (draft d3,
for ours30: `../outbox/VD-docs/`).

## Inputs

A reading i declares its data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$,
an error covariance $\Sigma_i$ (positive definite, any correlation), a candidate subspace
$N_i$ ($d\times k_i$, orthonormal columns; $I_d$ when nothing is known about the push's
direction), and, since 10 Oct, a window over which both channels are averaged
(`type.md`). Declared once for the series: the starting law $\pi$, half-Cauchy with centre $c$,
$\pi(m)=c/(c^2+m^2)$; by default $c=s=\sigma_F/\sigma_a$ (L3′, not yet ruled).

## Per reading

For every candidate slope on the whole projective line, $m=c\tan\phi$,
$\phi\in(-\pi/2,\pi/2]$ (negative and infinite m included):

$$
H_i(\phi)=\mathbf b^\top G^{-1}\mathbf b,\qquad
C(\phi)=\begin{pmatrix}c\sin\phi\,N_i\\ \cos\phi\,N_i\end{pmatrix},\quad
G=C^\top\Sigma_i^{-1}C,\quad \mathbf b=C^\top\Sigma_i^{-1}\mathbf y_i,
$$

the explained energy. $E_i$ is the whitened energy of the data inside the span of all
candidates, and

$$
Q_{+,i}=E_i-\max_{0\le\phi\le\pi/2}H_i(\phi),\qquad \kappa_i=k_i-\tfrac12,\qquad
g_i=\min\!\Big(1,\frac{\kappa_i}{Q_{+,i}}\Big),\qquad
\ell_i(\phi)=\tfrac12\,g_i\,H_i(\phi).
$$

$Q_{+,i}$ is the reading's misfit at its best **positive** mass ($m=0$ and $m=\infty$ included as
limits; T32's $Q_+$): the smallest error that reconciles the reading with Newton II as assumed,
$m>0$ included (R-2026-09-07-1). It counts the $k_i-1$ codirectionality conditions the reading's
declaration imposes (L8) and L1's sign condition $\mathbf F^\ast\cdot\mathbf a^\ast>0$, which holds for
1D projections too. $\kappa_i$ is the most misfit the declared noise produces in expectation over
every allowed signal size (O12's meaning, read without a size by L10): $k_i-1$ for the equalities
and $\tfrac12$ for the sign, which noise decides half the time when one channel carries no signal
(O13, 29 Oct; numerical for isotropic channels, proof: brief B02). $g_i$ is the misfit weight
(L21, T30): a reading is evidence, just less; $g=1$ whenever $Q_{+,i}\le\kappa_i$, so on an aligned
1D reading; never 0. The whole-line $Q_i=E_i-\max_\phi H_i$ with $\nu_i=k_i-1$ was ours29's (it
cannot see the sign: an anti-parallel pair has $Q=0$).

## Series

At $\hat\phi=\arg\max_\phi\sum_i\ell_i$ on the whole line (unchanged from ours29; checked against
R-2026-09-07-1 on 29 Oct in T31: the whole line is forced here by Z and the stationary point),

$$
\mathcal H=-\sum_i\ell_i''(\hat\phi),\qquad J=\sum_i\ell_i'(\hat\phi)^2,\qquad
\eta=\min\!\Big(1,\frac{\mathcal H}{J}\Big)\quad(\eta=1\ \text{if}\ J=0).
$$

$\eta$ is the Bartlett power (T31): one per series, from the readings' disagreement about m
alone, never from their sizes (L19). An exact zero reading has $\ell_i'=\ell_i''=0$ and
leaves $\eta$ alone (Z).

## Law

$$
p(m)\;\propto\;\pi(m)\,\exp\!\Big(\eta\sum_i\ell_i(m)\Big)\qquad\text{on } m>0 .
$$

With one noise ratio per series (isotropic noise) everything is closed form in
$C_i=(Q_i-P_i)+2\mathrm iD_i$, where $P_i=\lvert\mathbf x_i\rvert^2$, $Q_i=\lvert\mathbf y_i\rvert^2$,
$D_i=\mathbf x_i\cdot\mathbf y_i$ in noise units: the best line is $\arg(\sum_i g_iC_i)/2$,
$\mathcal H=\lvert\sum_i g_iC_i\rvert$, $\ell_i'=-\tfrac12 g_i\,\mathrm{Im}(\bar C_ie^{2\mathrm i\hat\theta})$
(tester README, "ours29" and "ours30"), and the positive-mass misfit is
$Q_{+,i}=\tfrac12\big(P_i+Q_i-\lvert C_i\rvert\big)$ when $D_i>0$ and $\min(P_i,Q_i)$ otherwise. Note the
two uses of the letter Q: the reading's misfit above, and the acceleration energy in this closed
form; the records use both.

## Point (28 Oct, provisional)

The definition's point is the law's **median**: $\widehat m=F_p^{-1}(\tfrac12)$. It is a functional
of the law alone (L19: readings are joined only through the mass), it is the only point that
does not depend on the coordinate the law is read in (L23, proposed; T35, with the swap), and
composition holds at the level of the law (L5, L14 restated, R-2026-10-28-1), not of the point.
It meets Z exactly (T24), returns $c$ at the all-zero experiment for every $c$, and respects the
swap and units. On a CDF flat to rounding it is the plateau's midpoint in $\ln m$.
In `met/full.py`: `ours30(...)["point"]`, the same number as `["median"]`.

## Candidates and the cloud readout (kept for comparison)


Each reading keeps its own flat candidate cloud, $\alpha_i\mid m\sim N(G^{-1}\mathbf b,\,G^{-1})$,
$\mathbf a^\ast_i=N_i\alpha_i$, $\mathbf F^\ast_i=m\,\mathbf a^\ast_i$. The readout is the
ratio of means on the stacked length (T4, T9, T10):

$$
\widehat m=\frac{\int m\,A(m)\,p(m)\,dm}{\int A(m)\,p(m)\,dm},\qquad
A(m)=E\big[\lVert\mathbf a^\ast_{\rm stack}\rVert\;\big|\;m\big],
$$

equivalently $E[\lVert\mathbf F^\ast\rVert]/E[\lVert\mathbf a^\ast\rVert]$. Neither $g_i$ nor
$\eta$ touches $A$ (shrinking the candidates by $\eta$ would infer signal size from the
series, against L19). **This readout fails Z (T29) and gives c at the all-zero experiment
only when c = s.** It joins readings through their sizes (against L19) and on null readings its
weight is set by the instrument alone (O2, 27–28 Oct). Since 28 Oct it is not the point.

## Output

The point $\widehat m$ (the law's median) and the 95% interval (the quoted width), the raw
law's interval ($\eta=1$) beside it, and $\eta$ itself (ruling 3 of 30 Sept). The misfit is kept
with the result. The ordinary mean of m under $p$ is infinite (T19), so "mass ± SD" is not
available.

## Exact inputs

$\sigma=0$ does not enter this function. It goes to the exact branch (T1): compute
$m_1=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and $m_2=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$;
one mass iff $m_1=m_2\in(0,\infty)$; no mass if they disagree or agree outside $(0,\infty)$;
every mass at exact 0/0, the only place "undefined" survives.

## Declared, still free

The shape of $\pi$ (O6) and, with several noise pairs, its centre (O11). The middle of $g$
was declared until 30 Sept and is now forced (T30, O12 closed). $\kappa$'s value under a general
covariance (cross-channel correlation) is unknown yet; $k-\tfrac12$ is shown for isotropic channels.

## In the tester

`law: {reference: cartesian2, nuisance: radius}` with
`combine: {rule: likelihood_product, prior: [{uniform_angle: {center: first_reading}}], misfit: kink, positivity: true, calibrate: bartlett}`
and readout `median` (the point): the batch pipeline, isotropic and per-axis noise. Without
`positivity` it is ours29.

The full input type above (any Σᵢ, a declared Nᵢ, the window carried) is `met/full.py`:
`ours30(readings, c)` on a list of `Reading(force, acceleration, cov, subspace=None, window=None)`,
with `c` required (π is an input). It agrees with the batch pipeline on isotropic noise, and on the
weights with per-axis noise (`tests/test_ours30.py`); the weight's properties are in
`tests/test_sign_numerator.py`. `ours29(readings, c)` is kept, and agrees with
`../reference/noise_vs_signal/ours29_reference.py` at its grid (`tests/test_full.py`). Presets and
the world do not generate correlated noise yet.

## History

- ours21 (21 Sept): flat in the two magnitudes; weighting $df\,d\alpha\,d\Omega$.
- ours24 (24 Sept): flat in the push vector (the $r^{d-2}$ change).
- ours25b (25 Sept): stacked-length readout for many readings.
- ours26 (26 Sept): per-axis noise. ours27 (27 Sept): full covariance per reading.
- ours28 (28 Sept): a per-reading size law; withdrawn 30 Sept.
- ours29 (30 Sept): misfit weight and Bartlett power.
- 10 Oct: the window joins the type; the exact branch restated in pairing form.
- 2026-10-24: this file created.
- 2026-10-28: the point is the law's median (provisional: O2, T35, L23 proposed, R-2026-10-28-1); the
  cloud readout kept for comparison. No recorded number moved; `met/full.py` gains `point`.
- **ours30 (29 Oct, R-2026-10-29-1): the misfit weight reads the misfit at the best positive mass,
  numerator k − ½. Current.** Medians on the robot data move by under 1e-4; acute readings are
  discounted less, by (k − ½)/(k − 1); anti-parallel readings get bounded weight in 1D–3D.

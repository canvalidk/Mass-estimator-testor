# The equation: ours29

Stated once, from `../docs/noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md` §3
and lemma list d4 (T8, T10, T30, T31). Edit in place when the equation changes; git holds the
history. The equation note in VD-docs is cut from this file when Can asks for one.

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
Q_i=E_i-\max_\phi H_i,\qquad \nu_i=k_i-1,\qquad
g_i=\min\!\Big(1,\frac{\nu_i}{Q_i}\Big),\qquad
\ell_i(\phi)=\tfrac12\,g_i\,H_i(\phi).
$$

$Q_i$ is the reading's off-line misfit at its own best line, the violation of the $k_i-1$
codirectionality conditions its declaration imposes (L8). $g_i$ is the misfit weight (L21,
T30): a reading is evidence, just less; $g=1$ in 1D and whenever $Q_i\le\nu_i$.

## Series

At $\hat\phi=\arg\max_\phi\sum_i\ell_i$ on the whole line,

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
(tester README, "ours29"). Note the two uses of the letter Q: the reading's misfit above,
and the acceleration energy in this closed form; the records use both.

## Candidates and readout

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
only when c = s; replacing it is O2** (`../problems/O2.md`).

## Output

The law's median and 95% interval (the quoted width), the raw law's interval ($\eta=1$)
beside it, the point $\widehat m$, and $\eta$ itself (ruling 3 of 30 Sept). The misfit is kept
with the result. The ordinary mean of m under $p$ is infinite (T19), so "mass ± SD" is not
available.

## Exact inputs

$\sigma=0$ does not enter this function. It goes to the exact branch (T1): compute
$m_1=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and $m_2=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$;
one mass iff $m_1=m_2\in(0,\infty)$; no mass if they disagree or agree outside $(0,\infty)$;
every mass at exact 0/0, the only place "undefined" survives.

## Declared, still free

The shape of $\pi$ (O6) and, with several noise pairs, its centre (O11). The middle of $g$
was declared until 30 Sept and is now forced (T30, O12 closed).

## In the tester

`law: {reference: cartesian2, nuisance: radius}` with
`combine: {rule: likelihood_product, prior: [{uniform_angle: {center: first_reading}}], misfit: kink, calibrate: bartlett}`
and readout `ratio_of_means_rss`. `met/` takes isotropic and per-axis noise only; the full
input type above is implemented in `../reference/noise_vs_signal/ours29_reference.py`.

## History

- ours21 (21 Sept): flat in the two magnitudes; weighting $df\,d\alpha\,d\Omega$.
- ours24 (24 Sept): flat in the push vector (the $r^{d-2}$ change).
- ours25b (25 Sept): stacked-length readout for many readings.
- ours26 (26 Sept): per-axis noise. ours27 (27 Sept): full covariance per reading.
- ours28 (28 Sept): a per-reading size law; withdrawn 30 Sept.
- **ours29 (30 Sept): misfit weight and Bartlett power. Current.**
- 10 Oct: the window joins the type; the exact branch restated in pairing form.
- 2026-10-24: this file created.

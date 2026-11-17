# The equation: ours33, the frontier probe

**Status (17 Nov 2026).** Not the definition. This page states the probe the work is at (R-2026-11-17-2): ours33,
assembled on 4 Nov (`../problems/O18.md`). Its reading level, its raw law and its point stand, and they are what
R-2026-11-17-1 asks for: computed from what is declared. It breaks at the report, the quoted width ("Series and
report", below). ours30, which this page stated until 17 Nov, is history: its statement is this file at commit
`688cd8f`. Edit in place when the equation changes; git holds the history.

## Inputs

A reading i declares its data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$,
an error covariance $\Sigma_i$ (positive definite, any correlation), how well that covariance is known,
$\nu_i=\nu_{{\rm eff},i}\in(0,\infty]$ (the effective degrees of freedom of the stated uncertainty; $\infty$ declares it
exact), a candidate subspace $N_i$ ($d\times k_i$, orthonormal columns; $I_d$ when nothing is known about the
push's direction), and a window over which both channels are averaged (`type.md`). Declared once for the series:
the starting law $\pi$, half-Cauchy with centre $c$, $\pi(m)=c/(c^2+m^2)$. The centre is an input; $c=s=\sigma_F/\sigma_a$
is the caller's convention for one noise ratio (O11).

## Per reading

**The frame (T36).** A reading lives in the stacked space of its force and acceleration
components; the only declared metric on it is $\Sigma_i^{-1}$ (the axes carry N and m s⁻², so space
supplies none). A mass $m$ is the $k_i$-dimensional subspace $\{(mN_i\alpha,N_i\alpha)\}$, the column
space of $C(\phi)$ below; $H_i$ is the squared length of the reading's projection onto it, and
$K_i=E_i-H_i$ the residual: the squared distance, in the noise metric, from the reading to the candidate.

For every candidate slope on the whole projective line, $m=c\tan\phi$,
$\phi\in(-\pi/2,\pi/2]$ (negative and infinite m included):

$$
H_i(\phi)=\mathbf b^\top G^{-1}\mathbf b,\qquad
C(\phi)=\begin{pmatrix}c\sin\phi\,N_i\\ \cos\phi\,N_i\end{pmatrix},\quad
G=C^\top\Sigma_i^{-1}C,\quad \mathbf b=C^\top\Sigma_i^{-1}\mathbf y_i,\qquad
K_i(\phi)=E_i-H_i(\phi),
$$

with $E_i$ the whitened energy of the data inside the span of all candidates. The reading's factor is

$$
\ell_i(\phi)=-\frac{\nu_i+k_i}{2}\,\log\!\Big(1+\frac{K_i(\phi)}{\nu_i}\Big)\quad(\nu_i<\infty),\qquad
\ell_i(\phi)=-\tfrac12K_i(\phi)\quad(\nu_i=\infty).
$$

This is the declared noise law's probability density of the reading's residual at $m$, up to a constant free of
$m$ (T48). Given its noise size $v$ the reading's error is $N(0,v\Sigma_i)$, with $v\sim\mathrm{IG}(\nu_i/2,\nu_i/2)$:
how well $\Sigma_i$ is known, as declared (T37, G20). The push is integrated flat in noise units, and
independently of $v$, since the declared noise-size law holds whatever the push (`../problems/O18.md`, 17 Nov).

- The law applies no weight for a reading's own misfit. There is no $g$, no $\kappa$, no cap.
- The sign of $\mathbf F\cdot\mathbf a$ enters only through the domain $m>0$.
- An exact-zero reading has $K_i\equiv0$ and contributes nothing (Z).
- With $\nu_i<\infty$ the evidence from a reading no positive mass fits stays bounded as its noise shrinks, and its
  direction is kept. With $\nu_i=\infty$ it grows without bound. Each is the consequence of its declaration
  (R-2026-11-17-1).

## Raw law

$$
p_0(m)\;\propto\;\pi(m)\,\exp\!\Big(\sum_i\ell_i(m)\Big)\qquad\text{on } m>0 .
$$

Proper: every factor is in $(0,1]$. On the all-zero series it is $\pi$.

## Point

The raw law's **median**, $\widehat m=F_0^{-1}(\tfrac12)$: a functional of the law alone (L19), respecting units and
the swap, $c$ on the all-zero series. Provisional: units and the swap leave a family of location functionals of
$\ln m$, and the median is the member in the code (G6). Nothing in the report moves it.

## Series and report (where it breaks)

At $\hat\phi=\arg\max_\phi\sum_i\ell_i$ on the whole line,

$$
\mathcal H=-\sum_i\ell_i''(\hat\phi),\qquad J=\sum_i\ell_i'(\hat\phi)^2,\qquad
\eta=\min\!\Big(1,\frac{\mathcal H}{J}\Big)\quad(\eta=1\ \text{if}\ J=0\ \text{or}\ \mathcal H\le0),
$$

one number per series, from the readings' disagreement about $m$ alone (T31). The quoted law is the normal-score
warp of the raw law's distribution function (T47),

$$
F_\eta=\Phi\big(\sqrt\eta\,\Phi^{-1}(F_0)\big),
$$

so the quoted 95% interval is the raw law's quantiles at ranks $\Phi(\mp z_{0.975}/\sqrt\eta)$. The raw law is never
raised to a power: $\eta$ widens the report and touches neither the law nor the point (R-2026-09-30-3: the law does
not carry its own uncertainty).

**Broken here (O18, 4 Nov).**

- **The width.** At small $\eta$ the warp asks for ranks that lie in $\pi$'s tails, which the data never touch. On two
  readings that both sit at 2 (T39's pair) the quoted interval is $[2\times10^{-36},\,2\times10^{42}]$ while the raw one
  is $[1.985,\,2.015]$. The quoted width has to be a statement about the evidence, not a quantile of a law whose
  tails are the declared $\pi$'s; and the evidence alone has no finite quantile interval (B16). What the width is, is
  not written (`proof.md` G28; B24).
- **How $\eta$ is computed.** The fitted $J$ gives $\eta\to1/99$ on that same pair, fails under one strong reading
  among nulls (T43), and no estimate reaches the target $\Delta/(\Delta+M)$ at small $N$ (O14; B15).
- **Composition.** A warp that commutes with composition is the identity (T51); under ours33 the order matters
  (B18: 87.75% against 91.75%).
- **The weak-signal centre.** The warp widens around the raw median, and the raw median carries the flat push
  measure's tilt whatever $\eta$ is (O15).

## Output

The point, the quoted 95% interval, the raw law's interval beside it, and $\eta$. The ordinary mean of m under $p_0$
is infinite (T19), so "mass ± SD" is not available.

## Exact inputs

$\sigma=0$ does not enter this function. It goes to the exact branch (T1): compute
$m_1=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and $m_2=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$;
one mass iff $m_1=m_2\in(0,\infty)$; no mass if they disagree or agree outside $(0,\infty)$;
every mass at exact 0/0, the only place "undefined" survives. A $\Sigma_i$ declared exactly known
($\nu_i=\infty$) is not this: it is a measured input with $\sigma>0$.

## What stands, and what is open

- **Stands:** the type with $\nu_{\rm eff}$; the per-reading factor, forced for readings with the full subspace once
  $\nu_{\rm eff}$ is declared (`proof.md` G19, 17 Nov; provisional); the product; the raw law.
- **Open at the reading level:** a declared subspace with $k_i<d$ (ours33 keeps the residual inside the candidates'
  span; keeping the whole reading gives exponent $(\nu+2d-k)/2$, T48); the normalisation of the noise-size law
  (G20); the old gaps in the flat push (G2, G17).
- **Open beside it:** the shape of $\pi$, an input (O6); which location functional the point is (G6); consistency
  at the rate as $N$ grows, unknown for this raw law (O19); everything under "Series and report".

## In the tester

`met/ours33.py`: `ours33(readings, c, nu_eff)` on a list of
`Reading(force, acceleration, cov, subspace=None, window=None)`; `nu_eff` one value or one per reading,
`math.inf` for exact; `calibrate=False` returns the raw law's interval as the quoted one. Returns `point`, `median`,
`interval`, `raw_interval`, `eta`, `law` on `lm`. Its premise report is `tests/test_ours33.py`; coverage,
`tools/ours33_coverage.py`. ours30 and ours29 remain callable in `met/full.py`, and ours31, ours32 beside it, as
history.

## History

- ours21 (21 Sept): flat in the two magnitudes; weighting $df\,d\alpha\,d\Omega$.
- ours24 (24 Sept): flat in the push vector (the $r^{d-2}$ change).
- ours25b (25 Sept): stacked-length readout for many readings.
- ours26 (26 Sept): per-axis noise. ours27 (27 Sept): full covariance per reading.
- ours28 (28 Sept): a per-reading size law; withdrawn 30 Sept.
- ours29 (30 Sept): misfit weight and Bartlett power.
- 10 Oct: the window joins the type; the exact branch restated in pairing form.
- 2026-10-24: this file created.
- 2026-10-28: the point is the law's median (provisional); the cloud readout kept for comparison.
- ours30 (29 Oct, R-2026-10-29-1): the misfit weight reads the misfit at the best positive mass, numerator k − ½.
  Stated here until 17 Nov. Its weight was retired on principle by R-2026-11-17-1.
- ours31 (30 Oct): the exact cut factor per reading, $\nu_{\rm eff}$ declared, no η. ours32 (1 Nov): η put back, as a
  power. Both in `../problems/O18.md`.
- ours33 (4 Nov): the declared noise's density of the residual (T48) and η as a warp of the report (T47).
- **2026-11-17: this page states ours33** (R-2026-11-17-2); $\nu_{\rm eff}=\infty$ accepted (R-2026-11-17-1).

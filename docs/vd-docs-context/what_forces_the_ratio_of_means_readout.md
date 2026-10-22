# What Forces the Ratio-of-Means Readout

## Reciprocal invariance does not select $\widehat m=\overline f/\overline\alpha$; the composition laws of Newtonian mass do

**Status:** Claude-authored estimator review, 2026-09-02. Kind unset; filing is the author's call.
**Project:** VD Newton / inertial-mass estimation
**Purpose:** answer one question about `inertial_mass_estimator_full_working_spec_d1` — is the central readout $\widehat m=\overline f/\overline\alpha$ the right estimator, and can its form be derived rather than chosen? — and record what the check turned up on the way: which arguments in the spec do the selecting, which only appear to, one place where the single-trial construction does not extend, and the calibration of the gate.
**Source basis:** `inertial_mass_estimator_full_working_spec_d1.md` (Parts III, IV, VI, IX), `median_based_estimator_inertial_mass.md`, `disagreement_rule_mass_estimator_trace.md` §7, `why_the_estimator_is_the_definition.md` §§1–3. All numbers below were computed with a grid/Monte-Carlo harness against the spec's own test cases (Appendix).

---

## 1. Thesis

The form is right. The reason given for it is not the reason it is right.

The spec presents reciprocal invariance (§11) as the hard requirement that the compatible-pair mean readout satisfies and that competitors fail. Reciprocal invariance is satisfied by every reasonable candidate — the induced-mass median, every quantile, the geometric mean, the RMS ratio, the profile minimiser, and any readout of the form $T(f)/T(\alpha)$ whatsoever. It excludes nothing that anyone would have proposed. What actually singles out $\overline f/\overline\alpha$ is a stronger property the spec never states: **the readout is linear in the latent magnitudes**, and linearity is forced by requiring the estimator to compose the way Newtonian mass composes —

$$
\text{mass adds under lumping at common acceleration,}
\qquad
\text{inverse mass adds under common force.}
$$

Each of these two composition laws kills a different class of competitor. The ratio of means is the unique readout that survives both. Once linearity is in hand, dimensional covariance and exact consistency fix the form, and reciprocal invariance follows as a corollary rather than being imposed.

Three further results came out of the check:

1. $\widehat m$ has a clean identity: it is the **leverage-weighted posterior mean of the mass**, $\widehat m=E[\alpha\,m]/E[\alpha]$, and $\widehat\mu=E[f\,\mu]/E[f]$; the two are reciprocal *because* $f=m\alpha$ is the law. Among power weights $\alpha^k$, $k=1$ is the unique reciprocal choice; $k=2$ (Fisher-information weighting, equivalently posterior least squares) reproduces regression dilution at the posterior level.
2. The flat latent-pair measure is symmetric between $m$ and $\mu$ **only for a single trial**. For $N$ trials with a common mass the two flat measures differ by $m^{2N-2}$. Level 3 cannot inherit §11 from measure symmetry.
3. The gate statistic $Q_{\min}$ is $\chi^2_{d-1}$ under the true model at interior fits. $Q_{\mathrm{crit}}=9$ is a $0.3\%$ false-rejection threshold in 2D and $1.1\%$ in 3D.

---

## 2. Reciprocal invariance is generic

Let $P$ be the normalised support over compatible latent pairs $(f,\alpha,\widehat{\mathbf u})$ constructed in spec §9, and let $m=f/\alpha$ be the induced mass variable.

**Any separable readout is reciprocal-invariant.** If $T$ is any functional of a positive random variable and one defines $\widehat m_T=T(f)/T(\alpha)$, then $\widehat\mu_T=T(\alpha)/T(f)=1/\widehat m_T$ identically. The mean is one $T$; the median, the geometric mean, the RMS, every quantile, are others.

**Any equivariant summary of the induced ratio is reciprocal-invariant.** $\mu=1/m$ is a monotone map, so $\operatorname{median}(\mu)=1/\operatorname{median}(m)$, every quantile transforms, and $\exp E[\log\mu]=1/\exp E[\log m]$.

**The profile minimiser is reciprocal-invariant.** $\arg\min$ commutes with reparameterisation: $\arg\min_\mu Q_{\min}(1/\mu)=1/\arg\min_m Q_{\min}(m)$. The old best-point estimator's defect is not asymmetry; it is that the minimum sits at $m\to0^+$ or $m\to\infty$ for every zero reading (§17, §18) and at a *maximum* for anti-aligned data (§7). It is boundary-seeking, not asymmetric. The spec's Part III → Part VI ordering invites the wrong diagnosis.

So §11 is a necessary condition that every candidate passes. In particular:

> The **median of the induced mass distribution under the spec's own $\rho_\Pi$** is reciprocal-invariant, needs no excitation model, needs no prior on $m$, and needs no parameterisation choice for such a prior. §23's three reasons for demoting the median apply to the *earlier* median proposal (which carried an excitation model and a mass prior), not to this one. The spec already uses this induced distribution for the interval (§12, route 1). It has the median at hand and does not use it, and gives no reason that applies.

The test cases cannot decide between the candidates either. Numerically, on the spec's own cases (Appendix, Table A):

| case | $\overline f/\overline\alpha$ | median$(m)$ | $\exp E[\log m]$ | RMS ratio |
|---|---|---|---|---|
| strong aligned, $50/10$ | 5.000 | 5.000 | 5.000 | 5.000 |
| §17: $\widetilde F=0$, $\widetilde a=10$ | 0.160 | 0.135 | 0.107 | 0.200 |
| §18: $\widetilde F=50$, $\widetilde a=0$ | 125 | 148 | 187 | 100 |
| §21: zero–zero | 4.00 | 4.00 | 4.00 | 4.00 |
| §20: anti-aligned (conditional) | 1.00 | 1.00 | 1.00 | 1.00 |

They coincide wherever the trial identifies the mass and differ by $O(1)$ factors exactly where the readout is labelled apparatus-dominated. Behaviour on the cases is not what selects the estimator. Only a structural argument can.

---

## 3. What the ratio of means is

Because $f=m\alpha$ on the support,

$$
\widehat m=\frac{E[f]}{E[\alpha]}=\frac{E[\alpha\cdot m]}{E[\alpha]},
\qquad
\widehat\mu=\frac{E[\alpha]}{E[f]}=\frac{E[f\cdot\mu]}{E[f]}.
$$

So $\widehat m$ **is** a posterior mean of the mass — under the posterior reweighted by the latent acceleration magnitude. Pairs in which the object barely accelerates carry little information about $m$ and are down-weighted in proportion; the weight $\alpha$ is exactly what cancels the $1/\alpha$ singularity that makes the plain posterior mean $E[m]$ diverge. Dually $\widehat\mu$ is the force-weighted mean of the inverse mass. The two reweighted means are reciprocal not by construction of the readout but because $\alpha\,m=f$ and $f\,\mu=\alpha$ — the law itself makes the two leverages each other's numerators.

**Uniqueness among power weights.** Take $\widehat m_k=E[\alpha^k m]/E[\alpha^k]$ and its dual $\widehat\mu_k=E[f^k\mu]/E[f^k]$. Then $\widehat m_k\widehat\mu_k=1$ for every posterior iff

$$
E[\alpha^{k-1}f]\,E[f^{k-1}\alpha]=E[\alpha^k]\,E[f^k],
$$

which holds identically only at $k=1$. At $k=2$ the condition is $E[f\alpha]^2=E[f^2]E[\alpha^2]$ — Cauchy–Schwarz, with equality only when $m$ is constant on the support.

$k=2$ is not an arbitrary competitor. It is the Fisher-information weight (information about $m$ from a pair scales as $\alpha^2/\sigma_F^2$) and it is what posterior least squares returns: $\arg\min_m E|\mathbf F^*-m\mathbf a^*|^2=E[f\alpha]/E[\alpha^2]$, while $\arg\min_\mu E|\mathbf a^*-\mu\mathbf F^*|^2$ gives $1/\widehat\mu=E[f^2]/E[f\alpha]$. These bracket $\overline f/\overline\alpha$ and fail reciprocity badly (Appendix, Table A):

| case | $E[f\alpha]/E[\alpha^2]$ | $\overline f/\overline\alpha$ | $E[f^2]/E[f\alpha]$ |
|---|---|---|---|
| §18 | 80 | 125 | 126 |
| §21 zero–zero | 2.53 | 4.00 | 6.32 |
| §20 anti-aligned | 0.07 | 1.00 | 14.2 |

This is the Deming/errors-in-variables asymmetry — regress $F$ on $a$ or $a$ on $F$ — reappearing after the latent pair has been integrated out. The corpus already rejected that asymmetry at the best-point level; $\overline f/\overline\alpha$ is the classical ratio estimator, which is the standard resolution of it for regression through the origin. That is a good pedigree, but it is still an identity, not a derivation. The derivation is next.

---

## 4. The forcing argument: composition laws

A readout $R$ maps a compatible-pair support to a positive number. Newtonian mass composes in two ways, and each is a definitional operation elsewhere in the corpus.

**(L) Lumping.** If material objects $A_1$, $A_2$ are lumped into $A$ and move with a common acceleration, then $\mathbf F^*_A=\mathbf F^*_1+\mathbf F^*_2$ and $m_A=m_1+m_2$. The compatible-pair support of $A$ is the support of $(f_1+f_2,\ \alpha)$ under the joint law of the parts. Requirement:

$$
R(f_1+f_2,\ \alpha)=R(f_1,\alpha)+R(f_2,\alpha)
\quad\text{for every coupling with common }\alpha.
$$

**(R) Common force.** If one force acts on $A_1$ and on $A_2$ (the same acting object; or the mutual force of an interaction pair, equal in magnitude), relative acceleration is $\mathbf a^*_1+\mathbf a^*_2=(\mu_1+\mu_2)\mathbf F^*$: inverse masses add. Requirement:

$$
R(f,\ \alpha_1+\alpha_2)^{-1}=R(f,\alpha_1)^{-1}+R(f,\alpha_2)^{-1}
\quad\text{for every coupling with common }f.
$$

**Proposition.** Among readouts that are (i) determined by the law of the latent pair, (ii) covariant under change of force units and of acceleration units, (iii) exact on point supports ($R(\delta_{(f_0,\alpha_0)})=f_0/\alpha_0$), and (iv) continuous and monotone in the usual sense, the only one satisfying both (L) and (R) is

$$
R=\frac{E[f]}{E[\alpha]}.
$$

*Sketch.* Fix the acceleration side at a point mass $\delta_{\alpha_0}$. (L) says $R(\cdot,\delta_{\alpha_0})$ is an additive, law-determined, monotone functional of the force law. An additive law-determined functional is determined on indicators: $T(\mathbf 1_A)=h(P(A))$ with $h$ additive on disjoint events and monotone, hence $h(p)=cp$; extending through simple functions and continuity, $T=c\,E[\,\cdot\,]$. Units covariance and (iii) give $c=1/\alpha_0$. So $R(f,\delta_{\alpha_0})=E[f]/\alpha_0$. The same argument on (R) with the force side at a point gives $R(\delta_{f_0},\alpha)=f_0/E[\alpha]$. For general joint laws, (L) forces $R=E[f\,c(\alpha)]$ for some coefficient $c$ that may depend on the acceleration law, and (R) forces $R^{-1}=E[\alpha\,c'(f)]$; evaluating both on product laws forces $E[\alpha]\,E[c(\alpha)]$ and $E[f]\,E[c'(f)]$ to be constants across all laws, which (within the regular family $c=\varphi(\alpha)/E[\alpha\varphi(\alpha)]$) forces $\varphi$ constant, i.e. $c=1/E[\alpha]$. $\square$

The regularity assumptions in (iv) and in "the regular family" are where a sceptic should push; I do not claim more than a sketch. What is not in doubt is the **separating power** of the two axioms, which is checkable and was checked (Appendix, Table B):

| readout | (L) lumping | (R) common force |
|---|---|---|
| $\overline f/\overline\alpha$ | exact | exact |
| $E[f\alpha]/E[\alpha^2]$ (m-regression) | exact | fails ($-27\%$) |
| $E[f^2]/E[f\alpha]$ (μ-regression) | fails ($-30\%$) | exact |
| median$(m)$ | fails ($+20\%$) | fails ($+16\%$) |
| $\exp E[\log m]$ | fails ($+30\%$) | fails ($+26\%$) |
| RMS ratio | fails ($-15\%$) | fails ($-14\%$) |
| profile minimiser (old best point) | fails (by type — a mode; not run) | fails (by type; not run) |

Each axiom kills exactly what the other lets through. (L) alone admits the whole family $E[f\varphi(\alpha)]/E[\alpha\varphi(\alpha)]$ — including the m-regression — and rejects every nonlinear summary of $f$. (R) alone admits the dual family and rejects every nonlinear summary of $\alpha$. Their intersection is the ratio of means.

**Corollary.** Reciprocal invariance is a theorem: $R=E[f]/E[\alpha]$ and its (R)-dual $E[\alpha]/E[f]$ are reciprocals by inspection. §11 should be demoted from requirement to consequence, and §32 item 8 ("the point readout is the ratio of mean compatible magnitudes") moves from *chosen* to *forced, given (L) and (R)*. Item 7 (integrate rather than take the stationary root) is likewise forced, since the profile minimiser fails both axioms.

**On the premise.** A statistician will object that the estimate of a sum need not be the sum of the estimates, and that the Bayes-optimal readout under the natural reciprocal-symmetric loss $(\log\widehat m-\log m)^2$ is the geometric mean, which fails both axioms. The VD reply is the one `why_the_estimator_is_the_definition` already makes: on the empirical branch the estimator *is* the definition of mass. Mass, as defined elsewhere in the same dictionary, adds under lumping and its inverse adds under common force. An estimator that returns a non-additive quantity is defining a different word from the one the lumping and interaction entries use; the readout must compose as the definiendum composes, or E10 and the lumping entry are not about the same quantity. That is a VD-internal argument, not a statistical theorem, and should be recorded as such. It is also the strongest argument available: it is the only one that distinguishes the candidates at all.

Two cautions on (L) as tested. First, additivity holds for the *induced* law of $f_1+f_2$ propagated from the parts; if the lumped body is instead re-estimated from a fresh total-force measurement with a fresh flat measure, the two posteriors differ, and additivity is a statement about the readout functional, not about that re-estimation. Second, in $d$ dimensions the common acceleration fixes the common direction $\widehat{\mathbf u}$, so the force marginals given $\widehat{\mathbf u}$ add and linearity of expectation carries through the direction mixture; the spec's choice to integrate scalar magnitudes rather than take magnitudes of vector means is what makes this clean (§5 below).

---

## 5. What remains genuinely policy

**The latent-pair measure $\rho_\Pi$.** The flat measure on $(f,\alpha)$ is not the only natural one. The geometrically honest choice is the Riemannian volume of the compatible cone under the noise metric $\Sigma^{-1}$, which for isotropic errors in $d$ dimensions is

$$
\rho_{\mathrm{cone}}\propto\left(\frac{f^2}{\sigma_F^2}+\frac{\alpha^2}{\sigma_a^2}\right)^{(d-1)/2}df\,d\alpha\,d\Omega,
$$

flat only in $d=1$. In $(m,\alpha)$ coordinates it carries a factor $(m^2\sigma_a^2+\sigma_F^2)^{(d-1)/2}$ — a soft, symmetric-in-form weighting with the apparatus scale $\sigma_F/\sigma_a$ built in. Numerically it changes the readouts by under $5\%$ in every spec test case, and by nothing in the identifying regime (Appendix, Table C). So the flat measure is a policy choice, as §33 says, but not a load-bearing one for any case the corpus has argued.

One measure *is* excluded: the scale-invariant $df\,d\alpha/(f\alpha)$. It puts non-normalisable mass at $f=0$ and $\alpha=0$, so every zero reading collapses onto the exact boundary and the §17/§18 repairs disappear. The flat measure's real content is "a latent magnitude of exactly zero is not special", which is precisely the disagreement rule's repair of $F/0$. The measure and the rule are the same commitment stated twice.

**Mean of magnitudes vs magnitude of means.** The spec integrates the scalar $f$ over all directions. The alternative $|E[f\widehat{\mathbf u}]|/|E[\alpha\widehat{\mathbf u}]|$ is also reciprocal-invariant and exact-consistent, but the two vector means need not be codirectional, and direction uncertainty would shrink a scalar readout that has nothing to do with direction. The spec's choice is right and should be recorded as a choice.

**Point and interval come from different functionals.** $\widehat m$ is a ratio of means; the interval (§12 route 1) is quantiles of the induced $m$. Nothing guarantees the point lies in the interval, though in every test case it does. Since the induced median is already computed for the interval, report it beside $\widehat m$: the two coincide exactly when the trial identifies the mass and separate as it stops doing so (5.000/5.000; 0.160/0.135; 125/148; 13.1/15.9 on the perpendicular case). The ratio $\widehat m/\operatorname{median}(m)$ is a candidate criterion for the `weakly identified` / `apparatus-dominated` labels §33 leaves open.

---

## 6. Where the single-trial construction does not extend: Level 3

§11 holds for the flat pair measure because in $(m,\alpha)$ coordinates it is $\alpha\,dm\,d\alpha$ and in $(\mu,f)$ coordinates it is $f\,d\mu\,df$ — the same measure, flat in $m$ and flat in $\mu$ at once. That is a property of two dimensions. For $N$ trials sharing one mass the latent space is $(m,\alpha_1,\ldots,\alpha_N)$, and the two candidate "flat" measures

$$
\rho_A=\prod_i\alpha_i\,d\alpha_i\,dm,
\qquad
\rho_B=\prod_i f_i\,df_i\,d\mu
$$

pull back to each other with ratio

$$
\frac{\rho_B}{\rho_A}=m^{2N-2},
$$

equal only at $N=1$. A common-mass posterior therefore requires a measure choice that *is* a prior on $m$ versus $\mu$; the symmetric option (the $\Sigma^{-1}$-volume of the constraint manifold, or the geometric mean $m^{N-1}\rho_A$) exists but is not flat in either.

The clean way out is the one the linearity principle already dictates: do not build a common-mass posterior for the readout at all. Pool first moments,

$$
\widehat m_{1..N}=\frac{\sum_i\overline f_i}{\sum_i\overline\alpha_i},
$$

which is reciprocal-invariant, is the $N$-trial version of $E[\sum_i(f_i-m\alpha_i)]=0$, and weights each trial by its mean latent acceleration — its leverage on $m$ — which is the right weighting. But it must then be said plainly that the Level-3 *test* (§27, the profile $Q_{\mathrm{joint}}(m)$, which genuinely imposes one mass) and the Level-3 *estimate* (pooled marginals, which do not) are different objects. That mismatch is acceptable; leaving it unstated is not.

---

## 7. Gate calibration

$Q_{\min}$ is the squared Mahalanobis distance from the $2d$-component data vector to the compatible cone, a $(d+1)$-dimensional set ($f$, $\alpha$, and $d-1$ angles). At an interior fit, under the true model and small noise,

$$
Q_{\min}\sim\chi^2_{d-1}.
$$

Simulation with Newton II true (Appendix, Table D): in 2D, mean $0.94$ against $1$, $95\%$ quantile $3.80$ against $3.84$; in 3D, mean $2.25$ against $2$, $95\%$ quantile $6.23$ against $5.99$. The playground's $Q_{\mathrm{crit}}=9$ is therefore a $0.27\%$ false-rejection threshold in 2D and $1.1\%$ in 3D. The anti-aligned case's $Q_{\min}=25$ on one degree of freedom has $p\approx6\times10^{-7}$.

Two consequences the spec should carry. In one dimension there are **zero** interior degrees of freedom: $Q_{\min}=0$ whenever the signs agree, and the gate reduces to a sign test with $Q_{\min}=\min(\widetilde F^2/\sigma_F^2,\ \widetilde a^2/\sigma_a^2)$ when they do not. And boundary fits ($f\to0$ or $\alpha\to0$, i.e. exactly the zero-reading cases) have fewer effective degrees of freedom, so a $\chi^2_{d-1}$ threshold is conservative there — which is the correct direction for a gate whose job is not to reject weak trials.

---

## 8. Corrections owed to the corpus

1. **Spec §11.** Reciprocal invariance is a consequence, not a selecting requirement; every candidate passes it. Replace with the composition requirement of §4 above, or at least stop presenting §11 as what demotes the median.
2. **Spec §23.** The three reasons given against the median target the earlier excitation-model median. The induced-mass median under the spec's own $\rho_\Pi$ needs none of the three inputs listed. The reason it is not the readout is (L)/(R), and §23 does not say so.
3. **Spec §32.** Items 7 and 8 are forced given (L) and (R), which are themselves definitional commitments of the corpus, not estimator choices. The list can be split as previously suggested, but the split lands differently than a "forced by the law / chosen by the estimator" reading would put it: the composition laws are part of what "the law" fixes.
4. **`why_the_estimator_is_the_definition` §2, Result 1.** "The estimator always returns exactly one positive mass" is correct as a root count and wrong as a claim about the estimator: for anti-aligned data the positive root of the Deming quadratic is the profile's *maximum* (spec §7; $Q_{\mathrm{prof}}(m)=25(1+m)^2/(1+m^2)$ on the §20 data, maximal at $m=1$). The argument form of that record survives; the estimator it identifies is superseded by this one, and the record's "200,000 null realisations gave a positive root" was counting roots, not minima.
5. **Spec §7 / Part VI framing.** The old best-point estimator is reciprocal-invariant. Its defect is boundary-seeking (and failing (L)/(R)), not asymmetry.

---

## 9. Verdict

$\widehat m=\overline f/\overline\alpha$ is the right readout. It is the unique estimator on the compatible-pair support that composes as Newtonian mass composes under both lumping and common force; it is the leverage-weighted posterior mean of the mass; its reciprocal symmetry is a theorem of the law rather than a design constraint; and its behaviour on every test case in the corpus is reproduced to the precision of the spec's closed forms. The guesswork the spec admits to was in the *justification*, not in the form — the form was being selected by the composition laws all along, and the reciprocal requirement was the visible shadow of that selection.

What is still open is exactly what §33 says is open, minus one item (the gate's calibration, now given) and plus one (the Level-3 measure asymmetry, now identified).

---

## Appendix — numerical checks

Harness: grid quadrature in $(\theta,f,\alpha)$ for the 2D posterior, $720\times600\times600$ points (midpoint grids for log moments); $Q_{\min}$ by direction scan with closed-form $(f,\alpha)$ optimum; Monte Carlo ($3\times10^6$ samples, skewed positive test laws) for the composition tests; $2000$ (2D) and $400$ (3D) simulated true-model trials for the gate.

**Table A — spec test cases, all readouts.** $\widehat m\widehat\mu=1$ to $10^{-4}$ in every row.

| case ($\sigma_F,\sigma_a$) | $Q_{\min}$ | $\overline f/\overline\alpha$ | median | geo | RMS | m-reg | μ-reg$^{-1}$ | $[q_{16},q_{84}]$ |
|---|---|---|---|---|---|---|---|---|
| $(50,0),(10,0)$; $1,0.2$ | 0.00 | 5.000 | 5.000 | 5.000 | 5.000 | 4.998 | 5.002 | $[4.86,5.14]$ |
| §17 $(0,0),(10,0)$; $2,0.5$ | 0.00 | 0.160 | 0.135 | 0.107 | 0.200 | 0.159 | 0.251 | $[0.04,0.28]$ |
| §18 $(50,0),(0,0)$; $2,0.5$ | 0.00 | 125.2 | 148.0 | 187.4 | 100.0 | 79.7 | 126.1 | $[71,509]$ |
| §21 $(0,0),(0,0)$; $2,0.5$ | 0.00 | 4.000 | 4.000 | 4.000 | 4.000 | 2.533 | 6.317 | $[1,16]$ |
| §20 $(-10,0),(10,0)$; $2,2$ | 25.00 | 1.000 | 1.000 | 1.000 | 1.000 | 0.070 | 14.20 | $[0.01,75]$ |
| perp $(10,0),(0,2)$; $2,0.5$ | 16.00 | 13.10 | 15.92 | 18.78 | 10.78 | 7.83 | 14.93 | $[6.7,56]$ |
| $30^\circ$, $10$ vs $2$; $1,0.2$ | 13.40 | 5.000 | 5.000 | 5.000 | 5.000 | 4.943 | 5.058 | $[4.29,5.82]$ |

Spec closed forms: §17 $\sigma_F\sqrt{2/\pi}/10=0.1596$; §18 $50/(\sigma_a\sqrt{2/\pi})=125.3$; §21 $\sigma_F/\sigma_a=4$; analytic half-normal median for §18 $=148.3$, geometric mean $=188.7$.

**Table B — composition tests** (relative deviation from additivity; skewed exponential/lognormal test laws with genuine coupling).

| readout | (L) $R(f_1{+}f_2,\alpha)/[R_1{+}R_2]-1$ | (R) $R(f,\alpha_1{+}\alpha_2)^{-1}/[R_1^{-1}{+}R_2^{-1}]-1$ |
|---|---|---|
| ratio of means | $-0.0000$ | $-0.0000$ |
| median$(f/\alpha)$ | $+0.198$ | $+0.159$ |
| $\exp E[\log f/\alpha]$ | $+0.300$ | $+0.257$ |
| RMS ratio | $-0.153$ | $-0.139$ |
| m-regression | $-0.0000$ | $-0.268$ |
| μ-regression | $-0.301$ | $+0.0000$ |

**Table C — measure sensitivity**, $\widehat m$ under flat vs cone-volume $\rho$ ($d=2$, $d=3$): strong $5.000/5.000/5.000$; §17 $0.1589/0.1587/0.1585$; §18 $125.9/126.0/126.1$; §21 $4.000/4.000/4.000$; §20 $1.000/1.000/1.000$; perpendicular $13.16/13.50/13.81$.

**Table D — gate under the true model.** $d=2$: mean $Q_{\min}=0.942$, $P(Q_{\min}>9)=0.0005$ (χ²₁: $0.0027$), $95\%$ quantile $3.80$ (χ²₁: $3.84$). $d=3$: mean $2.246$, $P(Q_{\min}>9)=0.0125$ (χ²₂: $0.0111$), $95\%$ quantile $6.23$ (χ²₂: $5.99$).

**Symbolic check** (sympy): pull-back of $\prod_i f_i\,df_i\,d\mu$ to $(m,\alpha_i)$ coordinates divided by $\prod_i\alpha_i\,d\alpha_i\,dm$ equals $1,\ m^2,\ m^4$ for $N=1,2,3$.

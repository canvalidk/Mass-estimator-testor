# Where the mass estimator's parts come from: credit and originality

**2026-09-26.** Record, filed at the user's request ("Will you do a record on
where the parts come from. So an originality honest account of credit and
originality"), written by Claude from the literature check of the morning of
26 September. It builds on `mass_estimator_originality_research_2026-09-05.md`
(the 5 September deep-research assessment) and extends it to the current
equation (ours25b). It edits nothing and binds nothing.

**What "not found" means here.** The sources below were read in the ways
listed in §7: some in full from `papers/`, some only through a secondary
source, some from general knowledge. About fifteen web searches were run on 26
September, on top of the 5 September search. "Not found" means not found in
that bounded search, not absent from the literature.

**In one paragraph.** The estimator has two stages: a **law** (a curve over
candidate masses, built from the readings) and a **readout** (one number taken
from that curve). The law is not the project's. Its fit measure is textbook
errors-in-variables. Its weighting of the true points, which makes the law
come out clean, is Beale's remark in Lindley's 1966 discussion of Sprent. Its
prior shape matches known default priors (Liseo's reference posterior in 1D;
uniform-angle priors in Lindley & El-Sayyad, Creasy and Hogg). The readout is
built from a known device, a ratio of posterior means, which Little uses as an
approximation and which is a size-biased mean in Arratia et al.'s sense. What
was not found anywhere is **the selection principle**. Two things are
unlocated: that the composition laws of Newtonian mass force exactly that
readout; and that axis splitting forces both the weighting of the true points
and, with many readings, the stacked-length pooling. The reading of the whole
construction as the definition of mass is the project's own position. The
honest summary is: the models are borrowed; the reasons for choosing them, and
the readout's exact selection, are the project's.

---

## 1. The two stages, in plain terms

- **The law**, $p(m)\propto\dfrac{s}{s^2+m^2}\exp\!\Big[-\dfrac{\sum_i|\tilde{\mathbf F}_i-m\tilde{\mathbf a}_i|^2}{2(\sigma_F^2+m^2\sigma_a^2)}\Big]$,
  is a curve saying how believable each candidate mass is. Intervals and the
  median come from it alone.
- **The readout**, $\widehat m=E\big[|\mathbf F^\ast|\big]/E\big[|\mathbf a^\ast|\big]$
  over the candidate true pairs behind that curve (stacked over readings), turns
  the curve into the one number the mass entry returns.

The peak of the same curve is the Deming / total-least-squares slope. The
curve's plain mean of $m$ is infinite. Its median is another option. Which
single number to report is where the project parts from the literature.

## 2. Part by part

| # | Part | Earliest source found | Read how | What is the project's |
|---|---|---|---|---|
| 1 | The fit measure $\sum_i(F_i-m a_i)^2/(\sigma_F^2+m^2\sigma_a^2)$ | Kummell (1879), Deming (1943); Lindley (1947) for functional relationships | general knowledge; Lindley 1947 only as cited by Sprent | nothing |
| 2 | Integrating over the true points instead of fitting them | Lindley (1966, discussion of Sprent); Lindley & El-Sayyad (1968) | read | the argument from totality (the best fit runs off at F/0 and on anti-aligned data), which L&ES make in general form against maximum likelihood |
| 3 | The weighting of the true points that makes the extra factor unity: data scaled to unit error variances, distance from the origin independent of the slope | Beale, reported in Lindley's contribution (1966, p. 290) | read | its derivation: from the 24 Sept premise (force and acceleration symmetric, either can take any value) and, decisively, from axis splitting (25 Sept) |
| 4 | A flat prior on each true point adds a per-point factor pulling the slope to 0 that does not fade with n | Lindley (1966, p. 289) | read | rediscovered as the drift of eq. 18 (contribution 12; 25 Sept convergence argument); credit is Lindley's |
| 5 | A uniform prior on the angle of the line | L&ES (1968, eq. 13, as an example, data units); Creasy (1956, slope as tan φ, equal known variances); Hogg, Bovy & Lang (2010, §7, data units) | L&ES and Hogg read; Creasy only via Schneeweiss (1982) | the angle measured in noise units, so the prior is centred on $\sigma_F/\sigma_a$; that centre forced by the zero-reading argument; positivity $m>0$. The shape itself is still an open step |
| 6 | The 1D single-reading law as a whole | Liseo (2003) eq. (14): the reference-prior posterior (Liseo 1993; Berger & Bernardo 1992), flat in the nuisance mean | read | nothing beyond the restriction to $m>0$ and scaling to noise units. (The 21 Sept law matched Liseo's Jeffreys posterior, eq. 12, instead.) |
| 7 | Positive magnitudes from zero readings; the zero–zero answer $\sigma_F/\sigma_a$; the half-Cauchy | Katz (1961); Roe–Woodroofe; normal-ratio theory; Liseo's null posterior (5 Sept record) | 5 Sept record | the 25 Sept argument that *any* always-defined estimator of this type must return $\sigma_F/\sigma_a$ at zero readings (dimensional analysis plus the swap symmetry); not located, and modest |
| 8 | A ratio of posterior means as a readout | Little (survey inference, as an approximation dropping $O(1/n)$ terms); size-biased mean, Arratia, Goldstein & Kochman §2.3 | Little only via the 5 Sept record; Arratia cited in note d1 | **using it exactly**, including at weak signal where the mean of the ratio diverges |
| 9 | The composition characterisation: lumping plus force/acceleration symmetry (any two of lumping, common force, symmetry) force $E[f]/E[\alpha]$ | **not found** (5 Sept search and this one) | — | **candidate contribution 1**: the 2 Sept argument and the 6 Sept theorem (C1–C5) |
| 10 | With many readings, the size of the series is the length of the stacked vector | **not found** | — | **candidate contribution 2** (small): forced by axis splitting given one joining rule (25 Sept) |
| 11 | Axis splitting as a premise that forces the weighting and the pooling | treating vector components as separate points is standard practice; the premise-and-consequence form **not found** | — | **candidate contribution 3**: the 25 Sept records; Herschel–Maxwell (1850s–60) supplies one step |
| 12 | The four-number state $(C,T,M)$ | the law's part, $C$, is a function of the classical sufficient statistics $S_{xx},S_{yy},S_{xy}$ of the functional relationship; the readout's extra $(T,M)$ **not found** | general knowledge | the readout part, as a consequence of #10 |
| 13 | The estimator read as the definition of mass on measured data | metrology: a result is a value plus uncertainty and information (VIM 2.9, 5 Sept record); the constitutive claim itself not found | 5 Sept record | the VD position (`why_the_estimator_is_the_definition`, 1 Aug; steering, 7 Sept). A philosophical position, not a theorem |
| 14 | No always-finite interval with guaranteed coverage | Gleser & Hwang (1987); Koschat (1987); Fieller (1954) | read | nothing |

**Mechanics and metrology.** Mass estimation from force and acceleration
already exists.

- Butler (2004) estimates a lithography stage's mass from controller force and
  acceleration, and discusses the constant-velocity phase where both are zero
  (5 Sept record).
- Robotics identifies payload inertial parameters by ordinary least squares,
  reports infeasible and negative masses under acceleration noise (Lambert
  2023), and repairs them with constraints (Wensing et al. 2017; not read).
- NIST calibrates force transducers by errors-in-variables regression
  (Bartel, Stoudt & Possolo). Their criterion, cut down to a line through the
  origin with the true forces minimised out, is the fit measure of #1, and
  they report its minimum.

None of these uses the law of #3–#6 or a composition-forced readout. The
project's potential difference is the positive, always-defined readout and
how the VD admits its result, as the 5 Sept record already said.

## 3. Where the project parts from the literature, and why that matters

- **The law: no disagreement.** The project's premises land on the Lindley–Beale
  treatment, which was chosen in 1966 for symmetry. Landing there from other
  premises is evidence the premises are sound, not a failure of originality.
- **The readout: a real disagreement.** Standard Bayesian practice reports the
  curve, or the number that minimises a declared loss. Under the natural
  symmetric loss that is the geometric mean, which fails both composition
  laws. The project instead requires the reported number to compose as mass
  composes, because on measured data the estimator is what "mass" means. That
  is a coherence requirement between dictionary entries, not an optimality
  claim. A statistician can reasonably reject it; the project's answer is the
  composition argument.
- **Scope of the composition claim.** It holds within one propagated joint law
  (parts of one object, or bodies whose uncertainty is carried together). It
  does not say that three separate experiments on A, B and A+B give estimates
  that add, and they need not.

## 4. Inside the project: who contributed what

From the records' own provenance lines.

| Piece | Source |
|---|---|
| The programme, the method (forced moves, "hand on the left wall"), the premises judged | the user |
| σ_F, σ_a > 0 as part of the type | the user, 25 Sept |
| The estimator is the definition | `why_the_estimator_is_the_definition` (1 Aug); steering record (the user's ruling, 7 Sept) |
| Ratio of mean magnitudes as the readout; reciprocal invariance | the user's spec d1 (30 Aug) |
| Composition, not reciprocity, forces the readout | Claude, 2 Sept, answering the user's ask |
| The C1–C5 uniqueness theorem | a VD-Newton Codex session, 6 Sept |
| The 21 Sept equation | VD-Newton (`mass_estimator_equation/`) |
| The symmetric rule and its comparison with eq. 18 / eq. 20 | contribution 12, VD-Newton lab (23 Sept) |
| The flat-push premise (P1, P2) | the user, 24 Sept; its vector reading and the Herschel–Maxwell route: Claude |
| Axis splitting | proposed by Claude, judged by the user ("too elementary to be wrong"), 25 Sept |
| Convergence proof, width formula | Claude, 25 Sept |
| Zero-reading σ_F/σ_a argument; λ as a fake reading; stacked readout; equal-energy result; four-number state | Claude, 25–26 Sept, in conversation with the user's rulings and reorderings |
| The literature positions in this record | Claude, 26 Sept; the 5 Sept deep research before it |

Many derivations were made by AI sessions under the user's direction. How that
is acknowledged in public work is the user's call.

## 5. What can honestly be claimed

> Starting from Newton's second law, positive mass, measured inputs with
> declared uncertainty, and two consistency requirements (the estimate must
> compose as mass composes; a vector reading and its components must agree),
> we derive a total estimator of inertial mass. Its law is the classical
> Lindley–Beale functional-relationship posterior, with a default prior pinned
> to the instruments' noise ratio; in one dimension it is the reference-prior
> posterior for Fieller's problem. Its readout, a ratio of posterior means of
> the true force and acceleration sizes, is selected exactly by the
> composition requirements, and the many-reading pooling by the component
> requirement. To our knowledge the selection results are new; the models are
> not.

Not to be claimed: a new statistical model; a better estimator in the
decision-theoretic sense; guaranteed coverage (Gleser–Hwang rules it out); that
separately measured masses must add.

## 6. What could change this assessment

- Little's use of a ratio of posterior expectations, read directly (only seen
  through the 5 Sept record).
- A targeted search, or a statistician's opinion, on composition-type or
  consistency-type characterisations of point estimators.
- The robotics paper not read (arXiv 2507.15604, rate-limited), Wensing et al.
  (2017), and Butler (2004) read directly.
- Creasy (1956) read directly, for whether its angle construction is in
  standardised units.

## 7. Sources and how each was read

- **Read in full from `papers/` (26 Sept):** Lindley, contribution to the
  discussion of Sprent (1966), JRSS B 28(2):288–297 (`jrsssb_28_2_288.pdf`,
  fetched 26 Sept), including Beale's remark on p. 290: with unit-variance
  scaling "the factor outside the exponential would have been unity"; Lindley &
  El-Sayyad (1968), JRSS B 30(1):190–202, pp. 192 and 201; Hogg, Bovy & Lang
  (2010), arXiv:1008.4686, §7 (fetched 26 Sept); Liseo (2003), METRON
  61(1):133–150, pp. 137–139; Berger, Liseo & Wolpert (1999), Statist. Sci.
  14(1):1–28, Example 5 and §4.3; Bartel, Stoudt & Possolo, NIST force
  calibration, §3.
- **Web pages read on 26 Sept:** Schneeweiss (1982), J. Multivariate Anal., on
  Creasy (1956); Yin & Ghosh (2001), abstract; Lambert (2023), MIT.
- **Through the 5 Sept record only:** Little; Butler (2004); Katz (1961) and
  Roe–Woodroofe; VIM 2.9; Arratia, Goldstein & Kochman.
- **General knowledge, not checked:** Kummell (1879); Deming (1943); the
  classical sufficient statistics of the functional relationship;
  Herschel–Maxwell.

# Old information in the prior slot

**2026-09-24.** Record, filed at the user's request. Written by Claude from
his idea stated in conversation that day: *use the exact combination equation
to separate old information from new, and the old information lands in the
places the parameters occupy.*

**What this record is.** An algebraic argument, stated once, and a
reproducible test of how much of it survives compression. The test is a
preset of the mass estimator tester (`canvalidk/Mass-estimator-testor`),
so anyone can rerun it with one command.

**What it scopes.** Contribution 11 (`exact_mass_curve_combination.md`,
17 Sept: exact combination, eqs. 17–18, 20, 28–29) and contribution 12
(`mass_estimator_regimes_priors_and_updating_2026-09-23.md`: λ as a prior
factor centred at the instrument ratio, the update-rule reading, the
symmetric rule). It binds nothing.

---

## 1. The claim

For readings that share only the mass, the exact law of all the readings
factorises into *old × new*. The old part is a function of mass alone, so
it enters the new reading's calculation exactly where a prior does. Two
numbers of it, the geometric mean $m_0$ and the log-spread $\sigma_{\ln m}$,
fill the slots of contribution 12's prior parameters. Carried that way:

- they reproduce the exact answer to within about 1% when both channels
  of the old readings are informative;
- they lose up to about 16% (median distance) when either channel is weak;
- a prior family that contains the zero-reading law holds weak old
  information much better than a log-normal does.

On this reading the tilt parameters are not fitted. They are state: they
start at $(m_0, \sigma_{\ln m}) = (s, \pi/2)$ and every reading updates them.

## 2. The algebra

### 2.1 Old × new

Write $u = \log m$. For reading $i$ the per-reading law gives a likelihood
of mass $L_i(m)$ (its own latent magnitudes and direction integrated out, under a
declared nuisance measure) and a conditional mean acceleration
$a_i(m) = E[\alpha_i \mid m]$. A combination rule that multiplies
likelihoods and states one prior $\pi$ gives

$$
p_{1:N+1}(m) \;\propto\; \pi(m)\prod_{i=1}^{N+1} L_i(m)
\;=\; \underbrace{\Big[\pi(m)\prod_{i=1}^{N} L_i(m)\Big]}_{\propto\; p_{1:N}(m)}\; L_{N+1}(m).
$$

The bracket is the old law. In the calculation for reading $N+1$ it stands
exactly where a prior stands. Nothing is approximated.

### 2.2 What the ratio-of-means readout needs as well

The readout is contribution 11's eq. (28),

$$
\widehat m \;=\; \frac{E\big[\sum_i f_i\big]}{E\big[\sum_i \alpha_i\big]}
\;=\; \frac{\int m\,A(m)\,p(m)\,dm}{\int A(m)\,p(m)\,dm},
\qquad A(m) = \sum_i a_i(m),
$$

and $A$ splits the same way: $A_{1:N+1} = A_{1:N} + a_{N+1}$. So the
exact old state is the pair $\big(p_{1:N}, A_{1:N}\big)$, which is
contribution 11's two-function output (eq. 29) read as *old information*.
A prior slot carries only $p$. Dropping $A_{1:N}$ leaves the median and
every interval exact, but moves the ratio of means.

### 2.3 Every rule splits this way; they differ in what each new reading brings

(This corrects a claim made in conversation, that only the count-the-prior-once
rule splits cleanly.) All three rules compared in contribution 12 factorise
as old × new. What differs is the factor the new reading contributes. With
$\theta$ defined by $m = s\tan\theta$ and a common instrument ratio $s$
(contribution 12, §7, verified in the tester's tests):

| rule | the new reading contributes | extra per reading, relative to symmetric |
|---|---|---|
| symmetric (radius nuisance, prior once) | $L_{N+1}$ | none |
| eq. (18), mass matching | $L_{N+1}\cos^2\theta$ | $\cos^2\theta$, which pulls towards $m \to 0$ |
| eq. (20), log matching | $L_{N+1}\,\tfrac12\sin 2\theta$ | the half-Cauchy again, which pulls towards $m = s$ |

Eq. (18) is also exactly the likelihood product with *acceleration* as each
reading's nuisance and a flat prior on $m$ stated once (tested). So the
choice that makes eq. (18) converge to the wrong mass is its nuisance split,
not a double-counted prior. Eq. (20) does re-count the reference with every
reading. Only the symmetric rule adds nothing per reading, so only there is
the old information *all* of the prior.

### 2.4 The slots: $(m_0, \sigma_{\ln m})$, not $\lambda$

Contribution 12's tilt, multiplied onto the flat law's half-Cauchy, is the
family

$$
\rho(z) \;\propto\; \operatorname{sech} z\; \exp\!\big(\lambda\,[\operatorname{sech} z - 1]\big),
\qquad z = \log(m/m_0).
$$

It is symmetric in $z$, so its mean in $\log m$ is $\log m_0$. Its spread
$\sigma_{\ln m}(\lambda)$ decreases monotonically in $\lambda$:

| $\lambda$ | $-4$ | $0$ | $0.35$ | $1$ | $4$ | $16$ | $64$ | $1000$ |
|---|---|---|---|---|---|---|---|---|
| $\sigma_{\ln m}$ | 2.853 | $\pi/2$ = 1.571 | 1.457 | 1.260 | 0.666 | 0.264 | 0.127 | 0.0316 |

Near $z = 0$, $\log\operatorname{sech} z \approx -z^2/2$ and
$\operatorname{sech} z - 1 \approx -z^2/2$, so for large $\lambda$ the law is
close to normal with $\sigma_{\ln m} \approx 1/\sqrt{1+\lambda}$. The "+1" is
the half-Cauchy's own curvature. The approximation is within 2% from
$\lambda \approx 64$.

So $(m_0, \lambda)$ is a reparametrisation of (old geometric mean, old
log-spread). Three readings follow:

- $\lambda = 0$, $m_0 = s$ is the zero-reading law: $\sigma_{\ln m} = \pi/2$
  is what *knowing nothing* looks like here.
- $\lambda < 0$ is a prior claiming less than the zero-reading law. That is
  what contribution 12's negative fits below the operating range were doing.
- Contribution 12's fitted $\lambda^* = 0.35$ is $\sigma_{\ln m} = 1.457$:
  about 7% narrower than knowing nothing. (A conversation estimate of
  "about 1.7" used the large-$\lambda$ formula far outside its range and
  is withdrawn.)

The spread is stated in $\log m$ because the ordinary SD of $m$ under the
flat law is infinite (equation README, 21 Sept). $\sigma_m$ would name
nothing.

## 3. The test

### 3.1 Design

Preset `presets/old_info_as_prior.yaml` in the tester. Its world:

- readings share only the mass (`new_excitation`, random direction per
  reading), so the old information is about the mass alone;
- 3D, isotropic unit noise in both channels, so $s = 1$;
- mass $\in \{s/4, s, 4s\}$, acceleration SNR per reading $\in \{1, 3\}$,
  old readings $\in \{1, 5, 20\}$;
- 18 cells, 1,000 series per cell.

In each series the last reading is new and the rest are old. Four estimators
see the same series, all with the symmetric rule for the old readings:

| estimator | old information carried |
|---|---|
| `exact` | $p_{1:N}$ and $A_{1:N}$ (equals the symmetric rule on all readings; tested) |
| `curve_only` | $p_{1:N}$ exactly, $A_{1:N}$ dropped |
| `lognormal_fit` | $(m_0, \sigma_{\ln m})$ as a normal law in $\log m$ |
| `tilt_fit` | $(m_0, \sigma_{\ln m})$ as the tilt family, $\lambda$ from the table above |

Each carry is scored against `exact` series by series as
$|\log(\text{carry}/\text{exact})|$, and against the true mass.

### 3.2 Results

Distance from the exact answer, median / 90th percentile over series. Intervals are
omitted here; the full summary has them.

<!-- generated from run.json, code 2c9da55, 1000 series per cell -->
| mass / s | acc. SNR | old readings | curve only: ratio of means | log-normal fit: median | tilt fit: median | 95% coverage: exact / log-normal / tilt |
|---|---|---|---|---|---|---|
| 0.25 | 1 | 1 | 3.7% / 11% | 12% / 21% | 1.6% / 9.7% | 93 / 94 / 95 |
| 0.25 | 1 | 5 | 3.7% / 11% | 15% / 24% | 8.4% / 20% | 84 / 86 / 87 |
| 0.25 | 1 | 20 | 2.3% / 7.4% | 15% / 26% | 12% / 24% | 81 / 82 / 83 |
| 0.25 | 3 | 1 | 1.6% / 4.6% | 16% / 24% | 6.9% / 21% | 96 / 96 / 96 |
| 0.25 | 3 | 5 | 1.1% / 3.2% | 15% / 27% | 13% / 25% | 96 / 95 / 94 |
| 0.25 | 3 | 20 | 0.41% / 1.1% | 4.9% / 17% | 5.1% / 16% | 92 / 94 / 94 |
| 1 | 1 | 1 | 3.2% / 10% | 11% / 20% | 1.7% / 9.7% | 99 / 99 / 99 |
| 1 | 1 | 5 | 2.6% / 8.3% | 10% / 22% | 5.7% / 17% | 92 / 97 / 95 |
| 1 | 1 | 20 | 1.2% / 4.2% | 2.7% / 15% | 2.3% / 13% | 81 / 94 / 92 |
| 1 | 3 | 1 | 0.63% / 2.1% | 5.6% / 16% | 2.5% / 9.8% | 94 / 96 / 94 |
| 1 | 3 | 5 | 0.28% / 0.76% | 0.39% / 1.2% | 0.31% / 1% | 93 / 95 / 95 |
| 1 | 3 | 20 | 0.087% / 0.22% | 0.041% / 0.11% | 0.039% / 0.11% | 94 / 94 / 95 |
| 4 | 1 | 1 | 0.82% / 2.8% | 17% / 26% | 10% / 24% | 97 / 97 / 96 |
| 4 | 1 | 5 | 0.53% / 1.6% | 11% / 26% | 10% / 22% | 95 / 95 / 95 |
| 4 | 1 | 20 | 0.2% / 0.53% | 2.5% / 6.5% | 2.5% / 6.5% | 93 / 95 / 95 |
| 4 | 3 | 1 | 0.057% / 0.15% | 4.4% / 12% | 4.1% / 11% | 95 / 94 / 94 |
| 4 | 3 | 5 | 0.022% / 0.06% | 0.99% / 1.7% | 0.98% / 1.7% | 95 / 95 / 95 |
| 4 | 3 | 20 | 0.0072% / 0.017% | 0.25% / 0.32% | 0.25% / 0.31% | 94 / 95 / 95 |

Checked identities (exact to numerical precision in every cell): the
`curve_only` median and 95% interval equal `exact`'s.

### 3.3 What the numbers say

- **Strong old information fits in two numbers.** Where both channels are
  informative (mass $s$ or $4s$ at acceleration SNR 3, so force SNR 3 or 12)
  with 5 or 20 old readings, both fits sit within 0.04–1.0% (median) of
  exact. The old law is then close to normal in $\log m$.
- **Weak old information does not.** Where either channel is weak (any cell
  at acceleration SNR 1, or mass $s/4$, whose force SNR is 0.25 or 0.75) the
  fits are 1.6–16% off (median), up to 27% at the 90th percentile. More
  weak readings help only slowly: at mass $s/4$, SNR 3, the fits are still
  5% off with 20 old readings.
- **Shape matters when information is weak.** With one old reading at
  SNR 1, the tilt fit is 1.6% off against the log-normal's 12% (mass $s/4$),
  and 1.7% against 11% (mass $s$). The tilt family contains the
  zero-reading law, and its heavy tails match a weakly informed old law.
  With more weak old readings the advantage shrinks.
- **$A_{1:N}$ matters only at weak signal.** Dropping it moves the ratio
  of means 0–4% (median), up to 11% at the 90th percentile in the weakest
  cells.
- **Point accuracy barely moves.** Across all 18 cells the within-×1.25
  success of every carry is within 6 points of exact (largest gaps: mass
  $4s$, SNR 1, 5 old readings, 42% against 36%; mass $4s$, SNR 3, 1 old
  reading, 65% against 60–61%).

## 4. What this does not show, and what it owes

- **The updating is textbook.** Turning yesterday's posterior into today's
  prior is Bayes' theorem applied twice. Compressing it by matching moments
  is assumed density filtering, the family expectation propagation grew out of
  (Opper; Minka). These should be cited, not claimed. What is specific here:
  the reading of contribution 12's fitted $\lambda$ as carried state; the
  measured advantage of a family containing the zero-reading law; §2.3's
  account of what each rule adds per reading.
- **The starting point is still the instrument.** The state starts at
  $(s, \pi/2)$. Carrying moves the instrument dependence into the first
  prior; it does not remove it.
- **What is carried is itself miscalibrated at weak signal.** In three weak
  cells the exact law undercovers: 95% coverage of 81–84% (mass $s/4$ at
  SNR 1 with 5 and 20 old readings; mass $s$ at SNR 1 with 20). This is the
  symmetric-rule undercoverage found on 24 Sept. Compressing
  it faithfully compresses that too.
- **An oddity, not a result.** In some weak cells the fits cover better than
  exact: 94% (log-normal) against 81% at mass $s$, SNR 1, 20 old readings.
  But the lightest weak cells barely move (81% against 82–83%). If it holds,
  it says the exact law's *shape* is what is miscalibrated there.
- **Same-pair data were not tested, and would behave differently.** When
  old and new readings share the latent pair, the old information is not a
  function of mass alone. It lands in the data slot (the pooled pair with
  $\Sigma/N$, contribution 11, eq. 5), not the prior.
- **Owed before this is presented as a result:** a second seed and a denser
  grid; readings from different instruments (different $s$ per reading);
  and real data.

## 5. Reproduction

Repository `canvalidk/Mass-estimator-testor`. The results above were
produced at commit `2c9da55`; `cfd6714` adds only documentation.

```sh
python -m pip install -r requirements.txt
python -m pytest -q tests/test_sequential.py     # the identities of §2 and §3.2
python -m met run presets/old_info_as_prior.yaml # the study of §3 (about 7 minutes)
```

Series $i$ of each cell is drawn from `SeedSequence(20260926, spawn_key=(cell_key, i))`,
so the run reproduces exactly. The run writes `summary.md` (every cell,
including intervals and coverage for all four estimators), `run.json` and
`values.npz` (every readout of every series).

Code: the sequential rule and its four carries are in `met/analyst.py`
(`_sequential_log_density`, `tilt_sd`, `tilt_lambda_for_sd`); the
agreement score is in `met/scores.py`.

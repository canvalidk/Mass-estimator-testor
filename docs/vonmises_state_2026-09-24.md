# The von Mises state: old readings as one doubled-angle vector

**2026-09-24.** Record, filed at the user's request. Written by Claude from
the conversation of that afternoon. It follows
`old_information_in_the_prior_slot_2026-09-24.md` (same day, same
discussion) and scopes it; it binds nothing.

**In one paragraph.** Each force/acceleration reading's likelihood of mass is,
in two dimensions exactly, a von Mises law in the doubled angle $2\theta$
(where $m = s\tan\theta$). Contribution 12's tilt $\lambda$ is a von Mises
concentration of the same kind. So readings and prior are the same kind of
object, and with a common instrument ratio they combine by vector addition:
all the old readings reduce to one vector (a direction and a length). In 2D
that vector carries the old information exactly, for the median and every
interval. In 3D it carries it approximately, and in most cells better than
the two-number fits of the morning. Its peak is the total-least-squares
axis. The estimator itself is unchanged; what this gives is an exact, compact
summary of repeated readings that share only the mass.

---

## 1. How it was found

1. **Morning.** The morning record showed that the exact law of $N+1$
   readings is the old readings' law times the new reading's likelihood, so the old
   information sits exactly in the prior slot. It then tried two-number
   summaries of the old law (a log-normal, and contribution 12's tilt shifted
   in $\log m$): good to about 1% when the old readings were strong, 1.6–16%
   off when they were weak.
2. **The naming step.** The user asked whether $\lambda$ is better called
   $\sigma_m$. Answer: it is a reparametrisation of the spread of $\log m$,
   $\sigma_{\ln m}$ (the ordinary SD of $m$ is infinite under the flat law).
3. **The user's question, 16:41:** *"$\lambda$ has some natural form; why can't
   we express the uncertainties of $f$ and $a$ also in that form? See if that
   simplifies the expression."*
4. **The rewrite.** The data enter each reading's likelihood only through
   $h(\theta)$, and $h^2$ rewritten in $2\theta$ has exactly the form of the
   tilt (§2). A numerical check on random readings: after removing the von Mises part,
   the log-likelihood's leftover was identically zero in 2D and up to about 2
   in 1D and 3D.
5. **The test.** A `vonmises_state` carry was added to the tester and run
   against the exact law, the whole old curve and the two morning fits, in 2D
   and 3D (§4).

The idea, and the conjecture that the measurement uncertainties share
$\lambda$'s form, are the user's. The algebra, code and study are Claude's.

## 2. The derivation

**Coordinates.** Standardise each reading, $\mathbf x = \mathbf F/\sigma_F$,
$\mathbf y = \mathbf a/\sigma_a$ (vectors in $d$ dimensions), and write
$P = |\mathbf x|^2$, $Q = |\mathbf y|^2$, $D = \mathbf x\cdot\mathbf y$,
$s = \sigma_F/\sigma_a$. The latent pair is
$(f/\sigma_F, \alpha/\sigma_a) = (r\sin\theta, r\cos\theta)$ along a common
unit direction, so $m = s\tan\theta$ with $\theta \in (0, \pi/2)$.

**The likelihood depends on the data through one function.** Integrating the
direction and the radius $r$ (the radius-nuisance split of the flat law) gives
$L(\theta) = K_d\big(h(\theta)\big)$ with

$$
h(\theta)^2 = \big|\mathbf x\sin\theta + \mathbf y\cos\theta\big|^2
= \frac{P+Q}{2} + \frac{Q-P}{2}\cos 2\theta + D\,\sin 2\theta .
$$

**The radial kernel.** $K_2(h) = e^{h^2/2}$ exactly;
$K_3(h) = \sqrt{\pi/2}\,e^{h^2/2}\operatorname{erf}(h/\sqrt2)/h$;
$K_1(h) = 1 + \sqrt{\pi/2}\,h\,e^{h^2/2}\operatorname{erf}(h/\sqrt2)$. For large $h$,
all three behave as $C\,h^{2-d}\,e^{h^2/2}$.

**So, in 2D exactly,**

$$
\log L(\theta) = \text{const} + \kappa\cos(2\theta - \varphi),
\qquad
\kappa\,e^{i\varphi} = c \equiv \frac{(Q-P) + 2iD}{4} = \frac14\sum_{k=1}^{d}(y_k + i\,x_k)^2 .
$$

Each reading is one complex number $c$: a vector in the doubled-angle plane.
It is a von Mises law in $2\theta$ with concentration $|c|$ and mean direction
$\arg c$.

**The prior is the same kind of object.** The zero-reading law (uniform in
$\theta$, the half-Cauchy at $s$) is the vector $0$. Contribution 12's tilt
centred at $s$ is $\exp(\lambda\operatorname{sech} z)$ with $z = \log(m/s)$, and
$\operatorname{sech}(\log\tan\theta) = \sin 2\theta$, so the tilt is the vector
$\lambda e^{i\pi/2}$: $\lambda$ is a von Mises concentration pointing at
$m = s$.

**Combination is addition.** With a common $s$, under the rule that counts
the prior once,

$$
p(\theta \mid \text{prior},\ 1..N) \;\propto\;
\exp\!\Big(\operatorname{Re}\big[(c_{\rm prior} + \textstyle\sum_{i} c_i)\,e^{-2i\theta}\big]\Big),
\qquad 2\theta \in (0, \pi),
$$

a von Mises law in $2\theta$ truncated to the half-circle that positive mass
allows. The whole old state is the vector
$C = c_{\rm prior} + \sum_i c_i$: two real numbers, updated by adding each
new reading's $c$. (Moving the centre is a rotation of $C$, not a shift in
$\log m$. The morning's tilt fit shifted in $\log m$, which is why it lost.)

**The peak is total least squares.** The mode $2\hat\theta = \arg C$ (when
it lies in $(0,\pi)$) maximises $\sum_i h_i^2$, the energy of the
standardised readings $(y, x)$ projected on the direction
$(\cos\theta, \sin\theta)$. That is the principal axis of their scatter
matrix: the orthogonal-regression (total-least-squares) fit through the origin
in standardised units, i.e. Deming regression with known error ratio.
So the law the estimator averages over is centred on TLS. The ratio of means
is a different summary of the same law. Anti-aligned data put $\arg C$ outside
$(0,\pi)$, and the mass then piles up at the boundary; this is how "no
positive mass fits" appears in these coordinates.

**What the vector does not carry.**

- **$A(m) = \sum_i E[\alpha_i \mid m]$.** Each term is
  $\sigma_a\cos\theta\,E[r \mid h_i(\theta)]$, nonlinear in each $h_i$, so
  the ratio-of-means readout still needs more than $C$ (as in the morning
  record).
- **In 1D and 3D, the factor $\approx h_i^{2-d}$ per reading**, exactly
  $\log\big(\sqrt{\pi/2}\operatorname{erf}(h/\sqrt2)/h\big)$ in 3D. It depends on each
  reading's own $h_i(\theta)$, so it does not fold into a fixed extra number.
- **Different instruments.** With different $s_i$ the readings' angles are
  different functions of $m$, and the vectors do not add.

## 3. Checks in the code

`tests/test_sequential.py` in the tester:

- in 2D the von Mises carry equals the whole old curve for the median, log SD
  and 95% interval, to $10^{-9}$, on random readings;
- in 1D and 3D it does not (the dropped factor is real);
- $\tan\hat\theta$ from $\arg\sum c$ equals the slope of the principal axis of the
  standardised scatter, to $10^{-10}$;
- the carry refuses settings it does not approximate (anything but the
  radius-nuisance likelihood product).

## 4. The study

Preset `presets/old_info_vonmises.yaml` (commit `9f27343`): readings share
only the mass (new excitation, random direction), unit noise so $s = 1$,
mass $\{s/4, s, 4s\}$, acceleration SNR $\{1, 3\}$, old readings
$\{1, 5, 20\}$, dimension $\{2, 3\}$; 36 cells, 500 series each. In each series the
last reading is new. Each carry's **median** is compared with the exact
answer's, $|\log(\text{carry}/\text{exact})|$, median / 90th percentile over
series.

<!-- generated from run.json, code 9f27343, 500 series per cell -->

**d = 2**

| mass / s | acc. SNR | old readings | log-normal fit | tilt fit (log-shifted) | von Mises state |
|---|---|---|---|---|---|
| 0.25 | 1 | 1 | 12% / 21% | 1.6% / 9.6% | 1e-14% / 4e-14% |
| 0.25 | 1 | 5 | 14% / 25% | 8.9% / 21% | 2e-14% / 6e-14% |
| 0.25 | 1 | 20 | 16% / 27% | 14% / 24% | 4e-14% / 1e-13% |
| 0.25 | 3 | 1 | 17% / 25% | 7.9% / 22% | 2e-14% / 4e-14% |
| 0.25 | 3 | 5 | 15% / 28% | 13% / 24% | 2e-14% / 7e-14% |
| 0.25 | 3 | 20 | 4.2% / 15% | 4.4% / 15% | 7e-14% / 2e-13% |
| 1 | 1 | 1 | 11% / 21% | 2% / 11% | 1e-14% / 4e-14% |
| 1 | 1 | 5 | 9.6% / 22% | 5.8% / 17% | 2e-14% / 4e-14% |
| 1 | 1 | 20 | 1.7% / 12% | 1.5% / 9.3% | 3e-14% / 9e-14% |
| 1 | 3 | 1 | 5.1% / 16% | 2.5% / 8.7% | 1e-14% / 3e-14% |
| 1 | 3 | 5 | 0.39% / 1.1% | 0.28% / 1% | 2e-14% / 4e-14% |
| 1 | 3 | 20 | 0.04% / 0.11% | 0.034% / 0.11% | 2e-14% / 4e-14% |
| 4 | 1 | 1 | 16% / 26% | 11% / 25% | 2e-14% / 5e-14% |
| 4 | 1 | 5 | 11% / 26% | 10% / 22% | 3e-14% / 9e-14% |
| 4 | 1 | 20 | 2.4% / 6% | 2.4% / 6.3% | 9e-14% / 2e-13% |
| 4 | 3 | 1 | 4.8% / 13% | 4.3% / 11% | 4e-14% / 1e-13% |
| 4 | 3 | 5 | 0.96% / 1.7% | 0.97% / 1.7% | 8e-14% / 2e-13% |
| 4 | 3 | 20 | 0.25% / 0.32% | 0.25% / 0.32% | 1e-13% / 3e-13% |

**d = 3**

| mass / s | acc. SNR | old readings | log-normal fit | tilt fit (log-shifted) | von Mises state |
|---|---|---|---|---|---|
| 0.25 | 1 | 1 | 11% / 21% | 1.5% / 9.6% | 6.6% / 20% |
| 0.25 | 1 | 5 | 14% / 25% | 7.1% / 21% | 7.7% / 24% |
| 0.25 | 1 | 20 | 15% / 26% | 12% / 23% | 6.3% / 20% |
| 0.25 | 3 | 1 | 17% / 24% | 6.6% / 21% | 2.1% / 7.1% |
| 0.25 | 3 | 5 | 15% / 28% | 13% / 25% | 1.7% / 5.4% |
| 0.25 | 3 | 20 | 5.2% / 17% | 5.6% / 16% | 0.94% / 2.6% |
| 1 | 1 | 1 | 12% / 21% | 1.6% / 10% | 4.6% / 16% |
| 1 | 1 | 5 | 12% / 22% | 5.3% / 17% | 3.8% / 14% |
| 1 | 1 | 20 | 2.2% / 16% | 2% / 12% | 2.3% / 7.3% |
| 1 | 3 | 1 | 5.1% / 16% | 2.8% / 9.8% | 0.55% / 2% |
| 1 | 3 | 5 | 0.4% / 1.3% | 0.29% / 1% | 0.36% / 1.1% |
| 1 | 3 | 20 | 0.042% / 0.11% | 0.039% / 0.11% | 0.19% / 0.42% |
| 4 | 1 | 1 | 16% / 26% | 9.9% / 24% | 1.1% / 4.1% |
| 4 | 1 | 5 | 11% / 26% | 9.7% / 23% | 0.8% / 2.6% |
| 4 | 1 | 20 | 2.4% / 6.1% | 2.4% / 5.7% | 0.41% / 1.1% |
| 4 | 3 | 1 | 4.4% / 13% | 4.2% / 11% | 0.06% / 0.15% |
| 4 | 3 | 5 | 0.92% / 1.6% | 0.93% / 1.6% | 0.023% / 0.061% |
| 4 | 3 | 20 | 0.25% / 0.33% | 0.25% / 0.33% | 0.0093% / 0.021% |

**What it shows.**

- **2D: exact.** The von Mises state matches the exact median and 95%
  interval to about $10^{-13}$ in all 18 cells, whatever the signal
  strength. The fits are 0.03–17% off.
- **3D: closest in 12 of 18 cells.** Mass $4s$: 0.009–1.1% against
  0.25–16% for the fits. It loses in six cells. Clearly in two: one weak
  old reading at mass $s/4$ and $s$, SNR 1, where the tilt fit is 1.5–1.6%
  off against 4.6–6.6%. Narrowly in four: by 0.3–0.6 points in two weak
  cells (mass $s/4$, 5 readings; mass $s$, 20 readings, both SNR 1), and at
  mass $s$, SNR 3 with 5 and 20 readings, where everything is already small
  (0.19% against 0.04% at 20 readings). This fits the dropped $1/h$
  factor, which matters most for a single weak reading and adds up over
  many readings.
- **Point accuracy against the true mass is unchanged.** This is about
  storing old information, not about estimating better.

## 5. Where to read further

Papers marked **[papers/]** are in `../papers/` but were not re-read against
these claims for this record. The others are cited from general knowledge
and have not been fetched or checked. Check each before citing it in public.

**The peak (TLS).**
- Golub & Van Loan (1980), *An analysis of the total least squares problem*
  **[papers/]**. The principal-axis solution of TLS; the identity in §2.
- Deming (1943), *Statistical Adjustment of Data*, and Adcock (1878), *A
  problem in least squares*, *The Analyst* 5. Orthogonal / known-error-ratio
  regression; the classical form of the peak.

**The law around the peak (Bayesian errors-in-variables).**
- Gull (1989), *Bayesian data analysis: straight-line fitting* **[papers/]**,
  and Lindley & El-Sayyad (1968) **[papers/]**. Straight lines with noise in
  both coordinates and the true values integrated out. Compare their slope
  posteriors with the von Mises law here, and check whether either already
  writes it in the doubled angle.
- Liseo (2003) **[papers/]**. The 1D flat law (his eq. 12). In 1D the
  dropped factor is $+\log h$.

**Doubled angles and the von Mises family (directional statistics).**
- Mardia & Jupp (2000), *Directional Statistics* (Wiley). The von Mises
  distribution, its conjugacy, and the standard treatment of *axial* data
  (directions without sign) by doubling angles. Our $2\theta$ is that move.
- Fisher (1993), *Statistical Analysis of Circular Data* (Cambridge).
  Practical estimation of $\kappa$ and $\varphi$, and truncated circular
  laws.
- Bingham (1974), *An antipodally symmetric distribution on the sphere*,
  *Ann. Statist.* 2, 1201–1225. Densities $\propto\exp(\mathbf v^\top A\mathbf v)$ on
  unit vectors. Our $e^{h^2/2}$ is exactly this form in
  $\mathbf v = (\sin\theta, \cos\theta)$; on the circle a Bingham law is a von
  Mises law in the doubled angle. This is probably the precise name for what
  was found.

**Why a vector state exists.**
- Diaconis & Ylvisaker (1979), *Conjugate priors for exponential families*,
  *Ann. Statist.* 7, 269–281. In 2D each likelihood is an exponential family
  in $(\cos 2\theta, \sin 2\theta)$, so its conjugate prior is von Mises and
  the natural parameters add. This is the general reason behind §2.

**The radial kernel.**
- Rice (1944–45), *Mathematical analysis of random noise*, *Bell System
  Technical Journal* 23–24. The Rice distribution; the 2D mean radius used for
  $E[\alpha \mid m]$ is the Rice mean.
- `../literature_review/` (the unsanctioned Codex search of 6–8 Sept; see
  that folder's README for its standing). It located the target as a
  Gaussian errors-in-variables likelihood with an improper prior
  $\|\mathbf v\|^{2-d}$. That power is the $h^{2-d}$ factor here: a pointer
  to follow, not a source.

**The 3D approximation.**
- Opper (1998), *A Bayesian approach to on-line learning*, and Minka (2001),
  *Expectation propagation for approximate Bayesian inference*. Assumed
  density filtering: projecting each update back into a family. In 3D the
  natural next test is to project into the von Mises family by matching
  trigonometric moments, rather than dropping the $h^{-1}$ factor.

## 6. Open

- **3D:** does ADF projection into the von Mises family beat the dropped-factor
  state? Is there a finite exact state in 3D (probably not, given §2)?
- **Mixed instruments:** readings with different $s_i$. The per-reading
  vectors are then in different angle coordinates.
- **The ratio of means:** can $A(m)$ be carried compactly, at least in 2D,
  where $E[r \mid h]$ is the Rice mean?
- **Literature:** is the doubled-angle / Bingham form of the Bayesian
  errors-in-variables posterior already written down? If so, this is a
  rediscovery with a new use; if not, it is new.

## 7. Reproduction

Repository `canvalidk/Mass-estimator-testor`, commit `9f27343`.

```sh
python -m pytest -q tests/test_sequential.py      # the identities of §3
python -m met run presets/old_info_vonmises.yaml  # the study of §4 (about 25 minutes)
```

Code: `met/analyst.py`, `_sequential_log_density` (carry `vonmises_state`).

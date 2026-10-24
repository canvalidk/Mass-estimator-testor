# The size measure: the Jeffreys check, and what calibrates it

**2026-09-28 (night).** Record, filed at the user's request ("Go ahead"),
written by Claude. It pays the first owed item of
`ours28_candidate_push_scale_as_nuisance_2026-09-28.md` §7: run the Jeffreys
route for the push scale and see whether it agrees with the size-angle rule.
It edits nothing and binds nothing. ours28 stays a candidate.

**In one paragraph.** Under ours28's own model the Fisher information about
the mass is $\omega^4/(1+\omega^2)$, independent of $\theta$, and the mass and
the push scale are orthogonal. So the reference prior on the mass is uniform in
$\theta$ (the half-Cauchy) for *any* treatment of the scale, and the fixed push
of the starting-law record §6 is no longer needed to reach it. For the scale
itself the reference prior is $\tan\varphi\,d\varphi$, not uniform, and the
joint Jeffreys prior is $\sin^3\!\varphi/\cos^2\!\varphi$. On coverage both do
worse than uniform $\varphi$, and the reference prior under-covers at weak
signal (81% at SNR ½ against 92%). What separates the rules is not their
shape but their symmetry: every prior symmetric under $\omega\leftrightarrow1/\omega$
(push and noise exchanging roles) calibrates, 94–98% in every cell for
$\mathrm{Beta}(a,a)$ with $a$ from ¼ to 2, and the three that fail (flat,
reference, joint Jeffreys) are the three tilted toward "the reading is all
push". So the size measure has the same structure as the starting law: one
symmetry premise carries the calibration, the shape is open, and the shape
matters only at weak signal. Jeffreys is not the operative principle for the
size, and was therefore not the reason the mass's rule works either.

---

## 1. The Fisher matrix of ours28's model

Per axis in noise units, a reading is $\mathbf z=v\mathbf u+\mathbf e$ with
$\mathbf u=(\sin\theta,\cos\theta)$, $v\sim N(0,\omega^2)$, $\mathbf e\sim N(0,I)$.
The marginal is $\mathbf z\sim N(0,\;I+\omega^2\mathbf u\mathbf u^\top)$ with
parameters $(\theta,\omega)$. For a zero-mean Gaussian
$I_{ab}=\tfrac12\operatorname{tr}(\Sigma^{-1}\partial_a\Sigma\,\Sigma^{-1}\partial_b\Sigma)$,
which gives (sympy, and by hand)

$$
I_{\theta\theta}=\frac{\omega^4}{1+\omega^2},\qquad
I_{\omega\omega}=\frac{2\omega^2}{(1+\omega^2)^2},\qquad
I_{\theta\omega}=0 .
$$

Two things follow.

- **The mass.** $I_{\theta\theta}$ does not depend on $\theta$, and $\theta$
  is orthogonal to $\omega$. So the Berger–Bernardo reference prior with
  $\theta$ the parameter of interest is uniform in $\theta$, whatever is done
  with $\omega$. The half-Cauchy on the mass is now derived from ours28's own
  model. The starting-law record §6 reached it by holding the push fixed; that
  assumption is no longer needed. (At large $\omega$, $I_{\theta\theta}\to\omega^2=E[v^2]$,
  the fixed-push value.)
- **The size.** In the size angle $\omega=\tan\varphi$:
  - reference prior ($\theta$ of interest, $\omega$ nuisance):
    $\sqrt{I_{\omega\omega}}\,d\omega\propto\tan\varphi\,d\varphi$;
  - joint Jeffreys: $\sqrt{I_{\theta\theta}I_{\omega\omega}}\,d\omega\propto\sin^3\!\varphi/\cos^2\!\varphi\,d\varphi$.

  Neither is uniform in $\varphi$. Both are improper at $\varphi\to\pi/2$
  (infinite push): the reference prior logarithmically, the joint one like
  $1/\cos\varphi$. Both give a proper per-reading evidence for $d\ge2$; the
  joint prior's evidence diverges in $d=1$.

## 2. The family, in closed form

With $t=\beta=\sin^2\varphi$ the per-reading evidence is
$(1-t)^{d/2}e^{t h^2/2}$ (ours28 record §3). A $\mathrm{Beta}(a,c)$ prior on
$\beta$ gives

$$
g_d(h^2)\;\propto\;\int_0^1 t^{a-1}(1-t)^{c-1+d/2}e^{t h^2/2}\,dt
\;\propto\;{}_1F_1\!\left(a;\;a+c+\tfrac d2;\;\tfrac{h^2}{2}\right),
$$

one Kummer function per rule. The rules on the table:

| rule | prior on $\beta$ | evidence | tempering at large $h^2$ |
|---|---|---|---|
| flat push (ours27) | point mass at $\beta=1$ | $e^{h^2/2}$ | none |
| uniform $\varphi$ (ours28 as filed) | $\mathrm{Beta}(\tfrac12,\tfrac12)$ | ${}_1F_1(\tfrac12;\tfrac d2+1;\tfrac{h^2}{2})$ | $h^{-(d+1)}$ |
| reference prior | $\mathrm{Beta}(1,0)$ | ${}_1F_1(1;\tfrac d2+1;\tfrac{h^2}{2})$ | $h^{-d}$ |
| joint Jeffreys | $\mathrm{Beta}(2,-\tfrac12)$ | ${}_1F_1(2;\tfrac d2+\tfrac32;\tfrac{h^2}{2})$ | $h^{-(d-1)}$ |

(${}_1F_1(a;b;z)\sim\Gamma(b)/\Gamma(a)\,z^{a-b}e^z$ at large $z$, so the
power $z^{a-b}$ is the tempering.) The quadrature used for the tests was
checked against the closed forms at four values of $h^2$: ratios 1 to $10^{-5}$.

Uniform $\varphi$ is the half-Cauchy on $\omega$: $\omega=\tan\varphi$ with
$\varphi$ uniform has density $(2/\pi)/(1+\omega^2)$. So ours28 as filed puts
on the push scale, in noise units, the same starting law it puts on the mass.

## 3. Coverage

Same protocol as the ours28 record §5: one noise pair, $d=3$ unless stated,
Gaussian pushes with per-component SD equal to the SNR, $\pi$ the half-Cauchy
at $c=s$, 600 series per cell ($\pm1$ point), 95% central interval of the law
of mass.

### 3.1 The three rules against flat

| cell | flat | uniform $\varphi$ | reference ($\tan\varphi$) | joint Jeffreys |
|---|---|---|---|---|
| $m=s$, SNR 1, $N=21$ | 81.5 | **93.3** | 88.5 | 85.2 |
| $m=s/4$, SNR 1, $N=21$ | 86.7 | **96.8** | 94.5 | 91.3 |
| $m=s$, SNR 1, $N=6$ | 85.2 | 97.5 | **94.3** | 89.7 |
| $m=s/4$, SNR 3, $N=21$ | 93.3 | 95.5 | **95.0** | 93.7 |
| $m=s$, SNR ½, $N=21$ | 65.7 | **92.3** | 80.7 | 73.0 |

Median 95% log-width, first cell: 1.09 / 1.64 / 1.37 / 1.21. The ordering is
the tempering column of §2: less tempering, less coverage.

### 3.2 Other worlds (uniform $\varphi$ vs reference; flat for scale)

| world | cell | flat | uniform $\varphi$ | reference |
|---|---|---|---|---|
| fixed-length pushes, $d=3$ | $m=s$, SNR 1, $N=21$ | 84.3 | 94.2 | 90.7 |
| fixed-length pushes, $d=3$ | $m=s$, SNR ½, $N=21$ | 65.5 | 91.2 | 81.8 |
| Gaussian pushes, $d=1$ | $m=s$, SNR 1, $N=21$ | 84.2 | 97.0 | 90.2 |
| Gaussian pushes, $d=1$ | $m=s/4$, SNR 1, $N=21$ | 91.2 / 1.26 | 97.8 / 1.52 | 94.5 / 1.34 |
| fixed-length pushes, $d=1$ | $m=s$, SNR 1, $N=21$ | 84.0 | 95.0 | 89.7 |
| mixed push sizes (half at SNR 1.5, half at 0.375), $d=3$ | $m=s$, $N=21$ | 87.3 | 95.2 | 92.5 |
| mixed push sizes, $d=1$ | $m=s$, $N=21$ | 86.3 | 96.2 | 90.2 |

The mixed world is the one where the pushes have no common scale, the case
the per-reading form is for. The 1D weak-mass row shows the price of
tempering where it is largest: median over truth 1.26 → 1.52 (uniform
$\varphi$) or 1.34 (reference).

### 3.3 Inside the symmetric family $\mathrm{Beta}(a,a)$

Priors on $\beta$ symmetric under $\beta\leftrightarrow1-\beta$, i.e.
$\omega\leftrightarrow1/\omega$. (Different seed from §3.1; the uniform-$\varphi$
column moves by about a point.)

| cell | $a=\tfrac14$ | $a=\tfrac12$ (uniform $\varphi$) | $a=1$ (uniform $\beta$) | $a=2$ |
|---|---|---|---|---|
| $d=3$, $m=s$, SNR 1, $N=21$ | 94.3 | 94.2 | 94.3 | 94.7 |
| $d=3$, $m=s$, SNR ½, $N=21$ | 96.7 | 94.7 | 93.2 | 91.7 |
| $d=3$, $m=s/4$, SNR 1, $N=21$ | 97.7 | 97.7 | 97.2 | 97.3 |
| $d=1$, $m=s$, SNR 1, $N=21$ | 97.7 | 97.7 | 97.5 | 97.7 |
| $d=3$, $m=s$, SNR 3, $N=21$ | 94.7 | 94.7 | 95.0 | 96.0 |

Median over truth is the same across $a$ to two decimals in every cell but
the weak-mass one (1.17 / 1.14 / 1.12 / 1.12).

**Reading.** Every symmetric prior calibrates. The three that fail in §3.1
are exactly the asymmetric ones, and all three are tilted toward $\beta=1$,
"the reading is all push". So the calibration comes from the symmetry, and
the shape is nearly free on this yardstick; it shows only at the weakest
signal, where more tempering (smaller $a$) covers more.

## 4. What this says

- **Jeffreys is not the principle.** For the size, the reference prior and
  the joint Jeffreys prior both under-cover, and the rule that works is not
  either of them. For the mass, Jeffreys and uniform-angle coincide (twice
  now: fixed push, and §1), so the agreement there was not evidence that
  Jeffreys is what step 8 rests on. The starting-law record §6's closing
  line ("why Jeffreys at all is itself a premise") stands, and the premise
  is not needed.
- **The size measure has the same shape of story as the starting law.**
  - Mass: L4 (force and acceleration on an equal footing; $m\leftrightarrow s^2/m$)
    is the load-bearing premise; the shape is any even $g(\ln(m/c))$; it
    matters only at weak signal (lemma record §3.1).
  - Size: **push and noise on an equal footing before data**
    ($\omega\leftrightarrow1/\omega$, $\beta\leftrightarrow1-\beta$) is the
    load-bearing premise; the shape is any $\mathrm{Beta}(a,a)$, or any even
    density in $\ln\omega$; it matters only at the weakest signal.

  L10 was the maximally asymmetric member: attribute every reading entirely
  to the push. T18 is the cost of that.
- **One rule, two instances, as a candidate.** Uniform $\varphi$ is the
  half-Cauchy on $\omega$ in noise units, the same law as on $m/s$. Taking
  "uniform in the angle the noise defines" for both is one shape choice with
  two uses rather than two choices. The yardstick does not single it out
  within the symmetric family; the parallel with the mass does.
- **The mass's rule is strengthened.** Uniform $\theta$ now follows from
  ours28's own Fisher information for any size prior (§1). It does not
  depend on holding the push fixed.
- **The remaining freedom.** Two shapes, both even, both nearly immaterial
  on coverage: $g$ for the mass, $a$ for the size. And C3. Nothing else.

**Proposed for the lemma record (not made).** Beside L4, a premise **L4′,
push–noise symmetry**: the starting measure on a reading's push scale is
symmetric under $\omega\leftrightarrow1/\omega$. Its status is the same as
L4's: strong, motivated by equal footing, and what carries the widths. T18
becomes "the width error of the asymmetric endpoint". §3.1's "open: its
shape" acquires a twin for the size.

## 5. Caveats

- The yardstick is coverage against a planted mass, in simulated push worlds.
  The folder already notes that this yardstick is not settled under "mass is
  its estimator". On the properties of mass every rule here ties (they are
  all prior-shape bends, which the perturbation matrix found free). So the
  decision between them needs a premise, and §4 names the candidate.
- One noise pair only. With several pairs the units fork reappears for
  $\omega_i$ (ours28 record §7).
- 600 series per cell; differences under about two points are noise.
- Not run: the readout point under any of these; widths of the point; the
  per-axis coverage cells of the sweep.

## 6. Reproduction

Scripts in Claude's workspace on 28 September, not kept. sympy for §1 (the
$2\times2$ covariance, its inverse, the three traces; the $\theta$-dependence
of $I_{\theta\theta}$ cancels by $\sin^6-3\sin^4+3\sin^2+\cos^6=1$). numpy and
scipy for §3: the law on a 400-point $\theta$ grid, the size prior on a
200-point $\varphi$ grid as a weight $w(\varphi)$ multiplying
$(\cos\varphi)^d\exp(\tfrac12\sin^2\!\varphi\,h_i^2)$; `scipy.special.hyp1f1`
for the closed-form check; `scipy.stats.beta` for §3.3.

## 7. How we got here

1. The ours28 record listed the Jeffreys check as owed and guessed
   $\tan\varphi$ from the one-parameter information. The guess was right;
   the two-parameter calculation added orthogonality and the mass result of
   §1.
2. The three rules were run in the five cells (§3.1), then in other worlds
   to rule out an artefact of Gaussian pushes or $d=3$ (§3.2).
3. Noticing that uniform $\varphi$ is $\omega\leftrightarrow1/\omega$
   symmetric and the reference prior is not, the symmetric family was run
   (§3.3), which is where the finding is.


# The starting law under ours26: what is forced, a no-go, and a candidate

**2026-09-27.** Record, filed at the user's request ("Go ahead", on Claude's
proposal to record the unrecorded 27 September results by thread), written by
Claude.

**Where the results come from.** §§2–5 come from a completion message VD-Newton
wrote on 27 September (not filed; its own checks were a scratch probe in
VD-Newton's `.tools/ours26_completion_20260927/`). Claude checked its algebra and
reproduced its numbers. §6 (the Jeffreys route and what it selects) is Claude's,
27 September.

**Relations.**
- **Corrects** `ours26_per_axis_noise_2026-09-26.md` §4, and only there: the
  tester's `first_reading` / `first_reading_geometric` centres depend on arrival
  order once readings come from different instruments (§4 below). Everything
  else in that record stands.
- Scopes the lemma record (`../mass_estimator_lemmas_2026-09-27.md`) §3.1 and
  supplies its T13 and T20. Binds nothing else.

---

## 1. Two separate things

- **The operation** that turns readings into a law and a readout:
  $p(m)\propto\pi(m)\,L(m)$ with $L$ the ours26 likelihood, and
  $\widehat m=\int mA\,p\,dm\big/\!\int A\,p\,dm$,
  $A(m)=E[\lVert\mathbf a^\ast\rVert\,\vert\,m]$.
- **The starting law** $\pi$ it is applied to.

The chain forces most of the operation. This record is about how much it forces
of $\pi$.

## 2. Any proper starting law gives a mass

$A(m)$ falls like $1/m$ as $m\to\infty$ and stays bounded as $m\to0$. So
$A$ and $mA$ are both bounded, and every proper $\pi$ gives a finite, strictly
positive $\widehat m$, even a $\pi$ with no finite mean.

- Check (Claude, a random 3D case with a full noise covariance): $mA(m)$ =
  6.1256 at $m=10^2$ and 6.1174 at $m=10^4$; $A(10^{-3})=4.04$.
- **Consequence:** totality (the output is always a mass) does not select
  $\pi$. It must not be cited as a reason for the half-Cauchy.

## 3. What the symmetries leave

With one instrument ratio $c$, every

$$
\pi_c(m)=\frac1m\,g\!\left(\ln\frac mc\right),\qquad g\ \text{an even probability density},
$$

is unit-covariant (L3), symmetric under the channel swap (L4), and has median
$c$ at zero readings. The half-Cauchy is the case $g(u)\propto1/\cosh u$.
They give different answers from the same informative reading. For
$\tilde F=1$, $\tilde a=0$, $\sigma_F=\sigma_a=c=1$ (1D), with the same law and
readout:

| starting law | readout |
|---|---|
| half-Cauchy | 1.33213 |
| lognormal, log-SD 0.3 | 1.03539 |
| lognormal, log-SD 1 | 1.23185 |
| lognormal, log-SD 2 | 1.48048 |

(VD-Newton's numbers; reproduced by Claude to every digit by quadrature.)

So symmetry, totality and exact combination together do not finish the
selection of $\pi$.

## 4. The mixed-instrument no-go

**A zero reading carries no evidence.** Its misfit is zero, so under ours26 it
contributes likelihood exactly 1: it leaves the law unchanged.

**Claim.** With two instruments whose own defaults differ ($c_1\neq c_2$), these
three cannot all hold:

1. the starting law is taken from the instruments that have reported so far
   (so a series from instrument $k$ alone starts from $\pi_{c_k}$);
2. a reading with no evidence leaves the law unchanged;
3. the order in which readings arrive does not matter.

**Proof.** Take one zero reading from each instrument. Instrument 1 first: the
law starts at $\pi_{c_1}$ (1), and the zero reading from instrument 2 leaves it
there (2). Instrument 2 first: $\pi_{c_2}$. By (3) these are equal, so
$c_1=c_2$. $\square$

A rule that picks $\pi$ from the whole set of instruments used, at the end, is
order-free but breaks (2): a zero reading from a new instrument would change the
law without carrying any evidence, and the old law could no longer serve as the
starting law for new readings.

**What gives.** (1). The starting law is **part of the specification**, fixed
before any data. Changing it is changing the specification, not updating.

**The correction to ours26 §4.** `first_reading` and `first_reading_geometric`
take the centre from the first reading. With readings, or combined series, from
different instruments, the answer then depends on which came first. The
tester's current worlds use one noise pair per axis, the same in every reading,
so every reading has the same ratios and nothing there depends on order. The
defect bites mixed instruments and sequential combination across instruments,
which ours26 allows and no tester world yet generates. The fix is to declare the
centre, or the instrument geometry, in the specification.

**The consequence for L3.** "Nothing else goes in" holds only with $\pi$ counted
as a declared input. The one mass scale allowed in is the declared $\pi$, which
transforms with the units.

## 5. A candidate: uniform arc length

ours26 gives each noise pair $k$ an angle $\theta_k(m)=\arctan(m/s_k)$. As $m$
runs from 0 to $\infty$ these trace one curve. Equal weight per unit of its
Euclidean arc length:

$$
\pi_{\rm arc}(m)\;\propto\;\Big(\sum_k\Big(\frac{s_k}{s_k^2+m^2}\Big)^2\Big)^{1/2}.
$$

- **Proper:** its total is at most $d\,\pi/2$ (3.0165 for ratios 5, 2, 0.8
  against the bound 4.712; Claude, quadrature).
- **Unit-covariant and swap-symmetric** (checked to rounding), and **the
  half-Cauchy when all ratios agree.**
- **No separately chosen centre.** The instruments fix it.
- **Not the only such rule.** Summed angular speeds $\sum_ks_k/(s_k^2+m^2)$,
  an equal mixture of half-Cauchys, passes the same checks and gives different
  answers; so does any other norm of the angular speeds.
- **Fixed geometry.** By §4 the instrument set behind $\pi_{\rm arc}$ must be
  declared up front. Recomputing it as instruments are added changes $\pi$.

## 6. The Jeffreys route, and what it selects (Claude)

**Setting.** Per axis, in noise units, the true pair is a push $v$ at angle
$\theta$: means $(v\sin\theta,\,v\cos\theta)$ with unit noise. The information
matrix in $(\theta,v)$ is diagonal: $v$ is **orthogonal** to $\theta$. The
efficient information about $m$ from axis $k$ is

$$
I_k(m)=r_k^2\Big(\frac{d\theta_k}{dm}\Big)^2=\frac{r_k^2\,s_k^2}{(s_k^2+m^2)^2},
$$

with $r_k$ the push in noise units. Informations add over axes, and the Jeffreys
prior is $\sqrt{\sum_kI_k(m)}$. Which prior this gives depends on what is held
equal as $m$ varies. Two readings:

| held equal on every axis | Jeffreys prior | one pair | L4 (swap) | proper |
|---|---|---|---|---|
| the push in noise units, $r_k$ | $\pi_{\rm arc}$ | $s/(s^2+m^2)$, the current half-Cauchy | yes | yes |
| the physical acceleration, $\alpha$ | $\big(\sum_k1/(\sigma_{a,k}^2(s_k^2+m^2))\big)^{1/2}$ | $\propto1/\sqrt{s^2+m^2}$ | **no** | **no** |

(Symbolic checks, Claude: the diagonal information matrix; the efficient
information with $\alpha$ as nuisance, $\alpha^2/(\sigma_a^2(m^2+s^2))$; swap
invariance and the integral of each one-pair prior.)

**What this settles.**

- **One noise pair:** L4 and properness both reject the physical reading. The
  current half-Cauchy is the Jeffreys prior with the orthogonal push held fixed.
  It now has two independent reasons: uniform in the angle the noise defines,
  and this.
- **Several noise pairs:** "equal $r_k$ on every axis" is isotropy in noise
  units. It gives $\pi_{\rm arc}$, the reason to take the Euclidean norm rather
  than summed speeds: informations add, so their square roots combine in
  quadrature. Whether pushes should be isotropic in noise units is the units
  fork (`README.md` in this folder), not settled here.

**What it does not settle.** Why Jeffreys at all is itself a premise, not a step
of the chain. Holding the orthogonal nuisance fixed is the standard
reference-prior move, but it is a choice.

## 7. What this leaves

- $\pi$ is declared in the specification (§4).
- One noise pair: the half-Cauchy, with the two reasons of §6. Its shape is
  still not forced by the lemma list alone (§3); the Jeffreys route is an extra
  premise.
- Several noise pairs: $\pi_{\rm arc}$ is the leading candidate, conditional on
  noise-unit isotropy and on the Jeffreys principle. It replaces the centre
  family of ours26 §4 if adopted.
- Owed: a tester option for $\pi_{\rm arc}$ and a declared-centre default; a
  world with mixed instruments, to show the order dependence and its fix.

## 8. Reproduction

Claude's checks were run in its sandbox on 27 September and not kept: the prior
table (1D quadrature of law and readout under each prior); $A$ and $mA$ at
extreme $m$ (Laplace-transform integral for $E\lVert N(\mathbf b,K)\rVert$);
$\pi_{\rm arc}$'s total, unit covariance, swap symmetry and one-pair reduction;
the information matrices and one-pair priors (sympy).

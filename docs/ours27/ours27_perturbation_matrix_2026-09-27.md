# Bending ours27: which moves are forced, and by what

**2026-09-27.** Record, filed at the user's request ("Make a folder called ours27
put it there. Make a record of it."), written by Claude from the work of the
night of 27 September. The user had asked for ours27 to be checked "against all
the algebraic moves. See if slightly changing the algebra changes what we can do
with it". It edits nothing and binds nothing.

**The equation.** ours27 (`ours27_full_covariance_2026-09-27.md`, in this
folder): law $\pi(m)\prod_ie^{-Q_i(m)/2}$, readout
$E\lVert\mathbf F^\ast\rVert/E\lVert\mathbf a^\ast\rVert$ over the stacked true
vectors, flat push in noise units, Gaussian noise, $\pi$ counted once.

**The method.** Take each move in the chain that leads to ours27. Bend it by a
parameter $\varepsilon$, so that $\varepsilon=0$ is ours27. Then test every
property of mass on the bent estimator. Three outcomes are possible:

- **a property breaks at first order in $\varepsilon$:** the move is forced, and
  that property is what forces it;
- **nothing breaks:** the bend is a live alternative, a rival definition of
  mass;
- **it breaks only at higher order:** forced, but weakly (none found).

**In one paragraph.** Nine bends were tested against eight properties. Six are
forced, each by at least one property. Three are free: the starting law's shape,
a direction weight in the readout (the counterexample C3 exists to exclude), and
a finite push scale $\omega$ in noise units. So ours27 is unique given C3 and
"no push scale", up to the starting law's shape. Two results matter beyond the
list:

- **Gaussian noise is forced, not declared.** With the flat push, a
  heavy-tailed noise model is biased even in its own world, and the bias does
  not shrink with more readings. The Gaussian model converges even when the
  real noise is heavy-tailed. Its criterion uses only the second moments the
  instrument declares.
- **Axis splitting carries the most weight.** It is the only property forcing
  the flat push size, and the only one forcing the stacked length.

---

## 1. The bends

Each bend is one parameter; $\varepsilon=0$ is ours27.

| bend | what changes | $\varepsilon$ |
|---|---|---|
| push size | push weight $\lVert\mathbf v_i\rVert^{\varepsilon}\,d^d\mathbf v_i$ per reading (ours21 is $\varepsilon=2-d$) | ±0.1 |
| push units | an extra $\prod_c(\cos\theta_c)^{\varepsilon}$ per reading, part-way from noise units to physical units | ±0.1 |
| prior count | $\pi$ counted again per reading, $(m\,\pi(m))^{\varepsilon}$, as a density in $\log m$ (the swap-symmetric version; eq. 20 is $\varepsilon=1$) | ±0.1 |
| noise law | Student-t per component with $\nu=1/\varepsilon$, scaled to unit variance | 0.1 |
| push scale | Gaussian push weight $e^{-\varepsilon\lVert\mathbf v\rVert^2/2}$, i.e. $\omega=\varepsilon^{-1/2}$ in noise units | 0.1 |
| prior shape | $\pi\to\pi\,e^{\varepsilon\,\mathrm{sech}\ln(m/c)}$, contribution 12's λ tilt | ±0.1 |
| joining rule | size of a series $(\sum_i\lVert\mathbf a^\ast_i\rVert^p)^{1/p}$, $p=2+\varepsilon$ | ±0.1 |
| readout, acceleration weight | $E[f\alpha^{\varepsilon}]/E[\alpha^{1+\varepsilon}]$ | ±0.1 |
| readout, direction weight | $E[hf]/E[h\alpha]$, $h=1+\varepsilon(\mathbf u\cdot\mathbf e)^2$, $\mathbf e$ the law's mean direction (the 6 Sept counterexample) | ±0.1 |

Every bend is built rotation-covariant, so rotation does not separate them. The
rotation-breaking variant, a diagonal-only noise model, is ours26, and is the
reason ours27 exists (`ours27_full_covariance_2026-09-27.md` §2).

## 2. The properties and how each was tested

- **Lumping; common force.** At the level of the readout: a law of 200,000
  candidates (two 3D readings each), with the mass or the acceleration split
  between two bodies by a candidate-dependent fraction. Residual of
  $R(X_1{+}X_2,Y)=R(X_1,Y)+R(X_2,Y)$ and of its inverse-mass mirror. Law bends
  do not change the readout's form, so they pass both by construction.
- **Swap.** Exchange the channels (readings, SDs, $c\to1/c$): readout product
  $-1$ and median product $-1$.
- **Units.** Force readings, force SDs and $c$ scaled by 3: readout ratio $-1$.
- **Axis splitting.** A series read as its readings, against the same series
  read component by component: readout and median ratio $-1$.
- **A zero reading leaves the law unchanged.** A 3D reading of zeros appended:
  median and 90th-percentile ratio $-1$.
- **Consistency.** Law bends only. True $m=5$, $s=c=2$, pushes of length 1.5 in
  noise units, $N$ from 100 to 51,200: does the median's shift from ours27 die
  away with $N$? Readout bends are consistent by algebra: within a candidate
  $\mathbf F^\ast=m\,\mathbf a^\ast$, so every readout of this family tends to
  $m$ as the law concentrates.
- **The noise law's limit.** Where the estimate converges, found exactly as
  the peak of the expected per-reading evidence, by quadrature over the noise
  (no sampling), in 1D.

Test data for the first five columns:

- **One pair:** two 3D readings, $\sigma_F=1$, $\sigma_a=0.5$, $c=2$; ours27
  gives 3.0474.
- **Per-axis:** the same readings with $\sigma_a=(0.2,0.5,1.25)$ per axis;
  ours27 gives 2.9146.

Candidates are drawn as ours27's conditional Gaussians, 12,000 per mass node on
481 nodes in $\log m$. Bends enter as reweighting, and the same random numbers
are used across $\varepsilon$, swaps and splits, so exact properties read as
zero to machine precision.

**Controls.** Unbent ours27 reads zero, to $10^{-16}$, in every column, in both
test problems. The prior-shape bend, known to be free, reads zero everywhere.

## 3. The matrix

✓ means the property holds exactly (zero to machine precision, or to the
simulation's noise floor where the algebra says it is exact). ✗ means it
breaks, with the size of the break at $\varepsilon=0.1$, one pair first and
per-axis second.

| bend | lumping | common force | swap | units | axis splitting: law / readout | zero reading | consistency | verdict |
|---|---|---|---|---|---|---|---|---|
| push size | ✓ | ✓ | ✓ | ✓ | ✗ $-1.7\!\times\!10^{-4}$, $+3.9\!\times\!10^{-4}$ / ✗ $-1.6\!\times\!10^{-4}$, $+2.3\!\times\!10^{-4}$ | ✓ | ✓ | forced by axis splitting alone |
| joining rule | ✓ | ✓ | ✓ | ✓ | ✓ / ✗ $+2.9\!\times\!10^{-5}$, $+1.6\!\times\!10^{-4}$ | ✓ | ✓ | forced by axis splitting alone |
| push units | ✓ | ✓ | ✗ $-6.0\%$, $-3.7\%$ | ✓ | ✓ / ✓ | ✗ $-1.9\%$, $-1.1\%$ | ✗ | forced three ways |
| prior count | ✓ | ✓ | ✓ | ✓ | ✗ $-1.6\%$, $-0.8\%$ / ✗ | ✗ $-0.4\%$, $-0.2\%$ | ✗ | forced three ways |
| noise law | ✓ | ✓ | ✓ | ✓ | ✓ / ✓ | ✗ $+0.5\%$, $+0.06\%$ | ✗ (§4) | forced by consistency and the zero reading |
| readout, acceleration weight | ✓ | ✗ $-5.1\!\times\!10^{-3}$ | ✗ $-1.1\%$, $-0.6\%$ | ✓ | – / ✓ | – | ✓ | forced by common force and the swap |
| **push scale** | ✓ | ✓ | ✓ | ✓ | ✓ / ✓ | ✓ | ✓ | **free** |
| **prior shape** | ✓ | ✓ | ✓ | ✓ | ✓ / ✓ | ✓ | ✓ | **free** |
| **readout, direction weight** | ✓ | ✓ | ✓ | ✓ | – / ✓ | – | ✓ | **free unless C3** |

**Consistency, the numbers.** Shift of the median in $\log m$ at
$\varepsilon=0.1$, against ours27 on the same data:

| bend | $d$ | $N=100$ | $N=800$ | $N=6400$ | $N=51{,}200$ |
|---|---|---|---|---|---|
| push units | 1 | −0.251 | −0.209 | −0.236 | |
| push units | 3 | −0.966 | −0.429 | −0.527 | |
| prior count | 1 | −0.207 | −0.168 | −0.190 | |
| prior count | 3 | −0.498 | −0.151 | −0.199 | |
| push size (mean of 4 seeds) | 1 | | $+7.5\!\times\!10^{-4}$ | $-6.3\!\times\!10^{-4}$ | $+1.8\!\times\!10^{-4}$ |
| push size | 3 | $-3.6\!\times\!10^{-3}$ | $+3.1\!\times\!10^{-4}$ | $+5.1\!\times\!10^{-5}$ | |
| prior shape | 1 | $-2.0\!\times\!10^{-3}$ | $-2.0\!\times\!10^{-4}$ | $-2.9\!\times\!10^{-5}$ | |
| prior shape | 3 | $-3.3\!\times\!10^{-3}$ | $-1.8\!\times\!10^{-4}$ | $-3.1\!\times\!10^{-5}$ | |

The push-units and prior-count shifts do not shrink: the estimate settles on a
wrong mass (about 4 against a true 5). The push-size and prior-shape shifts die
away. The push-scale bend is consistent by the tempering identity
(`../ours26/ours26_learned_push_scale_2026-09-27.md` §2): in noise units its
law is $\pi$ times the flat likelihood raised to a power, with the same peak.

**Read by column: what each property forces.**

| property | bends it kills |
|---|---|
| axis splitting | push size, joining rule, prior count |
| swap | push units, readout acceleration weight |
| common force | readout acceleration weight |
| zero reading leaves the law unchanged | push units, prior count, noise law |
| consistency | push units, prior count, noise law |
| lumping, units | none of these (both are already built into every bend) |

Two moves rest on a single property: the flat push size and the stacked length
each rest on axis splitting alone. Every other forced move is forced at least
twice.

## 4. The noise law is forced

The t bend breaks the zero-reading property, as found earlier tonight. A zero
reading's evidence under the flat push varies with the angle unless the noise
is Gaussian. The consistency test shows why that matters.

**Where each estimator converges.** 1D, $s=2$, true $m=5$, pushes $\pm r^\ast$.
The limit is the peak of the expected per-reading log-evidence, computed by
quadrature over the noise on a 0.1 grid; the grid floor is about $+0.001$ in
$\log m$.

| world | model | $r^\ast=1$ | $r^\ast=3$ |
|---|---|---|---|
| Gaussian | Gaussian | 5.0045 | 5.0045 |
| t₅ | **t₅** (the correct model) | **5.6389** | 5.0670 |
| Laplace | **Laplace** (the correct model) | **6.1403** | 5.1125 |
| t₅ | **Gaussian** (the wrong model) | **5.0038** | 5.0036 |

With the flat push, a heavy-tailed model is biased in its own world: 13% high
for t₅ and 23% high for Laplace at a push the size of the noise. The bias does
not depend on $N$. The Gaussian model converges in the t₅ world.

**Why.** The Gaussian model's expected per-reading criterion is

$$
E\big[(x\cos\theta-y\sin\theta)^2\big]=r^{\ast2}\sin^2(\theta-\theta_0)+\operatorname{Var}(e),
$$

which peaks at the true angle for any noise with the declared variance, equal
in the two channels in noise units. It uses only second moments, and those are
what the instrument declares. A non-Gaussian model with the flat push carries a
per-reading factor in the angle: the zero-reading evidence that varies with
$\theta$. That factor drags the answer, exactly as T7 predicts.

**Consequence for the lemma record.** L13 ("errors Gaussian …: declared, never
argued") becomes a theorem, given the flat push: consistency forces the Gaussian
model, and the Gaussian model needs only the declared SDs to be consistent. What
stays declared is that the SDs, and ours27's covariances, are right. The earlier
worry that real instruments are not Gaussian is answered: they need not be.

**Limits of this.** Tested against two heavy-tailed families (t₅ and Laplace),
in 1D, with pushes of fixed length. Not a proof that every non-Gaussian model
fails. Efficiency, meaning how wide the answer is when the noise is really
heavy-tailed, was not examined.

## 5. The three free directions

- **The starting law's shape.** Known to be free (lemma record §3.1). The bend
  here was contribution 12's tilt. A change of $\pi$'s centre was not tested:
  it is excluded by L3′ with one noise pair (the sweep's E1) and declared with
  several.
- **The readout direction weight.** It passes every property, including axis
  splitting, since the stacked direction survives regrouping. It was already
  known to pass composition, the swap and units (lemma record §4.2). Only C3
  (L6) excludes it, and C3's justification is still an argument, not a proof
  (§3.2 there). The size of its effect here: +0.03% at $\varepsilon=0.1$
  (per-axis).
- **The push scale $\omega$.** New as a clean result. In noise units $\omega$
  is dimensionless, so "no built-in scale" (L3) does not exclude it.
  - A Gaussian push weight is a product over axes, so axis splitting holds.
  - It is isotropic.
  - It is symmetric in the channels, because the push is shared.
  - It leaves a zero reading uninformative.
  - By the tempering identity it has the same peak as the flat weight.

  Only L10 ("no push scale", $\omega\to\infty$) excludes it. L10 is therefore
  an independent premise, not a consequence of the others. The many-reading
  width problem lives on exactly this free direction (`README.md` in
  `../ours26/`, failure (b); `ours26_learned_omega_is_a_hint_2026-09-27.md`).
  Effect on the point at $\omega^2=10$: −0.01% (one pair), −0.5% (per-axis).

**The uniqueness statement this supports.** Given C3 and "no push scale", the
properties force ours27 up to the starting law's shape. That covers composition,
the swap, units, rotation, axis splitting, a zero reading leaving the law
unchanged, and consistency. Each of the three residual freedoms is already on
the open list:

- the starting law's shape (lemma record §3.1);
- C3's justification (§3.2 there);
- the push scale (§3.3 there, and the width problem).

## 6. What this changes, proposed not made

For the next lemma record:

- **L13:** Gaussian noise moves from declared to forced, given the flat push
  and consistency (§4).
- **L10:** marked as the premise behind one of only three free directions, and
  not implied by L3.
- **L8, axis splitting:** marked as the sole support of the flat push size and
  of the stacked length. A record defending axis splitting would harden two
  moves at once.
- **L6, C3:** as before, the premise behind one free direction.

## 7. What was not tested

- Bends that break rotation. Only ours26 (diagonal-only noise) was covered, and
  that was earlier tonight.
- Correlation between the channels as a bend. ours27 allows it; no bend varied
  it.
- The centre of $\pi$ (see §5).
- More than one bend at a time. A pair of bends that together restore a broken
  property would be a free direction this matrix misses. The push-units bend
  combined with a compensating prior change is the obvious candidate to check.
- Coverage and widths under any bend.
- The size of the break for the push-size and joining-rule bends is small at
  these data ($10^{-4}$), though first order. The large version of the push-size
  bend is ours21 ($\varepsilon=2-d$), whose axis-splitting failure is measured
  in `axis_splitting_and_the_24_september_change_2026-09-25.md` §4, and the
  large version of the joining rule is ours24's plain sum ($p=1$), measured in
  `ours25b_the_many_reading_readout_2026-09-25.md` §5.

## 8. Reproduction

Scripts written in Claude's workspace on 27 Sept and handed to the user with
this record; they are not in this repo.

- `bend.py`: ours27 with diagonal noise, each bend as a reweighting of the
  conditional Gaussian candidates, and common random numbers.
- `matrix.py`: the swap, units, axis-splitting and zero-reading columns, for
  both test problems.
- `matrix2.py`: composition at the readout level, and consistency at growing
  $N$.
- `matrix3.py`: the 1D push-size consistency, with four seeds up to
  $N=51{,}200$.
- `limit2.py`: the noise-law limits of §4, by quadrature.

numpy and scipy only. The matrix run takes about 15 minutes, the limits about
6.

## 9. How we got here

1. ours27 was named the same evening.
2. The user asked what to do with it next, listing a new lemma record, a record
   of ours27, a check against data, and "check it against all the algebraic
   moves. See if slightly changing the algebra changes what we can do with it",
   which he thought might be best. Claude agreed and proposed the bends.
3. First run: two artefacts were fixed before any cell was read.
   - A half-grid-step bias in the median made every swap cell read −6.4%,
     including unbent ours27.
   - The prior-count bend first counted $\pi$ as a density in $m$, which breaks
     the swap for a reason unrelated to counting.
4. The noise-law consistency first ran as a simulation too noisy to read: the
   exactly consistent Gaussian control read −0.048. It was redone by quadrature.
5. The push-scale row was added last, when it became clear that $\omega$ in
   noise units is dimensionless and might pass everything. It did.
6. The user: "Make a folder called ours27 put it there. Make a record of it."

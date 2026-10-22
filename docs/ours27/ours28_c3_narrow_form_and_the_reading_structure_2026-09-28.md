# C3's narrow form: the surviving family, its size, and the premise that closes it

**2026-09-28 (late night).** Record, written by Claude in the session that
filed the ours28 reference record, on the user's instruction "Do the C3
thing." It edits nothing and binds nothing. It takes the item left open by
the no-borrowing record (§5) and the reference record (§5): C3's narrow form
under the ruling that the definition borrows nothing between readings.

**In one paragraph.** Without regrouping across readings, the ruled premises
leave exactly one family of readouts besides C3: ratios of means weighted by a
function $h$ of how the candidate's size is shared among the readings. On the
two ours28 examples, mild weights move the point by 0.1–0.8%. The premises put
no bound on $h$ (about ±15–20% at the extreme tried), though every value
stays inside the law's interval. One premise removes the whole family: the
ruling "same candidate law, same mass, whichever instruments produced it",
read so that **the reading structure** (which components were recorded
together as one reading) is one of the declared facts about the instruments,
like $\Sigma_i$ and $N_i$. The law uses it: a reading is one push. The score
does not see it. Then the one score has to respect the own frames of every
reading structure, those rotations generate every rotation of the stacked
candidate, and the C3 derivation runs to the end. Nothing is borrowed between
readings. The reading of the ruling is the user's to confirm (§6).

**Ruled the same night (§9):** the user confirmed it. C3 is derived.

---

## 1. What was open

`ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md` §1 derived C3
in three steps:

1. betweenness makes the readout the root of a fixed score,
   $E_P[\varphi(c,r)]=0$;
2. the bookkeeping moves make the score depend only on $(X,Y)$;
3. composition leaves the ratio of means.

Step 2 used regrouping axes across readings as a symmetry of the data. Under
ours28 it is not one. A reading's components share one push scale, so a
regrouped series has a different law: the six-reading example split into
eighteen 1D readings without the scale carried gives readout 4.0807 against
3.7328. The no-borrowing ruling keeps it that way, because a reading is one
push. So step 2 had only each reading's own rotations and signs, and the
permutations of readings. C3's narrow form stayed a premise, folded into L12
(`ours28_full_covariance_form_and_declared_subspace_2026-09-28.md` §4).

## 2. What survives, exactly

Premises used:

- betweenness (provisional);
- each reading's own frame and sign convention;
- permutations of readings;
- the swap and units;
- composition (lumping and common force, the 6 Sept theorem);
- same-law-same-mass, as applied to the declared noise in the C3 record §1.3.

The steps:

- **Own frames and permutations.** The score sees a candidate only through
  $m$ and the per-reading sizes $\lVert\mathbf a^\ast_i\rVert$, symmetrically
  (among readings of the same dimension).
- **Total and shape.** Write the sizes as the stacked size
  $Y=\lVert\mathbf a^\ast\rVert$ and the shape
  $u_i=\lVert\mathbf a^\ast_i\rVert/Y$. A split of mass or acceleration common
  to a candidate's readings scales every size alike, so lumping and common
  force act at fixed shape. On each shape the 6 Sept theorem runs as it
  stands and leaves $X-rY$, times a positive factor that may depend on the
  shape and on $r$.
- **Units.** Rescaling force multiplies $r$. For the root to move with it for
  every law, the factor's dependence on $r$ must separate from its dependence
  on the shape. So the factor is $h(\mathbf u)$.

The survivors are therefore exactly

$$
\widehat m_h=\frac{E\big[h(\mathbf u)\,\lVert\mathbf F^\ast\rVert\big]}{E\big[h(\mathbf u)\,\lVert\mathbf a^\ast\rVert\big]},\qquad h>0\ \text{and symmetric in the readings},
$$

and $h$ constant is C3.

Each $\widehat m_h$ passes betweenness (it is linear), own frames,
permutations, the swap (the shape of the forces $m\lVert\mathbf a^\ast_i\rVert$
is the shape of the accelerations), units, lumping and common force. It
changes nothing in the law, so the law's properties (the zero-reading
property, calibration) are untouched.

An example: $q=\sum_iu_i^4$. It is 1 when the whole candidate push is in one
reading and $1/N$ when it is shared evenly.

## 3. How much it matters

Candidates are drawn from ours28's own conditional law (reference record §1),
4000 per trial mass with common random numbers, on the two examples of the
reference record. The table gives the change of the readout relative to $h=1$
in the same draws:

| $h$ | one pair, six weak 3D readings | full covariance, three readings (one $k=1$) |
|---|---|---|
| $1$ (C3) | 3.7374 (exact 3.732789) | 2.2476 (exact 2.248458) |
| $1+q$ | +0.34% | +0.08% |
| $1+5q$ | +0.81% | +0.16% |
| $1-\tfrac12q$ | −0.30% | −0.08% |
| $q^4$ | +4.01% | +0.76% |
| $q^{-4}$ | −4.19% | −0.68% |
| law's 95% interval | [0.49, 64.7] | [0.75, 25.2] |

On the six weak readings, stronger weights give:

| $h$ | $q^{\pm8}$ | $q^{\pm16}$ | $q^{\pm32}$ |
|---|---|---|---|
| change | +5.7% / −6.6% | +10.2% / −9.5% | +20.1% / −15.6% |

- The effect is small against the law's width, and smaller with stronger
  data. That is expected: once the shape is pinned at each trial mass, $h$
  only reweights the law smoothly.
- The premises put no bound on it. So without a further premise, each $h$ is
  a rival definition, and they differ where the evidence is weak.

## 4. The premise that closes it

**The reading structure is declared.** A reading declares its data
$\mathbf y_i$, its covariance $\Sigma_i$ and its candidate subspace $N_i$.
Declared with them is which components belong to it: the reading structure of
the series.

The C3 record §1.3 already settled how declared instrument facts enter:

- **The law may use them.** ours28's $H_i$ uses $\Sigma_i$ and $N_i$, and the
  starting law may use the instruments.
- **The score may not.** A score that read the declared noise directions was
  excluded by "same candidate law, same mass, whichever instruments produced
  it".

**Apply the same to the reading structure.** The law uses it: a reading is one
push, with one scale, and nothing is borrowed. The score does not see it.
Then:

1. There is one score $\varphi$ on stacked candidates, whatever the reading
   structure.
2. Own frames must hold for whatever structure the data have. So $\varphi$ is
   unchanged by each reading's rotation, for every way the components could
   have been grouped into readings.
3. Any two components can lie in one reading, so every rotation in a
   coordinate plane of the stacked candidate is included. These generate all
   rotations (the argument of the C3 record §1.1, step 2).
4. So the score sees only $(X,Y)$, and step 3 of the C3 record gives the
   ratio of means.

**What it does not use.** Regrouping is not asked to be a symmetry of the
data, and the law is not asked to be the same under regrouping. The 4.0807
against 3.7328 stays, because different laws may give different masses. Only
the score is one function. Nothing learned from one reading is used on
another, so the no-borrowing ruling is untouched. ours28 already reads out
this way.

**Where the family breaks it.** ours27 is kept by the ruling. Take the
six-reading example (one pair, diagonal noise). Entered as six 3D readings,
or as the same eighteen numbers as eighteen 1D readings, it gives the same law
(to $5\times10^{-15}$) and the same candidates, draw for draw. But $q$
computed from the two structures differs, and so does the mass:

| $h$ | as six 3D readings | as eighteen 1D readings |
|---|---|---|
| $1$ (C3) | 4.1936 | 4.1936 |
| $1+q$ | 4.2015 | 4.2020 |
| $1+5q$ | 4.2163 | 4.2217 |
| $q^4$ | 4.4116 | 4.5716 |
| $q^{-4}$ | 4.0742 | 4.0120 |

These use 8000 draws per trial mass with common random numbers; the exact C3
readout is 4.193908. So every non-constant $h$ gives one candidate law two
masses.

## 5. Why the no-borrowing ruling cannot close it by itself

Two readings of the ruling were tried, and both fail.

- **"Prefer no sharing of push size among readings."** Applied to the score,
  it is exactly the rotation invariance of §4. But as a principle it binds
  the law too. A law with no preference about how push size is shared among
  readings has one scale for the series (Kingman 1972, C3 record §3.1), and
  that is borrowing. So it cannot be the principle. Only the law/score
  division of §4 binds the score and leaves the law free.
- **"The score factorises over readings."** It would also exclude the stacked
  size $Y$, which composition requires. It excludes C3 itself.

## 6. Where C3 stands

- **If the reading structure is among "which instruments produced it":** C3 in
  full follows from betweenness (provisional), own frames and signs,
  same-law-same-mass, and composition. It leaves the premise list, and the
  narrow-form clause folded into L12 becomes a theorem. The betweenness
  representation step is still owed in writing, as before.
- **If it is not:** C3's narrow form stays a premise. The family of §2 is the
  freedom it removes, of the size in §3.
- **Either way** the equation and the code are unchanged.

The question, in plain words: *six 3D readings, or the same eighteen numbers
recorded as eighteen 1D readings, under a law that cannot tell them apart:
one mass or two?* One mass is the premise of §4.

## 7. Reproduction

Scripts in Claude's workspace, not kept. They use numpy and
`ours28_reference.py` in this folder.

- **§3.** The law is on a 361-point grid in $\log m$. Per reading, at each
  trial mass, $\beta$ is drawn from its weights by inverse CDF and
  $\boldsymbol\alpha_i=\beta\boldsymbol\mu+\sqrt\beta\,G_m^{-1/2}\mathbf z$,
  with the same uniforms and normals at every trial mass. The examples are
  those of the reference checks (seed 7).
- **§4.** `ours28(..., flat=True)` gives the law. The candidates are drawn
  per component from the flat conditional
  $N\big(v(mF/\sigma_F^2+a/\sigma_a^2),\,v\big)$ with
  $v=(m^2/\sigma_F^2+1/\sigma_a^2)^{-1}$ (seed 1).

## 8. How we got here

The reference record listed C3's narrow form as not done, and the user said
"Do the C3 thing." Claude:

1. wrote down what the ruled premises leave without cross-reading regrouping
   (§2);
2. measured it on ours28 (§3);
3. looked for a premise already ruled that removes it, and rejected two
   readings of the no-borrowing ruling (§5);
4. found the reading of same-law-same-mass in §4. It is the reading the C3
   record §1.3 already used for the declared noise.

## 9. Ruling

The user, 28 September, 23:41, answering §6: "yes they can be a proprty of
the instrument". The reading structure is a declared fact about the
instruments, like $\Sigma_i$ and $N_i$. The law may use it and the score does
not see it.

Consequences:

- **C3 is derived**, in full, from betweenness (provisional), own frames and
  signs, same-law-same-mass and composition. It leaves the premise list.
- **The narrow-form clause folded into L12** by the full-covariance record §4
  is a theorem, not a premise.
- **The profile-weighted family of §2 is excluded.**
- **Still owed:** the betweenness representation step, written out.
- **Of the three freedoms of 27 September** (the starting law's shape, C3,
  the push scale), C3 is gone. The push scale's structure is fixed by the
  no-borrowing ruling (one scale per reading, from a declared law); its
  shape remains. What is left is two shapes, both currently set by one rule:
  uniform in the noise-defined angle.

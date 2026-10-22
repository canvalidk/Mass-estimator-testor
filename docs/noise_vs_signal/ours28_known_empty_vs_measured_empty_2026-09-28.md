# Known-empty against measured-empty under ours28: a declaration that now costs something

**2026-09-28 (night).** Record, filed at the user's request ("Sure go ahead"),
written by Claude from his story of the evening and the check made on it. It
edits nothing and binds nothing. It extends
`ours27_candidate_size_and_null_readings_2026-09-28.md` §5 to ours28
(`ours28_candidate_push_scale_as_nuisance_2026-09-28.md`).

**His story (28 Sept, his words lightly compressed).** Take two aligned
measured vectors with uncertainties. Rotate the axes so the signal sits on one
of them; the other axes then read 0/0. If we expect that reading to have the
same certainty as it would with the 0/0 axes absent, we miss the possibility
that the true vectors, though aligned, point in a different direction. Reading
in one dimension implicitly assumes the true direction is the one we rotated to.

---

## 1. Where the story needs one substitution

Direction by itself cannot touch the mass. Under L1 every candidate is an
aligned pair, $\mathbf F^\ast=m\,\mathbf a^\ast$, so if the true pair points
off the chosen axis both vectors gain the same off-axis factor and
$\lVert\mathbf F^\ast\rVert/\lVert\mathbf a^\ast\rVert$ is unchanged. This is
why ours27 gives a 3D reading with two null axes exactly the 1D law
($Q_k=0$ on a null axis) and why the candidate-size record §5 could resolve
the question with the point moving and the law standing still.

What a measured-empty axis carries is not direction but **scale**: it says
the push in that direction is at most noise-sized. Under ours27 that is inert,
because the flat push weight cannot express "the push might be small". Under
ours28 it is not inert: a reading's axes share one $\omega_i$, so two null
axes are evidence for a smaller $\omega_i$, which tempers the signal axis's
evidence. With that substitution the story's conclusion holds: recording in
1D is a declaration ("the push has no other components, by construction"),
and a declaration should buy certainty. Under ours27 it was free. Under
ours28 it is paid for.

## 2. The numbers, one noise pair

Aligned reading at $m=s$ ($\theta=45^\circ$), push of length $r$ in noise
units on axis 1, the other axes exactly 0/0, $\pi$ the half-Cauchy at $c=s$.
Median and 95% log-width of the law of mass, 4000-point grid.

| $r$ | ours27 (1D = 2D+null = 3D+nulls) | ours28, 1D declared | ours28, 2D + one null | ours28, 3D + two nulls |
|---|---|---|---|---|
| 1 | 6.17 | 6.37 | 6.39 | 6.41 |
| 2 | 5.26 | 5.90 | 6.05 | 6.13 |
| 4 | 2.61 | 3.04 | 3.32 | 3.65 |
| 8 | 1.04 | 1.06 | 1.08 | 1.09 |
| 16 | 0.50 | 0.50 | 0.50 | 0.50 |
| 32 | 0.25 | 0.25 | 0.25 | 0.25 |

The median is $s$ to four figures in every cell. With realised noise on the
null axes instead of exact zeros (2000 draws, $r=4$): ours27 2.62 against
its 1D 2.61; ours28 3D 3.53 against its 1D-declared 3.04.

**Reading.**

- Under ours28 the 1D declaration narrows the law where the signal is weak
  to moderate, by an amount graded by how many directions were left
  undeclared, and by nothing once the signal is strong. That is the shape
  the story asks for.
- The mechanism is the determinant factor: a reading's evidence is
  $(1-\beta_i)^{d/2}e^{\beta_i h_i^2/2}$ with $d$ the number of axes the
  reading *declares*. A 1D reading has $d=1$; the same numbers embedded in
  3D have $d=3$. The evidence functions are ${}_1F_1(\tfrac12;\tfrac32;\cdot)$
  and ${}_1F_1(\tfrac12;\tfrac52;\cdot)$ of the same $h_1^2$.
- This is a behavioural difference between ours27 and ours28 that the
  perturbation matrix has no column for: known-empty and measured-empty now
  differ in the **law**, not only in the point.
- It is consistent with axis splitting as ours28 reads it: splitting a
  reading into its axes carries the scale sharing along, so "axis 1 of a 3D
  reading" and "a 1D reading" are different objects, and the story says
  they should be.

## 3. What the story exposes: the units fork, sharpened

With one noise pair the above is complete. With per-axis noise the same
rotation finds the fork the lemma record §3.3 left open ("in which units
isotropy holds").

Give the empty axes very precise instruments. In their own noise units they
still read $N(0,1)$, so they still say "push $\lesssim1$ noise unit here".
With $\omega_i$ shared across the reading's axes **in noise units** (which is
what T7 and the push-units bend force), that discounts the signal axis by an
amount that does not shrink as the empty axes' instruments improve.
Physically that is backwards: a precise zero should establish the direction
and say nothing about the push along the track.

The alternatives each fail a forced move:

| where the shared scale lives | what goes wrong |
|---|---|
| noise units, shared across axes (ours28 as filed) | precise empty axes discount the signal axis; the discount is set by their noise units, not by what they measured |
| physical units, shared across axes | a $\theta$-dependent per-reading factor ($\omega_k=\omega_{\rm phys}/\sigma_{a,k}\cos\theta$): the push-units bend, forced out three ways (matrix §3) |
| one scale per component | isotropy and rotation invariance break (ours27 exists to keep them) |

**The candidate.** One $\omega$ per noise pair per reading. The sweep found
the same requirement for the shared case on other grounds (§5, "tempering
needs one $\omega$ per noise pair": one $\beta$ for all components leaves a
width ratio of 1.69, one per pair gives 0.97). A precise empty axis then
informs only its own pair's scale. What is unexamined is how "per noise
pair" behaves under a full covariance $\Sigma_i$, where the pairs are not
axis-aligned; presumably per principal axis of $\Sigma_i$, which is
rotation-covariant by construction, but that is a guess.

**Consequence for the owed list.** This moves the several-pair form of
$\omega_i$ from the bottom of the ours28 record's §9 to the top: ours28 as
filed is wrong in exactly the story's per-axis scenario, and the story is the
test any repair has to pass (precise empty axes → the 1D-declared law).

## 4. Reproduction

Script in Claude's workspace, not kept: the two laws on a 4000-point $\theta$
grid, ours27 as $e^{h^2/2}$ and ours28 as `hyp1f1(1/2, d/2+1, h^2/2)`, for
$d=1,2,3$; widths from the cumulative law; the realised-noise check draws
two $N(0,I_2)$ null axes and adds their $(\mathbf z\cdot\mathbf u)^2$ to $h^2$.

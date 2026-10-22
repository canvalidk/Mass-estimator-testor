# ours28: reference implementation, its readout, and its checks

**2026-09-28 (late night).** Record, filed at the user's request. He asked
"But were there fixes to ours28 you could do? Did you fix what you said could
be fixed?", was told three were fixable and not yet done, and said "Yes please".
Written by Claude in the session that filed the no-borrowing record. It edits
nothing and binds nothing.

**What it does.** It implements ours28 in the full-covariance form of
`ours28_full_covariance_form_and_declared_subspace_2026-09-28.md`, computes its
**readout** for the first time, and runs the property checks. The code is in
this folder: `ours28_reference.py`, `ours28_reference_checks.py` and
`ours28_reference_sim.py`. It closes three items of that record's owed list:

- the readout point and its interval (item 4);
- a reference implementation;
- the property checks, including whether Gaussian noise is still forced
  (item 2, in part).

---

## 1. The equation as implemented

A reading declares:

- data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$;
- a positive-definite error covariance $\Sigma_i$, with any correlation,
  including between force and acceleration;
- a candidate subspace $N_i\in\mathbb R^{d\times k}$ with orthonormal
  columns (default $N_i=I_d$, $k=d$).

The starting law $\pi$ is the half-Cauchy at $c$. The per-reading scale law
$G$ on $\beta=\omega^2/(1+\omega^2)$ is the arcsine law.

At trial mass $m$, with $C=B_mN_i$ and $B_m=(mI;I)$:

$$
G_m=C^\top\Sigma_i^{-1}C,\qquad \mathbf b=C^\top\Sigma_i^{-1}\mathbf y_i,\qquad
H_i(m)=\mathbf b^\top G_m^{-1}\mathbf b .
$$

**Law:**

$$
p(m)\propto\pi(m)\prod_i\int G(d\beta)\,(1-\beta)^{k_i/2}e^{\beta H_i(m)/2},
$$

which with the arcsine law is ${}_1F_1(\tfrac12;1+\tfrac{k_i}2;\tfrac{H_i}2)$
per reading.

**Candidates.** This is the new piece, needed for the readout. The push is
isotropic in the whitened candidate subspace with scale $\omega$, so the
whitened data are $N(0,I+\omega^2P_m)$. The whitened candidate given the
data is then $N(\beta P_mz,\beta P_m)$. In the subspace coordinates
$\boldsymbol\alpha$ (with $\mathbf a^\ast_i=N_i\boldsymbol\alpha_i$):

$$
\boldsymbol\alpha_i\mid m,\beta_i\sim N\big(\beta_iG_m^{-1}\mathbf b,\;\beta_iG_m^{-1}\big),
\qquad
\beta_i\mid m\;\propto\;G(d\beta)(1-\beta)^{k_i/2}e^{\beta H_i(m)/2},
$$

independently per reading. At $\beta=1$ this is ours27's conditional
($K_i=G_m^{-1}$ when $N_i=I$). Then $\mathbf F^\ast_i=m\,\mathbf a^\ast_i$.

**Readout.** Unchanged in form (T9, T10):
$\widehat m=\int mA(m)p(m)\,dm\big/\int A(m)p(m)\,dm$, with
$A(m)=E[\lVert\mathbf a^\ast\rVert\mid m]$ over the stacked candidates. $A$
comes from the Laplace transform of $\sum_i\lVert\boldsymbol\alpha_i\rVert^2$,
which factorises over readings given $m$, with each reading's $\beta$ mixed in.

**Numerics.**

- $\beta$ on a 481-point logit grid; the law on a 721-point grid in $\log m$
  (half-width 9).
- For the readout, the $\beta$ weights are coarse-grained into bins of 8 and
  the Laplace integral uses 241 points. On the six-reading example below this
  changes the readout by $4\times10^{-5}$ relative (3.732955 at full
  resolution, 3.732789 coarse-grained).

## 2. Checks

| check | data | result |
|---|---|---|
| per-reading evidence against ${}_1F_1(\tfrac12;1+\tfrac k2;\tfrac H2)$ | $H$ = 0.5, 7, 40, 200; $k$ = 3, 3, 1, 2 | agree to $4\times10^{-7}$ in the log |
| flat limit against ours27 | one pair, two 3D readings; and per-axis | readout 1.983766 = 1.983766; 1.917844 = 1.917844 |
| swap, $c\to1/c$ | one pair, six weak 3D readings | readout product 1.000000, median product 1.000000 |
| swap, units (force ×3) | full covariance with cross-channel correlation, one reading with a declared 1D subspace | products 1.000000; ratios 3.000000 |
| each reading rotated with its covariance and its subspace | same | readout $3\times10^{-13}$, median $4\times10^{-16}$ |
| exact zero reading appended: the law | both examples | unchanged ($8\times10^{-17}$, $4\times10^{-16}$) |
| exact zero reading appended: the readout | six weak readings | 3.732789 → 3.722316: the point moves, the law does not, as in ours27 |
| zero experiment alone | three 3D zero readings, $c=s=2$ | median 2.000000, readout 2.000000 |
| declared subspace | 3D reading $F=(4,0,0)$, $a=(2,0,0)$, unit noise | $k=3$: readout 1.8836, 95% [0.552, 16.5]; $k=1$ along $x$: 1.9174, [0.664, 14.6]; the same numbers as a 1D reading: 1.9174, [0.664, 14.6], identical |
| totality | perpendicular, anti-parallel, force only, acceleration only | readouts 1.3015, 1.3300, 2.2914, 0.4364, all finite; the last two multiply to 1.000 (swap) |

**Axis splitting needs the scale carried.** The six 3D readings, split into
eighteen separate 1D readings each with its own scale, give readout 4.0807
against 3.7328. A reading's axes share its scale. Splitting is the same law
only when that shared scale travels with the components (L8 as restated in
the full-covariance record §4). This is the no-borrowing ruling's "a reading
is one push" seen from the other side.

**The declared subspace is known-empty, exactly.** A 3D reading declared
$k=1$ along an axis gives the same law and readout as the 1D reading, to every
digit. So "recorded in 1D" and "recorded in 3D with the direction declared"
are one input.

**Gaussian noise is still forced**, now by the zero-reading property alone.
Under ours28's push, a zero reading's evidence must not depend on the mass
angle. Take one component in noise units, integrated over the arcsine scale
law, at mass angles 2.9°, 13.4°, 24°, 34.5° and 45°, relative to 45°:

- Gaussian noise: 1.0000 at every angle;
- Student-$t_5$ noise (unit variance): 1.0158, 1.0121, 1.0064, 1.0017,
  1.0000.

So heavy-tailed noise would let a zero reading tilt the law of mass. The
matrix's consistency route to L13 used the flat push and has not been redone
for ours28.

## 3. The readout under repetition

One noise pair, $\sigma_F=1$, $\sigma_a=0.5$ ($s=2$), true $m=2$, starting law
centred on it. 3D readings with Gaussian pushes of per-component SD $\omega$
in noise units: ours28's own world. "Law SD" is the median 95% log-width
divided by 3.92.

| readings | $\omega$ | series | law coverage | readout / true, median | law median / true | SD of log readout | law SD |
|---|---|---|---|---|---|---|---|
| 5 | 1 | 200 | 98.0 | 0.913 | 0.898 | 0.531 | 1.130 |
| 21 | 1 | 150 | 93.3 | 1.012 | 1.012 | 0.422 | 0.460 |
| 21 | 3 | 150 | 95.3 | 0.995 | 0.995 | 0.098 | 0.093 |
| 100 | 1 | 60 | 96.7 | 1.040 | 1.040 | 0.169 | 0.166 |
| 400 | 1 | 40 | 92.5 | 1.036 | 1.036 | 0.084 | 0.084 |
| 1600 | 1 | 6 | 6 of 6 | 0.972 | 0.972 | 0.039 | 0.043 |

- **The readout goes to the truth**, within about two standard errors of the
  median at the larger sizes, where the series counts are small.
- **From 21 readings up, the readout's own scatter over repeats matches the
  law's width** (0.42/0.46, 0.098/0.093, 0.169/0.166, 0.084/0.084,
  0.039/0.043). So in the regular regime the law's interval is also the
  readout's interval. That is the owed "interval of the point", answered for
  this world.
- **With 5 weak readings the law is far wider than the readout's scatter**,
  because the starting law dominates. The readout sits close to the law's
  median throughout (they agree to 0.1% from 21 readings up).

These worlds are ours28's own. Worlds whose pushes share an undeclared
pattern are covered by the no-borrowing record (for example, 84% when every
push is one noise unit).

## 4. Worked comparison

| data | ours27 | ours28 |
|---|---|---|
| one pair, six weak 3D readings, $c=2$ | 4.1939, median 4.276, [1.61, 31.4] | 3.7328, median 4.014, [0.49, 64.7] |
| full covariance with cross-channel correlation, three readings (one declared $k=1$), $c=2$ | 2.2066, median 2.238, [0.99, 12.4] | 2.2485, median 2.360, [0.75, 25.2] |

With weak data ours28 is honestly wider than ours27 and leans on the starting
law. With strong data the two converge.

## 5. Still not done

- C3's narrow form (no weight on the profile of per-reading sizes) is a
  premise; it has not been derived within the no-borrowing ruling.
- The betweenness representation step, written out.
- The full perturbation matrix on ours28 (bends of every move), and the
  consistency route to Gaussian noise under the tempered evidence.
- The tester build (`canvalidk/Mass-estimator-testor`): the reading type with
  $N_i$, the law, the readout, and presets for the worlds above.
- Per-axis coverage cells, and the physically-isotropic-push world under an
  anisotropic instrument (full-covariance record, owed item 3).

## 6. Reproduction

- `ours28_reference_checks.py` reproduces §2 in about 25 s.
- `ours28_reference_sim.py N omega series seed` reproduces a row of §3. Rows
  with many readings take tens of minutes.
- The Student-$t$ check used scipy quadrature over the push for each of 40
  arcsine-spaced scales; that script was not kept.
- numpy and scipy only.

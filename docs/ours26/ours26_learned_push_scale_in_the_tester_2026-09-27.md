# The learned push scale in the tester: a correction, the ours26 extension, and the study

**2026-09-27 (evening).** Record, filed at the user's request ("Great! Yes do those.":
record the learned push scale, then build it into the tester), written by Claude.

**Relations.**
- **Corrects** `ours26_learned_push_scale_2026-09-27.md` on three points (§1): its
  status line ("not yet the tester"), §6 ("which prior: open") and §9 ("the tester
  implementation is to follow"). Its algebra, the fork table and the standalone
  coverage table stand.
- Extends the tester's `hierarchical` rule (24 September) to ours26. Binds nothing
  else.

---

## 1. The correction: the tester already had a learned push scale

The tester's `hierarchical` rule, built on 24 September (commit `0c3235f`), is the
noise-unit learned-$\omega$ model of the earlier record:

- every standardised push gets $N(0,\omega^2I)$, with $\beta=\omega^2/(1+\omega^2)$;
- $\beta\sim\mathrm{Beta}(1,b)$ is integrated out exactly (a lower incomplete gamma
  function), so the prior on $\omega$ is **proper**. It uses the noise as its scale,
  as the earlier record's §6 said a proper prior must. Only $b$ is a choice;
- `excitation_weighting_premise_2026-09-24.md` §6 already showed it repairing the
  undercoverage with 50 readings: 80% intervals covering 81.4–87.6% against
  52.8–56.8% for the flat limit.

So the earlier record was wrong to call the prior open and the tester work unbuilt.
Claude wrote it without reading the tester. What that record adds that was not
already there:

- the units fork decided by consistency (its §4);
- the tempering identity stated as such, and the exact match
  $1/\hat\beta\approx1+M/R^{\ast2}$ with lemma T18 (its §§2–3);
- the observation that $d\omega/\omega$ fails (its §6), which explains why a proper
  $\beta$ prior is needed.

## 2. The extension: ours26 and the stacked readout

Before today the `hierarchical` rule refused per-axis noise and had only the plain-sum
readout. Now:

- **Per-axis noise.** Each component's push is standardised in its own noise pair's
  units, and all share one $\omega$ (isotropy in noise units, the reading the fork
  selects). At fixed $\beta$ the law is ours26's likelihood raised to the power
  $\beta$. The plain-sum readout stays refused, as for ours26.
- **The stacked readout** `ratio_of_means_rss`. Given $\theta$ and $\beta$ the stacked
  push is $N(\beta\mathbf W,\beta I)$, so the weight is $E\lVert\mathbf a^\ast\rVert$
  of that: $\sigma_a\cos\theta\sqrt\beta\,\mu_M(\sqrt\beta H)$ for one noise pair,
  ours26's weighted noncentral chi mean with both arguments scaled by $\beta$
  otherwise. It is averaged over $\beta\,\vert\,\theta$ by tanh-sinh nodes in
  $\beta$'s conditional quantile (22 nodes; about $10^{-11}$ against direct
  quadrature over $\beta$).
- **Tests** (`tests/test_hierarchical_ours26.py`, 19): the tempering identity
  ($10^{-10}$); the peak unmoved for $\beta$ = 0.05, 0.5, 0.99; the law against
  brute-force integration of every push and $\beta$ with per-axis noise ($10^{-6}$);
  the shrunk weight against $10^6$ Monte Carlo draws; the two implementations of the
  weight against each other ($10^{-9}$); the $\beta$ average against quadrature
  ($10^{-9}$); $\beta\to1$ gives ours26 ($10^{-8}$); a 3D series and its axes agree
  ($10^{-8}$); the swap; units; the refusals. One older test pinned the refusal this
  lifts and was updated. Full suite: 861 pass.

## 3. The study: `presets/ours26_learned_omega.yaml`

3D, a new push per reading, fixed length, random direction (not the Gaussian the rule
assumes). Masses 0.5, 2, 8 (with $\sigma_F=1$); acceleration SNR 1 and 3; $N$ = 2, 6,
21. Two noise settings: one noise pair for all axes ($s=1$), and ours26's per-axis SDs
(ratios 5, 2, 0.8, prior centre 2). 36 cells × 400 series. `ours26` against `learned`
($b=1$). Monte Carlo error: about 1.1 points near 95%, 2 near 80%.

**Coverage, summarised.**

| noise | interval | ours26: range, mean miss | learned: range, mean miss | cells within 2.5 points: ours26 / learned |
|---|---|---|---|---|
| one pair | 95% | 74.5–98.0, 3.4 | 94.2–99.8, 1.8 | 12 / 12 of 18 |
| one pair | 80% | 52.5–88.5, 5.7 | 77.8–97.5, 5.4 | 9 / 10 |
| per axis | 95% | 85.8–97.0, 2.4 | 92.5–99.2, 1.6 | 13 / 16 |
| per axis | 80% | 66.8–83.5, 3.2 | 78.2–92.8, 2.7 | 11 / 13 |

**The cells where it matters** (SNR 1; 95% / 80% coverage; median $\widehat m/m$ of
the readout):

| $m$ | $N$ | noise | ours26 | learned | median $\widehat m/m$, ours26 → learned |
|---|---|---|---|---|---|
| 0.5 | 21 | one pair | 74.5 / 52.5 | 99.0 / 90.8 | 0.96 → 1.27 |
| 0.5 | 6 | one pair | 79.0 / 58.8 | 99.0 / 94.0 | 1.48 → 1.69 |
| 0.5 | 21 | per axis | 85.8 / 66.8 | 92.8 / 78.2 | 0.99 → 1.05 |
| 2 | 21 | one pair | 90.8 / 74.2 | 94.8 / 84.2 | 0.99 → 0.99 |
| 2 | 21 | per axis | 87.5 / 72.0 | 92.5 / 79.5 | 0.96 → 0.95 |
| 0.5 | 2 | one pair | 94.2 / 81.2 | 99.8 / 97.5 | 1.68 → 1.85 |

At $m=8$, and at $m=2$ with SNR 3, the two agree to within about a point and the
points to under 1% (typically exactly). At $m=0.5$ with SNR 3 the learned intervals
are still up to 3.5 points wider in coverage, with the points unchanged from $N=6$.

**Reading it.**

- **The correctable failure is largely repaired.** The flat law's undercoverage
  (down to 74.5% at 95%, 52.5% at 80%) goes. The worst learned cell is 92.5% at 95%.
- **It overshoots where the data are weak.** With few readings or weak pushes,
  $\beta$'s law is broad and the intervals become too wide (99.8% at nominal 95%).
  $b=1$ is the obvious first knob.
- **The point moves toward the prior at weak signal.** The *likelihood's* peak is
  ours26's for every $\beta$. But a small learned $\beta$ flattens the likelihood, so
  the prior's centre pulls harder: at $m=0.5$ (centre $s=1$) the median ratio went
  from 0.96 to 1.27 with 21 readings. The point is still consistent (the tempered
  likelihood grows with $N$), but at $N=21$, SNR 1, the effective number of readings
  is about $\beta N\approx6$. Within ×1.25 still improved there (21.5% to 30.8%).
- **Per-axis noise calibrates only partly.** The world's pushes are isotropic in
  *physical* units, so in noise units they differ by axis (the acceleration SDs
  differ by a factor 6). One shared $\omega$ cannot match every axis. The earlier
  record's §8 predicted this: one $\omega$ calibrates exactly only if the pushes are
  equally large in each pair's noise units.

## 4. What this leaves

- **Keep:** the fork's answer for the *point* (noise-unit isotropy never moves the
  likelihood's peak). The learned scale as the in-chain fix for the many-reading
  widths.
- **Open:**
  - **One $\omega$ per noise pair.** This would calibrate the per-axis case. The peak
    then moves with the $\beta_g$, but the model contains both worlds of the fork: with
    pushes isotropic in physical units, each pair's noise-unit push scale is constant at
    the true angle. So it may be consistent in both. Untested; this is the natural next
    study.
  - **The hyperprior $b$**, for the over-coverage at weak signal.
  - **The pull toward the prior at weak signal:** accept it as what weak data mean, or
    count it against the learned rule. A ruling.
  - **A Gaussian-push world** for the same grid, where the rule's model is exactly
    right.

## 5. Reproduction

Tester `canvalidk/Mass-estimator-testor`: commit `d5af0ce` on top of `e288110`. The
changes were written by Claude into the user's clone on 27 September and committed
unchanged on 30 September. Files: `met/analyst.py`, `met/law.py`, `README.md`,
`tests/test_stacked.py`, and new `tests/test_hierarchical_ours26.py`,
`presets/ours26_learned_omega.yaml`. The study was run in Claude's sandbox on the same
files (about 70 minutes of computing).

```sh
python -m pytest -q tests/test_hierarchical_ours26.py
python -m met run presets/ours26_learned_omega.yaml
```

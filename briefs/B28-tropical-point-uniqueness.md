# B28: does tropical composition on E13's curve uniquely force $\arg\min K(m)$, and what happens at boundary minima?

**Status:** open
**Item:** `problems/O20.md` (finding 1, open question on E13's point, 22 Nov section), `problems/O2.md` (G6, O2);
T73, T74, T80 in `definition/theorems/`; L2, L3, L5, L7, L14, L19, L20 in `definition/premises.md`;
R-2026-10-28-1, R-2026-11-22-1, R-2026-11-22-2.
**Branch:** `codex/b28`, from `origin/claude/intake-2026-11-22` (or `origin/main` once `claude/intake-2026-11-22` is merged).
**Written:** 2026-11-22, Antigravity (at Can's request; provisional, pending Claude's review).
**Size (R-2026-11-07-2, merged as delivered):** one dated section in `problems/O20.md` (under ~100 lines), a short
cross-note in `problems/O2.md`, one test file `tests/test_tropical_point.py`, one `docs/STATE.md` log line. No change
to `met/`, `definition/premises.md` or `rulings.md`.

Read `problems/O20.md` (especially the 22 Nov section "The point on E13: two zero-padding endpoints, and tropical
composition of argmin K"), `definition/theorems/T73.md`, `T74.md`, `T80.md`, `problems/O2.md` (27–28 Oct),
`problems/O7.md` (B11 entry), and `docs/vd-docs-context/mass_readout_composition_uniqueness_2026-09-06.md`.

## Standard mathematics first

Infimal convolution (epigraphical addition) of cost functions is Rockafellar (1970), *Convex Analysis*, §5 (the
tropical / max-plus counterpart of probability-measure convolution; see also `T74.md`). The classical probability-law
composition uniqueness proof is `docs/vd-docs-context/mass_readout_composition_uniqueness_2026-09-06.md` (C1–C5
forcing $E[f]/E[\alpha]$); its finite-noise breakdown on probability laws of mass (Jensen between arithmetic and
harmonic means, plus `T19` infinite mean) is `problems/O2.md` (28 Oct). On the unnormalized curve $K(m)$ (`T73`),
`B16` proved $e^{-K/2}$ is non-integrable in $dm$ and $d\ln m$ on every finite series, ruling out integral functionals
without a starting law $\pi$.

## The sorted fork (provisional; do not re-sort it)

Under `R-2026-11-22-2`, E13 still owes a point on every input for type safety: where every positive mass is
compatible its value does not matter, and where the set is informative its value is open (`problems/O20.md`). Starting
directly at the $N_{\rm null}\to\infty$ zero-padding endpoint of a moving rule eliminates exact-zero drift:

- On the proper-law line (`T10`), that endpoint is the noise-angle circular mean $s\,E_p[\sin\theta]/E_p[\cos\theta]$
  (loses finite-noise `L5`/`L14` and retains an $s$ or $c$ anchor).
- On the unnormalized curve (`B11` / `T73`), that endpoint is $\arg\min_{m>0} K(m)$ (prior-free, exact `Z` across all
  instruments, and under `T74` infimal convolution satisfies **both** `L5` and `L14` simultaneously and exactly at
  finite noise wherever the minima are unique and interior).

The fork: **(1) Does tropical `L5` + `L14` on `T74` curves (with `L3`, `T5`, `L7`) uniquely force $R(K)=\arg\min_{m>0} K(m)$
on interior-minimum curves, or do other curve functionals pass too? And (2) what exact obstruction or admissible
extension exists when $\inf_{m>0} K(m)$ is non-unique or sits at the boundary $0$ or $\infty$ (especially on one-sided
compatibility sets `[L, ∞)` or `(0, U]`, which are not "every mass compatible")?**

## Deliverable

1. **Tropical composition uniqueness (interior-minimum domain).**
   Let $\mathcal K_{\rm int}$ be the class of finite-series `T73` curves $K:(0,\infty)\to[0,\infty)$ (and separately a
   natural broader epigraphical class closed under `T74` infimal convolutions $\square_+$ and $\square_\parallel$,
   explicitly stated) whose global minimum on $(0,\infty)$ is unique and interior. Require a readout
   $R:\mathcal K_{\rm int}\to(0,\infty)$ to satisfy:
   - `L5` under `T74`: $R(K_1 \square_+ K_2) = R(K_1) + R(K_2)$;
   - `L14` under `T74`: $1/R(K_1 \square_\parallel K_2) = 1/R(K_1) + 1/R(K_2)$;
   - `L3` (units covariance), `T5` (swap $K_{\rm swap}(m)=K(1/m)\implies R(K_{\rm swap})=1/R(K)$), `L20` (`Z`: shifting
     $K$ by a constant or adding $K_0\equiv 0$ changes nothing when defined), and exact-curve calibration (`L7`:
     $K(m)=0$ iff $m=m_0 \implies R(K)=m_0$).

   Either **prove** that $R(K)=\arg\min_{m>0} K(m)$ is the unique functional meeting these (naming the exact domain and
   any regularity condition needed — e.g. scale-invariance $R(\alpha K)=R(K)$, continuity, or closure under `T74`), **or
   construct an explicit counterexample functional** $R\neq\arg\min$ that passes all of them. Check explicitly whether
   level-set midpoints (arithmetic and geometric), tilted-minimum / Legendre-conjugate functionals, or multi-level
   averages survive or fail.

2. **Deterministic check of the two zero-padding endpoints.**
   In `tests/test_tropical_point.py`, verify deterministically:
   - `T10`'s circular-mean endpoint $\widehat m_\infty^{(s)}$: exact same-$s$ `Z`, failure of finite-noise `L5`
     (explicit two-point law counterexample), and anchor dependence on $s$ or $c$ (explicit `L17` counterexample for the
     $c$-anchored version on $\tfrac12\delta_1+\tfrac12\delta_4$).
   - `B11`'s endpoint $\arg\min_{m>0} K(m)$: exact `L5` and `L14` composition on `T74` curves (both Gaussian and genuine
     multi-reading vector `T73` curves with full correlated $\Sigma_i$ and subspaces $N_i$), exact `Z`, `L3`, `T5`, `L23`,
     and failure of level-set midpoints on `T74 (2)` Gaussian curves.

3. **Boundary minima and ties on the full `T73` reading type.**
   Characterize all valid measured inputs (`L2`: $\Sigma_i\succ 0$) where $\arg\min_{m>0} K(m)$ fails to be a unique
   interior positive mass:
   - all-zero (`0/0`, $K\equiv 0$, every mass compatible — covered by `R-2026-11-22-2`);
   - one-sided boundary infima (e.g. `F = 10 ± 1, a = 0 ± 1`, where $\inf_{m>0} K(m)=K(\infty)$ and the `T80` set is
     `[L, ∞)`, not all of `(0, ∞)`);
   - anti-aligned pairs (e.g. `F = (1,0), a = (-1,0)` with equal SD, where $K(0)=K(\infty)$ is a tie at both ends);
   - interior ties / multi-modal curves from conflicting readings (e.g. `{5 ± 0.1, 6 ± 0.1}` at small $\sigma$ or
     symmetric pairs).

   Check whether any total extension $R_{\rm ext}:\mathcal K_{\rm all}\to(0,\infty)$ agreeing with $\arg\min K$ on
   $\mathcal K_{\rm int}$ can preserve `L3`, `T5`, and `T74` composition (or continuity at the boundary), or prove an
   obstruction. Clearly separate what fails on `(0, ∞)` vs. on `[0, ∞]`.

4. **Write** items 1–3 into `problems/O20.md` as a dated section (with a short cross-note in `problems/O2.md`), add
   `tests/test_tropical_point.py`, and add one line to `docs/STATE.md`'s log.

## Do not

- Adopt $\arg\min K(m)$ or any boundary rule in `met/` or `definition/equation.md`.
- Introduce a starting law $\pi$ or a reference centre $c$ to rescue boundary minima.
- Change any premise row in `definition/premises.md` or edit `rulings.md`.
- Run large Monte Carlo coverage sweeps; deterministic algebraic/numerical checks and exact counterexamples are what
  this brief asks for.

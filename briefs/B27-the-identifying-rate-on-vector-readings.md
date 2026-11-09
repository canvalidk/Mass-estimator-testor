# B27: the identifying rate on vector readings

**Status:** open.
**Item:** G25 (L11's rate on the full type), `problems/O19.md`; T41, T42, T53; B19's item 3.
**Branch:** `codex/b27`, from `origin/main`.
**Written:** 2026-11-09, Claude.
**Size (R-2026-11-07-2, merged as delivered):** one section in `problems/O19.md` under about 120 lines, one test
file, one STATE log line. No change in `met/`.

Read O19's B19 section (why the matrix TLS conditions do not transfer to one scalar mass), T42's and T53's
statements and proofs, L11's row in `definition/premises.md`, and `definition/type.md`.

## Standard mathematics first

Minimax lower bounds from the second moment of a mixture likelihood ratio (Ingster's method; it is T41's and T53's,
and Mikusheva & Sun's Theorem 1, arXiv:2004.12445). For attainment: the leading eigenvector of a rank-one signal
plus Gaussian noise Gram matrix and its perturbation (Davis–Kahan; Yu, Wang & Samworth, as in T42).

## The sorted fork (Claude's, provisional; do not re-sort it)

L11 as restated on 7 Nov carries the literature's rate for scalar readings on one instrument,
$\Delta_N/\sqrt N\to\infty$. The type has readings with $k_i\ge1$ components in a declared subspace and its own
$\Sigma_i$. B19 found that the published matrix conditions are about an unrestricted matrix slope and say nothing
sharp for one scalar mass read through vectors. The fork: **what replaces $N$ in the rate.** Claimed, to be proved
or refuted: with one common isotropic noise ratio it is the number of noise dimensions in the candidate subspaces,
$D_N=\sum_ik_i$, because every component is a scalar reading of the same line; so the rate is
$\Delta_N/\sqrt{D_N}\to\infty$.

## Deliverable

1. **Common isotropic noise ratio, general $k_i$.** The lower bound for every estimator (T53's class and
   quantifiers, stated) and an estimator attaining it. Is the count $D_N$, $N$, or something else when the $k_i$
   differ or grow.
2. **Different noise ratios between readings.** The hypotheses are then not one line in a common plane. State the
   information quantity that replaces $\Delta_N$ (the efficient information for the mass summed over readings is the
   candidate) and the count that replaces $N$; prove the threshold or give a counterexample to the natural guess.
3. **Correlated channels and proper subspaces.** Say whether items 1–2 survive a general $\Sigma_i$ after whitening,
   and what a declared subspace $N_i$ with $k_i<d$ changes.
4. **Write** the result into O19 as a dated section, each statement with its class and quantifiers (one series or
   triangular arrays; the point or the law); tests for the finite identities used; one line in `docs/STATE.md`'s log.

## Do not

- Change L11's row, grade or wording, or propose the replacement as the premise.
- Declare a population of signals, or make the growth condition a caller input (TASK §5).
- Treat the null (O15, 9 Nov) or the quoted width (G28, B24): this brief is the rate only.

# B12: which temperings are admissible under L24 (and is η?)

**Status:** merged 2026-10-30: corrected Bernstein/proper-prior criterion; checker and 57-case report under both mean conventions; shared-noise-size obstruction for η.
**Item:** `definition/premises.md` L24 (candidate), `problems/O14.md`; theorems T36, T37. **Branch:** `codex/b12`.
**Written:** 2026-10-30, Claude.

Read `definition/theorems/T36.md` (homogeneity: tempering is a declared Σ in disguise), `T37.md` (the
per-reading noise size, the cut precision $g(Q)=E[1/v\mid Q]$), L24's row, and O14's three 30 Oct sections.

## The sorted fork (Claude's, provisional; do not re-sort it)

L24 (candidate): a tempering of a reading's evidence is admissible iff it is the cut posterior of a
noise-size law $P$ with $E_P[v]=1$. At strong signal $Q\mid v\sim v\chi^2_\nu$, $\nu=k-1$. Claude's sketch of
the test, **to be derived and checked, not assumed**: put $u=1/v$ and $d\mu(u)=u^{\nu/2}dP(u)$. Then
$g(Q)=-2L'(Q)/L(Q)$ with $L(Q)=\int e^{-Qu/2}d\mu(u)$, so

$$L(Q)=L(0)\exp\Big(-\tfrac12\int_0^Q g(q)\,dq\Big).$$

A law $P$ exists iff this $L$ is **completely monotone** (Bernstein's theorem). The mean-1 condition is
$\int u^{-1-\nu/2}d\mu\big/\int u^{-\nu/2}d\mu=1$, and $\int u^{-s}d\mu=\frac{2^{-s}}{\Gamma(s)}\int_0^\infty Q^{s-1}L(Q)\,dQ$
turns it into a condition on $L$'s Mellin moments (take care at $\nu=0$).

**Added 30 Oct, afternoon (G20, T37):** practice states σ so that the *precision* is right on average,
$E_P[1/v]=1$, not $E_P[v]=1$. Report both moments for every $P$ found; the verdict is per normalisation.

## Deliverable

1. **The derivation**, written into O14's attempts: the identity for $L$, the Bernstein criterion,
   the mean-1 condition in terms of $L$, with the $\nu=0$ case handled. If any step is wrong, say so.
2. **A checker** in a study module (for example `met/admissible.py`, imported by nothing in
   production): given $g$ as a function of $Q$ and $\nu$, report (a) whether $L$ is completely
   monotone, analytically where the form allows and otherwise by a stated numerical criterion with
   its limits; (b) $E_P[v]$ when $P$ exists (finite or not).
3. **Run it on**, for $k=1,2,3$, and report pass / fail / unknown for each, with the reason:
   - ours30's $g=\min(1,(k-\tfrac12)/Q)$ and ours29's $\min(1,(k-1)/Q)$ (the kink);
   - the smooth form $\nu/(\nu+Q)$ (`misfit_form="smooth"` in `met/full.py`);
   - the floor: $v\ge1$ with prior $v^{-a}$, its plug-in and its exact cut;
   - **controls that must pass:** T37's inverse-gamma family ($\alpha=1.5,3,10$), $g=(2\alpha+\nu)/(2\alpha-2+Q)$;
     and $g\equiv1$ (the fixed-size Gaussian: passes L24, fails T30).
   Claude's guesses, to be confirmed or refuted, not used: the kink fails Bernstein outright, because
   its $L$ is not smooth at the kink; the smooth form passes Bernstein but has $E[v]=\infty$.
4. **η.** State the series-level version of the test: η multiplies every reading's evidence, so it is
   the cut posterior of a noise size **shared by the series** (one draw for all readings), read from
   the readings' disagreement. Derive, in the linear-Gaussian limit of B04 (readings
   $x_i\sim N(\theta,1/c_i)$), whether any mean-1 law on a shared size has $\min(1,\mathcal H/J)$ as its
   cut posterior. Report what you find; a shared size is a correlation between readings, and say
   whether the type declares one.
5. Tests for 2–4 in `tests/test_admissible.py`; one line in `docs/STATE.md`'s log.

## Do not

- Change `met/full.py`, `ours30`, or any premise row. L24 stays a candidate.
- Propose a replacement for $g$ or η, pick T37's tail index α, or state any result as a ruling.
- Treat a numerical complete-monotonicity check as a proof: say which results are analytic.
- Touch the positive-mass misfit $Q_+$ or weak signal: strong signal, whole-line $Q$ only.

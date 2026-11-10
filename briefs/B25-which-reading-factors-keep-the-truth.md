# B25: which per-reading factors keep the truth as the population optimum, and what bounded confidence adds

**Status:** done: T72 proves a residual-increasing, sign-dependent Fisher-consistent counterexample, with deterministic quadrature and T30-failure checks; explicit counterexample stop before tail-family classification (Codex, 10 Nov 2026; not checked by Claude).
**Subsequent local continuation (10 Nov, authorized by Can):** T72 §4 and O18 record a family satisfying Fisher consistency and T30 together while retaining sign dependence. The fixed Gaussian near-line condition still fails; mean-one mixture representability remains unestablished. This extends the investigation after the original brief's stop; no factor adopted.
**Item:** G19 (the body of the per-reading factor), G22 (the sign in the weight against consistency);
T30, T36, T37, T39 (parts 1–3), T46, T48; `problems/O18.md` (4 Nov).
**Branch:** `codex/b25`, from `origin/main`.
**Written:** 2026-11-09, Claude.
**Size (R-2026-11-07-2, merged as delivered):** `definition/theorems/T72.md` under about 200 lines (T71 is B23's),
one test file, one paragraph in `problems/O18.md`, one STATE log line. No change in `met/`.

Read T36 (the frame), T46 and its test, T39 parts 1–3, T30, T37 (1), and G19's and G22's rows in
`definition/proof.md`.

## Standard mathematics first

M-estimation for orthogonal regression: Zamar (1989, *Biometrika* 76:149–160: orthogonal-regression M-estimates, a
consistency proof under equal error variances, robustness when the loss is bounded; 1992, *Ann. Statist.*
20:1875–1888); Cheng & Van Ness (1992, *Ann. Statist.* 20:385–397: Fisher consistency of generalized M-estimators in
the errors-in-variables model). Fisher consistency is the name of the property G22 is about. For the tail: scale
mixtures of normals (Andrews & Mallows 1974, *JRSS B* 36:99–102) and Bernstein's theorem, already used in B12.

## The sorted fork (Claude's, provisional; do not re-sort it)

On 4 Nov a candidate for what forces G19 was recorded: consistency (L11) together with bounded confidence (T30).
What is proved so far is one side. T46: among weights of the form $w(Q_+)$ on the quadratic $H_i$, every nonconstant
monotone $w$ moves the population optimum off the truth. T39: ours32's factor, which depends on the mass only
through the reading's residual $K_i(m)$ (criterion $-\sum_i(\nu_i+k_i-1)\log(\nu_i+K_i(m))/2$), keeps it there, and
its $Q_+$ cancels as a mass-free constant. Neither is a statement about the class. The fork: **is "a fixed function
of the reading's own residual $K_i(m)$" what consistency at every signal size forces, and within that class what
does T30 force.** If it is, G19 stops being a list of four candidates and becomes a family with a named freedom.

## Deliverable

1. **The class.** Per-reading factors $\exp(-\rho_i)$ allowed by the frame: a function of reading $i$ alone, invariant
   under its frame rotations, covariant under units and the swap. Show what is left: $\rho$ as a function of
   $K_i(m)$ and of mass-free statistics $s_i$ of the reading (its total energy $E_i$, the whole-line misfit $Q_i$, the
   positive misfit $Q_{+i}$).
2. **Fisher consistency.** Under the declared Gaussian noise, for every true mass and every fixed signal size, when
   is the truth a stationary point (and the optimum) of the population criterion $\sum_iE\rho(K_i(m);s_i)$. State
   and prove the necessary and sufficient condition on how $\rho$ may depend on $s_i$. Expected, to be confirmed or
   refuted: a fixed increasing $\rho(K)$ passes by the reflection symmetry of the residual about the true line;
   mass-free statistics that are blind to the sign of $F\cdot a$ may enter; one that sees the sign may enter only as
   a mass-free factor that cancels, which is T46's obstruction and T39's exception as one statement.
3. **Bounded confidence.** Within the class of item 2, which $\rho$ give T30's bounded, non-vanishing evidence
   contrast on a fixed pair no positive mass fits as its noise shrinks. Expected: a logarithmic tail,
   $\rho(K)\sim a\log K$ with $a>0$. Then say what is still free (the constant $a$, the shape away from the tail),
   and what two further conditions cut from it: agreement with the Gaussian factor to first order near the reading's
   own line (the declared σ's believed, R-2026-10-30-1), and being the cut posterior of a mean-1 noise-size law
   (L24's Bernstein test, B12). Is what remains the Student family with its index, or larger.
4. **Write** T72 with the statements proved, their scope, and the counterexamples; a test file checking the
   population stationarity identities by quadrature and the tail claims on T30's pairs; one paragraph in O18; one
   line in `docs/STATE.md`'s log.

## Do not

- Adopt a factor, name a probe, or change `met/`, any premise row or `rulings.md`.
- Re-sort G18 (the cap against the σ's believed) or treat the series level (η): this brief is the reading level only.
- Fix the tail constant or the index by a coverage study.
- If a Fisher-consistent factor outside the expected class exists, that is the result: write it into T72 and stop.

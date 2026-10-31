# B05: T23's consistency route under the tempered evidence (O4)

**Status:** merged 2026-10-31 (checked by Claude): both routes refuted as stated; ours30's L11 counterexample confirmed with production code (G22).
**Item:** `problems/O4.md`, `definition/theorems/T23.md`; `definition/proof.md` step 2, hole G1.
**Branch:** `codex/b05`. **Written:** 2026-10-29, Claude.

Read `definition/proof.md` §2.1 step 2, then `problems/O4.md` and T23.

## The sorted fork

None to sort: this is a proof owed. T23 claims the noise model is forced to be Gaussian by two
routes. The Z route holds under the flat weight (tested). The consistency route, that L11
($\widehat m\to m$ as the total true signal grows) forces the Gaussian among second-moment
families, is stated for the flat push and never written for the tempered law (misfit weight
$g_i$ on $Q_+$, power $\eta$).

## Deliverable

1. A derivation, into T23.md's "Proof or support": for one noise pair and a family of
   zero-mean error densities with the declared covariance (take the elliptical family, or any
   family the derivation can handle, stated), which members give a per-reading factor whose
   point converges to $m$ under L11 with ours30's $g$ and $\eta$ applied. Say exactly what is
   shown: the Gaussian is forced, or only a class containing it, or the route fails.
2. The Z route written out alongside (why only the Gaussian leaves an exact-zero reading's
   factor constant in $m$), since T23's file has only a test for it.
3. A test for whatever 1 establishes (a counterexample family that breaks L11, or the
   convergence under the Gaussian), in `tests/test_theorems.py` beside T23's test.
4. One line in `docs/STATE.md`'s log; O4's status updated.

## Do not

- Change `met/`, the equation, or any premise.
- Treat "the data need not be Gaussian" (R-2026-09-27-12) as contradicted: the claim is about
  the model, not the world.

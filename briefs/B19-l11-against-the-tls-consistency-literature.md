# B19: L11's growth condition against the errors-in-variables consistency literature

**Status:** merged 2026-11-07, as delivered (not checked by Claude): the literature's scalar rate is Δ_N/√N → ∞ (Gallo, Kukush–Van Huffel, Shklyar), equal to T42's and T53's threshold; T41's Δ_N ~ log N fails it; L11 unchanged (a rule, G25)
**Item:** `problems/O19.md`; L11 in `definition/premises.md`; G25; theorems T41, T42, T52, T53.
**Branch:** `codex/b19`, from the head of `codex/ours32-consistency` (O19 and T41–T53 live there).
**Written:** 2026-11-04, Claude.

Read O19, T41's statement, T42's and T53's thresholds, L11's row.

## The sorted fork (Claude's, provisional; do not re-sort it)

T41 proves literal L11 ("$\hat m\to m$ whenever the total true push diverges") is unattainable by any
estimator on the unrestricted type; T42 and T53 give $\Delta_N^2/(N+\Delta_N)\to\infty$ as the sharp threshold
for the flat Gaussian procedure and for a specified uniform Gaussian energy class. O19 asks for the
information-growth condition that a valid L11 needs. **This is a solved problem in standard terms and the
branch does not cite it.** The model is the functional errors-in-variables / total-least-squares model
with one slope; its consistency conditions are Gleser (1981, *Ann. Statist.* 9:24–44), Gallo (1982, *Comm.
Statist. A* 11:973–983), Kukush & Van Huffel (2004, *Metrika* 59:75–97: a condition of the form
$\lambda_{\min}^2(A_0^\top A_0)/\lambda_{\max}(A_0^\top A_0)\to\infty$ on the true-signal scatter), and Shklyar (2011,
*Theor. Probab. Math. Statist.* 83:175–190; 2018, *Mod. Stoch. Theory Appl.* 5:247–295), whose conditions are
weaker than Kukush–Van Huffel's. Standing rule 6: bring the mathematics in before working the fork.
The fork: **which of these conditions, stated in the type's terms, is T42/T53's threshold, and which is
weaker?** L11's restatement is a sort (Claude's) and a ruling (Can's); this brief supplies what the
literature already knows so the sort is not re-derived.

## Deliverable

1. **The conditions, in the type's terms.** For 1D readings with common unit noise (T42's setting), write
   each paper's sufficient condition for TLS consistency as a condition on $(N,\Delta_N,\{r_i\})$: Gleser's,
   Gallo's, Kukush–Van Huffel's eigenvalue condition, Shklyar's 2011 and 2018 conditions. For the
   one-slope model $A_0^\top A_0$ is a scalar, so say what each condition reduces to and whether it is
   $\Delta_N\to\infty$, $\Delta_N/\sqrt N\to\infty$ (which is T42's $\Delta^2/(N+\Delta)\to\infty$), or something else.
2. **Compare with T41, T42, T53.** State which literature condition T41's counterexample ($r_i=b/\sqrt i$:
   $\Delta_N\sim b^2\ln N$, $\Delta_N/\sqrt N\to0$) violates, and whether any of the published sufficient conditions
   is weaker than $\Delta_N/\sqrt N\to\infty$ (Shklyar's 2018 abstract says his conditions do not require
   Kukush–Van Huffel's eigenvalue ratio to diverge). If a published condition is both weaker than T42's
   and sufficient for the TLS point, T53's "sharp for the class" and that condition must be reconciled:
   say how (different classes, different estimators, or an error).
3. **Heteroscedastic and multi-dimensional readings.** Kukush & Van Huffel's elementwise-weighted TLS
   covers row-wise different covariances; say what their condition becomes with the type's $\Sigma_i$ and
   $k_i>1$, and whether the branch's T61 (changing covariance shapes) is the same statement.
4. **Write** the result into O19 as a section ("what the literature already knows"), with the conditions as
   a table, each with its citation and the page or theorem number; mark anything you could not verify
   against the paper as unverified. No code is required; a short check script for the one-slope reductions
   is welcome (`tools/`, seeded if random).
5. One line in `docs/STATE.md`'s log.

## Do not

- Change L11's row, grade or wording; do not propose the replacement condition as the premise.
- Make the growth condition a caller input (TASK §5: an unobserved excitation threshold is not an input).
- Re-derive what the papers prove; cite and translate.
- Extend the anchor/score-law constructions (T56–T65).

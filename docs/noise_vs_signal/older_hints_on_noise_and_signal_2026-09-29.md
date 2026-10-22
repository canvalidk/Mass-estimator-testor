# Older hints on noise and signal: where they already are

**2026-09-29 (night).** Record, written by Claude at the user's request
("Ok now do A and B lists"), list B. A pointer record: it quotes the one
sentence that matters from each earlier record that touched the
noise-versus-signal problem before it had a name, and says what that
sentence now reads as. The records stay where they are. Nothing here is
new mathematics.

**Relations.** Points at eleven records; binds nothing. Companion to the
four records written today in this folder, and to the records moved into
it the same night (list A, see the folder README).

---

| # | where | what it said | what it reads as now |
|---|---|---|---|
| 1 | `median_based_estimator_inertial_mass.md` §13 (27 Aug) | "although every mass possesses the exact latent solution $a^*=0$, the total probability associated with the observed zero-zero result varies with mass." And §12: "a $0/0$ reading does not have one universal interpretation." | The first statement that *preference* and *exclusion* come apart at $0/0$: the profile is flat, the integrated support is not. Under ours27/28 the data factor at exact $0/0$ is constant (T24), so the preference it found came from the flat latent measure of the time, not from the reading. |
| 2 | `mass_estimation_zero_zero_two_cases.md` Case I §8, Case II §9 (5 Sept) | Case I: "$\{m>0:Q(m)\leq q\}=(0,\infty)$ for any $q\geq0$ … Every positive mass passes." Case II: the five-line output — point, compatibility set, induced distribution, no finite SD, standing — and "If those intervals were reported as though the null trial itself had excluded a large mass, the Sun objection would remain unanswered." | Property N and the resolution as a third output, a month early. The five-line output is what the view record §5 makes a requirement. |
| 3 | `what_forces_the_ratio_of_means_readout.md` §7 (2 Sept) | "boundary fits … exactly the zero-reading cases … have fewer effective degrees of freedom, so a $\chi^2_{d-1}$ threshold is conservative there — which is the correct direction for a gate whose job is not to reject weak trials." | The residual $Q_{\min}=(T_i-\lvert C_i\rvert)/2$ is smaller under noise than under a push: the alignment record §3, seen from the gate's side. |
| 4 | same record §6 (2 Sept) | "For $N$ trials with a common mass the two flat measures differ by $m^{2N-2}$." | The first place the flat measure's choice of *what is uniform* compounded with $N$. The $\sqrt N$ tilt is the same kind of compounding, in the data rather than the measure. |
| 5 | `mass_estimator_originality_research_2026-09-05.md` (5 Sept) | "divergent latent-ratio moments are not confined to zero–zero. Under the default Gaussian/flat-magnitude policy, they occur at every finite empirical input." And: "Do not use divergent moments as the null detector." | The flat push in one sentence: every push infinitely large in noise units, so no input can be told apart as noise by its moments. The warning against a moment-based null detector is still right: the noise measure is $\lvert C\rvert$ against $2\sqrt N$, not a moment. |
| 6 | `ours26/README.md`, "Coverage: two different failures" (27 Sept) | "(b) Many readings: correctable … At SNR 1 the cell average went 95.8% → 91.1% → 86.3% for $N$ = 2, 6, 21, and the worst cell was 76.5% … the law is too narrow by $\sqrt{1+M/R^{\ast2}}$ (T18)." And: "The user's position: wait for the fix to drop out naturally rather than impose it." | The measured symptom. Failure (b) grows with $N$ as the null-series tilt does. Failure (a), Fieller/Gleser–Hwang at the weak-signal corner, is the $0/0$ limit and is accepted, not fixed. |
| 7 | `ours26/ours26_learned_push_scale_2026-09-27.md` (27 Sept) | The tempering identity: a finite push scale raises the likelihood to the power $\beta=\omega^2/(1+\omega^2)$; the peak does not move and the law widens by $1/\beta$. Scores: too wide at weak signal (99.8% at nominal 95%), the point pulled toward the prior when data are thin. | The mechanism any collective fix would use. Its failure mode (too wide at weak signal) is what a noise measure that tempers everything by one global $\omega$ does; the null-series record §5 asks for tempering by $\lvert C\rvert/2\sqrt N$ instead. Stays in `ours26/`; its ruling record moved here. |
| 8 | lemma record T18 | $\mathrm{Var}_{\rm true}/\mathrm{Var}_{\rm law}=1+\sum_c\theta_c'^2/\sum_c\theta_c'^2r^{\ast2}_c$; "the width error of L10, the asymmetric endpoint. Absent under ours28 to Monte Carlo precision." | The width half of the problem, closed by ours28 for signal-bearing series. The null-series record shows the preference half is not. |
| 9 | lemma record T14, T15 | "readout at zero readings is $s$"; "$P=Q\Rightarrow\widehat m=\text{median}=s$". | Noise units leaking into the answer at zero and equal signal. The $c$ record: this is a declaration ($c=s$), not something the noise supplies. |
| 10 | lemma record T26 | "no per-reading size rule has the right width at every push size; none is never overconfident (a pure-noise reading adds score noise and no curvature)." | The per-reading half of the dichotomy, stated on 28 Sept. The alignment record §6 and the null-series record §5 extend it: no per-reading rule removes the $\sqrt N$ either. |
| 11 | `ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` (28 Sept) | L19, his ruling: nothing inferred from one reading about pushes is used on another; readings may come from any regime in any combination ("F = 1, then F = 50000000, then F = 0.00001"). | The other horn of the dichotomy. The narrowest weakening (share one scalar, the series' signal energy against chance) is open; whether it survives his example is parked item 1 of the view record. Moved here with `ours28b/`. |

## What the list shows

Every piece of today's picture was in the corpus before today, spread over
five weeks: preference versus exclusion (1, 2), the resolution as an output
(2), the residual pointing the wrong way (3), the flat push as the root (5),
the symptom (6), the mechanism (7), the width half closed (8), the
per-reading limit (10), and the ruling that blocks the collective fix (11).
What was missing was the object that ties them: $\lvert C\rvert$ against
$2\sqrt N$, a function of the carried state.

## Sources

The eleven records in the table, at their current paths.

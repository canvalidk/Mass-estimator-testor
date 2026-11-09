# B26: is the starting law's family closed under lumping, and is that B18's weak-signal under-coverage

**Status:** open.
**Item:** `problems/O6.md` (the family, B08), G27 (the joint composition protocol), `problems/O2.md` (B18's section);
L4, L5, L14 in `definition/premises.md`; R-2026-10-28-1.
**Branch:** `codex/b26`, from `origin/main`.
**Written:** 2026-11-09, Claude.
**Size (R-2026-11-07-2, merged as delivered):** one section in `problems/O6.md` under about 80 lines, one table in
`problems/O2.md`, one test file, one script in `tools/` (seeded), one STATE log line. No change to the prior in
`met/full.py`.

Read O6's family derivation (B08), R-2026-10-28-1's three clauses, and B18's section in O2 with its three presets.

## Standard mathematics first

Closure of a scale family under convolution is stability: the strictly stable laws on $(0,\infty)$ have index
$\alpha\in(0,1)$ (Feller, *An Introduction to Probability Theory and Its Applications* II, ch. VI and XIII;
Zolotarev, *One-dimensional Stable Distributions*). The map $m\mapsto1/m$ exchanges the sum with the parallel sum, so
L14's operation is L5's on the reciprocal. Propagation of laws through a measurement function is JCGM 101:2008,
as B18 already cites.

## The sorted fork (Claude's, provisional; do not re-sort it)

R-2026-10-28-1 puts composition at the level of the law: the composite's law is computed from the parts' laws, and a
direct measurement of the composite must agree with it as coverage. With no information the parts' laws are their
starting laws, so the computed composite is the convolution $\pi_{c_1}*\pi_{c_2}$ (L5), or its parallel-sum analogue
(L14), while the direct measurement starts from a member $\pi_c$ of the family $m^{-1}g(\ln(m/c))$. If no member of
the family is closed under these two operations, the two laws differ before any reading is taken, and at weak
signal the difference is still there. B18 found weak centred sums under-cover on every probe. The fork: **is that
under-coverage the starting law's, not the probes'.**

## Deliverable

1. **Closure.** For the family of O6 (L3, L4, properness): is there a member with $\pi_{c_1}*\pi_{c_2}$ in the family
   for all $c_1,c_2$, and likewise for the parallel sum; and the weaker question for equal centres only. Expected, to
   be confirmed or refuted: none. Closure for all centres makes the member strictly stable on the half-line, whose
   upper tail is a power law while its reciprocal has every moment, against L4. Prove it or give the member.
2. **The half-Cauchy case, in numbers.** The median and central 95% interval of $\pi_c*\pi_c$ and of the parallel
   sum, against $\pi_c$'s and against $\pi_{2c}$'s and $\pi_{c/2}$'s (analytic or by quadrature).
3. **Against B18.** For the weak cells of the three B18 presets (true acceleration 0.1), report beside each recorded
   coverage the coverage the starting laws alone give at zero signal for the same truths and centres. Then, by
   reweighting each saved or recomputed composite law on its grid by the ratio of starting laws (no change to
   `met/full.py`'s prior), the coverage when the direct measurement starts from the convolution instead of $\pi_c$.
   Say how much of each weak-cell shortfall the starting law accounts for.
4. **Write** items 1–2 into O6 as a dated section and item 3's table into O2 under B18's section; tests for the
   closure identities and the quadrature; one line in `docs/STATE.md`'s log.

## Do not

- Propose or adopt a starting law, a default centre for a composite, a joint prior or a dependence model.
- Change B18's worlds, seeds or protocol beyond the reweighting of item 3.
- Change any premise row or `rulings.md`. What a non-closure means for L4 or for π as an input is Claude's sort.

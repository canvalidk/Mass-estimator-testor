# B24: the quoted width and η against the weak-instrument literature

**Status:** done: corrected canonical mapping and threshold quantifiers; AR/LM/CLR premise and exact null-shape audit, with T45 identified; Bekker/HHN/MS variance targets distinguished in O14/O15 and deterministic checks supplied (Codex, 10 Nov 2026; unverified original statements marked).
**Item:** G28 (what the quoted width promises on noise), G3 and G26 (η's finite-$N$ form and its variance target);
`problems/O15.md` (9 Nov), `problems/O14.md` (4 Nov: η's target); T42, T45, T53; B21 uses the same model.
**Branch:** `codex/b24`, from `origin/main`.
**Written:** 2026-11-09, Claude.
**Size (R-2026-11-07-2, merged as delivered):** one section in `problems/O15.md` and one in `problems/O14.md`, each
under about 100 lines, a table in each; an optional check script in `tools/`; one STATE log line. No change in `met/`.

Read O15's 30 Oct, B13 and 9 Nov sections, T45, T42's and T53's statements, O14's 4 Nov section, B19's section in
`problems/O19.md` (this brief is its counterpart for the width and for η), and B21.

## Standard mathematics first

Linear instrumental variables with one endogenous regressor, in its canonical form with known reduced-form
covariance (Moreira 2003, *Econometrica* 71:1027–1048; 2009, *J. Econometrics* 152:131–140); many weak instruments
(Bekker 1994, *Econometrica* 62:657–681; Chao & Swanson 2005, *Econometrica* 73:1673–1692; Hansen, Hausman & Newey
2008, *JBES* 26:398–422; Mikusheva & Sun 2022, *Rev. Econ. Stud.* 89:2663–2686, arXiv:2004.12445);
identification-robust confidence sets (Anderson & Rubin 1949; Kleibergen 2002, *Econometrica* 70:1781–1803;
Moreira's conditional likelihood ratio); and why an almost surely bounded set cannot have positive confidence level
when identification can fail (Dufour 1997, *Econometrica* 65:1365–1387; Gleser & Hwang 1987, already in B13).
The volume and page numbers here are from Claude's memory; correct them when citing. B21's Moreira (2009) is a
different paper (*Ann. Statist.* 37, the incidental-parameter likelihood).

## The sorted fork (Claude's, provisional; do not re-sort it)

**The claimed mapping (Claude's, from memory; item 1 checks it).** A 1D reading with known noise and its own unknown
push is one coordinate of the canonical model: reading $i$ ↔ instrument $i$ (so $K=N$), the push $r_i$ ↔ the
first-stage coefficient, $\Delta_N$ ↔ the concentration parameter, the mass ↔ the structural coefficient, the declared
$\Sigma$ ↔ the known reduced-form covariance. The mapping is to the canonical sufficient statistics, not to an
$N$-observation regression with one dummy instrument per observation: there the projection matrix is the identity and
every leave-one-out statistic is empty.

Under it, two things main lacks are standard. The threshold: LIML is consistent when $\sqrt{K_n}/r_n\to0$ (Chao &
Swanson), and no test, hence no estimator, is consistent while the concentration parameter over $\sqrt K$ stays
bounded, even with known covariance (Mikusheva & Sun, Theorem 1); that is T42's and T53's rate. The width: that field
met Gleser–Hwang's wall (Dufour) and answered it by quoting sets obtained by inverting tests, which may be unbounded,
instead of always-finite intervals. T45 on main is such a set and was set aside because step 18 defines the width as
an interval of the law; G28 reopens exactly that. And η's job, the variance under many weak instruments, has named
estimators there.

The fork: **which of that field's sets and variance corrections, stated in the type's terms, meet the premises a
width or an η must meet, and which do not.** The sort of G28 is Claude's; this brief supplies what the literature
already knows.

## Deliverable

1. **The mapping, verified or corrected.** Each object of the type's 1D slice against the canonical model, with the
   place it is defined in Moreira (2003) or Mikusheva & Sun. Say where it stops: readings with $k_i>1$, different
   $\Sigma_i$ per reading, the restriction $m>0$, and that the type declares no order.
2. **The threshold.** Chao & Swanson's condition and Mikusheva & Sun's Theorem 1 in the type's terms, beside T42 and
   T53: the same class and quantifiers, or not (as B19's item 2 did for the TLS papers).
3. **The sets.** The Anderson–Rubin set, Kleibergen's score set and the conditional likelihood ratio set, written
   for known $\Sigma$ in the type's terms. For each, passes / fails / unknown, with the argument or a counterexample:
   an exact-zero reading leaves it alone (Z); the order of readings does not matter; units; the swap of force and
   acceleration; restriction to $m>0$; nothing used across readings but the mass (no borrowing, no size law). When
   is it the whole half-line, empty or disconnected, and with what probability on pure noise. Is T45's set one of
   them. What each loses as $N$ grows at weak signal (the many-instrument power loss of Anderson–Rubin).
4. **η.** The Bartlett η and its target $\Delta/(\Delta+M)$ (O14, 4 Nov) beside Bekker's many-instrument variance and
   the variance estimators of Hansen–Hausman–Newey and Mikusheva–Sun: is the correction the same quantity; what the
   estimators become when the covariance is known; whether any is a function of the readings that is Z-safe and uses
   no count of readings.
5. **Write** items 1–3 into `problems/O15.md` and item 4 into `problems/O14.md` as dated sections, each claim with
   its citation and theorem or page number; mark anything not verified against the paper as unverified. One line in
   `docs/STATE.md`'s log.

## Do not

- Adopt a set as the quoted width, change the output type, or change any premise row or `rulings.md`.
- Read an unbounded set as "undefined": the law and the point are still returned on every measured input (L2).
- Re-derive what the papers prove; cite and translate.
- Extend the mapping to $k_i>1$ beyond saying where it stops (B27 is the rate there).
- Run a coverage study. If one is needed to settle an item, say which cells and stop.

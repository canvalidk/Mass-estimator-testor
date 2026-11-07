# B22: ours30 with η from the series' Gram (η_G) in place of the Bartlett η

**Status:** open.
**Item:** G3 (η's finite-$N$ form); `problems/O14.md` (7 Nov: η from the series' real density matrix); T70 (7).
**Branch:** `codex/b22`, from `origin/main`.
**Written:** 2026-11-07, Claude, from a conversation with Can (misalignment as imaginarity; T70).
**Size (R-2026-11-07-2, merged as delivered):** keep to what is asked: one option in `met/full.py`, one test
file, one script, one section in `problems/O14.md` under about 80 lines, one STATE log line.

Read T70, O14's 4 Nov and 7 Nov sections, `met.full.gram_eta`, `met.full.bartlett_eta`, and T30, L11, L20 (Z) in
`definition/premises.md`.

## Standard mathematics first

$\eta_G=1-Q/H(\hat\phi)$ is the empirical-Bayes shrinkage of the structural errors-in-variables model
(Lindley & El-Sayyad 1968; James–Stein form $1-M/H^2$, O14 4 Nov) with the declared latent count $M$ replaced by
the series' whole-line misfit $Q$, which is $\lambda_{\min}$ of a noncentral 2×2 Wishart on one instrument
(Anderson 1946). Its small-$N$ bias is the extreme-eigenvalue bias of a small Wishart; the 2×2 case has closed
forms (Muirhead 1982, *Aspects of Multivariate Statistical Theory*): cite rather than re-derive if
item 3 needs them.

## The sorted fork (Claude's, provisional; do not re-sort it)

T70 (6): g stays per reading (T30 with Z need a per-reading null count); T70 (7): the series-level tempering is a
function of the series' Gram. So the probe keeps ours30's per-reading $g_i$ and replaces only η, computed on the
$g$-weighted series: $\eta_G=$ `gram_eta(readings, c, g)[0]`. Nothing else changes: π, the misfit weight, the
median point. It is a probe: presumed to break; what matters is which premise the break names.

## Deliverable

1. **The option.** `ours30(..., eta_form="bartlett")` (default, unchanged) and `eta_form="gram"`, which uses
   `gram_eta` on the $g$-weighted series in place of `bartlett_eta`; the result dict reports which. Every existing
   test and golden output must be unchanged at the default.
2. **Checked against the premises**, passes / fails / unknown each, in `tests/test_eta_gram.py` (seeded): Z (an
   exact-zero reading leaves law and point alone), L3 and the swap (as `tests/test_short_law_proofs.py` does for
   ours30), frame rotation, T30 (the pairs of `tests/test_sign_fork.py` and
   `tests/test_ours31.py::test_t30_bounded_as_the_declared_noise_shrinks`), and L11 at fixed per-reading signal
   (`tests/test_weight_consistency.py`'s cells: does the population optimum move from ours30's?).
3. **Coverage beside ours30**, one script `tools/ours30_eta_gram_coverage.py` (seeded, base seed 2211): the cells of
   `tools/ours30_coverage_vs_n.py` (3D, signal 1 and 2, $N$ = 5, 20, 80, 200) and the gaussian-world cells of
   `tools/ours32_coverage.py` (k = 1, 3; signal 0, 0.5, 1, 2, 5; N = 5, 20), at least 150 series a cell. Report
   95% coverage of the truth and median width for ours30 (Bartlett η), ours30 with $\eta_G$, and ours30 with
   η = 1. Name the cells where $\eta_G$ under- or over-covers by more than 3 points.
4. **Write** the result into `problems/O14.md` as a dated section: the premise table from 2, the coverage table
   from 3, and which premise (if any) the probe breaks on. One line in `docs/STATE.md`'s log.

## Do not

- Change the default η, π, g, the point, any premise row, or `rulings.md`. The option is opt-in and gets no
  `oursNN` name; if it survives, naming it is Claude's call.
- Remove $g$ or apply $\eta_G$ per reading (T70 (6) sorts that).
- Add a cap, a floor or a bias correction to $\eta_G$ in the option. If item 3's numbers suggest one, say so in
  O14 with the cell and stop.
- Run it on ours32 or ours33: their reading level is a cut factor, not $g H$, and what $\eta_G$ is there is not
  sorted.
- Re-sort G18, G19 or G22.

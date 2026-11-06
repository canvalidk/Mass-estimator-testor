# B21: Moreira's maximal invariant likelihood for the latent pushes

**Status:** open.
**Item:** `problems/O15.md` (the 30 Oct table of non-hierarchical cures, closed by B14 for the per-reading
constructions); G3 (η's finite-$N$ form); `problems/O14.md` (4 Nov: η's target is ours28b's shrinkage $\Delta/(\Delta+M)$).
**Branch:** `codex/b21`, from `origin/main` (O15's table and B14's row live on `main`; nothing on
`codex/ours32-consistency` is needed).
**Written:** 2026-11-06, Claude, from a conversation with Can (the covector recasting of the mass corner).

Read O15's 30 Oct (later) section and its table, B14's section in O15, O13 §1 of the 3 Nov positivity section
(the Gram tilt $e^{-L^2/2}\psi_k(L\sqrt{v^\top Gv})$), O14's 4 Nov section, T36, and L9, L10, L16, L19, L20 in
`definition/premises.md`.

## Standard mathematics first

Moreira (2009), *A maximum likelihood method for the incidental parameter problem*, Ann. Statist. 37(6A):3660–3696,
doi:10.1214/09-AOS688 (open: arXiv:0909.0613). Seek a group acting on the sample space that preserves the model
and the structural parameter, so that its action on the parameters is $(m,\eta)\mapsto(m,A_3(g,\eta))$ (his
Definition 2.1); the maximal invariant statistic's law depends on the parameters only through the parameter
space's maximal invariant; when that has **fixed dimension**, the incidental-parameter problem is gone and the
invariant likelihood is an ordinary likelihood. His IV example is the case that matters here: the invariant is a
noncentral Wishart, the maximum invariant likelihood estimator is LIML, and it is consistent under many weak
instruments. LIML and the functional errors-in-variables maximum likelihood (TLS with known error covariance) are
the same estimator (Anderson 1976, JRSS B 38:1–36); the noncentral Wishart is Anderson (1946), Ann. Math. Statist.
17:409–431. So the mass corner on one instrument should be his case almost verbatim.

O15's table checked every per-reading cure (Cox–Reid, flat integrated, conditional, marginal, Barndorff-Nielsen).
Moreira's is the one cure in the family not in the table, and it is the only one that **acts across readings**.

## The sorted fork (Claude's, provisional; do not re-sort it)

**Per reading:** the largest group of affine isometries of one reading's whitened space $\mathbb R^{2k}$ that maps
the model $\{N(U_m t, I_{2k})\}$ into itself and fixes $m$ is the frame group $(F,a)\mapsto(RF,Ra)$, $R\in O(k)$
(L16): translations along the candidate plane depend on $m$, so they are not an action on the sample space, and
nothing rescales the push without rescaling the declared noise. $O(k)$ is not transitive on the latent $\alpha_i$:
its orbits are spheres $\lvert\alpha_i\rvert=r_i$, so the push size stays incidental and the per-reading construction
brings back O15's question. Nothing new there; that is why the table closed.

**Across readings, on the one-instrument slice:** if every reading in the window has the same declared isotropic
covariance $(\sigma_F^2I_k,\sigma_a^2I_k)$, stack the series' $Nk$ components as rows of a whitened $Nk\times2$ matrix
$Z$, row $j=(F_j/\sigma_F,\ a_j/\sigma_a)$. Then $Z=\alpha v_m^\top+E$ with $v_m\propto(m/c,1)$ a unit vector,
$\alpha\in\mathbb R^{Nk}$ incidental, $E$ iid $N(0,I_2)$ rows. **$O(Nk)$ acting on the rows** preserves the model and
fixes $m$; the sample maximal invariant is the $2\times2$ Gram $G=Z^\top Z$; the parameter maximal invariant is
$(m,\Delta)$ with $\Delta=\lvert\alpha\rvert^2$, the series' total true push in noise units. **$Nk$ incidental parameters
collapse to one.** The law of $G$ is noncentral Wishart $W_2(Nk, I, \Delta v_mv_m^\top)$; by O13 §1's argument with
$k\to Nk$, its $m$-dependence is the tilt $e^{-\Delta/2}\psi_{Nk}\big(\sqrt{\Delta\,v_m^\top Gv_m}\big)$. Since $\psi$ is
increasing, the maximiser in $m$ is the top eigenvector of $G$ restricted to $m>0$ for every $\Delta$: the point is
positive-restricted TLS, which is consistent. Off the slice (different $\Sigma_i$, or anisotropic channels) whitening
gives each reading its own line direction and the row group is gone; the per-reading result above is all there is.

**Why it matters (Claude's guess, to be checked, not a claim):** the flat law is $\exp(v_m^\top Gv_m/2)$ on the
slice. The invariant likelihood replaces it by a function of the same quadratic whose curvature is tempered by
an amount that depends on $\Delta$ against the $Nk$ noise dimensions. If that tempering tends to η's target
$\Delta/(\Delta+M)$ (O14, 4 Nov), then on the slice **η's job is forced by invariance**, with no declared size law and
no quasi-likelihood argument, and ours28b's shrinkage is recovered as a theorem on one instrument. If it does not,
the row reads "a different law, checked against the premises", and that is a result too.

**What the sort leaves to Can, and the brief only reports:** whether a likelihood that depends on the readings
only through the summed Gram, with one shared nuisance $\Delta$, is borrowing under L19 (ours28b was struck for a
shared size *scale*; $\Delta$ is a parameter, not a declared population, but the brief does not decide this); and
whether treating $\Delta$ by integration is a size law under L10 (profiling is not).

## Deliverable

1. **The groups**, into a theorem file `definition/theorems/T66.md` (T66 is the next number free on both `main`
   and `codex/ours32-consistency`; T39–T65 live on that branch), with its test named:
   (a) per reading, the largest group of affine isometries of $\mathbb R^{2k}$ preserving the model and fixing $m$ is
   $\{R\oplus R: R\in O(k)\}$; prove it or give the larger group. Say whether non-affine model-preserving maps can do
   more (for a Gaussian location family with fixed known covariance, argue or cite).
   (b) on the slice, $O(Nk)$ on the rows preserves the model and fixes $m$; the maximal invariants are $G$ and
   $(m,\Delta)$; the slice's exact conditions (same $\Sigma$, isotropic, same $k$, one window) stated. Say what the
   largest group is when the slice holds only within classes of readings (several instruments): a product of
   $O(N_gk)$, one $\Delta_g$ per class. O14's top note ("one $\eta_g$ per noise pair") is the comparison.
2. **The invariant law in closed form** on the slice: the noncentral Wishart density of $G$ as a function of
   $(m,\Delta)$ (the $_0F_1$ / $\psi_{Nk}$ form), the MILE, and two laws of mass from it, each times π:
   (i) $\Delta$ profiled; (ii) $\Delta$ integrated against one stated measure (say which and why; flag that this is
   the L10 question). Compare each, as a function of $m$, with the flat law $\exp(v_m^\top Gv_m/2)$ and with
   ours30 on the same series.
3. **The tempering, derived:** the curvature of $\log$ (i) at its maximum divided by the flat law's, as $N\to\infty$
   at fixed per-reading signal and at fixed $N$ as signal grows. Is it η's target $\Delta/(\Delta+M)$, and what is $M$
   here? Exact or asymptotic; say which.
4. **Checked against the premises**, passes / fails / unknown each, in `tests/test_maximal_invariant.py` (seeded):
   L2 (every measured input on the slice, all-zero included), L3, L4 (the swap with $c=\sigma_F/\sigma_a$), Z (L20:
   add an all-zero measured reading with the same $\Sigma$: $G$ is unchanged but $Nk$ grows by $k$; does the law
   move? O14 found $\hat\beta$ fails Z this way), T30 (bounded confidence on a fixed misaligned pair as σ → 0; the sign enters through
   $G_{12}=\sum F\cdot a$; the pairs of `tests/test_ours31.py::test_t30_bounded_as_the_declared_noise_shrinks` and
   `tests/test_sign_fork.py`'s anti-parallel pair), L11 at fixed per-reading signal as $N$ grows
   (`tests/test_weight_consistency.py`'s cells: truth 2, signal 1 and 2, 1D and 3D: where is the population optimum?),
   and property N: coverage against the truth on O15's cells (3D; $N=5,20,80,200$; signal 0.5, 1, 2; 200 series a cell,
   seeded), beside ours30 and the flat law.
5. **Write** the result into `problems/O15.md`: a new row in the 30 Oct table ("Moreira (2009), maximal invariant
   likelihood") and a dated section with the derivation's summary, the premise table from 4, and what it means for
   G3. One line in `docs/STATE.md`'s log.

## Do not

- Change `met/`, η, π, any premise row, or `rulings.md`. A law found here is reported, not adopted, and gets no
  `oursNN` name.
- Decide L19 or L10 for the shared $\Delta$: report what the law depends on, and stop.
- Extend the row group off the slice by whitening readings with different $\Sigma_i$ into one matrix (the line
  directions differ; the group does not exist there). If the slice is all the construction reaches, that is the result.
- Put a size law on the per-reading pushes (L10), or estimate anything about one reading's push from another (L19).
- Count the swap or units as group elements: they move $m$, so they are equivariance checks (item 4), not part of $G$.
- Re-sort G22, B17's latent-measure exponent, or the O15 per-reading rows.

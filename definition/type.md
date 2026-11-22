# The type

What the definition of mass takes and returns. From lemma list d4 §1 (10 Oct 2026). Edit in
place; git holds the history.

## Input

A series of readings $i=1,\dots,N$, each declaring:

- data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$, any values
  including zero, in the reading's own frame;
- an error covariance $\Sigma_i\in\mathbb R^{2d\times2d}$, positive definite, any correlation,
  including between the channels;
- how well that covariance is known: the effective degrees of freedom $\nu_{{\rm eff},i}\in(0,\infty]$ of the
  stated uncertainty, as a calibration certificate gives it (GUM, JCGM 100:2008, Annex G); $\infty$ declares it
  exact, and the caller owns what follows (added 17 Nov, R-2026-11-17-1; the form is G20's sort, provisional).
  ours31–ours33 take it; ours30 does not. **17 Nov, afternoon (R-2026-11-17-3):** this is not an input of the
  definition. It is the input of the case study in which the discount's parameter is handed over; "the caller
  owns what follows" is withdrawn;
- a candidate subspace $N_i\in\mathbb R^{d\times k_i}$ with orthonormal columns (default $I_d$:
  nothing known about the push's direction);
- membership: which components belong to which reading (the reading structure);
- the reading's **window** $\Delta t_i$, over which both its force and its acceleration are
  averaged: one window per reading (T33; added 10 Oct). Part of the reading structure (L18).

Also declared for the series:

- the starting law $\pi$: its shape and its centre $c$.

No size law is declared (L10, 30 Sept): the definition says nothing about how large a
signal is in noise units, per reading or per series.

**The ">0" is part of the type** (L2, ruled 25 Sept). We cannot promise that the
acceleration is nonzero; we can promise that each instrument has noise. That promise is what
keeps the function from ever returning "undefined" on measured input (ruled 29 Sept).

## Output

A law of mass on $(0,\infty)$ and a point $\widehat m\in(0,\infty)$; since 30 Sept also
$\eta$ and the raw law's interval beside the quoted one (ruling 3 of 30 Sept). The misfit is
kept with the result (an output clause, not a premise the mathematics uses). The output is
never "undefined".

**22 Nov (T80; R-2026-11-22-1, R-2026-11-22-2; provisional).** The paragraph above is the output of the
law-and-point line, and it stands as that line's. It is a true statement relative to a declared
starting law. Read as a true statement for every true mass, the quoted uncertainty has to be a set of
masses that can be unbounded (T80), and that is the kind of output the E13 direction returns: the
evidence curve (T73), sets of candidate masses with a kept promise (T75), the leftover beside them. Where
every positive mass is kept, no point is owed. Still never "undefined": the same kind of answer on
every measured input. Which set is open (`../problems/O20.md`).

## Not this function

Exact inputs ($\sigma=0$) go to the exact branch (T1), which has three outcomes: one mass,
no mass fits, every mass fits. "Undefined" survives only there, at exact 0/0.

## History

- 2026-11-22: the kind of output T80 forces added beside the law-and-point output, which is kept (R-2026-11-22-2).
- 2026-11-17: how well each covariance is known added to the input (R-2026-11-17-1).
- 2026-11-17, afternoon: that input is the case study's, not the definition's (R-2026-11-17-3).
- 2026-10-24: created from lemma list d4 §1.

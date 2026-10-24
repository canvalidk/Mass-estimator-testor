# The type

What the definition of mass takes and returns. From lemma list d4 §1 (10 Oct 2026). Edit in
place; git holds the history.

## Input

A series of readings $i=1,\dots,N$, each declaring:

- data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$, any values
  including zero, in the reading's own frame;
- an error covariance $\Sigma_i\in\mathbb R^{2d\times2d}$, positive definite, any correlation,
  including between the channels;
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

## Not this function

Exact inputs ($\sigma=0$) go to the exact branch (T1), which has three outcomes: one mass,
no mass fits, every mass fits. "Undefined" survives only there, at exact 0/0.

## History

- 2026-10-24: created from lemma list d4 §1.

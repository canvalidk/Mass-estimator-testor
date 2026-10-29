# B06: write out T16, T17, T18, T19, T24, T25 from the law

**Status:** open
**Item:** `definition/theorems/T16.md, T17.md, T18.md, T19.md, T24.md, T25.md`; `definition/proof.md`
hole G15. **Branch:** `codex/b06`. **Written:** 2026-10-29, Claude.

**Frame (added 29 Oct):** `definition/theorems/T36.md` states the law's geometry ($H_i$ the projection onto the
mass subspace in the noise metric, $E_i-H_i$ the residual). Use it where it shortens a proof; do not change T36.

## The sorted fork

None: these are theorems about ours30's law whose files say "carried from the records; not yet
written out here". Each is algebra on the law of `definition/equation.md` and each already has a
test. The deliverable is the proofs, in the files.

## Deliverable

For each file, replace "Carried from the records named in the status line; not yet written out
here" with the derivation, from the law as stated in `equation.md` (ours30: $g_i$ on $Q_+$, $\eta$):

- **T19**: $E_p[m]=\infty$ for every series (the half-Cauchy tail in $\ln m$).
- **T17**: $\widehat m\to m$ as the total true signal grows, with $g_i\to1$ and $\eta\to1$ shown
  on the way (or state the conditions under which they do).
- **T16**: one reading with $A=B$ and $C=2\lambda$ gives the raw law a tilt $\lambda\sin2\theta$ at any
  centre.
- **T18**: the flat factor's width ratio $1+\sum_c\theta_c'^2/\sum_c\theta_c'^2r_c^{\ast2}$ for one noise
  pair, and that $1/\eta$ estimates it in the long-series limit.
- **T24**: an exact-zero reading has $H_i\equiv0$ for any $\Sigma_i$ and $k_i$ (two lines).
- **T25**: measured zeros on some axes leave $H_i$ as without those axes; a $d$-dimensional
  reading declared $k=1$ is the 1D reading.

Keep each file's Statement unchanged; if the derivation shows the statement needs a condition
the file does not state, add the condition in the proof and say so in the file's History and in
one line of `docs/STATE.md`'s log. The existing tests must keep passing; add one where the
derivation gives a sharper check.

## Do not

- Change `met/`, `equation.md`, `premises.md`, or any statement.
- Write proofs for T6, T30, T31: those are hole G2, G13 and G3 in `proof.md`, Claude's.

# B23: What makes mass real: isotropy in 3D, parity in 2D

**Status:** open.
**Item:** G10 (L22 as a clause of L1), `problems/O17.md`; T70 (3) (the sign is not a phase; imaginarity is what a
non-scalar inertia would explain); the phase-two sort (T70, "Phase-two sort").
**Branch:** `codex/b23`, from `origin/main`.
**Written:** 2026-11-07, Claude, from a conversation with Can (misalignment as imaginarity).
**Size (intake, 7 Nov):** merged as delivered: `definition/theorems/T71.md` under about 150 lines, one test file,
one paragraph in `problems/O17.md`, one STATE log line.

Read T70, L16 and L22 in `definition/premises.md`, T36 (2), and O17's "three metrics".

## Standard mathematics first

Schur's lemma for real representations. The linear maps $\mathbb R^k\to\mathbb R^k$ commuting with a group $G$
of frame changes form the commutant of $G$'s representation. For $k=3$ and $G=SO(3)$ the representation is
irreducible of real type, so the commutant is $\mathbb R$ (real scalars). For $k=2$ and $G=SO(2)$ it is
$\{aI+bJ\}\cong\mathbb C$ (rotation-scalings); adding one reflection ($G=O(2)$) cuts it to $\mathbb R$, because a
reflection conjugates $J\mapsto -J$. For $k=1$ every map is a real scalar. (Numerically checked by Claude, 7 Nov:
commutant dimensions 2, 1, 1.) In the 2D qubit picture of T70 a complex scalar $m\,e^{i\theta}$ is exactly a
rotation-scaling, so *every* 2D pair $(F,a)$, $a\ne0$, is explained by one complex $m$: the whole-line misfit over
$\mathbb C$ is zero.

## The sorted fork (Claude's, provisional; do not re-sort it)

L22 says mass is a scalar; L16 says frames, rotations and signs are bookkeeping. Read as "inertia commutes with
every frame change", these give one condition in 3D (rotations alone force a real scalar) and two in 2D (rotations
leave a phase; the sign convention, i.e. parity, removes it). So **in 2D the reality of mass is parity, not
isotropy**, and the codirectionality conditions of L8 ($k-1$ of them) are, in 2D, the one condition that
parity supplies. Sort: L22 is Newton's bin (isotropy of inertia) *and* L16's "signs" is doing work in 2D that
L22 alone does not do. The deliverable is the theorem with its tests, not a restatement of L22 (that is
Claude's, after).

## Deliverable

1. **T71**, `definition/theorems/T71.md`, with its tests named: (a) the commutant of $SO(k)$ and of $O(k)$ on
   $\mathbb R^k$ for $k=1,2,3$ (cite Schur; state the real-type fact for $SO(3)$); (b) in 2D the set of inertias
   commuting with $SO(2)$ is $\mathbb C$ and the misfit of any pair over it is zero; with $O(2)$ it is $\mathbb R$
   and the misfit over it is T70's $Q$ (distance to the real states); (c) in 3D the $SO(3)$-commutant is already
   $\mathbb R$, and the maps commuting with rotations *about the push axis only* (the stabiliser $SO(2)\subset SO(3)$)
   are the 2D case plus a free scalar along the axis: say what hypothesis set that is and what its misfit would be,
   since that is the set a single 3D reading can actually test (its $k-1=2$ conditions).
2. **Tests**, `tests/test_what_makes_mass_real.py` (seeded): the commutant dimensions by linear algebra (as in the
   7 Nov check); zero complex misfit on random 2D pairs; the real misfit equals `reading_state(...)["Q"]`;
   in 3D, random pairs have positive misfit over real scalars and zero over the stabiliser set of (c).
3. **One paragraph in `problems/O17.md`**: what this says about O17's anisotropic-inertia question (a symmetric
   positive-definite $g$ is the $O(k)$-breaking case; the complex case is the $SO(2)$-but-not-$O(2)$ case; these
   are different alternative mechanics), and which bin each premise goes in. One STATE log line.

## Do not

- Restate L22 or L16 in `premises.md`; report, and Claude restates.
- Build an equation or a law for the complex-mass alternative (that is exploration E7).
- Touch T70, η, g, or any ruling.

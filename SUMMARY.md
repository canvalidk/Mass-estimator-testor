# E13: Unnormalized Evidence ("Bell with σ = ∞")

Branch: `explore/e13-unnormalized-evidence` (from `main` `d2a3352`, 21 Nov 2026). Exploration only; nothing in `definition/`, `problems/`, `rulings.md`, or `met/` changed.

**Core idea:** Drop area-1 normalization of the output law. Let `0/0` return the flat compatibility curve $K(m) \equiv 0$ ($L(m) \equiv 1$, literally a bell curve with $\sigma = \infty$) and a prior-free 95% compatibility set allowed to return $(0, \infty)$.

## Five numbered insights (for Claude to harvest)

1. **`0/0` forces the identity of combination ($L \equiv 1$), not a proper prior $\pi_c$** (`E5`, `B16`, `O6`, `O11`, `O2`, `G28`): area-1 normalization on `T36`'s infinite-tail likelihood is what created `O6`, `O11`, `G6`, and the Gleser–Hwang null floor. Dropping it lets `0/0` combine cleanly with `5 ± 0.5 kg` with zero drift and zero $\pi$.
2. **Exact Gaussian `L13` + believed $\sigma$ (`L18`) vs. discordance widening (`R-2026-11-17-3` / `T30`) is a Kac–Bernstein contradiction** (`O14`, `G13`, `G18`, `G19`, `G22`): splits cleanly into **Door A** (exact Gaussian evidence + `T45`/CLR compatibility set + separate leftover surprise $K_{\min}$) vs. **Door B** (declared heavy-tailed error model, evaluated unnormalized and without $\eta$).
3. **`Z` as `y = 0` at finite $\sigma$ vs. $\sigma \to \infty$, and `T45` vs. Moreira CLR** (`L20`, `O7`, `T45`, `T69`, `B24`): both `T45` and Moreira CLR hit ~95% coverage on pure noise ($m = 20$), weak-signal 3D $N=20$, `T43` leverage, and strong signal with zero $\pi$ and zero $\eta$, and return `(0, inf)` on right-angle Scene C (`G29`). On anti-parallel `(2,-1), (1,-2)`, `T45` shows a hilltop ghost (`[0.954, 1.048]`) while Moreira CLR returns `EMPTY` ($K_{\min} = 500$), at the sole cost of counting noise dimensions $M = \sum k_i$.
4. **Prior-free composition via infimal convolution** (`L5`, `L14`, `B26`): $K_{\rm lump}(M) = \inf_{m_1+m_2=M}[K_1(m_1)+K_2(m_2)]$ (and reciprocal dual for common force) is exact on noiseless inputs and Gaussian variances, immune to `B26`, and maps `5 ± 0.5 kg` glued to `0/0` to `[4.05, inf)`. Weak-signal calibration of composed curves remains open.
5. **Unnormalized Door B (`T48` Student-$t$, $\nu=10$, no $\pi$, no $\eta$)** (`O18`, `O14`): removes `ours33`'s 79-decade double-discount blowup (`[4.915, 6.09]` at 10σ) and its 50σ crash (splits into `[4.968, 5.041] U [5.963, 6.035]`), and returns `(0, inf)` on `0/0`.

## Files & verification

- Writeup & tables: `explore/e13-unnormalized-evidence.md`
- Pure-stdlib prototype & killer gallery: `python3 explore/e13_unnormalized_evidence.py`
- Unit & coverage checks: `python3 -m unittest -v tests/test_e13_unnormalized_evidence.py`

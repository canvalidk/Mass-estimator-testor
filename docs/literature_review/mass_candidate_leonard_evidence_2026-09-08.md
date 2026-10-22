# Leonard: first comparison evidence — 2026-09-08

Status: local active/unsubmitted research note. Saving here does not submit it to VD-docs. This bounded audit consumes the primary publication below and the comparison contract supplied in the current task; it consumes no managed VD-docs source. The coordinating review records the managed-context snapshot separately. No historical-priority conclusion is established.

Verdict: **demonstrated limitation for the published diffuse posterior**. A modified positive-scale, fixed-origin, known-covariance model remains **unresolved**; it is a different candidate requiring its own specification and proof.

## Source facts and equation anchors

[David Leonard, *Estimating a bivariate linear relationship*](https://arxiv.org/pdf/1202.0957), originally published in 2011, DOI [10.1214/11-BA627](https://doi.org/10.1214/11-BA627). Page numbers below are one-based PDF page numbers.

Equations (1)–(4), pp. 3–4, model independent observation pairs, normally distributed latent coordinates, an estimated intercept, an unrestricted slope, and unknown independent Gaussian channel-error variances. These are not a supplied general joint measurement-error covariance.

The diffuse posterior, equations (13)–(21), pp. 5–7, depends on centered sample covariance through

$$
r=\frac{S_{12}}{\sqrt{S_{11}S_{22}}},\qquad
l=\sqrt{\frac{S_{22}}{S_{11}}},\qquad \widetilde\beta=\beta/l.
$$

Equation (19) gives a scale family; equation (14) retains correlation sign through $r\operatorname{sign}(\beta)$. Page 7 and Appendix 2, p. 23, state the integral prescription is well-defined for $n>2$. Equations (16)–(18) also require nondegenerate correlation. Equation (51), p. 23, gives reciprocal invariance. Equation (12) retains prior covariance parameters before the diffuse limit.

## Deductions for the comparison contract

The contract permits any finite measured vectors with a known positive-definite error covariance. Its true-vector relation has zero intercept and positive scale. The following are deductions from the cited equations, not assertions attributed to Leonard.

### 1. Both-zero input is outside the displayed diffuse prescription

Treat the three spatial components as three observation pairs; thus this argument does not reject regression merely because it has multiple observations. Set

$$
\mathbf a=\mathbf F=(0,0,0),\qquad C=I_6.
$$

This is an allowed measured input with strictly positive error variances. Nevertheless, $S=0$, making both $r$ and $l$ undefined. The supplied error covariance does not substitute for $S$: one is measurement uncertainty, the other is centered scatter of measured components. The same problem arises with nonzero constant-component vectors such as $\mathbf a=(1,1,1),\mathbf F=(2,2,2)$.

### 2. No unique continuous extension repairs zero–zero

For any $c>0$, take the following path with the same $C=I_6$:

$$
\mathbf a_\epsilon=\epsilon(-1,0,1),\qquad
\mathbf F_\epsilon=\frac{c\epsilon}{\sqrt3}(1,-2,1),\qquad\epsilon>0.
$$

Both component means vanish. Direct calculation gives

$$
S_{11}=\epsilon^2,\quad S_{22}=c^2\epsilon^2,\quad
S_{12}=0,\quad r=0,\quad l=c.
$$

Every path approaches the same zero–zero measurements while keeping nondegenerate scatter for every positive $\epsilon$. Equation (19) therefore gives $p_c(\beta)=c^{-1}p_1(\beta/c)$, independent of $\epsilon$ but dependent on arbitrary $c$. These are distinct proper distributions, so the diffuse posterior has no unique continuous limit at zero–zero.

Even conditioning on $\beta>0$ does not remove this problem. Reciprocal invariance at $l=1$ implies equal conditional probability below and above 1: substituting $u=1/\beta$ exchanges the two intervals. Thus the positive-conditioned posterior median along this path is $c$. Different paths can approach zero with medians 1, 2, or any other positive value. The unconditioned native posterior is symmetric when $r=0$, giving median zero instead.

### 3. Centering loses information required by a fixed-origin vector problem

Consider $\mathbf a=(1,2,3)$, $\mathbf F=(2,5,6)$, and $\mathbf F'=(-8,-5,-4)$, again with $C=I_6$. Both pairs have

$$
S_{11}=1,\quad S_{12}=2,\quad S_{22}=13/3,
$$

so the published posterior is identical. However, $\mathbf a\cdot\mathbf F=30$ and $\mathbf a\cdot\mathbf F'=-30$. The raw vector relative-direction sign changes. This is appropriate translation invariance for an unknown-intercept regression, but it does not encode the supplied fixed-origin directional relation.

## What this does and does not settle

The displayed method is not a complete answer to the stipulated inputs: the zero–zero counterexample alone settles that bounded verdict. It also cannot directly consume arbitrary correlated and unequal directional errors. Conversely, it does use signed centered correlation and supplies a posterior on its native domain; it should not be described as direction-blind or as merely a point estimator.

Fixing the intercept, specifying positive support, supplying the full covariance, or retaining an informative prior may yield another legitimate method. These changes require a fresh derivation and coverage analysis. The present counterexample does not rule out such adaptations or prove novelty of the VD equation. No numerical sweep is needed for this verdict because the failure is algebraic.

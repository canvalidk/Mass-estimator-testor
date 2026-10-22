# von Luxburg-Franz / Fieller geometry: second bounded evidence pass

Date: 2026-09-08.
Status: local active research record; awaiting transfer to VD-docs. Saving here does not submit it.
Scope: one source, its exact pivot and concrete scalar boundary cases. Managed-context provenance belongs to the parent comparison record. No broad history search or implementation was performed.

## Source and access

[Ulrike von Luxburg and Volker H. Franz, *A Geometric Approach to Confidence Sets for Ratios: Fieller's Theorem, Generalizations, and Bootstrap*](https://arxiv.org/pdf/0711.0198), arXiv:0711.0198v1, 1 November 2007. The full 24-page primary PDF was readable through the web tool. Relevant equation and theorem text was inspected, with the pivot's algebra checked below.

Source anchors: section 1.1, PDF page 3, defines the native point as the ratio of sample means and allows correlated scalar observations. Construction 1 and Theorem 1, pages 5-7, establish geometric confidence sets. Definition 2, Theorem 3 and equation (2.1), pages 8-9, state the Fieller set and pivotal proof. Theorem 4, page 10, equates the two constructions. Construction 2/Theorem 5, page 12, generalize confidence-set inversion to other projected distributions. The native Gaussian presentation estimates covariance from repeated pairs.

## Published pivot, translated

Identify the numerator reading with F and denominator with a. Their equation (2.1) is the signed residual divided by its standard error:

$$
T_m=\frac{F-ma}{\sqrt{C_{FF}-2mC_{Fa}+m^2C_{aa}}}.
$$

In the native sample formulation C denotes the estimated covariance of the sample means, with a Student-t critical value. For our fixed known covariance input, the same projection argument gives a standard-normal pivot at the true m and therefore critical value q=Phi^(-1)(1-alpha/2). This known-covariance specialization is derived here; one should not feed a known covariance into the paper's unknown-covariance Student-t formula without changing the reference distribution.

For the scalar relation F_*=m a_*, define

$$
\mathcal C(F,a)=\{m\in\mathbb R:
(F-ma)^2\le q^2(C_{FF}-2mC_{Fa}+m^2C_{aa})\}.
$$

The covariance quadratic is positive for every finite m because C is positive definite. This single quadratic inequality handles all readings without division by the measured a. The cross-term explicitly uses correlated errors; changing the relative reading signs changes the set. Solving it yields bounded intervals, two unbounded pieces, or the entire line. The set can change topology continuously through its coefficients; the explicit root formula's zero leading coefficient is not a failure of the defining inequality.

## Our decisive cases

For transparent examples set C=I and q=1.96, a rounded 95% normal critical value.

| Readings (F,a) | Derived positive-restricted confidence set C intersect (0,infinity) | Native point F/a |
|---|---|---|
| (0,0) | All positive values | Undefined |
| (0,10) | (0, q/sqrt(100-q^2)] | Zero |
| (10,0) | [sqrt(100/q^2-1), infinity) | Undefined |
| (10,10) | A bounded interval around 1 | 1 |
| (-10,10) | Empty | -1 |

For the last row, every m>0 has (F-ma)^2=100(1+m)^2 >=100(1+m^2)>q^2(1+m^2); no positive ratio is accepted. This is an ordinary finite Gaussian reading with positive-definite covariance. It is a concrete limitation against a contract demanding a positive estimate and a reconciled positive distribution for every reading.

The empty restricted set does not make confidence-set inference invalid. If the true m is positive, intersecting a confidence set with (0,infinity) leaves its coverage event unchanged. It can still miss every admissible value for some data. Similarly, unbounded sets express genuine lack of constraint and should not be counted as uncertainty failure merely for lacking finite endpoints.

These statements concern measured zeros. The original ratio parameter E(Y)/E(X) presupposes a nonzero true denominator. If instead the parameter is an asserted relation indexed by m even when both true means are zero, the residual-pivot argument still works, although m is then unidentifiable; that is a relational extension of the ratio wording.

## Uncertainty interpretation and verdict

A confidence set is not a probability density over m. Normalizing its width is impossible when unbounded and would not convert frequentist confidence into posterior probability. A positive-conditioned Gaussian-ratio posterior can have finite interior quantiles on the same weak data because it uses a prior/base measure and a different probability statement. Neither result can be substituted for the other without saying what changed.

**Demonstrated limitation of the native point-estimator-plus-set method against the full requested contract:** F/a fails zeros and positivity. **Strong coverage result for scalar uncertainty under its own frequentist meaning:** the covariance-aware pivot defines confidence sets throughout the measured-data domain. Positive restriction preserves that coverage but can produce an empty set, so it does not alone supply the requested all-readings positive point or posterior distribution. The paper does not establish our full vector estimator or its exact identity.

The main incremental learning is a fair comparison rule: distinguish estimator existence, posterior propriety and repeated-sampling confidence coverage. This source adds a covariance-aware uncertainty mechanism, rather than another candidate positive-scale point formula. Further generic Fieller geometry is likely to have diminishing returns for the current vector-estimator question unless it supplies a new estimator or a specifically relevant vector extension.

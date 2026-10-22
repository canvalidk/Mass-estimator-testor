# Zhou, Kou, Li and Fang (2016): full-covariance structured EIV

Date: 2026-09-08. Status: local evidence note, awaiting transfer; NOT submitted to VD-docs. Scope: one new publication and its direct scalar-scale specialization. The combined review records the managed target snapshot; no managed source is modified here.

Source: Yongjun Zhou, Xinjian Kou, Jonathan Li and Xing Fang, *Comparison of Structured and Weighted Total Least-Squares Adjustment Methods for Linearly Structured Errors-in-Variables Models*, Journal of Surveying Engineering (2016), DOI 10.1061/(ASCE)SU.1943-5428.0000190. [Author-hosted full publication](https://uwaterloo.ca/geospatial-intelligence/sites/default/files/uploads/files/2016-zhou-kou-li-fang-jse.pdf). Page references below use the article suffix, e.g. page 04016019-4 is PDF page 4.

## Published prescription

Equations (3)-(6), p. 2, model $y=A(x)\xi$, with a matrix linear in noisy input coordinates, Gaussian input/output errors, and their complete block cofactor matrix, including cross-covariance. Equation (18), p. 4, minimizes $e^TP_le$ subject to the corrected exact relationship, with $P_l=Q_l^{-1}$. Equation (23) eliminates the measurement corrections. Equation (32), p. 5, supplies a first-order parameter covariance approximation:

$$
Q_{\widehat\xi}=
\left[\widetilde A^T(\widehat H Q_l\widehat H^T)^{-1}
\widetilde A\right]^{-1}.
$$

The formulation takes $\xi\in\mathbb R^p$; positivity is not specified. It assumes full column rank and notes that parameter constraints can be imposed. Its Newton algorithm has no guaranteed global convergence. The covariance formula is explicitly an approximation, not a full parameter posterior.

## Our direct input mapping

Set $p=1$, $x=a$, $y=F$, $A(x)=x$ as a one-column matrix, and $\xi=m$. This meets the published linear structure without introducing a fitted rotation, intercept, or new observations. The measurement vector is exactly the two vectors in common fixed coordinates. For $d\ge2$, the published overdetermined size condition is satisfied; the counterexample below also has a nonzero input column.

Set its covariance scale to the known value and its block matrix to the supplied $\Sigma$. This is a straightforward known-variance specialization. Full covariance is a strength of this candidate; the audit does not fault it for independent or isotropic errors.

Write $\Sigma$ in acceleration/force block order. Eliminating the corrections in Eq. (18), equivalently specializing Eq. (23), gives

$$
\chi^2(m)=(F-ma)^TV_m^{-1}(F-ma),
\qquad
V_m=\Sigma_{FF}+m^2\Sigma_{aa}
-m(\Sigma_{aF}+\Sigma_{Fa}).
$$

This is our notation for the published profiling objective. $V_m$ is SPD for every finite $m$, because it equals $[-mI\ I]\Sigma[-mI\ I]^T$. Its finite inverse does not ensure a finite positive optimum. The residual retains relative signs and all relevant cross-covariance terms. There is no marginalization determinant factor in this profiling objective.

## Our counterexample with correlated errors

Take $a=(1,0)$, $F=(-1,0)$, and

$$
\Sigma=\begin{bmatrix}I_2&\rho I_2\\\rho I_2&I_2\end{bmatrix},
\qquad \rho=0.4.
$$

This covariance has eigenvalues $1\pm\rho$, all positive. The unconstrained fit $m=-1$ has zero residual. If strict positivity is added as an adaptation, then

$$
\chi^2(m)=\frac{(1+m)^2}{1+m^2-2\rho m}
=1+\frac{2(1+\rho)m}{1+m^2-2\rho m}>1
\quad (0<m<\infty).
$$

Its infimum 1 is approached at $m\downarrow0$ or $m\to\infty$, never attained at a finite positive value. The positive-scale model remains statistically possible under Gaussian error; this is not an exact, uncertainty-free contradiction. Closed constraints permitting zero would attain a nonpositive boundary answer. Bounds strictly away from zero could force a positive estimate, but would add a parameter restriction whose inferential effect needs a separate audit.

At both measured vectors zero, the residual is zero for every $m$, so the profiled criterion does not identify mass. Extending the published full-rank setting to this case also makes the information term in Eq. (32) vanish. Thus the stated first-order covariance does not provide the requested uncertainty law at this input. This claim concerns the formula's boundary, not a general impossibility of uncertain mass inference.

## Verdict and contribution

**Demonstrated limitation for the published optimization prescription under the full positive-scale contract.** The method already handles the covariance and coordinate structure needed by the task; it fails positivity without constraints and may lack an attained finite solution with a strict positive constraint. Its local covariance approximation is not a globally defined positive-scale distribution. The specific prior measure, pushforward mass law and point statistic of the target are not identified with this method.

What this adds: it closes a weakness in comparisons limited to white-noise TLS. Arbitrary cross-covariance is already within an earlier EIV formulation. The limiting issue is the objective and parameter domain, not absence of directional information or a crude covariance model. The core boundary mechanism repeats the earlier TLS result, so further sweeps of this prescription would have low additional value.

Coverage: full primary formulation through Eq. (32), pp. 1-5, and conclusions inspected. Proof above is an independent derivation from the specified specialization. No code run or numerical convergence claim is needed to establish its counterexample. Further Bayesian regularizations and other cited methods were not audited in this bounded task.

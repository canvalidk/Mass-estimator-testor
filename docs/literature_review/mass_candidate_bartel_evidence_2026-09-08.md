# Bartel, Stoudt and Possolo (2016): bounded mass-candidate audit

Date: 2026-09-08. Status: local evidence note, awaiting transfer; NOT submitted to VD-docs. Scope: one primary publication and the supplied comparison contract. No managed VD-docs source is interpreted here; the combined review records the managed target snapshot.

Source: Thomas Bartel, Sara Stoudt and Antonio Possolo, *Force Calibrations using Errors-in-Variables Regression and Monte Carlo Uncertainty Evaluations*, Metrologia 53, 965-980 (2016). [NIST publication record](https://www.nist.gov/publications/force-calibrations-using-errors-variables-regression-and-monte-carlo-uncertainty); [NIST-hosted full manuscript](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=919758). Page locators below use manuscript page numbers; add two for the PDF viewer page number.

## Published prescription

Section 3.2.2 (p. 7) specifies noisy force and response pairs, $F_j=\phi_j+\delta_j$, $R_j=C(\phi_j)+\epsilon_j$, initially with independent Gaussian errors and known component standard uncertainties. Section 3.2.6, Eq. (3), p. 11, minimizes

$$
S^*(A,\phi)=\sum_j\left[
\frac{(R_j-C_A(\phi_j))^2}{u(R_j)^2}
+\frac{(F_j-\phi_j)^2}{u(F_j)^2}
\right],\qquad
C_A(t)=A_0+A_1t+A_2t^2+A_3t^3.
$$

Degree selection is described in section 3.2.3. Section 5.3 (pp. 20-22) simulates input errors, refits the EIV curve and constructs coverage bands. It includes a shared electrical-calibration error and within-run common-sign perturbations. Thus uncertainty propagation does include particular correlations; the fitting objective is diagonal and the simulation is not a general full joint Gaussian covariance prescription. No strictly positive slope constraint is stated. Section 3.2.3 (p. 9) also describes bounded latent-force searches during repeated fits.

## Our equation-level specialization

The many calibration pairs are not, by themselves, a reason to exclude this method. Associate index $j$ with the fixed coordinates of the measured vectors, map the paper's $F_j$ to our $a_j$, and its $R_j$ to our $F_j$. Select degree one and fix the intercept to zero, so $C(t)=mt$. Fixing the intercept is an explicit specialization for the physical relationship; it is not the native free-intercept polynomial fit.

For joint measurement covariance $I_{2d}$, the specialized objective becomes

$$
Q(m,\alpha)=\|a-\alpha\|^2+\|F-m\alpha\|^2,
\qquad
\alpha_m=\frac{a+mF}{1+m^2},
$$

and our elimination of the latent vector gives

$$
Q_{\rm prof}(m)=\frac{\|F-ma\|^2}{1+m^2}.
$$

This objective retains the signed cross term $-2m\,a^TF$. It does not inherently discard the observed relative signs. The problem is what happens when the admissible scalar must remain positive.

For $a=(1,0)$, $F=(-1,0)$, covariance $I_4$: allowing any real slope gives a zero-residual fit at $m=-1$. Adding $m>0$ as an adaptation gives

$$
Q_{\rm prof}(m)=\frac{(1+m)^2}{1+m^2}
=1+\frac{2m}{1+m^2}>1.
$$

The infimum 1 is approached at zero or infinity and is not attained at a finite positive $m$. At both measured vectors zero, $\alpha=0$ yields objective zero for every slope. Thus minimization does not select a unique mass there. These inputs remain legitimate noisy observations; unit covariance has no zero uncertainty directions.

Resampling and refitting an optimizer is not by itself a rule for selecting among tied optima, attaining a nonexistent positive optimum, or ensuring a positive scale distribution. Bounds or tie-breaking could produce an implementation output, but their inferential effect would need a separate specification and audit. Arbitrary correlated covariance would likewise require replacing or extending the native criterion. No claim is made that these extensions are impossible.

## Verdict and contribution

**Demonstrated limitation for coverage of the supplied positive-scale contract by the published prescription and the explicit linear specialization above.** The native artifact is a fitted calibration function and propagated uncertainty; it is not established as the requested universal positive-scale estimator. Exact identity with the workspace equation is not established.

What this adds: a real uncertainty-propagation procedure that does account for some dependence, and the distinction between propagating uncertainty through fitted values and defining a globally valid positive-scale estimator. The anti-aligned optimization obstruction repeats the TLS boundary mechanism; further simulations of that obstruction would be diminishing returns.

Access: NIST-hosted full manuscript, relevant text of sections 3.1-3.3 and 5.2-5.4, including Eq. (3) and Monte Carlo steps UC-1 through UC-7, inspected. Web screenshot calls did not deliver a viewable image; formulas above were checked against the extracted equation and its surrounding text. This note neither claims an exhaustive audit of EIV literature nor treats specific simulated correlations as support for arbitrary covariance matrices.

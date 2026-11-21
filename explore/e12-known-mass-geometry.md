# Adding a known mass in the same shared-acceleration record

Consider an unknown mass \(m>0\), a known exact mass \(b>0\), and independent records \(y_i=(F_i,a_i)\). The added force channel in this protocol is calculated as \(b a_i\) from the **same measured acceleration**. The observation for total mass \(M=m+b\) is \(y_i'=(F_i+b a_i,a_i)\). The declared medium, spatial subspace and sampling window stay the same. This protocol is one invertible change of the existing observations. A separately observed force for the known mass has its own noise and requires a joint acquisition model. [Physical composition protocols](e12-joining-bodies.md)

## Spatial pairing, defects and noise

Where \(a\ne0\), the reciprocal acceleration covector \(a^\vee=g(a,\cdot)/g(a,a)\) obeys
\[
 a^\vee(F+b a)=a^\vee(F)+b,\qquad
 F+b a-a^\vee(F+b a)a=F-a^\vee(F)a .
\]
The exterior product \(F\wedge a\) is also unchanged. Under the full isotropic linear inertia premise the noiseless inertia map becomes \(mI+bI\). The measured action on a single spatial line, \((F+b a)\otimes a^\vee\), gains \(b\,a\otimes a^\vee\), the projector onto that measured line times \(b\). This distinguishes what the measurement establishes from what isotropy supplies on other spatial directions. At a measured \(a=0\) the residual identities below continue to hold.

In the stacked channel space set \(S_b=\left(\begin{smallmatrix}I&bI\\0&I\end{smallmatrix}\right)\). Apply \(y'=S_b y,\ \Sigma'=S_b\Sigma S_b^\top\). In particular \(\Sigma'_{Fa}=\Sigma_{Fa}+b\Sigma_{aa}\), even for initially independent channels. With \(D_m=[I,-mI]\),
\[
 D_{m+b} S_b=D_m,\qquad
 D_{m+b}\Sigma' D_{m+b}^\top=D_m\Sigma D_m^\top .
\]
The measured defect, its covariance and its squared noise length coincide at the corresponding masses. For example, \((F,a)=(2,1)\) with unit independent errors, \(m=1,\ b=1\) gives original cost \(1/2\). The total data are \((3,1)\) with covariance \(\left(\begin{smallmatrix}2&1\\1&1\end{smallmatrix}\right)\). At total mass 2 the cost is again \(1/2\). Replacing that covariance by its diagonal gives \(1/6\).

## Projector and tangent transport for arbitrary instruments

Let \(N_i\) contain the columns of the \(k_i\)-dimensional declared spatial subspace, \(G_i(m)=\binom{mN_i}{N_i}\), and \(L_iL_i^\top=\Sigma_i,\ L'_iL_i'^\top=\Sigma'_i\). The matrix \(O_i=L_i'^{-1}S_b L_i\) is orthogonal. The physical graph and its whitened projector satisfy
\[
 S_bG_i(m)=G_i(m+b),\qquad P'_i(m+b)=O_iP_i(m)O_i^\top ,
\]
where a prime on \(P_i\) labels the total-mass experiment. The same conjugation holds for projection onto the span of both copies of \(N_i\). The T36 candidate cost \(K_i\) within that declared span, its fitted acceleration in physical units, and the entire residual curve therefore transport. For the fitted acceleration this also follows directly from equality in \(\alpha\) of the two quadratic forms \(\|y_i-G_i(m)\alpha\|_{\Sigma_i^{-1}}^2\) and \(\|y'_i-G_i(m+b)\alpha\|_{\Sigma_i'^{-1}}^2\). Each has a unique minimizer. The identity holds across nonidentical instruments with general positive-definite cross-channel covariance. [T36](../definition/theorems/T36.md)

Differentiate the conjugation in the **physical mass coordinate**. The fitted transverse tangent \(t_i=(\partial_m P_i)P_iL_i^{-1}y_i\) and residual \(r_i\) both transform by \(O_i\). Thus the reciprocal correction in mass units \(t_i^\vee(r_i)\), its weight \(\|t_i\|^2\), the score \(\langle t_i,r_i\rangle\), and the zero-tangent condition all transport. The same is true of their multi-reading contractions and of the exact within, between and common squared energies. The implementation uses the angle \(\phi=\arctan m\) in fixed mass units. Comparing numerical scores in those coordinates introduces the factor
\[
 \lambda=\frac{1+(m+b)^2}{1+m^2}.
\]
At the total candidate, the angle-coordinate score is multiplied by \(\lambda\), the angle weight by \(\lambda^2\), and any defined angle displacement by \(1/\lambda\). The signed standardized score and its test probability agree. Hence the exact conditional-score acceptance set on \(M>b\) is the original set shifted by \(b\), including the degenerate acceptance convention. Inverting on the same parameter domain is essential. For instance the reading \((-1,1)\) with unit errors has least cost 1 on the closed positive mass arc. Adding \(b=2\) gives a reading \((1,1)\) with cost zero at total mass 1, which corresponds to an excluded original mass. On \(M\ge2\) the least cost is again 1. [T45](../definition/theorems/T45.md) · [Reciprocal tangent decomposition](e12-tangent-pairings.md)

**Consequence for the complete law.** If a correction is a functional of the transported per-reading \(K_i\) curves and their invariant tangent energies, evaluating and minimizing it on the corresponding mass domain keeps that calculation unchanged. A likelihood factor \(L'(m+b)=c(y)L(m)\), paired with the pushforward \(\pi'\) of the starting law, gives the translated normalized law by change of variables. Pushing the starting law with the experiment sends its support to \(M>b\). For example, the shifted unit half-Cauchy density is \(2/[\pi(1+(M-b)^2)]\) above \(b\). A newly declared unshifted half-Cauchy on all positive total masses describes a different starting law. A rule which separately uses channel-axis angles, a precision cap or the instrument ratio must also specify how those quantities are carried. Invariance of the residual yields no such continuation automatically. [E12 finite readouts](e12-joining-residuals.md)

## Two existing project discounts both satisfy the comparison

The existing series discounts already illustrate the freedom after the geometric tests. Fix corresponding per-reading weights \(g_i\) and put \(\ell_i=g_iH_i/2\), where \(H_i=\|P_iL_i^{-1}y_i\|^2\). Each \(H_i\) and each total declared-span energy are invariant. The whole-projective maximizer shifts when unique, with infinity fixed. At a regular maximizer the total first derivative vanishes. By the chain rule both the minus second derivative of \(\sum_i\ell_i\) and \(\sum_i(\partial_\phi\ell_i)^2\) acquire the *same* squared Jacobian when changing angle charts. Thus the repo's Bartlett ratio and its cap at one commute with the shear. The repo's alternative \(1-Q/H\) with weighted whole-line misfit \(Q\) and weighted explained energy \(H\) commutes as well, directly by projection. Any positivity-dependent per-reading weights must use corresponding mass cones on both sides.

For scalar readings \((4,1),(1,4)\) with unit independent errors both rows admit a positive fit, giving unit per-reading weights. The whole-line best mass is 1. The curvature and score-square sum in the \(\arctan m\) coordinate are 16 and \(225/2\), hence the Bartlett discount is \(32/225\). The whole-line costs are \(Q=9,H=25\), giving the other discount \(16/25\). After adding any known positive mass the corresponding discounts are unchanged while the best mass moves to \(1+b\). The different finite relative evidence survives, using previously implemented formulas. These ratios do not inherit the full E12 conflict or repeated-data guarantees from this comparison. [Series calculations in met.full](../met/full.py)

## Known-mass addition also preserves two distinct Liseo specializations

For a scalar reading with equal errors, Liseo's equations (12) and (14) give different ratio posteriors. The following extensions through the known shear can be checked directly. Start with the positive restrictions of both laws. For (12), draw the two signed means \((X,Y)\) from their joint Gaussian posterior under a flat Cartesian prior. Translating the numerator mean to \(X+bY\) and applying the full covariance congruence gives exactly the pushforward of that Gaussian, since \(S_b\) has determinant one. For \(Y\ne0\) the ratio becomes \(X/Y+b\), and the event \(X/Y>0\) maps to \((X+bY)/Y>b\).

For the reference posterior in (14), its equal-error form is uniform projective angle times the profile factor \(e^{-K(m)/2}\). A direct covariant extension uses uniform angle in the whitened channel plane. With \(v(m)=(1,-m)\Sigma(1,-m)^\top\), its volume element in mass coordinates is \(\sqrt{\det\Sigma}\,dm/v(m)\). Indeed the whitened graph direction is \(L^{-1}(m,1)^\top\); its angular speed has magnitude \(|\det L^{-1}|/\|L^{-1}(m,1)^\top\|^2\). Both \(\det\Sigma\) and \(v(m)\) are exactly preserved by the matched shear. The volume and the whole profile curve therefore translate, including their normalization on \(M>b\). This calculation derives the correlated transport of (14) from its equal-error expression. [Liseo (2003), pp. 137–138](https://www.dss.uniroma1.it/RePec/mtn/articoli/2003-1-10.pdf) · [Earlier equation comparison](e12-joining-correction.md#a-boundary-obstruction-even-before-joining-two-readings)

The laws already differ for \(F=2,\ a=0\) in unit noise. Conditional on positive mass, (12) puts \(2\Phi(\sqrt2)\Phi(-\sqrt2)\) below \(m=1\). For (14) the probability is the integral of \(\exp[-2\cos^2\theta]\) over \(0<\theta<\pi/4\) divided by the same integral from 0 to \(\pi/2\). Numerically the two probabilities are (0.144928) and (0.219508), respectively. Adding a known mass \(b\) preserves the corresponding two probabilities below \(M=1+b\). Thus this shared-record operation by itself also respects historical prescriptions which demonstrably disagree. [Executable checks](../tests/test_e12_known_mass_geometry.py)

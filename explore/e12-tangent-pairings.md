# Joining the reciprocal tangent pairings at one candidate mass

**Claim.** A candidate mass and the declared noise metric determine a signed local mass correction for each nondegenerate reading. Enforcing one common local correction splits the total residual into three orthogonal pieces. Their squared lengths have independent, exactly known distributions at every fixed true mass in the declared Gaussian experiment. This gives an exact finite-sample meaning to within-reading residual, between-reading disagreement and common mass score. It also restricts how disagreement may widen an exact confidence set. [T36, reading geometry](../definition/theorems/T36.md) · [T45, conditional score](../definition/theorems/T45.md)

## The reciprocal of the candidate tangent

Use the independent, positive-definite Gaussian reading type of T45 with deterministic true accelerations in their declared subspaces. The algebra before the distributional statement requires no Gaussian sampling assumption. Whiten an observation to \(y_i\), and let \(W_i\) be the span of the whitened force and acceleration copies of its declared \(k_i\)-dimensional spatial subspace. In \(W_i\), let \(P_i\) project orthogonally onto the candidate graph of \(F=m a\). Let \(U_i=P_i y_i\), \(r_i=P_{W_i}y_i-U_i\). Thus \(r_i\) lies in a \(k_i\)-dimensional residual space, and \(K_i=\|r_i\|^2\) is the T36 candidate cost. Discarding \(W_i^\perp\) removes a mass-independent residual of dimension \(2(d_i-k_i)\). It is orthogonal to all graph tangents.

Differentiate the graph in a smooth regular mass coordinate \(x\), and write
\[
 t_i=P_i'(x)U_i,\qquad h_i=\langle t_i,t_i\rangle,\qquad s_i=\langle t_i,r_i\rangle .
\]
Differentiating \(P_i^2=P_i\) shows \(t_i\) is perpendicular to the graph, and it remains in \(W_i\). It is the transverse first-order motion of the fitted signal. At \(t_i\ne0\), its reciprocal covector in the noise metric is \(t_i^\vee(v)=\langle t_i,v\rangle/h_i\). Applying it to the observed residual gives \(b_i=t_i^\vee(r_i)=s_i/h_i\). This is the unique signed displacement minimizing \(\|r_i-bt_i\|^2\). It is a tangent displacement at the supplied candidate. A finite correction on the curved family need not equal the mass that fits this reading best. At zero tangent there is no reciprocal and \(s_i=0\).

With full spatial subspace there is a direct tensor expression in the original coordinates, using \(x=m\). Write \(d=F-ma\), \(D=[I,-mI]\), \(V=D\Sigma D^\top\), \(G=V^{-1}\), \(C=\Sigma_{aF}-m\Sigma_{aa}\) and \(\hat a=a-CGd\). This is the fitted acceleration on the graph. The projected reading is \((m\hat a,\hat a)\) since the correction to \((F,a)\) is \(\Sigma D^\top Gd\). Applying \(D\) to the transverse motion of a graph point gives \(\hat a\), and the squared norm of the shortest correction with image \(z\) under \(D\) is \(z^\top Gz\). Hence
\[
 h_i=\hat a^\top G\hat a,\qquad s_i=\hat a^\top Gd=-\tfrac12\,\partial_m K_i,\qquad b_i=\frac{\hat a^\top Gd}{\hat a^\top G\hat a}
\]
when \(h_i>0\). The metric and the fitted acceleration both depend on the candidate. At truth the covariance of \(\hat a\) with \(d\) vanishes by construction, so they are independent under the joint Gaussian law. This addresses the correlation that using the measured \(a\) in this pairing would introduce.

Under a reparameterization \(x=x(z)\) with nonzero derivative, \(t_z=x't_x,\ s_z=x's_x,\ h_z=(x')^2h_x,\ b_z=b_x/x'\). Thus the numerical displacement has the units of its coordinate and all three energies below are invariant. T45 and the implementation use \(x=\arctan m\) in the supplied mass units. Then \(s_x=(1+m^2)s_m\) and \(b_x=b_m/(1+m^2)\).

## Three nested subspaces determine the join

Form the direct sum of the residual spaces. Let \(L_{\rm sep}\) be the orthogonal sum of the individual tangent lines (the zero space for an inactive reading), and let \(L_{\rm com}=\operatorname{span}\{(t_i)_i\}\). These are the transverse motions allowed if each reading can move independently in mass and if every reading must move by the same amount. Set \(H=\sum_i h_i\), \(S=\sum_i s_i\). For \(H>0\) the reciprocal of the stacked tangent applied to the stacked residual gives \(\bar b=S/H\). It is the unique minimizing common displacement. The residual splits over \(L_{\rm sep}^\perp\), \(L_{\rm sep}\cap L_{\rm com}^\perp\) and \(L_{\rm com}\), giving
\[
 \sum_i K_i
  = \underbrace{\sum_i\|r_i-b_i t_i\|^2}_{W}
    +\underbrace{\sum_{h_i>0} h_i(b_i-\bar b)^2}_{J}
    +\underbrace{S^2/H}_{T}.
\]
For an inactive reading use zero projected vector in \(W\). At \(H=0\) both tangent spaces vanish, \(J=T=0\) and \(W=\sum K_i\). The three energies are unique given these subspaces and the metric by uniqueness of orthogonal projection. Only within-reading inner products occur, so separate spatial frames need no alignment. At an interior stationary point of the pooled Gaussian profile \(S=0\), the deterministic identity still holds. Evaluating at that fitted point generally changes the sampling distribution.

Let \(A,B\) be disjoint blocks at the same candidate and in the same coordinate. The summaries \(S,H,W\) add. The between-reading energy obeys
\[
 J_{A\cup B}=J_A+J_B+
 \frac{H_AH_B}{H_A+H_B}(\bar b_A-\bar b_B)^2.
\]
The last term is zero if either block has zero weight. This follows by expanding \(b_i-\bar b_{A\cup B}\) about its block mean and using zero weighted deviation in each block. Hence regrouping the readings gives exactly the same summary. The last term is the additional squared residual from requiring the two block displacements to coincide. This is an exact same-object composition law for the candidate statistic. Physical mass addition for different bodies requires a joint law of their masses.

For example, at \(m=1\) in unit channel noise, two scalar readings \((F,a)=(4,1),(1,4)\) have fitted accelerations \(5/2\), mass-coordinate weights \(25/8\) and signed mass corrections \(+6/5,-6/5\), so the common score vanishes while \(J=9\) and \(W=0\). If those same components are declared as one planar reading, the common score still vanishes, but now \(W=9\) and \(J=0\). In either declaration the total candidate residual is 9. Requiring one mass already within a vector reading changes which tangent directions count as individual motions. The two scalar readings \((2,2),(2,2)\) give all three energies zero at the same candidate. These are exact distinctions relative pooled Gaussian evidence alone erases. [Curve collisions](e12-pooled-information.md)

## One instrument gives a candidate-independent disagreement tensor

Suppose the full declared spatial coordinates of every reading have independent isotropic channel errors with the same noise ratio \(s\). Individual channel standard deviations may vary with the reading at that ratio. Form the noise-unit vectors \(f_i,a_i\), which live in a separate Euclidean space \(V_i\) for each reading. In their orthogonal direct sum let \(f=(f_i)_i,\ a=(a_i)_i\), and let \(B\) be their two-by-two channel Gram. The induced metric on the exterior square gives Lagrange's identity \(\det B=\|f\wedge a\|^2\). The orthogonal decomposition of the exterior square of the direct sum gives the more informative identity
\[
 \det B=
 \sum_i\|f_i\wedge a_i\|^2+
 \sum_{i<j}\|f_i\otimes a_j-a_i\otimes f_j\|^2.
\]
To prove it, use the wedge products of an orthonormal basis in each summand. Basis pairs within one reading give its exterior square. Pairs from distinct summands give an orthonormal basis of \(V_i\otimes V_j\) and the displayed difference of tensor products. Equivalently its squared length is \(\|f_i\|^2\|a_j\|^2+\|a_i\|^2\|f_j\|^2-2(f_i\cdot a_i)(f_j\cdot a_j)\); no cross-reading identification of directions occurs. For locally compatible nonzero pairs \(f_i=\mu_i a_i\) the cross term is \((\mu_i-\mu_j)^2\|a_i\|^2\|a_j\|^2\), with \(\mu_i=m_i/s\). It measures disagreement of the two local divisions exactly even when their signal magnitudes differ. On misaligned vector readings the cross term still has its tensor meaning but need not be a difference of two scalar masses.

For \(L=\lambda_{\max}(B)>0\), dividing this identity by \(L\) gives T32's whole-projective minimum \(Q=\lambda_{\min}(B)\) as a sum of within- and cross-reading contributions. The minimum on the closed positive-mass arc is instead \(\min(B_{11},B_{22})\) when \(B_{12}<0\). Zero exterior product expresses linear dependence, including opposite signs or axis data. Compatibility with a *finite positive* mass requires the sign and nonzero-channel conditions too, with the all-zero list compatible with every mass. [T32, domain of the misfit](../definition/theorems/T32.md)

There is an exact link to the preceding candidate-wise construction. In the noise-angle coordinate \(x=\arctan(m/s)\), write \(u=(\sin x,\cos x)\), \(v=(\cos x,-\sin x)\). Then \(H=u^\top B u,\ S=v^\top B u,\ \sum_iK_i=v^\top B v\). Taking the determinant in the rotated basis gives \(\det B=H(W+J)\) for \(H>0\). Its left side is independent of the candidate. For the conflicting example above \(\det B=225,\ L=25,\ Q=9\). Two scalar readings put the 225 entirely in the cross-reading tensor; one vector reading puts it entirely in its within-reading exterior square. This computes an exact magnitude of incompatibility before choosing a probability law.

## Exact Gaussian noise benchmarks and a widening obstruction

At a fixed true interior mass, condition on all fitted graph projections \(U_i\). The transverse \(r_i\) have independent standard Gaussian coordinates in their respective residual spaces, because their true means lie on the graph. The three nested-space projectors are now fixed. Gaussian independence on orthogonal subspaces implies that \(W,J,T\) are mutually independent chi-squares with respective degrees of freedom
\[
 \sum_i k_i-n_*,\qquad n_*-\mathbf{1}_{H>0},\qquad \mathbf{1}_{H>0},
\]
where \(n_*=\#\{i:h_i>0\}\) and zero degrees means identically zero. This is an application of Gaussian orthogonal projection (Cochran's theorem). T45 shows each tangent is nonzero almost surely in a nonsingular Gaussian experiment with \(k_i\ge1\). The displayed ranks also specify the observed projectors on exceptional deterministic inputs. Our numerical convention drops all-measured-zero rows before the rank bookkeeping. Such padding changes no energy, and has probability zero at every fixed candidate of a nonempty nonsingular sampling experiment. The conditional calibration is almost sure in \(U\). [Cochran (1934), distribution of Gaussian quadratic forms](https://doi.org/10.1017/S0305004100016595)

At a different fixed true mass \(m_0\), the same conditional orthogonality gives independent **noncentral** chi-squares. Write \(\mu_i\) for the whitened true mean and \(\gamma_i=(P_{W_i}-P_i)\mu_i\). Conditional on the fitted graph projections at the candidate, the residuals have means \(\gamma_i\) and identity covariance on their residual spaces. The three noncentrality parameters are exactly the squared lengths of the three orthogonal projections of \((\gamma_i)_i\). For a full spatial reading the mean of \(F-ma\) is \((m_0-m)a_i^*\). The projected mean generally uses the unknown true push, whereas the fitted tangent uses the observed projection. Thus an exact noise distribution at truth does not supply a nuisance-free distribution across candidate masses.

In particular \(Z=S/\sqrt H\) is conditionally standard normal and independent of **both full orthogonal residual vectors** contributing \(W\) and \(J\). Let \(g\) be a finite widening factor determined by \(U\) and those two residuals. A level-\(p\) candidate test accepting when \(|Z|\le z_pg\), where \(z_p=\Phi^{-1}((1+p)/2)\), has conditional coverage \(2\Phi(z_pg)-1\) given these data. If \(g\ge1\) and widens with positive probability, its unconditional coverage is strictly greater than \(p\). Requiring nominal conditional coverage given the disagreement forces \(g=1\) almost surely. Inversion of the inequality is a candidate confidence set under this sampling model. A quoted interval cut out of a separately defined mass law has additional requirements. The argument permits conservative widening and leaves open a different exact test that changes how type-I error is allocated among disagreement levels.

T45 supplied the common statistic \(T\) and its conditional normal score. T50 treats the higher-dimensional common-operator tangent family. B22 compares disagreement and sandwich variance at a *fitted* mass. The added step here is the intervening subspace of independently permitted scalar graph displacements, its uniquely measured distance to the common displacement, its associative join, and the resulting conditional restriction on widening. These statements hold candidate by candidate, and do not yield a nuisance-free likelihood across candidates or select a prior, a discount, or a map from a confidence set to a proper law. [T45](../definition/theorems/T45.md) · [T50](../definition/theorems/T50.md) · [B22](../problems/O14.md)

[Exploration implementation](e12_tangent_pairings.py) · [Checks](../tests/test_e12_tangent_pairings.py)

# A force sensor on the known body

Let the unknown mass be \(m>0\) and let an exact standard \(b>0\) share its spatial acceleration \(a\). On one push measure the unknown body's force \(F\), the acceleration \(A\), and the standard's force \(G\). Suppose the observed triple is Gaussian with known positive definite covariance and mean \((ma,a,ba)\). Each channel has \(d\) components. Independent pushes may have arbitrary fixed accelerations and distinct declared covariances. This gives three physical sensor channels. Using an accelerometer to *calculate* \(bA\) gives a different record. [Shared-drive experiment](e12-joining-bodies.md)

## The exact reduction

The standard supplies the vector defect \(r=G-bA\). Its mean is zero for every parameter. Set \(Z=(F,A)\), \(V=\operatorname{Var}(r)>0\), and define
\[
 B=\operatorname{Cov}(Z,r)V^{-1},\qquad Z_* = Z-Br,\qquad S=\operatorname{Var}(Z)-BVB^\top.
\]
Then \(Z_*\sim N((ma,a),S)\) independently of \(r\sim N(0,V)\), since their cross covariance vanishes and they are jointly Gaussian. The transformation from the original triple to \((Z_*,r)\) is invertible. Thus \(Z_*\) alone retains the full likelihood for the unknown mass and drive. Conversely the full sampling law can be recreated from \(Z_*\) by drawing a parameter-free independent \(r\) and applying the inverse. This is ordinary Gaussian regression and the Schur complement. [Adityanand Guntuboyina, *Statistics 201A*, chapter 5](https://www.stat.berkeley.edu/~aditya/resources/FullLectureNotes201AFall2019.pdf)

Now form the measured lump \(Z_T=(F+G,A)\). Write \(H_b(F,A)=(F+bA,A)\) and \(Jr=(r,0)\). Since \(Z_T=H_b Z+Jr\), its regression coefficient on \(r\) is \(B_T=H_bB+J\). The same correction on the total route gives
\[
 Z_{T,*}=Z_T-B_T r=H_b Z_*,\qquad S_T=H_b S H_b^\top .
\]
The equality holds observation by observation, before any prior or choice of mass law. The mass coordinate is \(M=m+b\), and the starting law on this route is the translated starting law on \(M>b\). Any shear-equivariant readout defined on these effective pair experiments passes the same-record comparison. The calibration sensor reduces uncertainty through the subtraction before the shear. Under physical spatial rotations these vector covariances transform by congruence and both effective records rotate in their respective frames. The reduction applies unchanged to parameter means restricted to a declared subspace.

## Which information does a bare force sum retain?

Fix either raw pair \(Z'\) (part or lump) and write \(B'=\operatorname{Cov}(Z',r)V^{-1}\), \(C'=\operatorname{Var}(Z')\). The conditional law of \(r\) given \(Z'\) has fixed covariance and mean
\(\operatorname{Cov}(r,Z')C'^{-1}(Z'-E Z')\).
Thus \(Z'\) is sufficient for the full triple's family indexed by the mass and the unknown drives exactly when \(\operatorname{Cov}(r,Z')C'^{-1}\) vanishes on the span of the differences of possible means. With unrestricted spatial drives, this span is \(\mathbb R^{2d}\): choosing two distinct masses and the same arbitrary drive generates every force direction by subtraction, then one mass generates the acceleration directions. The criterion becomes \(B'=0\). Since \(B_T=H_bB+J\) with \(J\ne0\), the two bare pairs can never both be sufficient on this full family.

At a single candidate \(x\) for the corresponding pair let \(D_x=[I,-xI]\). The raw transverse covariance exceeds the effective transverse covariance by
\[
 (D_x B')V(D_x B')^\top .
\]
At that candidate the residual vector and its covariance agree on both records for every observation exactly if \(D_x B'=0\). Agreement at every candidate requires \(B'=0\). These statements concern the full spatial force and acceleration channels.

**Example.** On one axis take \(b=1\) and three independent unit-variance errors. The effective part record is \((F,(A+G)/2)\) with covariance \(\operatorname{diag}(1,1/2)\). At \(m=1\) the transverse variance is \(3/2\). The raw measured lump has variance \(6\) in the transverse direction at \(M=2\). Both triples \((F,A,G)=(1,1,1),(2,1,0)\) give the raw lump record \((F+G,A)=(2,1)\). The effective part records fit masses \(1\) and \(4\), respectively. Both source records occur in the support of the same nonsingular Gaussian experiment. [Numerical checks](../tests/test_e12_known_mass_sensor.py)

**No-go for equality of the two bare readouts.** Consider the stronger requirement that a law computed solely from \((F+G,A)\) be the translation of a law computed solely from \((F,A)\), for each observed triple with its fixed declarations. At a fixed \(A\), the map \((F,G)\mapsto(F,F+G)\) is onto. Hence both laws must be independent of their force observations at that \(A\). Even almost-sure equality under a single nonsingular Gaussian distribution gives this conclusion almost everywhere. Fubini's theorem applied at fixed \(A\) to the probabilities of rational endpoint intervals suffices, since those intervals determine the mass law. The same argument applies to each full list of independent triples. At fixed repeated nonzero true drive, the accelerometer list has the same sampling distribution at every true mass, so such a procedure cannot concentrate at two distinct true masses. This rules out pathwise composition of these two incomplete records as a universal learning premise. Equality on equivalent effective records makes a different demand.

**The two candidate mass laws on a measured calibration experiment.** Keep \(b=1\). Give the unknown-force sensor variance 1, and the accelerometer and standard-force sensor variance 2 each, all independent. Take raw triples \((F,A,G)=(3,0,2),(1,1.5,2.5)\). The effective part readings are exactly the unit-noise pairs \((3,1),(1,2)\). On the total route the corrected records are \((4,1),(3,2)\) with covariance \(\left(\begin{smallmatrix}2&1\\1&1\end{smallmatrix}\right)\). The calibration defects are nonzero on both pushes. Starting the part from the unit half-Cauchy law and the lump from its translate gives total medians \(2.528611\) and \(2.510203\) for the translated E12 rules at \(\tau=0,64\), respectively. Both rules respect the samplewise equality of the two effective routes on this three-channel experiment. Experimental equivalence leaves the finite separation intact. [Original readouts and prior](e12-joining-residuals.md) · [Numerical integration](../tests/test_e12_known_mass_sensor.py)

## Where a nonlinear discount can be applied

With full spatial drives, minimizing the triple's Gaussian squared distance over the nuisance drive at each candidate gives exactly
\[
 K_{\rm triple}(m)=K_*(m)+c,\qquad c=r^\top V^{-1}r.
\]
Here \(K_*\) is the ordinary T36 cost of the effective reading. The identity follows by adding independent squared distances in the effective coordinates and observing that \(r\) has mean zero at every mass and drive. Thus \(c\) can test the calibration apparatus, and Gaussian relative likelihoods for mass cancel it. This is the geometric distinction between the residual in the span of all possible means and a residual orthogonal to that span. [T36](../definition/theorems/T36.md)

Suppose instead one uses a continuous scalar \(g(K_{\rm triple}(m))\) as log evidence. Require the resulting normalized law, for a fixed full-support prior, to equal the law with log evidence \(g(K_*(m))\) for every data triple. Then \(g\) is affine. Indeed \(c\ge0\) can be varied arbitrarily while holding the effective reading fixed. Choose a strictly positive scalar effective reading that fits one mass exactly, and scale it to attain any nominated \(x>0\) at another mass. Equality of the two laws implies equality of log factors up to a common constant, everywhere by continuity. Comparing those two masses gives \(g(x+c)-g(x)=g(c)-g(0)\). The continuous additive function \(g-g(0)\) is linear on the nonnegative reals. Conversely every affine \(g\) obeys the equality. For the bounded transform \(h(x)=c_k(1-e^{-x/2})\), the corresponding difference is \(h(x+c)-h(c)=e^{-c/2}h(x)\): the realized calibration error would rescale the mass contrast. Projecting away the parameter-free residual *before* this transform retains the experiment's likelihood equivalence. The condition constrains which residual a discount sees, while leaving the choice among the surviving E12 readouts open. [E12 disagreement construction](e12-joining-residuals.md)

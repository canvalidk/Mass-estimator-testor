# The aperture threshold for the mean residual

**Result.** For a full \(k\)-dimensional reading with separable Gaussian noise \(C\otimes S\), the supremum of the expected residual at the best positive mass is \(k-\tfrac12\) for every channel correlation \(\rho\ge-1/2\). In dimension one that threshold is sharp. For every \(\rho<-1/2\), some finite positive true mass and nonzero push give a mean strictly greater than \(1/2\). The corresponding threshold for the *entire* scalar upper-tail bound is different, at \(\rho=0\). This evaluates when T38's original mean benchmark carries across to a correlated instrument. It does not assume that this benchmark is a weight in the mass law. [Aperture and tail theorem](e12-positive-aperture.md) · [T38](../definition/theorems/T38.md)

## A reflection bound for every spatial dimension

In noise coordinates the positive channel rays occupy the closed projective arc \(I_\beta=[0,\beta]\) of width \(\beta=\arccos(-\rho)\). Using the closure evaluates the infimum over positive masses. Write \(d_\pi\) for distance on the circle of circumference \(\pi\). Let the observed channel Gram matrix have eigenvalues \(l\le L\) and whole-line principal angle \(\psi\). The residual minimized over the arc is
\[
 Q_\beta(\psi)=l+(L-l)\sin^2 d_\pi(\psi,I_\beta).
\]
At each fixed true strength, the largest mean over the true angle occurs at an endpoint of the closed arc, by symmetric decreasing rearrangement of the conditional principal-angle density. The [aperture proof](e12-positive-aperture.md) gives the rearrangement for each superlevel event.

Put the true line at 0. In an orthonormal channel frame let \(X=r+\xi\) lie along this line and \(Y=\zeta\) across it. The noises are independent standard \(k\)-vectors. Reflection \(Y\mapsto-Y\) preserves the sampling distribution, keeps the two eigenvalues and sends the principal line to \(-\psi\). At \(\beta=\pi/3\) the following pointwise inequality holds:
\[
 \sin^2d_\pi(\psi,I_\beta)+\sin^2d_\pi(-\psi,I_\beta)\le \sin^2\psi .
\]
By reflection it suffices to take \(0\le\psi\le\pi/2\). On \(0\le\psi\le\pi/3\) the first term vanishes, and 0 is available to the second. On \(\pi/3\le\psi\le\pi/2\) the sum is \(\sin^2(\psi-\pi/3)+\sin^2(\psi+\pi/3)=1+\cos(2\psi)/2\), which is at most \((1-\cos(2\psi))/2\) since \(\cos(2\psi)\le-1/2\). Every wider arc containing this arc can only reduce the left side.

The residual at the edge itself is \(K(0)=\|Y\|^2=l+(L-l)\sin^2\psi\). The preceding inequality therefore gives
\[
 Q_\beta(\psi)+Q_\beta(-\psi)\le l+K(0).
\]
Here \(E K(0)=k\). Fitting an unrestricted scalar regression of \(Y\) on the observed \(X\) while keeping \(X\) fixed gives \(l\le\|Y\|^2-(e^\top Y)^2\) with \(e=X/\|X\|\), almost surely. Conditional on \(X\), the right side has expectation \(k-1\). Thus \(2 E Q_\beta\le2k-1\) at the edge, and hence at every allowed truth. At any fixed aperture \(0<\beta<\pi\), a diverging signal on an interior sequence approaching an edge faster than inverse strength attains the limiting mean \(k-\tfrac12\) by the tangent limit and uniform integrability in the [aperture proof](e12-positive-aperture.md). This proves the sharp supremum throughout \(\pi/3\le\beta<\pi\).

## A strictly larger scalar mean for every narrower aperture

Let \(k=1,\ \rho=-h\), with \(1/2<h<1\). Work in the original standardized channel coordinates, so the two observed scalars \(X,Y\) have unit variances and correlation \(-h\). Start at the boundary mean \((R,0)\), \(R>2\). The compatible readings fill the same-sign quadrants. Their least squared residual is zero there, and in the opposite-sign quadrants it is \(\min(X^2,Y^2)\) by the two distances to the axes. Compare \(Q_\beta\) to \(Y^2\,1_{\{Y<0\}}\), whose expectation is \(1/2\).

There is a gain on \(X\in[-2/R,-1/R]\), \(Y>1\) of at least \(1/R^2\). Throughout this strip \(X\)'s marginal density is at least \(\varphi(R+2/R)\). The conditional mean of \(Y\) is \(h(R-X)\ge h(R+1/R)\) and its conditional variance is \(1-h^2\). Thus the expected gain is at least
\[
 G_R=R^{-3}\varphi(R+2/R)
     \Phi\!\left(\frac{h(R+1/R)-1}{\sqrt{1-h^2}}\right).
\]
Every positive loss occurs on \(Y<0\) and \(T=X+Y\le0\). Indeed, if \(X\ge |Y|\) the baseline and the fitted residual coincide. The loss is therefore bounded above by \(E[Y^2 1_{\{T\le0\}}]\). With \(\tau^2=2(1+\rho)\), \(v=R/\tau\), one has \(T\sim N(R,\tau^2)\), conditional mean \(E[Y\mid T]=(T-R)/2\), and conditional variance \((1-\rho)/2\). Truncated Gaussian integration gives the bound
\[
 L_R=\Phi(-v)+\frac{1+\rho}{2}\,v\varphi(v).
\]
The gain has exponential rate \(\exp(-R^2/2)\) times a polynomial factor. The loss bound has exponential rate \(\exp[-R^2/(4(1+\rho))]\) times a polynomial factor. Since \(\rho<-1/2\), \(G_R>L_R\) for all sufficiently large finite \(R\). In particular \(E Q_\beta>\tfrac12\) there. For example \(\rho=-3/4,\ R=7\) already gives \(G_R>L_R\) by more than four orders of magnitude. This is a strict finite inequality, even though the overshoot may be small.

The cost in these coordinates is a continuous function with \(0\le Q_\beta\le X^2+Y^2\). Fix such an \(R\) and couple perturbations of the mean by adding the same noise. Dominated convergence as the second mean increases from zero to \(b>0\) transfers strictness to sufficiently small \(b\), where the true positive mass is finite. This establishes failure for every narrower scalar aperture without assuming an exactly zero or infinite true mass.

Consequently the mean claim at \(k=1\) holds exactly for \(\rho\ge-1/2\), whereas the right-angle chi-square mixture is an upper bound at **every threshold and signal** exactly for \(\rho\ge0\). For \(-1/2\le\rho<0\) the mean bound survives the failure of the tail bound. For narrower arcs this result gives a strict scalar counterexample and leaves the numerical value of the largest mean and the higher-dimensional classification open. [Finite identities, reflection checks and tail-integral checks](../tests/test_e12_aperture_mean.py)

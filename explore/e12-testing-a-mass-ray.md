# How the readings can test whether both masses share one scale

The common-scale theorem gives an observable equality. There are two ways to test it. The first assumes scalar weighing channels for the parts and total on the same preparation law. The second uses the original shared-acceleration instrument with separate force channels. It allows arbitrary unrecorded drive sizes and does not require moments of either mass. [Composition theorem](e12-mass-rays.md) · [Shared-acceleration experiment](e12-mass-observables.md)

## Three scalar predictions can certify one fixed allocation

Choose the proposed scales \(a,b>0\) in fixed units and \(t>0\). Write \(U=X/a,\ V=Y/b,\ w=a/(a+b)\). Three separate scalar predictions about the same underlying pair law give the quantity
\[
 D_t=w E e^{-tX/a}+(1-w) E e^{-tY/b}
             - E e^{-t(X+Y)/(a+b)}.
\]
For **arbitrary** positive marginals, \(D_t\ge0\), and \(D_t=0\) if and only if \(X/a=Y/b\) almost surely. This is pointwise strict Jensen, followed by averaging. It requires neither a shape assumption nor access to paired observations when the three population expectations are supplied. At equal scales the discrepancy is exactly half the expected squared difference of the two transformed scales, and in general it bounds that difference with coefficient \(\min(w,1-w)\). The same experiment in reciprocal mass uses scales \(1/a,1/b\) and the reciprocal of the reduced mass as its sum channel. This explains how three output laws can certify a particular joint coupling even though finitely many scalar predictions do not determine a general joint law. [General identification theorem](e12-mass-projections.md)

With centered Gaussian error independent of the latent quantity, of known variance \(\sigma^2\) in the chosen scalar coordinate, a reading \(\widetilde U=U+\epsilon\) gives \(E\exp[-t\widetilde U-t^2\sigma^2/2]=E e^{-tU}\). Apply the identity separately to the three channels using their own variances and scales. These corrected variables have finite second moments bounded by \(e^{t^2\sigma^2}\), since \(U>0\). This statement concerns an instrument that actually observes a scalar mass coordinate with additive Gaussian error. A noisy ratio of force and acceleration generally has a different error distribution.

## One bounded diagnostic on the shared force readings

Consider two labelled bodies pushed at a shared acceleration. On scalar push \(i\), let \(f_i=(F_{1i},F_{2i})^\top\) be their measured force components with conditionally centered Gaussian noise of known positive definite covariance \(\Sigma_{F,i}\), independent between pushes. The accelerometer need not be exact or independent of the forces. For a proposed allocation \(r\in(0,1)\), the force covector \(c=(1-r,-r)^\top\) annihilates the proposed allocation direction \((r,1-r)^\top\). Its noise variance is \(v_i=c^\top\Sigma_{F,i}c\). Given the true masses and drive, the standardized reading is
\[
 Z_i=\frac{c^\top f_i}{\sqrt{v_i}}\sim N(\delta_i,1),\qquad
 \delta_i=\frac{a_i[(1-r)X_i-rY_i]}{\sqrt{v_i}}.
\]
No latent-drive density enters. For a positive total, \(\delta_i=0\) means that either the drive vanishes or the allocation is precisely \(r\). The inputs can be fixed physical pairs and drives or conditioned latent values with the declared independent Gaussian errors.

**Theorem.** One bounded score has exactly calibrated mean under every such input, including weak and zero drives:
\[
 B_i=1-\sqrt2\,e^{-Z_i^2/2},\qquad
 E B_i=1-e^{-\delta_i^2/4}.
\]
Indeed integrating the product of the two Gaussian exponentials completes the square at \(\delta_i/2\). Hence the average \(\overline B\) is unbiased for \(D=n^{-1}\sum_i(1-e^{-\delta_i^2/4})\). In a mixture, the same identity applies after averaging over the preparation. A zero population discrepancy certifies the allocation on every active push up to sets of probability zero. It permits arbitrary allocations on zero-drive preparations. At an individual noisy observation the score may be negative.

The interval of possible \(B_i\) has length \(\sqrt2\). Hoeffding's inequality on independent noise therefore gives \(D\le U_\alpha=\min(1,\max(0,\overline B+\sqrt{\log(1/\alpha)/n}))\) with probability at least \(1-\alpha\), for fixed latent inputs. Consequently for any threshold \(h>0\) the fraction of these pushes with \(|\delta_i|\ge h\) is at most \(\min\{1,U_\alpha/(1-e^{-h^2/4})\}\) with the same confidence. The same proof bounds an average of unconditional discrepancies for independent draws from specified preparations. Under a prespecified common ray, \(\sum_iZ_i^2\sim\chi_n^2\) exactly at arbitrary fixed drives. Fitting \(r\) to these same data changes this null calibration. [Executable diagnostic](e12_mass_ray_diagnostic.py) · [Gaussian and calibration checks](../tests/test_e12_mass_ray_diagnostic.py) · [Hoeffding](https://doi.org/10.1080/01621459.1963.10500830)

For vector forces in \(k\) spatial dimensions, use \(C=[(1-r)I,-rI]\) and \(K=(Cf)^\top(C\Sigma_F C^\top)^{-1}(Cf)\). At the candidate ray this is the squared norm of whitened noise. At any pair and shared vector acceleration let \(\lambda\) be the squared noise-metric norm of the corresponding mean contrast. Applying the scalar Gaussian integral to the independent whitened coordinates gives
\[
 E\!\left[1-2^{k/2}e^{-K/2}\right]=1-e^{-\lambda/4}.
\]
Simultaneously rotating the force vectors and covariance leaves both \(K\) and \(\lambda\) unchanged. The finite-sample score measures *force separation in the actual drive*, so arbitrarily weak drives cannot furnish a uniform bound on the raw allocation error. This gives physical evidence about a ray even though the original measurement declarations do not furnish posterior odds between rays.

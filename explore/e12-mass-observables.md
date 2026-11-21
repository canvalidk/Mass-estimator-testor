# Which two-body measurements identify the mass allocation?

Write \(T=m_1+m_2>0,\ r=m_1/T\in(0,1)\). Fix a one-dimensional push with arbitrary deterministic signed true acceleration \(a\). With known positive definite error covariance \(\Sigma\), the Gaussian reading \(Y=(F_1,F_2,A)^\top\) has mean \(a(Tr,T(1-r),1)^\top\). Independent pushes may have different true accelerations and noise covariances. This is the same shared accelerometer on the two bodies. The questions here concern the *global* experiment and binary likelihood ratios at finite signal. [Local reciprocal geometry](e12-joining-bodies.md)

## How the instrument distinguishes directions

More generally, suppose push \(i\) measures the force covectors in the rows of a known matrix \(H_i\) on body index space, together with its acceleration. Its mean is \(a_i(H_i m,1)^\top\). Two choices of masses and drives yield identical data laws exactly when their mean vectors match, since the covariances are declared. The acceleration coordinates force equality of their \(a_i\); equality of the force coordinates then requires \(H_i(m-m')=0\) on every push with \(a_i\ne0\). The set of indistinguishable masses is
\[
 (m+\bigcap_{i:a_i\ne0}\ker H_i)\cap(0,\infty)^2 .
\]
In particular summed force and acceleration alone have exactly the same sampling law for all allocations at a fixed \(T\) on every push. With two separately measured forces and one active push the mass pair is identifiable as a parameter of the sampling law. Identification from the law leaves finite noisy observations uncertain.

There is a quantitative version. Set \(g_i(m)=\Sigma_i^{-1/2}(H_i m,1)^\top\) in any whitening, and let \(P_i(m')\) project onto the span of \(g_i(m')\). Against an alternative mass pair with all its own signed drives allowed, the minimum of twice the Gaussian relative entropy from a fixed true experiment is
\[
 D=\sum_i a_i^2\|(I-P_i(m'))g_i(m)\|^2.
\]
Indeed twice the relative entropy of equal-covariance Gaussians is the squared whitened mean distance, and minimizing each distance fits the alternative drive by orthogonal projection. The closest simple alternative has optimal equal-prior binary error \(\Phi(-\sqrt D/2)\) by thresholding the whitened mean difference at its midpoint. This gives a lower bound for worst-case error over the two parameter families if both contain this pair. At \(m=(1,3),m'=(2,2)\) with one push \(a=1\), independent unit noise and separate force channels, the closest alternative has drive 1 and \(D=2\). A total-force-only sensor has \(D=0\). Both conclusions use the declared acquisition protocol.

## Two residuals and the shared noise between them

Let \(S=F_1+F_2,\ B=F_1\). At candidates \(T,r\) form \(U=S-TA,\ V=B-rS\). At the true total, \(U\) is a centered Gaussian regardless of the allocation or drive. At the true allocation, \(V\) is centered regardless of the total or drive. Thus \(V^2/\operatorname{Var}V\) is an exact \(\chi^2_1\) candidate statistic. It is also the T36 ratio residual for the two force channels at \(m_1/m_2=r/(1-r)\), since \((1-r)F_1-rF_2=V\). Full spatial readings use the residual vector and its noise-metric squared norm. This identifies a genuine way to examine the allocation without deciding which mass interval to put on the total. [T36](../definition/theorems/T36.md) · [Fieller's ratio tests](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x)

Denote the variances of \(U,V\) by \(p,q\), and their covariance by \(c\). The two contrast covectors annihilate the candidate graph and are linearly independent. Therefore the minimum squared Gaussian distance to that graph, profiled over the unknown drive, is exactly
\[
 K(T,r)=\frac{U^2}{p}+\frac{(V-cU/p)^2}{q-c^2/p}.
\]
Indeed the shortest correction with the two required residual components has squared norm \((U,V)\operatorname{Cov}(U,V)^{-1}(U,V)^\top\), and completing the square gives the result. The two standardized innovations are independent normals at joint truth. Reversing the order of the conditioning gives the same \(K\). Simply adding \(U^2/p+V^2/q\) gives the correct cost on all data iff \(c=0\). In these coordinates
\[
 c=\operatorname{Cov}(S,B)-r\operatorname{Var}(S)-T\operatorname{Cov}(A,B)+Tr\operatorname{Cov}(A,S).
\]
Positive definiteness makes \(\operatorname{Var}(S)>0\), so \(c\) cannot vanish on the whole candidate rectangle for any such instrument. Even with three independent unit-variance errors, \(c=1-2r\). Each marginal question can have an exact candidate test while counting their separate mass evidence as independent gives a wrong joint geometry.

## Does an unknown drive allow a unique Bayesian update between two rays?

Even give the total mass \(T\) exactly. For two candidate fractions \(r,r'\), write \(g=(Tr,T(1-r),1)^\top,\ g'=(Tr',T(1-r'),1)^\top,\ P=\Sigma^{-1}\), and
\[
 u=(g-g')^\top P y,\qquad v=g^\top P g-g'^\top P g'.
\]
With the drive fixed at \(a\), the log likelihood ratio at observation \(y\) favoring \(r\) over \(r'\) is exactly \(a u-a^2 v/2\). If the common drive is assigned any proper law \(\rho\) independently of the choice of ray, the marginal likelihood ratio is the ratio of integrals of the two Gaussian likelihoods against \(\rho\). Equivalently this averages \(\exp(a u-a^2v/2)\) against the conditional drive law under \(r'\).

**Criterion.** A ratio invariant under every common drive law supported on an open drive interval exists at that record iff \(u=v=0\). Its value is then one. Indeed point masses test the likelihood ratio at every drive, and the displayed polynomial is constant on an open interval only when both coefficients vanish. Smooth proper densities concentrated near each test point give the same necessity by continuity. When every positive drive is allowed, some common drive laws favor \(r\) and others favor \(r'\) iff \(uv>0\), since this is exactly the condition for a strict sign change on \((0,\infty)\). With both drive signs allowed a preference reversal occurs whenever \(u\ne0\).

For \(T=2,\ r=1/4,\ r'=1/2,\ y=(1,2,1)\) and three independent unit errors, \(u=v=1/2\). Equal prior odds on the two fractions become odds \(e^{1/4}\) at known drive 1 and \(e^{-3/4}\) at known drive 3. The observed accelerometer value is identical. Proper smooth common drive distributions can preserve either strict ordering. The calculation shows exactly what is missing for a *drive-prior-independent Bayesian* finite update. The project's unrestricted deterministic drives supply exact residual tests without declaring such a population; profiling or another physically justified construction is a different inference premise. [No-borrowing and no-size-law premises](../definition/premises.md)

## Sharp odds bounds as the drive becomes known

The preceding ambiguity has sharp bounds when the apparatus specifies a drive range \(\mathcal A\). For the fixed observation, let \(R(a)=\exp(au-a^2v/2)\). The lower and upper bounds on the marginal Bayes factor over common proper nuisance laws are precisely \(\inf_{\mathcal A}R\) and \(\sup_{\mathcal A}R\). The averaging identity proves containment. Every value strictly between those bounds is attainable: choose two drives with ratios on either side and form their convex average with weights \(w_1,w_2\) giving that value. Choose the corresponding prior masses proportional to \(w_j/L_{r'}(y\mid a_j)\). They are finite and positive and yield exactly those averaging weights. Endpoints are attained precisely when the corresponding extrema are attained by allowed drives. In the numerical example with \(a>0\) the sharp interval is \((0,e^{1/4}]\); the maximum occurs at \(a=1\). As a calibrated drive range narrows to a known \(a_0\), both bounds tend to the same ordinary Gaussian factor \(R(a_0)\). These are sensitivity bounds on a Bayesian completion when the nuisance population is unspecified.

An independent Gaussian accelerometer with variance tending to zero gives the same pointwise limit under any *fixed* proper drive density continuous and positive at the observed \(A\). The likelihood integral for either hypothesis is the convolution of that density times the continuous bounded force likelihood with the shrinking Gaussian noise kernel, evaluated at \(A\). It tends to the product of their values at \(A\); the common drive density cancels in the ratio. This step does not give a bound uniform over drive laws that may change with the instrument precision.

## A force-controlled experiment that supplies the kernel

Suppose the apparatus applies a **known total force** \(s_i\) to the common pair on each scalar push. Then \(a_i=s_i/T\). In \((r,\nu)\), \(\nu=1/T\), the three observed means are \((s_i r,s_i(1-r),s_i\nu)^\top\). This is an ordinary linear Gaussian experiment on \((0,1)\times(0,\infty)\). Set \(e=(1,-1,0)^\top,\ f=(0,0,1)^\top,\ D_i=s_i[e,f],\ y_i^0=y_i-(0,s_i,0)^\top\). The exact likelihood kernel for \(\beta=(r,\nu)^\top\) is
\[
 L(\beta)\propto \exp\{\beta^\top z-\beta^\top J\beta/2\},\quad
 z=\sum_i D_i^\top\Sigma_i^{-1} y_i^0,\quad J=\sum_iD_i^\top\Sigma_i^{-1}D_i.
\]
Every nonzero known load gives \(J>0\). Under Bayesian conditioning, any declared prior on the two masses pushes forward to a prior on \((r,\nu)\); multiplying it by this kernel and normalizing uniquely determines the finite law. The likelihood factors into separate functions of \(r\) and \(\nu\) exactly when \(J_{12}=\sum_i s_i^2 e^\top\Sigma_i^{-1}f=0\), by the mixed derivative of its logarithm. An independent accelerometer and force subsystem meet this criterion, even with correlations between the forces. A coupled prior can still couple the two posterior coordinates. A total force *reading* with positive error does not declare the physical setting \(s_i\); substituting that reading for the setting changes the experiment. Here the physical control, the Gaussian law and the explicit conditioning premise jointly supply the update. The rule for the original uncontrolled pushes, given the declared mass prior, remains open. [Sensor calibration comparison](e12-known-mass-sensor.md) · [Finite checks](../tests/test_e12_mass_observables.py)

For example take two equally likely labelled pairs \((1,3),(3,1)\), delivered force 4, and three independent unit errors. The means are \((1,3,1)\) and \((3,1,1)\). A record at the first mean changes its probability from \(1/2\) to \(1/(1+e^{-4})\), about \(0.9820\). Summing the force channels before acquisition makes both predicted joint means of summed force and acceleration identical, so the same two candidates retain equal odds on that instrument. More generally the log likelihood ratio under the first of two candidates has law \(N(d^2/2,d^2)\) where \(d^2\) is their squared joint whitened mean distance. Expanding the equal-covariance Gaussian likelihoods proves the statement. With positive initial weights on a finite candidate set, sufficiently many independent pushes concentrate the posterior on the true atom if its squared distance from each competitor tends to infinity. [Finite calculator](e12_mass_measurement.py) · [Numerical checks](../tests/test_e12_controlled_force_update.py)

# Can the E12 update be a measurement whose unread outcomes average to the prior?

The proposed projection analogy has a test even before introducing quantum phases. If a measurement reads a persistent mass and its outcome is then forgotten, the posterior averaged over the experiment's predictive outcome probabilities returns the starting law. This is the tower property of conditional probability. The question here is whether one mass-only update can obey even this **unconditional mean** identity when the unknown pushes in the Gaussian instrument are fixed and arbitrary. No posterior reliability conditional on the reported answer is assumed below. Finite-sample coverage of a confidence set is a different requirement. [Physical reading model](../definition/theorems/T36.md) · [What calibrated joint conditioning requires](e12-conditioning-nuisance.md)

## A sharp obstruction using the actual half-Cauchy priors

Use one fixed declared Gaussian instrument with positive definite noise. For a list of \(n\) readings let \(Q_n\) be its observation distribution when all true pushes vanish. It is identical at every mass. It has a strictly positive density on the same data space as every fixed-signal experiment at that \(n\). The measured observations are random and almost surely nonzero. L20's requirement on the observed all-zero list does not say what their average posterior must be.

Suppose a rule has a positive finite data factor \(L_y(m)\) multiplying the supplied prior, with a finite normalizer. The factor may depend on the whole list and all declared geometry, and is held fixed when only the starting law changes. Consider the additional calibration condition
\[
             E_{Q_n} P_\pi(B\mid Y)=\pi(B)
             \quad\text{for all Borel mass events }B .
\tag{1}
\]
This would hold for a Bayesian readout required to average correctly in each separately specified fixed-drive preparation, including zero drive, **using the same readout even when the preparation is hidden**. It can fail for a correctly constructed Bayesian posterior averaged under a single unannounced submodel instead of its own predictive distribution. The stronger requirement is being tested as a candidate measurement premise.

**Two-prior lemma.** Let \(d\pi'=h\,d\pi\), where \(h>0,\ \int h\,d\pi=1\). If (1) holds for both priors under the same outcome distribution, then
\[
              Z_y=\int h(m)\,P_\pi(dm\mid y)=1
              \quad Q_n\text{-almost surely}.
\]
Indeed, changing the prior gives \(dP_{\pi'}=h\,dP_\pi/Z_y\). Applying (1) to \(h\) under the first prior gives \(E Z_Y=1\). Applying it to \(1/h\) under the second gives \(E(1/Z_Y)=1\), since \(\int (1/h)\,dP_{\pi'}=1/Z_y\). For unbounded functions use monotone truncation. Now
\[
       E[(Z_Y-1)^2/Z_Y]=E Z_Y+E(1/Z_Y)-2=0 .
\]
If the two mean identities hold to absolute errors \(\delta,\delta'\), the left side is at most \(\delta+\delta'\). For a fixed whitened Gaussian signal \(\mu\), suppose \(Z\ge1+t\) on an event of probability at least \(p\) under that signal, with \(t>0\). Its probability under \(Q_n\) is at least \(p^2e^{-\|\mu\|^2}\), by Cauchy–Schwarz and \(E_{Q_n}(dP_\mu/dQ_n)^2=e^{\|\mu\|^2}\). Hence the two errors sum to at least \(p^2 e^{-\|\mu\|^2}t^2/(1+t)\). The finite-signal lower bound weakens with separation from zero.

For the project's positive half-Cauchy starting laws with distinct centres \(c,s\),
\[
              h(m)=\frac{s}{c}\frac{c^2+m^2}{s^2+m^2}.
\]
It is bounded, continuous and strictly monotone on positive mass, with \(h(m)=1\) exactly at \(m=\sqrt{cs}\). The almost-sure constraint transfers from \(Q_n\) to *every* fixed-signal Gaussian experiment at this list length by absolute continuity. If (1) holds at each length of an identifying sequence and the returned law concentrates in probability at \(m_0\), then \(\int h\,dP_\pi\to h(m_0)\) in probability. It follows that **the only possible concentration point is \(\sqrt{cs}\)**. This is a necessary condition, with no claim that learning even there is sufficient. With centres 1 and 4, for example, \(h(3)=1.6\) while the proposed axiom requires its posterior expectation to be exactly 1 on almost every dataset at every signal strength. This contradicts concentration at true mass 3. No independence of the individual evidence factors, asymptotic expansion, special choice of discount, or two-point prior is used. The starting priors can differ even when the instrument is unchanged because independent old information is a declared input. [Prior family](../problems/O6.md) · [Tests](../tests/test_e12_forgetting.py)

Two priors can permit genuine changes on finite data. Inversion \(m\mapsto cs/m\) interchanges the half-Cauchy laws of centres \(c\) and \(s\). The invariant band \(A=[a,cs/a]\), \(0<a<\sqrt{cs}\), therefore has the same probability \(p\) under both. Let a fair independent sign select \(L_\pm=1\pm\epsilon({\bf1}_A-p)\) with \(0<\epsilon<1\). Both normalizers equal 1 and both priors satisfy (1), even though either signed update changes the probability of \(A\).

**Three centres force a lower bound on uncertainty on almost every dataset.** Take \(c=1\) and the two comparison centres 2 and 4. Under the scale-1 posterior \(\mu\), put \(U=1/(4+m^2)\). The two required ratio means give \(E_\mu U=1/6\) and \(E_\mu[U/(1+12U)]=1/20\). Define the finite measure \(d\nu=d\mu/(1+12U)\) and \(v_j=\int U^j\,d\nu\). Thus \(v_1=1/20\) and \(v_2=(1/6-v_1)/12=7/720\). Cauchy–Schwarz gives \(v_3\ge v_2^2/v_1\). Consequently
\[
 \operatorname{Var}_\mu(4U)=16(v_2+12v_3-1/36)\ \ge\ 2/27.
\]
The corresponding moment bound is attained on the compactified mass space by weights \(6/7\) at \(m^2=8/7\) and \(1/7\) at infinity. This is a lower bound under the added mean calibration requirement, independent of how strongly an observed block favors one interior mass.

If (1) holds for a countable set of scales with an accumulation point in \((0,\infty)\), the same argument forces the *entire* posterior to equal the prior almost surely. Fix one centre \(c\) in the family and intersect the full-measure sets over the other centres \(s\). From the displayed \(h\) identity, the posterior and prior of \(X=m^2\) have the same values of \(\int(X+s^2)^{-1}dP\). These Stieltjes transforms are analytic in positive argument and hence agree everywhere there. For a direct uniqueness argument differentiate at argument 1 to recover every positive integer moment of \(1/(1+X)\). Both are probability laws on \([0,1]\) in this coordinate; polynomial approximation identifies the measures. If the data factor is continuous in the data as a law in the weak topology, the equality extends to every observation by full support.

## One starting prior if each recorded update can itself be repeated

There is a separate obstruction under the stronger operational requirement that the normalized law alone can be used as the prior for an identical new reading. Start with one \(\pi\). At a true-null preparation draw two independent outcomes \(Y,Z\sim Q_1\). Write \(f_y=dP_\pi(\cdot\mid y)/d\pi\), where the normalized densities are positive and jointly measurable. Suppose resetting the prior after seeing \(y\) uses the same new-reading factor \(L_z\). Require mean preservation on forgetting the new outcome, both at the initial prior and at almost every reached prior \(f_y\pi\). The two-prior lemma applied conditionally for fixed \(y\) gives
\[
                  \int f_y f_z\,d\pi=1
                    \quad (Q_1\otimes Q_1)\text{-almost surely}.
\]
For any outcome event \(A\) define \(g_A(m)=\int_A f_y(m)Q_1(dy)\). Integration of nonnegative functions in either order gives \(\int g_A\,d\pi=Q_1(A)\). Integrating the displayed identity over \(A\times A\) gives \(\int g_A^2d\pi=Q_1(A)^2\). Thus \(g_A=Q_1(A)\) almost surely under \(\pi\). For any mass event \(B\), integrating over \(B\) yields \(E[{\bf1}_A P_\pi(B\mid Y)]=Q_1(A)\pi(B)\) for every \(A\). It follows that \(P_\pi(B\mid Y)=\pi(B)\) almost surely. Use a countable determining family of mass intervals to obtain equality of the entire laws. Absolute continuity transfers this equality to every fixed push. Reusing the unchanged law yields no learning on independent repetitions. The assumption here is calibration at the reached priors as well as the initial one. T12's whole-block recomputation can keep other state. [Combination rule](../definition/theorems/T12.md)

A binary calibration fixture exhibits the distinction at one step. Start with equal probabilities of two masses, send a fair outcome sign, and for the positive sign give the first mass likelihood ratio 3 versus the second, using \(1/3\) for the negative sign. The posterior weights \(3/4,1/4\) average to \(1/2\). Starting the same update from the reached weight \(3/4\), the new weights \(9/10,1/2\) average to \(7/10\). Forgetting this second, equally uninformative sign changes the mean by \(1/20\). The first mean identity alone permits informative-looking weights.

## What the obstruction selects

On a declared calibrated experiment the predictive law genuinely depends on which preparation was made. Averaging at that predictive law gives ordinary conditional probability, including an exact projection for a sharp measurement. On the original force–acceleration instrument the rule receives covariance and arbitrary readings, but no distribution of latent push sizes. The extra demand of (1) at a common hidden zero drive cannot force a useful universal likelihood. Under two ordinary starting priors it defeats identifying consistency. The same conclusion applies if zero drive is only approached **at each finite list length**, the readout is held fixed, and (1) is required along a sequence of predictive preparations converging in total variation to \(Q_n\); pass to the limit on each bounded mass event. Approximate calibration on identified large-signal data entails no such demand.

This gives a difference between updating from a fully calibrated sensor, learning from a working mass-evidence factor, neutral treatment of a measured zero, and exact confidence-test calibration. Each is a property of a specified experiment. [Controlled-force example](e12-mass-observables.md) · [Concrete obstruction to representing an E12 factor by a proper Gaussian push mixture](e12-conditioning-nuisance.md)

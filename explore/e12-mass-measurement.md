# What a persistent measurement can do to a superposition of mass laws

On a shared nonzero exact acceleration \(\mathbf a\), let \(\mathbf a^\vee=\langle\mathbf a,\cdot\rangle/\|\mathbf a\|^2\) in the spatial metric. The tensor \(\mathbf a\otimes\mathbf a^\vee\) projects onto its direction. Contracting both exact force vectors with this **one** covector gives the labelled pair \(x=(X,Y)\). Contraction commutes with force addition and with a simultaneous spatial rotation of forces, acceleration and metric. It supplies the linear mass coordinates used below. When the measured acceleration and forces are noisy, substituting the observations into this formula does not give an outcome probability. [Noise-metric reciprocal for shared acceleration](e12-joining-bodies.md)

Every joint uncertainty law decomposes into components on fixed allocation rays. The two physical operations close on each ray, and there is a converse for any common uncertainty shape. A measurement now asks which components to retain and in what proportions. [Common-shape theorem](e12-mass-rays.md)

## From a persistent body to the conditional law

Let \(\mathcal X=(0,\infty)^2\) be the space of labelled masses and \(\mathcal Y\) a standard Borel outcome space. Consider an instrument with a probability kernel \(J(dy,dx'\mid x)\) giving its outcome and outgoing state for an incoming state \(x\). A preparation law \(\pi\) enters by integration against \(J\); this includes ordinary mixing of preparations. Write \(\kappa(dy\mid x)=J(dy,\mathcal X\mid x)\) for the response. Call an observable \(g\) **preserved** if, after forgetting the outcome, its outgoing distribution equals its incoming distribution for every preparation, including an exact pair. In a loaded-standard experiment \(x\) describes the original two bodies, with the known attachment accounted for at readout.

**Theorem.** Given the response, the possible classical instruments preserving finitely many real observables \(g_i\) are precisely the kernels which, conditional on \(x,y\), can redistribute mass pairs only among \(x'\) with every \(g_i(x')=g_i(x)\). If the combined observables separate mass pairs, the instrument is uniquely
\[
             J(dy,dx'\mid x)=\kappa(dy\mid x)\delta_x(dx').
\]

**Proof.** Apply preservation to an exact incoming \(x\). The outgoing law of \(g_i\) is \(\delta_{g_i(x)}\). Each equality has probability one, and their finite intersection does too. Disintegrating over the outcome yields the stated support restriction for \(\kappa(\cdot\mid x)\)-almost all outcomes. Conversely such a kernel preserves each value pointwise and hence its law for every preparation. A separating family leaves only \(x'=x\). This also specifies the outgoing mass law conditioned on any outcome event of positive probability,
\[
        P_\pi(B\mid y\in A)
          =\frac{\int_B\kappa(A\mid x)\pi(dx)}
                 {\int_{\mathcal X}\kappa(A\mid x)\pi(dx)} .
\]
An outcome density gives the analogous pointwise likelihood formula at predictive-almost-every observation. Preserving only the ensemble law at a single diffuse prior would allow, for example, redrawing the pair from that prior.

**Which quantities separate?** A total \(T=X+Y\) and a reduced mass \(H=XY/T\) determine the roots of \(z^2-Tz+TH=0\), including multiplicity. Both are symmetric, so a swap of the two labels is allowed. Either a labelled mass and the total, or two distinct loaded reciprocal readouts, separate all pairs. In the latter experiment use
\[
 W_b=\frac1{X+b}+\frac1Y,\qquad
 W_b-W_d=\frac{d-b}{(X+b)(X+d)},\quad 0\le b<d.
\]
The difference strictly decreases with \(X\), and then \(W_b\) fixes \(Y\). These are two settings at which the original masses can be checked. Preservation at every exact state is the premise. Two unpaired histograms across an ensemble do not supply an individual pair of readings. [Exact identification and its limits](e12-mass-projections.md)

**Sharp and noisy readings.** For the event \(g(x)\in A\), the response is the indicator of its inverse image. On square-root densities \(\psi=\sqrt{d\pi/d\nu}\), updating means multiplying by that indicator and normalizing in \(L^2(\nu)\). A sharp continuous value uses a regular conditional law, which can be singular relative to \(\nu\). For an ordinary noisy detector, a response density \(\kappa(y\mid x)\) gives the multiplier \(\sqrt{\kappa(y\mid x)}\). A total-only response preserves the original conditional allocation given \(T\); marginal allocation weights may move because \(r\) and \(T\) may be dependent. Independent detector outcomes conditional on the unchanged pair give commuting products of likelihoods. A common unobserved drive can correlate the readings even after conditioning on the masses. In that case the joint response is needed.

## Could amplitudes over different mass allocations interfere?

Test a precise additional hypothesis. Represent states on a finite labelled mass set, or on \(L^2(\nu)\) for a probability reference on the positive quadrant with one complex coordinate per mass pair. Allow a quantum instrument with finite or countably many Kraus operators \(K_{yj}\) and outcomes, satisfying \(\sum_{y,j} K_{yj}^*K_{yj}=I\). The probability of \(y\) at a unit vector \(\psi\) is \(\sum_j\|K_{yj}\psi\|^2\). Require the unread instrument to preserve the **whole distributions**, in every input state, of the two labelled standards above.

**Conclusion.** Each \(K_{yj}\) is multiplication by some \(k_{yj}(x)\). All the measurable outcome probabilities and conditional mass laws then depend only on the starting mass probabilities. This holds through any finite sequence of such instruments, including adaptive choice of standard based on earlier outcomes.

**Proof.** For any spectral projection \(P\) of either protected observable, preservation in all states gives \(\sum K^*PK=P\). For \(v\) in the kernel of \(P\), the sum of \(\|PKv\|^2\) is zero, so each \(K\) preserves that kernel. Apply the same argument to \(I-P\) for the range. Thus every \(K\) commutes with \(P\). Two separating observables with the explicit measurable inverse above generate the full mass-event sigma algebra. Commutation extends from their spectral projections to intersections, complements and countable disjoint unions by strong limits. In the finite case this makes every \(K\) diagonal. In the continuous case put \(k=K1\). Commutation gives \(K{\bf1}_B=k{\bf1}_B\) for every event \(B\). Boundedness forces \(|k|\le\|K\|\) almost everywhere by testing on its larger level sets. Simple-function approximation gives \(K\psi=k\psi\). Consequently \(\kappa_y(x)=\sum_j|k_{yj}(x)|^2\), and the outgoing mass-event probabilities at \(y\) are integrals of \(\kappa_y(x)|\psi(x)|^2\), divided by their total. The same argument applies to mixed states by linearity, and to sequences by conditioning at each outcome. Phase changes of an outgoing vector may still occur, but all subsequent instruments satisfying the same preservation condition have the displayed mass probabilities. This specializes the established sharp-observable nondisturbance argument to the project's reciprocal standards. [Heinosaari and Wolf](https://arxiv.org/abs/1005.5659)

**What would reveal a relative phase?** Take \( |0\rangle=|1,3\rangle,\ |1\rangle=|3,1\rangle\). They have identical \(T=4,H=3/4\). On this two-dimensional subspace a measurement in the basis \(|\pm\rangle=(|0\rangle\pm|1\rangle)/\sqrt2\) preserves those two symmetric quantities. Its probability for \(+\) is 1 for the plus state, 0 for the minus state, and \(1/2\) for the equal ordinary mixture. Yet \(W_1\) takes values \(5/6\) and \(5/4\) on the two labelled pairs. Starting from either known pair, ignoring the sign changes its subsequent \(W_1\) probabilities to \(1/2,1/2\). Thus a specific hypothetical measurement distinguishes relative phase, with a specific physical check it must disturb. The Newtonian reading model supplies no such sign detector. T70's phase already has a distinct operational meaning: it is the spatial angle between force and acceleration *inside an observed reading*. [T70 reading qubit](../definition/theorems/T70.md)

## A sharp bound if the labels are only approximately preserved

Use any quantum instrument with outcomes \(y\), followed by a binary ideal tag test \(P_0,P_1\) on its outgoing state, with \(P_0+P_1=I\). Label 0 has the tag error \(e_0=\Pr(1\mid0)\), and label 1 the error \(e_1=\Pr(0\mid1)\). Let \(p_\pm\) be the instrument's outcome laws for the two balanced, opposite-phase preparations. Then
\[
 \|p_+-p_-\|_{\rm TV}
 \le \sqrt{(1-e_0)e_1}+\sqrt{e_0(1-e_1)}
 \le 2\sqrt{\epsilon(1-\epsilon)}
 \quad\text{if }e_0,e_1\le\epsilon\le\tfrac12.
\]
Indeed the total variation is \(\sum_y|\operatorname{Re}\sum_j\langle K_{yj}0,K_{yj}1\rangle|\), where the two kets are abbreviated by their labels. Insert \(P_0+P_1\) into each inner product. The triangle inequality followed by Cauchy–Schwarz over \(y,j\) on each tag gives the first bound. For the second put \(e_i=\sin^2\alpha_i\); the sum is \(\sin(\alpha_0+\alpha_1)\) in the required range.

The first bound is attained for every pair of errors by the two real Kraus matrices
\[
 K_\pm=2^{-1/2}
 \begin{pmatrix}\sqrt{1-e_0}&\pm\sqrt{e_1}\\
                 \sqrt{e_0}&\pm\sqrt{1-e_1}\end{pmatrix}
\]
with the labelled output basis as tag. Their signed cross terms cancel in \(\sum K^*K\); inserting either balanced input gives the stated variation. The tag errors measure physical change of a labelled state under the hypothetical operation. Converting a Gaussian readout error into such an error would require an additional experiment. The related which-way/visibility tradeoff is established in interferometry. [Englert](https://doi.org/10.1103/PhysRevLett.77.2154) · [Checks of both measurement models](../tests/test_e12_mass_measurement.py)

## What determines the missing response

For two bodies given a known *delivered* total force \(s\) at shared acceleration, Newton's law fixes the three noiseless channels to \((sr,s(1-r),s/T)\). A declared Gaussian sensor covariance then determines \(\kappa\); persistence fixes the posterior given a starting joint law. An arbitrary unrecorded drive gives a family of sampling models. At fixed total \(T=2\), the same three readings can favor fraction \(1/4\) over \(1/2\) at one allowed drive and reverse those odds at another. A calibrated covariance alone does not select which evidence to use. [Gaussian experiment and exact odds](e12-mass-observables.md) · [Discrete implementation](e12_mass_measurement.py)

There is a direct test on the original one-body instrument. Any proper distribution for the unrestricted Gaussian push, even a different one for each candidate mass, forces an additive separation identity for four log odds. The two E12 candidates violate it with four positive readings. Requiring a mass-only update to average back to the starting law at each hidden preparation has a separate obstruction already at two ordinary half-Cauchy starting scales. Thus the fixed-ray representation and persistence classify the update **once an outcome response is specified**. For an arbitrary push the construction still owes a justified mass-evidence rule and its interpretation. [Four-reading calculation](e12-conditioning-nuisance.md) · [Two-prior proof](e12-forgetting-a-reading.md)

A bounded contrast of the two labelled forces tests a proposed allocation at arbitrary fixed drives, with exact mean and a finite-sample error bound. There is also an identifying experiment on the original Gaussian channels if independent draws from the same preparation are available: the full predictive observation law determines all active mass allocations, using the reciprocal of the deconvolved acceleration, and isolates the precise ambiguity at zero drive. [Testing a mass ray](e12-testing-a-mass-ray.md) · [Unmixing a repeatable preparation](e12-unmixing-readings.md)

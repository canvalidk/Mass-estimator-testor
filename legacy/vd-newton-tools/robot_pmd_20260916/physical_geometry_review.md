# Physical geometry in robot inertial identification

16 September 2026. Bounded research contribution; scratch only. No novelty or
performance claim for the VD estimator follows from this review.

## Sources and coverage

**A. Wensing, Kim, Slotine, 2018**, *Linear Matrix Inequalities for
Physically-Consistent Inertial Parameter Identification: A Statistical
Perspective on the Mass Distribution*.
[Author arXiv v3, eight pages](https://arxiv.org/pdf/1701.04395v3);
[DOI](https://doi.org/10.1109/LRA.2017.2729659).
Read the complete main text, Sections I–VII; no independent reproduction of
the experiment or supplement. PDF page numbers below are preprint pages.

### Factual synopsis A — under 180 words

The paper characterizes realizable rigid-body inertial parameters through a
positive pseudo-inertia matrix. Positive mass and positive rotational inertia
alone are insufficient; the principal inertias must satisfy triangle
inequalities. Its covariance is the spatial, density-weighted second central
moment of the body's mass, not sensor-error covariance. Moment realizability
also supplies convex support constraints for bounding ellipsoids and their
unions. The robot example fits link and rotor inertias plus friction with
least-squares residuals, a small CAD-centered quadratic regularizer, and
physical constraints. The conclusion explicitly distinguishes these
deterministic realizability results from stochastic inference using uncertain
measurements. Locators: physical consistency, p.3 Eq.(9); spatial covariance,
p.4 Eqs.(15)–(20); pseudo-inertia, p.5 Definition 7/Theorem 3; shape constraints,
p.6 Theorem 4/Eq.(23); identification, pp.6–7 Eq.(25); distinction from
measurement stochasticity, p.8 conclusion.

**B. Lee, Wensing, Park, 2020**, *Geometric Robot Dynamic Identification:
A Convex Programming Approach*.
[DOI](https://doi.org/10.1109/TRO.2019.2926491);
[author-uploaded full preprint](https://www.researchgate.net/publication/337753338_Geometric_Robot_Dynamic_Identification_A_Convex_Programming_Approach);
[SNU author/laboratory abstract](https://sites.google.com/robotics.snu.ac.kr/fcp/publications/articles/abstract/2011-2020/2019_geometric-robot-dynamic);
[author-linked reference implementation](https://github.com/alex07143/Geometric-Robot-DynID).
The full-text page identifies Patrick Wensing as uploader on 23 January 2020;
the document is marked author preprint, nineteen pages. The publisher page was
bot-blocked. Read from the actual full preprint, focusing on Sections II–IV,
VI–VII and Appendix A; experimental descriptions were inspected rather than
reproduced. Locators below use preprint numbering, not journal pagination.

### Factual synopsis B — under 180 words

The fit minimizes a covariance-weighted inverse-dynamics residual plus a
nominal-parameter regularizer (p.2 Eq.(3)). Residual covariance comes from an
unregularized fit; a non-diagonal covariance can be diagonalized without
changing the objective (p.7 footnote 4). Section II-C, p.4, explains a Gaussian
likelihood conditional on the observed regressor and a MAP interpretation.
The geometric contribution is a coordinate-invariant LogDet divergence
between pseudo-inertias (p.5 Eqs.(19)–(22)), yielding convex regularization
(p.7 Eqs.(27)–(30)). It penalizes approach to degenerate inertias. The spatial
second moment in the pseudo-inertia remains distinct from residual covariance.
At unit mass, Appendix A Proposition 3, p.16, identifies the divergence as
twice Gaussian KL divergence for spatial mass distributions. Section IV,
p.8, compares regularizers at matched training residuals rather than a shared
numerical regularization coefficient. The formulation is parameter
optimization; it is not integration over uncertain force and acceleration
channels followed by a ratio of magnitude means.

## Our comparison and deductions

The following applies the source structures to our problem. It is our analysis,
not a claimed theorem or proposed VD extension from either paper.

### 1. Three mathematically different objects must remain separate

| Object | Meaning in a possible extension |
|---|---|
| Spatial mass density `rho(r)` | Where the object's material is located |
| Measurement-error law | How force, torque, position, velocity, and acceleration readings differ from their true values |
| Probability law over candidate inertias and latent motion | Our uncertainty about the object and experiment after imposing mechanics |

Normalizing a spatial mass density by its total mass creates a probability
distribution over locations. That mathematical device does not create a
probability distribution over unknown masses or sensor errors. A covariance
matrix appearing in a mechanics paper therefore needs its random variable
specified before comparison with our input covariance.

In our existing scalar estimator, the declared covariance belongs to the joint
measurement vector `(F_observed,a_observed)`, and the conditional law is over
latent compatible pairs. Merely substituting the spatial covariance matrix for
that measurement covariance would be a category error.

### 2. The cone supplies necessary extra structure for a tensor extension

Write the object's first moment as `h=m c`, its second raw spatial moment as
`S`, and define

\[
J=\begin{pmatrix}S&h\\h^T&m\end{pmatrix},\qquad
I_O=\operatorname{tr}(S)I_3-S.
\]

The Schur complement of the mass entry is

\[
S-hh^T/m=\int(r-c)(r-c)^T\rho(r)\,dr.
\]

This directly explains the constraint: a negative quadratic form on the left
would assert a negative integral of squared spatial distances. For a
nondegenerate three-dimensional mass distribution, use `J>0`; allowing ideal
plates, rods, or point masses requires the appropriate semidefinite boundary
with `m>0`. Do not silently change between these domains.

An elementary counterexample shows why enforcing only positive parameters
would fail. With `m=1`, center of mass at the origin, and
`I_C=diag(1,1,3)`, the associated central second moment is
`diag(1.5,1.5,-.5)`. The inertia is positive definite but cannot describe a
nonnegative spatial mass distribution.

Therefore these papers supply a useful admissible domain for rotational
inertia and center-of-mass extensions. Our scalar condition `m>0` is enough
for its present scalar role; it is not the full physical domain of a rigid
body's ten inertial parameters.

### 3. What remains missing for our probabilistic construction

A proposed extension still needs all of the following, stated independently:

1. A joint measurement model for the actually recorded channels, including
   uncertainty in kinematics and correlations created by differentiation,
   filtering, or common sensors.
2. A mechanics constraint linking latent quantities and one shared object
   inertia across all trials.
3. A justified reference measure or prior on the admissible domain and latent
   quantities. A feasible cone alone assigns no relative weights inside it.
4. A definition and justification of the reported inertia and its uncertainty.

Conditioning on a measured regressor in `p(b|A,Phi)` and using a covariance for
its output residual can be useful, but it does not automatically propagate
uncertainty in `A` itself. A latent-channel model could sometimes reduce to
that likelihood after approximation or integration; that reduction would
need to be derived for the actual experiment. Neither the feasibility cone
nor the LogDet penalty selects our flat `df d alpha dOmega` measure.

### 4. A constructive compatibility observation for a future tensor readout

There is a simple extension candidate worth retaining as mathematics, without
claiming it is selected or optimal. Suppose a joint law already exists over
physically admissible `J`, total force magnitude `f`, and center-of-mass
acceleration magnitude `alpha`, with `f=m alpha` almost surely. Assume
`0<E[alpha]<infinity` and the weighted entries below are integrable. Define

\[
\widehat J=\frac{E[\alpha J]}{E[\alpha]}.
\]

Its mass entry is exactly

\[
\widehat J_{44}=\frac{E[\alpha m]}{E[\alpha]}
=\frac{E[f]}{E[\alpha]}.
\]

Thus it recovers the original mass readout. It also preserves positive
semidefiniteness: for every vector `v`,
`v^T Jhat v=E[alpha v^T Jv]/E[alpha]>=0`. It is positive definite if `J` is
positive definite almost surely on positive-weight samples. Fixed linear
support constraints are preserved by the same positive averaging argument.

With every component expressed in a common frame and the *same* weight law,
this averaging also commutes with addition of pseudo-inertias. That property
does not hold automatically if each component is fitted separately with a
different acceleration weight. This observation gives a possible physically
consistent readout architecture; it does not justify the joint law or choose
acceleration weighting for all rotational tasks.

In particular, rotational dynamics involve an inertia tensor and gyroscopic
terms. Torque and angular acceleration need not be parallel. Replacing scalar
mass by a matrix in the existing compatible-vector parameterization would
therefore be an inadequate extension.

### 5. Regularization would be a deliberate additional assumption

For a fixed shape and center of mass, write `J=m B` with fixed positive-definite
`B`, and let `J_0=m_0 B`. Direct substitution into the LogDet expression gives

\[
D(J\Vert J_0)=4\left(m/m_0-\log(m/m_0)-1\right).
\]

This explicitly prefers a supplied nominal mass and discourages zero mass.
It is not flat weighting and has an adjustable strength when used as a
regularizer. Its ability to stabilize a weakly excited inverse problem would
come with information from that nominal model, which must be credited in any
comparison with our estimator.

### 6. Investigation consequence

The most useful connection is a separation of responsibilities: physical
moment geometry defines admissible inertias; a measurement model describes
evidence; a measure/prior and readout define inference. These sources give
strong machinery for the first responsibility and practical, explicit choices
for parts of the latter two. They neither invalidate nor establish the
original ratio-of-means inference. They do show that a broader robotics
application should compare against physically constrained, geometrically
regularized identification using the same data and prior information.

## Coverage limitations

- This is a two-paper focused review, not an exhaustive robot-identification
  search or novelty determination.
- No claim that Lee's author preprint is text-identical to the final journal
  version; published metadata was checked against SNU and the author pages.
- No independent reproduction of either paper's experimental performance.
- The proposed weighted pseudo-inertia readout above is an algebraic
  feasibility observation. Calibration, identifiability, optimality, and
  computational tractability remain to be investigated.
- No production files or managed VD documents were changed.

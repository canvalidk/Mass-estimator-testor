# Internal continuation: shared size and direction, 2026-09-28

Scratch only. No new equation adopted, no managed writes, no production edits.
Continue ../ours27_bridge_2026_09_28/WORKING_STATE.md. User authorized pursuing
the most promising route after identifying that rotating aligned observed
vectors does not make transverse candidate directions known.

## Source boundary

Managed snapshot remains canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
Read through connector: catalog, mass_estimator_lemmas_2026-09-27.md,
ours27_full_covariance_2026-09-27.md, candidate_size_and_null_readings,
ours26 learned-scale and hint records, and the September25 axis-splitting
record. Agent premise audit also checked the ratified August1 shared-quantity
source. All paths and distinctions are in premises/direction_and_null_premises.md.
No use of unsanctioned expansions as an authority. No requested record.

## Main logical advance

The September25 axis-splitting record defines combination as integrating each
axis's candidate separately before multiplying mass evidence. Its section4(b)
then infers product candidate weighting from axis splitting, and product plus
isotropy/no size scale leads to flat weighting. This is a valid conditional
argument, but the independent nuisance elimination is ADDITIONAL to elementary
representation invariance. A shared size state or tensor can be carried through
an axis recoding. Integral(product likelihoods) is then correct; product of
separately integrated likelihoods resets the shared state and changes the model.

Thus preserving the elementary axis bookkeeping does not force flatness. This
is a proposed premise correction, not a silently adopted change of definition.
Measured-null neutrality is also a property of ours27, not a theorem of the
ratified shared-quantity argument or of rotation symmetry.

## Exact geometric mechanism

In the common-noise experiment put A=x sin(theta)+y cos(theta), H=||A||^2,
B=x cos(theta)-y sin(theta), Q=||B||^2, T=H+Q. For a fixed isotropic candidate
measure with finite damped moments, evidence=e^(-T/2) F_M(H). Define
kappa(H)=2 d log F_M/dH. Exact identities:

    E[v|A] = kappa A
    Cov(v|A) = kappa I + 2 kappa'(H) A A'
    E||v||^2 = H kappa^2 + M kappa + 2H kappa'
    at a profile mode: -ell'' = kappa(h)(h-q)

Flat dv has kappa=1 and candidate squared-size mean H+M. In a specified fixed
candidate experiment E[H]=R^2+M, so that mean averages R^2+2M. This is a
diagnostic of the high-dimensional size weighting, NOT a claim that posterior
variance is itself erroneous or that an ontic mass must exist.

For radial density r^(gamma M+o(M))dr with an interior regular saddle and H=hM,
h kappa^2+(1-h)kappa-gamma=0. Flat Cartesian volume has gamma=1, kappa=1;
dimension-subexponential radial weighting has gamma=0, kappa=1-1/h for h>1.
This explains the common-noise regular width correction without selecting
finite-data radial shape. Exact direction identities and numeric checks are
in geometry/size_direction_identities.py + JSON.

## Reopened scalar reference-size candidate

This was already considered in the earlier geometry_reference_probe; it is
not a novel newly discovered prior. It is reopened because its measured-null
feedback is no longer automatically a disqualification under the corrected
premise interpretation. Let t=1+omega^2>=1, conditional candidate covariance
(t-1)I. Fixed-m Fisher information I_tt=M/(2t^2) yields dt/t as conditional
reference weighting. This needs an extra reference-model premise; it is not
the unique law forced by Newton II. The half-line is not the full scale group,
so avoid an unqualified Haar claim.

    p(m,t|data) proportional pi(m) t^(-M/2-1)
                       exp[-Q_m/2-H_m/(2t)]

With beta=1-1/t, candidates are N(beta b_m,beta K_m). Integrating t gives
pi(m) exp(-Q_m/2) gamma(M/2,H_m/2)/(H_m/2)^(M/2), with the H=0 limit.
Posterior proper for every finite datum and M>0. Readout remains ratio of mean
stacked candidate sizes. The full-covariance version preserves profile ties
and quartic peaks, unlike the earlier fitted-local-power construction.

An all-zero experiment keeps p(m)=pi but updates the size state: beta is
Beta(1,M/2), E beta=2/(M+2). Appending k measured zero components multiplies
the shared-state law by (1-beta)^(k/2). This indirectly changes the mass law
when other data are active. Keeping only the unchanged mass marginal loses
information needed for subsequent measurements. Joint-state sequential and
one-shot updates agree, checked in shared_size_probe.py.

For aligned F=3,a=1, unit active-axis noise, point is2.2064082220 in1D;
adding two zero unit-noise axes gives2.0184357974. As their physical noise
tends to0 in the LOCAL-NOISE model, point tends2.0696736242, not the1D result.
Scalar-law quadrature refinement200->400 differed by <4e-14 in density.

## Zero conditioning versus omission

The epsilon discontinuity is real under the old comparison. For added
candidate b=delta sqrt(t-1) xi and zero measurement noise epsilon, the size
factor is [1+(delta/epsilon)^2(t-1)]^(-1/2). Structural delta/epsilon->0
removes the coordinate continuously. Local-noise delta=epsilon leaves
t^(-1/2) for every epsilon>0. A precise observation of an unknown coordinate
is conditioning information, not automatically omission of a structural dummy.

In shared local-noise Gaussian mixtures, demanding epsilon-to-fresh-default
omission equality forces a single fixed scale (not necessarily infinity).
Physical isotropy changes the limiting conditioning measure rather than
automatically restoring equality. See measure/null_conditioning_audit.py.
An explicit design metric with standardized extra covariance going to0
restores omission by dominated convergence, but adds a nonunique premise.
Do not dismiss these effects as changed pi: same channel ratios can leave
the angle starting law exactly unchanged.

## Why a scalar size state is still insufficient

Regular replicated counterexample: m=1, instrument ratios(1,10), squared
standardized candidates(.1,9.9). Per replicate,
J=.12204930889128517, G=.2598029604940692. A scalar state learns beta=5/6,
while the calibration multiplier is J/(J+G)=.31962441676133263. Hence
sampling/law variance tends2.607226762514809, despite information growing
without bound. This is separate from Fieller and must not be called a
vanishing-signal failure. geometry/shared_scale_heterogeneous_width.py proves it.

More generally a single radial evidence function of H cannot distinguish
two sequences with equal total energy but different allocations to the
instrument's mass-sensitive directions. The full candidate shape matters.

## Tensor route reopens, with actual remaining limitations

The earlier full-tensor route was rejected partly by comparing a shared
2D covariance model with new independent scalar models after splitting.
Keeping the SAME tensor, covariance, prior and sharing metadata through
the recoding removes that representation objection. It does not decide which
trials share shape, how frames are related across different instruments, or
whether readings repeat one candidate or only one mass.

For a replicated covariance group, the earlier tensor derivation remains:
in whitened candidate coordinates v_i~N(0,R), observed covariance I+U_m R U_m',
Beta=R(I+R)^(-1), and conditional candidates N(Beta A_i,Beta). This captures
direction-specific regular information without requiring Gaussian generating
pushes. Reference tensor volume |I+R|^(-(d+1)/2)dR has a proper posterior
only for N>d-1. Existing proper projective alternatives leave shape choices.

We additionally checked whether a canonical boundary extension cures small N.
For d3,N1, truncating eigenvalues R<L yields a different limiting mass/candidate
law than the equally rotation-invariant cutoff |I+R|<L. At zero projected
energy their mean shrinkage matrices are respectively (11/15)I and(2/3)I.
The generalized singular-Wishart Beta extension gives(4/5)I instead. The
different slopes give different finite-data mass laws; rotation does not pick
one. Normalizing a boundary posterior before adding another reading also
loses covariance directions needed for ordinary sequential Bayesian update.
See geometry/tensor_boundary* and measure/tensor_boundary* for final agent
derivations. Primary singular-Beta reference located through the search index:
Uhlig1994, https://home.uchicago.edu/~huhlig/papers/uhlig.aos.1994.pdf.
The indexed abstract was readable; opening the full PDF returned404. The
boundary comparisons above are our own derivations, not claims read in that
paper. measure/tensor_boundary_audit.py also directly checks mass laws:
for x=(2,0,0),y=0, uniform angle, limiting P(m<1) is .2705582155 for the
eigenvalue ceiling, .2859423924 for the determinant ceiling, and .2557286651
for the singular-Wishart extension. Thus the ambiguity reaches mass itself.

## Next justified work

Pursue one joint candidate size/direction law, preserving its actual dependence
under coordinate recoding. Do not choose another scalar power or call a finite
scale itself the solution. The specific unsupported passage is elementary
axis splitting -> independent nuisance integration -> flat measure. Fixing that
passage opens the tensor route, but neither its reference measure nor its
all-input continuation has been uniquely forced. No ours28 is yet justified.

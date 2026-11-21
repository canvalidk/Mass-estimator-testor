# What associativity forces for reciprocal disagreement when blocks are joined

The reciprocal tangent proof supplies a block weight \(h\), mean local displacement \(b\), and accumulated between-reading squared residual \(J\) at each candidate mass. On joining positive-weight blocks the extra residual from requiring a common displacement is
\[
 d_{AB}=\frac{h_Ah_B}{h_A+h_B}(b_A-b_B)^2,\qquad
 J_{AB}=J_A+J_B+d_{AB}.
\]
A singleton has \(J_i=0\). Both \(d\) and \(J\) are invariant under changing the local mass coordinate with nonzero derivative. We study how this exact geometry could alter **relative evidence** when new readings arrive. Multiplication of ordinary singleton factors, with the prior counted once, supplies the baseline. [Reciprocal tangent and boundary conventions](e12-tangent-pairings.md) · [T12](../definition/theorems/T12.md)

## Classification on realizable scalar readings

The increase in within-group squared residual is the classical least-squares merging cost used by Ward. Here the weights and local displacements come from the measurement tensor and vary with the candidate mass. The question concerns the factor that can enter the returned law. [Ward, *Hierarchical Grouping to Optimize an Objective Function* (1963)](https://doi.org/10.1080/01621459.1963.10500845)

Assume the extra log factor for joining two blocks, in addition to adding their log evidences, is a continuous symmetric function \(g(J_A,J_B,d_{AB})\). The same function is used at every candidate, and the same prior of positive continuous density on the open positive angle arc is counted once. Different choices of grouping must yield the same *normalized* law. The function uses the two accumulated discrepancies and the extra discrepancy. Dependence on other quantities, including precision separately from discrepancy, is a different class. Subtract \(g(0,0,0)\) without loss, since it changes every merge by a mass-independent constant.

**Theorem.** On nondegenerate common-instrument scalar data, necessarily
\[
 \ell_D(\theta)=\sum_{i\in D}\ell_i(\theta)+f(J_D(\theta))
                  \quad \pmod{\text{constants in }\theta}
\]
for a continuous \(f\) with \(f(0)=0\). Conversely **every** such \(f\) gives an associative law through the explicit merge
\[
 g_f(u,v,d)=f(u+v+d)-f(u)-f(v).
\]
Here and below a finite joint law requires positive singleton factors and a finite normalizing integral. The proof provides the stated identities of normalized laws; no physical rule in the project requires that this is the whole state of a join. If only the *new* discrepancy \(d\) can affect the extra factor, the sharper conclusion is \(f(J)=cJ\) for a universal constant \(c\). Thus even that smaller class leaves an undetermined constant.

**Why equality of normalized laws suffices.** At a given angle \(\theta_*\) any singleton tangent weight \(h_i>0\) and displacement \(b_i\in\mathbb R\) are realizable in the one-instrument experiment. In the angle chart choose the whitened observed row
\(y_i=\sqrt{h_i}(u_*+b_i v_*)\), where \(u_\theta=(\sin\theta,\cos\theta)\), \(v_\theta=(\cos\theta,-\sin\theta)\). At any \(\theta\), write \(x_i=u_\theta^\top y_i,\ z_i=v_\theta^\top y_i\). Then \(h_i=x_i^2,\ s_i=x_i z_i,\ b_i=s_i/h_i\) when \(h_i>0\). Shifting all the chosen \(b_i\) by one constant at \(\theta_*\) leaves every joining discrepancy there unchanged. Choose this shift large enough that the **middle** row of a three-singleton join has negative slope in the physical channel coordinates. It then has \(h_2=0\) at an interior angle \(\theta'\). For distinct displacements the other two rows are not parallel to that row, so \(h_1,h_3>0\) there.

The log-factor difference between the two grouping trees is constant almost everywhere in \(\theta\) by probability-law equality. The positive singleton factors cancel. Off the finitely many zero-tangent angles the difference is continuous, so this same constant applies throughout. Near \(\theta'\), \(s_2\to0\) and \(s_2^2/h_2\to\|y_2\|^2\). Therefore the first discrepancies \(d_{12},d_{23}\) have the same limit \(\|y_2\|^2\), and the second discrepancies \(d_{(12)3},d_{1(23)}\) both tend to \(d_{13}\). Symmetry and continuity of \(g\) make the tree difference tend to zero. Thus the constant was zero, and *pointwise* equality of the merge sums holds at \(\theta_*\) for every choice of the positive weights and local displacements. Continuity in these displacements extends the equality to coincident ones. This reasoning uses actual negative-slope **observations** allowed by positive-mass Gaussian measurements. It assigns no positive prior mass to an exceptional candidate value.

**What three readings force.** Put \(f(x)=g(0,0,x)\). Choose weights \((1,t,1)\) with \(t>0\) and take \(b_2=b_3\). If the total discrepancy is \(p+q\), the left tree joins with increments \(p,q\), and the right tree with increments \(0,p+q\), where \(p/(p+q)=1-(1+t)^{-2}\). Vary the scale of \(b_1-b_2\) and \(t\). This realizes every \(p,q>0\). Equality of the trees gives
\[
 f(p)+g(p,0,q)=f(p+q).
\]
Continuity covers zero. Sequentially append singletons and use the exact \(J\)-merge identity. The extra log factors telescope to \(f(J_D)\). Associativity of the law ensures the same outcome for every grouping. The converse follows by telescoping for arbitrary block merges. The canonical merge also works with zero-weight blocks under the reciprocal theorem's conventions.

If \(g\) depends only on \(d\), the last displayed equation becomes \(f(p)+f(q)=f(p+q)\). Additivity first fixes \(f(q)=q f(1)\) for nonnegative rationals; continuity extends this to every nonnegative real. If the additional factor is assumed to discount discrepancy, then \(c\le0\). Converting an observed one-dimensional standardized Gaussian contrast into its **density relative to its value at zero** gives \(e^{-d/2}\) and \(c=-1/2\). The density interpretation is an extra premise. At a fixed true candidate the conditional difference of the two independent corrections has precisely this standardized Gaussian form. It is a diagnostic calculated with the candidate held fixed. [Local contrast experiment](e12-tangent-diagnostic.md)

## Can the two individual precisions choose the coefficient?

There is a stronger rigidity statement when the merge costs themselves are required to compose **pointwise**, before normalization of the mass law. Allow the new log correction \(p(h,k,s)\) to use both precisions as well as the signed difference \(s=b_A-b_B\). Assume it is symmetric in \(h,k\), even in \(s\), continuous, and zero when \(s=0\). Old internal \(J\)'s do not enter this version. Addition of these corrections must be independent of every binary grouping for all positive weights and real displacements. Then necessarily
\[
 p(h,k,s)=c\,\frac{hk}{h+k}s^2 .
\]
No preselected functional dependence on \(d_{AB}\), or Gaussian sampling hypothesis, is used. Dropping the equal-displacement condition admits pure weight terms, for example \(a(h+k)-a(h)-a(k)\), which can affect the law because the precisions are candidate-dependent.

Write \(q_h(s)=p(h,h,s)\). Join \(n\) copies each of the \(h\)-weight observations at \(x\) and at \(y\). Grouping identical displacements first gives \(q_{nh}(x-y)\), whereas joining as \(n\) identical mixed pairs gives \(nq_h(x-y)\). Thus \(q_{nh}=n q_h\). Rational approximation and continuity in the weights now give \(q_h=h q_1\). With four equal weights and displacements \(a+b,a-b,-a+b,-a-b\), compare pairing by common \(a\) with pairing opposite displacements. The result is
\[
 2q_h(2b)+q_{2h}(2a)=q_h(2a+2b)+q_h(2a-2b).
\]
Consequently \(q=q_1\) satisfies the parallelogram equation. Evenness, \(q(0)=0\), and induction using \(q((n+1)s)+q((n-1)s)=2q(ns)+2q(s)\) give \(q(ns)=n^2q(s)\); continuity then gives \(q(s)=q(1)s^2\) for every real \(s\).

Finally take two \(h\)-weight blocks at \(x,y\) and two \(k\)-weight blocks at \(x,y\). Pair equal weights first, or pair each \(h\)-block with the \(k\)-block at the opposite displacement. With \(s=x-y\) this yields
\[
 2p(h,k,s)+q_{h+k}\!\left(\frac{h-k}{h+k}s\right)=q_h(s)+q_k(s).
\]
Substitution gives the theorem with \(c=2q(1)\). Conversely this multiple of the between-block squared residual composes by the previous \(J\) identity. If only normalized posteriors are to associate, individually weight-dependent penalties require a separate argument to control the mass-independent difference between trees. The earlier interior anchor establishes that control for the stated \(J_A,J_B,d\) class.

## Limits of using only this disagreement to change a Gaussian law

For scalar readings with positive channel components let \(B=\sum y_i y_i^\top\), \(H=u_\theta^\top B u_\theta\), \(K=v_\theta^\top Bv_\theta\), \(S=v_\theta^\top Bu_\theta\). Every nonzero singleton is active throughout the open positive arc and the within-reading residual vanishes. Thus \(J=K-S^2/H=\det B/H\). Starting with the face-value Gaussian singleton factors gives
\[
 \ell_D = H/2+ f(\det B/H)\quad\pmod{\text{constants}} .
\]

**Weak noncollapse.** Let \(f(x+d)-f(x)\le M d\) at all sufficiently large \(x\) and all \(d>0\), for a finite \(M\ge0\). This includes *every* nonincreasing penalty, every linear response of either sign, and many increasing rewards. Fix the incompatible positive rows \((r,1),(1,r)\) with \(r>1\). Their reference spectrum is \(L=(r+1)^2,\ l=(r-1)^2\), so on the positive arc \((L+l)/2\le H_0\le L\), with unique maximizer \(\pi/4\). Choose \(r\) close enough to 1 that \(4M Ll/(L+l)^2<1/2\). Shrinking the errors gives \(B=tB_0\) and log factor \(tH_0/2+f(tLl/H_0)\). For large \(t\), increasing \(H_0\) by \(v>0\) increases the log factor by at least
\[
 t v\left\{\tfrac12-\frac{4M Ll}{(L+l)^2}\right\}.
\]
Outside any fixed neighborhood of the unique maximizing line there is a positive gap in \(H_0\) compared to a smaller fixed neighborhood inside. The latter has positive prior mass. Hence the law concentrates on mass 1 and its angular variance tends to zero, despite the fixed conflict. For a nonincreasing \(f\), **every** \(r>1\) works. Penalizing a large \(J\) at each candidate reinforces preference for large \(H\) when \(\det B\) is held fixed.

**Information at a balanced conflict.** The fixed scalar readings \((1,0),(0,1)\) give \(B=t\,\mathrm{Id}\) under common shrinkage. On the open positive arc both rows are active, \(H=K=J=t\), and \(S=0\). For *every* \(f\) the joined factor is constant in mass, so the returned law equals the prior. Consequently none of these Gaussian-singleton constructions satisfies the literal requirement of **nonconstant** limiting evidence on every incompatible shape, even with unrestricted \(f\). Informative spread on a balanced shape is a stronger demand than spread alone. [Boundary premise](e12-product-obstruction.md)

**Bounded contrasts.** Even an unrestricted continuous \(f\), without a growth or sign assumption, cannot keep every relative log contrast bounded on a *single* such fixed shape. Write \(P=\det B_0>0\), \(j=P/H_0\), which sweeps a nontrivial interval. Pick three attained levels \(j_0,kj_0,k^2j_0\) with \(k>1\). At scale \(t\) the log law relative to the prior is \(tP/(2j)+f(tj)\) up to one additive constant. Compare the contrast for \((j_0,kj_0)\) at scale \(kt\) to that for \((kj_0,k^2j_0)\) at scale \(t\). Their \(f\)-contributions are identical, while the difference between their Gaussian contributions is
\[
 \frac{tP}{2j_0}(1-1/k)(k-1/k).
\]
At least one of these contrasts becomes unbounded along a sequence. This stronger requirement differs from requiring only persistent angular spread. The unrestricted nonlinear case under only the weaker requirement is not settled by this comparison.

**Bounded contrasts with arbitrary singleton factors in the new-discrepancy-only class.** Let \(x=(1,1),\,y=(2,1),\,z=(1,2)\). Every one of \(xxy,xy,xxz,xz\) is incompatible. In the combination
\[
 (\ell_{xxy}-\ell_{xy})-(\ell_{xxz}-\ell_{xz}),
\]
the two contributions of singleton \(x\) cancel and grouping leaves \(c(d_{xy,x}-d_{xz,x})\) up to a constant. At reference scale and angle zero the two increments are \(1/6,\ 2/15\), respectively, by direct substitution. At angle \(\pi/2\) they swap. These are continuous endpoint limits of the positive-domain data factors in the difference after cancellation. Shrinking the errors multiplies their difference by \(t\). Thus bounded log contrasts on all four fixed conflicts force \(c=0\), regardless of the positive singleton factors. Strict product composition then precludes compatible singleton law concentration by the earlier duplicate argument. [Duplicate obstruction](e12-product-obstruction.md)

These theorems identify where more state or more physics must enter if the law is to retain informative spread on conflict and learn from compatible data. Retaining further instrument invariants, or defining a different singleton law, changes the assumptions as indicated in each result. Same-object joining gives no physical joint law for distinct bodies by itself. [Finite checks](../tests/test_e12_joining_correction.py)

## A boundary obstruction even before joining two readings

**Which Liseo equation is closest?** For one scalar reading with equal channel errors and the unit half-Cauchy starting law, the ordinary Gaussian working factor gives exactly the positive restriction of Liseo's *reference* posterior, equation (14). His mass density is proportional to \((1+m^2)^{-1}\exp[-(F-ma)^2/(2\sigma^2(1+m^2))]\); the first factor is the Jacobian of the noise angle. For any strictly positive singleton the E12 family has \(Q=0\). Wherever the channel cap is inactive it gives that very same law. At \((F,a)=(2,2)\) in unit noise the precision is 8 in both. At \((2,1)\) Liseo's (14) gives precision 5 and E12 gives 2. Thus there is an exact agreement region and an explicit finite difference. Liseo's equation (12) uses a different prior on the two means. [Liseo (2003), equations (10)–(14)](https://www.dss.uniroma1.it/RePec/mtn/articoli/2003-1-10.pdf)

There is a direct comparison here with the earlier scalar law in the repository. Liseo’s equation (12), restricted to positive mass, is the ratio of two independent posterior normal means in the equal-variance flat-means experiment. On \((F,a)=(f,0)\), \(f>0\), the means have laws \(X\sim N(f,\sigma^2)\) and \(Y\sim N(0,\sigma^2)\). Conditional on \(X/Y>0\), the probability that \(m\le1\) is \(\Pr(|X|\le|Y|)=2\Phi(u)\Phi(-u)\), \(u=f/(\sqrt2\sigma)\). Here conditioning on the relative sign leaves the magnitudes unchanged, and \(Y-X,Y+X\) are independent. This probability tends to zero as \(\sigma\to0\). More generally \(X\to f,\ Y\to0\) in probability, so every fixed upper mass bound has vanishing probability after positive conditioning. The E12 construction on the same axis with the same zero-data half-Cauchy law retains a positive probability below 1. This separates the new law from that specific scalar posterior. Both behaviors have a mathematical derivation, and the project's boundary premise is what distinguishes them here. [Liseo (2003), pp. 137–138](https://www.dss.uniroma1.it/RePec/mtn/articoli/2003-1-10.pdf) · [Repository equation comparison](../docs/literature_review/mass_candidate_liseo_evidence_2026-09-08.md)

**Lemma.** Fix a proper angular prior with support throughout \(0<\theta<\pi/2\). Suppose that at every positive finite error scale the final law is weakly continuous in a scalar reading, and for every \(f,a>0\) its normalized singleton law has exactly the flat Gaussian working factor \(e^{(f\sin\theta+a\cos\theta)^2/2}\), with \(f,a\) in noise units. Then as both errors shrink on a fixed nonzero axis reading, the angular variance tends to zero.

Indeed, hold a whitened force \(b>0\) fixed and take positive accelerations down to zero. The Gaussian factors converge uniformly in angle. Uniqueness of the weak limit forces the same normalized law on \((b,0)\), with factor \(e^{b^2\sin^2\theta/2}\). Now put \(b=\sqrt{t}\,b_0\), \(b_0>0\), and let \(t\to\infty\). Outside any fixed neighborhood of \(\pi/2\), \(\sin^2\theta\) has a strictly smaller supremum than its infimum on some smaller interval nearer \(\pi/2\) having positive prior mass. Its probability outside the neighborhood therefore tends to zero. Angular variance tends to zero by boundedness. Swapping channels gives concentration at zero for an acceleration-only reading.

Each axis has positive leftover at every *strictly positive finite* mass, although the infimum over the closed arc is zero. Thus the November 17 whole-leftover reading of T30, **if it includes axes**, is incompatible with exact flat final laws for all strictly positive singletons and continuity. L10 explicitly calls the flat factor a working likelihood. This lemma concerns the stronger premise that the **final normalized singleton law** must equal that working likelihood times the prior. The two limit operations in the proof only use continuity at each fixed error scale. [L10 and L19](../definition/premises.md) · [T30 and its provisional domain](../definition/theorems/T30.md)

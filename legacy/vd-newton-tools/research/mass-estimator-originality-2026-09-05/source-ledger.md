# Claim-to-source ledger

Access date for all entries: 2026-09-05. This is an internal provenance record. No quoted prose is used in the deliverable. Mathematical derivations explicitly labeled as our specializations are evaluated against the pinned VD definition rather than represented as verbatim results of cited authors.

| ID | Author / title / date | URL and inspected locator | Supported claim and limitation |
|---|---|---|---|
| S01 | Morris W. Katz, Admissible and Minimax Estimates of Parameters in Truncated Spaces, 1961 | https://doi.org/10.1214/aoms/1177705146 | Historical origin of positive normal-mean estimator. Original full text inaccessible; mathematical verification through S02. |
| S02 | Yuan-Tsung Chang, Nobuo Shinozaki, William E. Strawderman, Simultaneous Estimation of p Positive Normal Means with Common Unknown Variance, manuscript 2017-10-26 | https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf ; p2 Eq1.4, p3 Lemma2.1 | Exact h_sigma formula and generalized Bayes interpretation. Full PDF available; decisive passage independently verified by root, turn27view0. Not a force/acceleration ratio construction. |
| S03 | Byron P. Roe and Michael B. Woodroofe, Setting Confidence Belts, 2000 preprint | https://arxiv.org/pdf/hep-ex/0007048 ; continuous example pp3–4 Eqs1–2 | Flat-positive-prior normal posterior. Verified by root turn28view5. Physical measurement precedent, not VD mass estimator. |
| S04 | George Marsaglia, Ratios of Normal Variables, JSS16(4), 2006 | https://www.jstatsoft.org/article/view/v016i04 ; full PDF https://pdfs.semanticscholar.org/9ead/12ca2b04f00533587ecfe90df7ea334d2f58.pdf ; sections2–4 | Gaussian ratio scale, Cauchy null, missing moments beyond null. Worker full read, root independent PDF read turn30view0. General normal ratios rather than sampling law of ratio-of-positive-means estimator. |
| S05 | Hea-Jung Kim, On the Ratio of Two Folded Normal Distributions, 2006, online2007 | https://doi.org/10.1080/03610920600672229 | Abstract explicitly includes half-standard Cauchy. Worker publisher abstract inspected; root direct open failed. Full theorem not inspected. Folded differs from truncated away from zero. |
| S06 | Brunero Liseo, Bayesian and conditional frequentist analyses of the Fieller's problem. A critical review, Metron61, 2003, pp133–150 | https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review ; Eq5 p136 and Eqs10,12,13 pp137–138 | Flat profile and null Cauchy posterior, transformed original-means prior. Full article text read by worker; root independently checked key equations turn27view1. Publisher PDF https://www.dss.uniroma1.it/RePec/mtn/articoli/2003-1-10.pdf failed. Half-positive unequal-scale adaptation is our derivation. |
| S07 | Ulrike von Luxburg and Volker H. Franz, A Geometric Approach to Confidence Sets for Ratios: Fieller's Theorem, Generalizations, and Bootstrap, 2007 preprint / 2009 publication | https://arxiv.org/pdf/0711.0198 ; Construction1, Theorem1 p6, Figure2.2 p7 | Bounded/unbounded/all-real sets and origin geometry. Root independently verified turn25view2/turn26view5. Positive-mass restriction is our specialization; threshold calibration not established for full VD cone. |
| S08 | Hans Butler, Controlling a position of a mass, especially in a lithographic apparatus, EP1455231A2 published2004-09-08 | https://patents.google.com/patent/EP1455231A2/en ; sections2.3–2.4 | Explicit noisy force0/acceleration0 mass-estimation issue; force offset at a0; adaptation suspension/bounding. Root read turn27view2,turn28view4. Technical disclosure only, not legal novelty analysis or identical estimator. |
| S09 | Thomas W. Bartel, Sara Stoudt, Antonio M. Possolo, Force Calibrations using Errors-in-Variables Regression and Monte Carlo Uncertainty Evaluations, 2016 | https://www.nist.gov/publications/force-calibrations-using-errors-variables-regression-and-monte-carlo-uncertainty | EIV plus uncertainty methods in force calibration. Official abstract inspected turn25view0. Does not support exact mass readout. |
| S10 | Roderick J. Little, Calibrated Bayes, an Alternative Inferential Paradigm for Official Statistics, JOS28(3)309–334, 2012 | https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/calibrated-bayes-an-alternative-inferential-paradigm-for-official-statistics.pdf ; p319 | Explicit ratio of posterior expectations as O(1/n) approximation. Worker read and root verified turn29view1. Not an exact justification in low-signal divergent-mean regime. |
| S11 | Richard Arratia, Larry Goldstein, Fred Kochman, Size bias for one and all, first posted2013, inspected revision2018 | https://arxiv.org/pdf/1308.2729 ; section2.3 p11 Eqs19–20 | Coordinate size-biasing, applied to acceleration-tilted ratio mean. Worker read, root opened and located coordinate section turn28view0/turn29view2. Reciprocal uses different tilted measure. |
| S12 | JCGM, International Vocabulary of Metrology, entry2.9, third-edition online text | https://jcgm.bipm.org/vim/en/2.9.html | Measurement value with information and uncertainty. Root verified turn28view3. Does not establish estimator-quantity identity. |
| S13 | JCGM100:2008, Guide to the expression of uncertainty in measurement, section5.2 | https://www.iso.org/sites/JCGM/GUM/JCGM100/C045315e-html/C045315e_FILES/MAIN_C045315e/05_e.html | Correlated input covariance propagation. Root inspected official section turn28view2. Vector sum formula is componentwise application. |

## Deduction ledger

| ID | Deduction | Verification and conditions |
|---|---|---|
| D01 | Positive normal means at x0 and resolved one-zero ratios | Substitute phi0/Phi0 in S02 formula; fixed-direction scalar approximation; checked numerically. |
| D02 | Default null half-Cauchy, median/readout coincidence, quantiles | Direct Jacobian integration of independent half-normal magnitudes in pinned model; CDF and quantiles checked. |
| D03 | Divergence at every finite input for default full vector support | Continuous strictly positive Gaussian density on a compact f interval and acceleration boundary gives lower bound times integral da/a; second moment follows. Covariance nondegenerate, flat magnitude measure. Analytic proof, not inferred from finite integration. |
| D04 | Finite sampling moments possible for scalar ratio-of-means rule | Mills-ratio asymptotic h(x)~sigma²/abs(x) as x→−infinity; positive tail linear; product-polynomial bound under Gaussian readings. Not a full vector sampling proof. |
| D05 | Alternate symmetric chi3 measure preserves point but gives finite variance | Independent chi moments: Echi3=2sqrt(2/pi), E1/chi3=sqrt(2/pi), Echi3²=3, E1/chi3²=1. Numerical quadrature checks. |
| D06 | One-zero compatibility thresholds | Solve scalar Q(m)≤k with nonnegative readings; direct substitution checks; not unconditional frequentist coverage claim. |
| D07 | Reciprocal quantile correction and composition limitation | Monotone inverse swaps p with1−p. Linear functional ratios satisfy both identities; expectation selection needs additional law assumptions. No uniqueness theorem claimed. |

## Gap matrix at stopping point

| Question | State | Reason to stop / residual limitation |
|---|---|---|
| Prior positive-zero mechanism | Answered | Exact classical formula located and independently verified. |
| Bayesian null posterior vs profile | Answered | Close published construction with equations located and verified. |
| Infinite-moment inference | Answered/corrected | Primary ratio literature plus direct proof for actual VD support. |
| Direct zero mass-estimation discussion | Answered | Explicit 2004 force0/acceleration0 disclosure read. |
| Exact complete VD estimator | Open | Focused follow-ups returned adjacent methods, no complete match. No absence claim warranted. |
| Two-composition uniqueness | Open | Current managed proof explicitly a sketch; no exact primary theorem located. Requires mathematical work beyond literature novelty assessment. |
| Universal operational definition originality | Open | Broad historical philosophy not exhaustively searched; standards establish only reporting precedent. |

Discovery, independent primary verification, disconfirming follow-up, and a final exact-phrase gap wave are complete. Additional undirected searching is unlikely to change the component assessment; residual claims are explicitly bounded in the report.

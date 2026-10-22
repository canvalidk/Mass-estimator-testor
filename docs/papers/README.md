# papers/ — published literature for the mass estimator

Published papers that the estimator work cites or should cite, collected on the
user's instruction on 2026-09-23 ("These are standard citations anyways").

- These are **references, not records.** They are not members of the
  discussion, and nothing in its argument is bound to them. Catalog role:
  `reference`.
- **Publisher filenames are kept**, as in `Theory/machine matching problem/`.
- To add a paper: drop the PDF in, catalogue it, and move its line below from
  "Still to fetch" to "In this folder".
- 15 of the PDFs were fetched on 2026-09-23 by Claude through the built-in
  browser on the user's computer (his Imperial sign-in for the OUP papers).
  Downloaded files keep descriptive names; publisher filenames are kept where
  the publisher supplied one (`jrsssb_*`, `biometrics_*`, `SIM-*`).
- Not to be confused with `../literature_review/`, which is the unsanctioned
  2026-09-06..08 Codex search *about* the literature (audits, not papers).

## In this folder

| File | Citation | Why it is here |
|---|---|---|
| `SIM-38-2074.pdf` | Morris, White & Crowther (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine* 38:2074–2102. doi:10.1002/sim.8086 | ADEMP: the study-design discipline the VD-Newton lab's simulations follow. |
| `biometrics_77_4_1355.pdf` | Wang, Xu, Wang, Wu, Fung, Gao, Liang & Liu (2021). Penalized Fieller's confidence interval for the ratio of bivariate normal means. *Biometrics* 77(4):1355–1368. doi:10.1111/biom.13363 | Makes Fieller's interval always bounded with a penalty λ·log\|μ₂\| on the denominator, a power-law weight \|μ₂\|^λ. Pays for it with coverage (about 84% at weak denominator power, Table 1). |
| `Golub_VanLoan_1980_TLS_Cornell_TR80-411.pdf` | Golub & Van Loan (1980). An analysis of the total least squares problem. *SIAM J. Numer. Anal.* 17(6):883–893. (Cornell CS tech report TR80-411 version.) | Total least squares: the standard 'best compatible pair' fit when both variables are noisy. |
| `Bartel_Stoudt_Possolo_NIST_force_calibration_EIV.pdf` | Bartel, Stoudt & Possolo. Force calibrations using errors-in-variables regression and Monte Carlo uncertainty evaluations (NIST; published in *Metrologia*). | The metrology precedent closest to the problem: force measurement treated as errors-in-variables. |
| `jrsssb_30_1_190.pdf` | Lindley & El-Sayyad (1968). The Bayesian estimation of a linear functional relationship. *JRSS B* 30(1):190–202. doi:10.1111/j.2517-6161.1968.tb01519.x | Bayesian errors-in-variables: integrate over latent true values instead of fitting one pair. |
| `Gull_1989_straight_line_fitting.pdf` | Gull (1989). Bayesian data analysis: straight-line fitting. In *Maximum Entropy and Bayesian Methods*, 511–518. | Bayesian straight-line fitting with noise in both coordinates, nuisance values integrated out. |
| `jrsssb_16_2_175.pdf` | Fieller (1954). Some problems in interval estimation. *JRSS B* 16(2):175–185. doi:10.1111/j.2517-6161.1954.tb00159.x | The exact confidence interval for a ratio of means; unbounded when the denominator is not clearly nonzero. |
| `Koschat_1987_AnnStat_15_462_Fieller_characterization.pdf` | Koschat (1987). A characterization of the Fieller solution. *Annals of Statistics* 15(1):462–468. doi:10.1214/aos/1176350282 | Fieller's is the only exact method. (Often mis-dated 1981, including in Wang et al. 2021.) |
| `Liseo_2003_Fieller_problem_review_METRON.pdf` | Liseo (2003). Bayesian and conditional frequentist analyses of the Fieller's problem. A critical review. *METRON* LXI(1):133–150. | His eq. (12) is the one-dimensional version of the estimator's flat law (see ../literature_review/mass_candidate_liseo_evidence_2026-09-08.md). |
| `vonLuxburg_Franz_2009_ratios_arXiv_0711.0198.pdf` | von Luxburg & Franz (2009). A geometric approach to confidence sets for ratios: Fieller's theorem, generalizations and bootstrap. *Statistica Sinica* 19(3):1095–1117. arXiv:0711.0198 | Fieller's theorem made geometric; generalisations. |
| `Franz_2007_ratios_arXiv_0710.2024.pdf` | Franz (2007). Ratios: a short guide to confidence limits and proper use. arXiv:0710.2024 | A readable introduction to the ratio problem. |
| `Katz_1961_AnnMathStat_truncated_spaces.pdf` | Katz (1961). Admissible and minimax estimates of parameters in truncated spaces. *Annals of Mathematical Statistics* 32(1):136–142. doi:10.1214/aoms/1177705146 | Estimating a parameter known to be positive; the known-direction quotient. |
| `Arratia_Goldstein_Kochman_size_bias_arXiv_1308.2729.pdf` | Arratia, Goldstein & Kochman. Size bias for one and all. arXiv:1308.2729 | The point E[f]/E[alpha] is an acceleration-weighted (size-biased) mean of the mass. |
| `Sun_Berger_2008_sequential_arXiv_0805.3064.pdf` | Sun & Berger (2008). Objective Bayesian analysis under sequential experimentation. *IMS Collections* 3:19–32. arXiv:0805.3064 | Objective priors such as Jeffreys' depend on the experimental design, which clashes with updating. Cited by the lab's contribution 3. |
| `Berger_Liseo_Wolpert_1999_StatSci_integrated_likelihood.pdf` | Berger, Liseo & Wolpert (1999). Integrated likelihood methods for eliminating nuisance parameters. *Statistical Science* 14(1):1–28. doi:10.1214/ss/1009211804 | Default priors when nuisance parameters are integrated out. |
| `Gleser_Hwang_1987_AnnStat_15_1351.pdf` | Gleser & Hwang (1987). The nonexistence of 100(1−α)% confidence sets of finite expected diameter in errors-in-variables and related models. *Annals of Statistics* 15(4):1351–1362. doi:10.1214/aos/1176350597 | No always-finite interval can guarantee its coverage. The lab's weak-excitation impossibility result looks like a case of it. |
| `Feldman_Cousins_1998_arXiv_physics_9711021.pdf` | Feldman & Cousins (1998). A unified approach to the classical statistical analysis of small signals. *Phys. Rev. D* 57:3873. arXiv:physics/9711021 | How physicists handle bounded parameters and coverage; a reference a physics audience already knows. |
| `jrsssb_28_2_288.pdf` | Discussion on Mr Sprent's paper (1966). *JRSS B* 28(2):288–297. doi:10.1111/j.2517-6161.1966.tb00642.x | Lindley's opening contribution (pp. 288–290) and Beale's remark (p. 290): with known error variances, scaled to unit variance, and the distance from the origin independent of the slope, the factor outside the exponential is unity. That is the ours24 law. A flat prior on each true x instead pulls the slope to 0 and does not diminish with n (the eq. 18 drift). Fetched 2026-09-26. |
| `Hogg_Bovy_Lang_2010_fitting_a_model_to_data_arXiv_1008.4686.pdf` | Hogg, Bovy & Lang (2010). Data analysis recipes: fitting a model to data. arXiv:1008.4686 | §7: the standard astronomy recipe for uncertainties in both coordinates. Its projected likelihood through the origin is the ours24 exponent, with a prior flat in the angle (data units). It says the model is unspecified along the line, which axis splitting supplies. Fetched 2026-09-26. |

## Still to fetch

- **Dellaportas & Stephens (1995).** Bayesian analysis of errors-in-variables regression models. *Biometrics* 51(3). JSTOR 2533007
  — Bayesian errors-in-variables. JSTOR put up a bot check, so it needs fetching by hand.
- **Dufour (1997).** Some impossibility theorems in econometrics with applications to structural and dynamic models. *Econometrica* 65(6):1365–1388.
  — The general form of the coverage limit: near parameters the data cannot identify, a valid confidence set must be unbounded with high probability; one that never is has zero confidence level. JSTOR / Wiley.
- **Jeffreys (1961).** *Theory of Probability*, 3rd ed. Oxford University Press. A book: library.
  — The log-uniform prior for scale parameters, and Jeffreys' rule.
- **Marchand & Strawderman (2004).** Estimation in restricted parameter spaces: a review. *IMS Lecture Notes* 45. doi:10.1214/lnms/1196285377 (book chapter; the Project Euclid download did not return a PDF).
- **Wang, Wang & Carroll (2015).** The direct integral method for confidence intervals for the ratio of two location parameters. *Biometrics* 71:704–713.

### Held back: they bear on the 2026-09-23 parameter findings, which come later

- **Hartigan (1964).** Invariant prior distributions. *Ann. Math. Statist.* 35(2):836–845. doi:10.1214/aoms/1177703583 — relatively invariant priors are power laws.
- **Jaynes (1968).** Prior probabilities. *IEEE Trans. Systems Science and Cybernetics* 4(3):227–241 — priors from transformation groups.
- **Sclove, Morris & Radhakrishnan (1972).** Non-optimality of preliminary-test estimators for the mean of a multivariate normal distribution. *Ann. Math. Statist.* 43(5):1481–1490. doi:10.1214/aoms/1177692380
- **Leeb & Pötscher (2005).** Model selection and inference: facts and fiction. *Econometric Theory* 21(1):21–59.
- **Mach (1883).** *Die Mechanik in ihrer Entwickelung* (*The Science of Mechanics*) — mass ratios from mutual accelerations.
- **Wensing, Kim & Slotine (2017).** Linear matrix inequalities for physically-consistent inertial parameter identification. arXiv:1701.04395 — robotics payload identification.

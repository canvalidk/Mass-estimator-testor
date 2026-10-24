"""Write the predeclared next studies before any of their simulations run."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
studies = root / '03_trace_and_evaluator/mass_estimator/investigation/studies'
pilot = json.loads((studies / 'mass_excitation_pilot.json').read_text())
pilot_checkpoint = root / '.tools/mass_estimator_studies/assessment-pilot-20260920/study_checkpoint.json'
pilot_hash = hashlib.sha256(pilot_checkpoint.read_bytes()).hexdigest()
metrics = deepcopy(pilot['primary_metrics'])
scenario = {'type': 'mass_excitation', 'masses': [.25, 1, 4], 'acceleration_snrs': [0, .5, 2, 8]}

main = {
    'schema_version': 1,
    'name': 'Mass and excitation: fixed-size main assessment',
    'phase': 'main_assessment',
    'ademp': {
        'aims': [
            'Measure cell-specific accuracy, systematic log error, invalid outputs and repeated-experiment interval coverage for the mass estimator and existing comparators.',
            'Determine how true mass, excitation and stable repeated readings jointly affect performance; distinguish no information from statistical and numerical error.',
            'Use fresh random numbers and a fixed 10000 replicates per cell to improve Monte Carlo precision after the exploratory pilot. No estimator or loss is tuned from the main results.'
        ],
        'data_generating_mechanisms': [
            'Newton II is assumed correct: the fixed positive physical mass satisfies true force = mass times true acceleration, with a common unknown direction [1,0,0]. The simulator uses that direction; inference does not receive it.',
            'Cross masses 0.25, 1 and 4 with per-reading acceleration SNR 0, 0.5, 2 and 8. Each force and acceleration coordinate has independent Gaussian noise with known SD1. The force SNR equals mass times acceleration SNR in these units.',
            'Separate experiments use n=1 and n=4 independent readings of the same latent pair. Their means have SD1/sqrt(n). R=10000 independently simulated experiments per condition is distinct from n readings in an experiment.',
            'Training noise is shared across conditions within an experiment to support matched comparisons. Independent seeds are used across n experiments and for held-out observations. Neither noise realization nor data from the pilot is reused.'
        ],
        'estimands': [
            'Fixed positive physical mass for each condition, rather than a distribution averaging mass values or conditions.',
            'Repeated-experiment expected bounded squared log loss and probability of a factor-two error, including invalid positive-mass estimates.',
            'Repeated-experiment coverage, lower and upper noncoverage, and mean log-width of each flat and tube law central interval. 95% coverage is primary; 50% and other saved levels diagnose calibration.',
            'Null-excitation diagnostic: re-evaluate the same saved null estimates and interval endpoints at stipulated masses [1/64, 1/16, 1/4, 1, 4, 16, 64]. The null observation law is exactly independent of mass, so these are valid additional fixed-truth assessments without fitting again.'
        ],
        'methods': [
            'Retain all methods supplied by the existing comparison engine and their assumptions. Focus interpretation on flat_joint, flat_geometric, flat_median, norm_ratio, norm_floor and the positive Gaussian profile/TLS comparator. The n=4 engine includes covariance-aware controls, so inventories may differ.',
            'Keep the flat-joint point E[f]/E[alpha] unchanged. Compare it with other point summaries of the same law separately from comparing flat and tube probability laws.',
            'Apply the existing quadrature refinement checks. They check selected trials and detect tested discrepancies; passing them is not a numerical error bound for every replicate.',
            'No law-rejection gate, output clipping, failed-trial deletion or selection of favorable conditions is introduced.'
        ],
        'performance_measures': [
            'Primary point loss: min(log(m_hat/m)^2, log(2)^2), assigning the cap to invalid estimates. Primary failure: m_hat outside [m/2, 2m] including invalid estimates. Factor two is an explicit exploratory multiplicative tolerance, not a domain-specific requirement or a guarantee of practical usefulness.',
            'Primary uncertainty metric: fixed-truth coverage of the 95% central interval. Report interval log-width and both miss directions alongside coverage; broad intervals and the no-information boundary must not be treated as successful identification.',
            'Use paired flat-minus-comparator differences and directly computed paired MCSEs within each cell. Interpret n=1 versus n=4 differences as independent Monte Carlo comparisons because the experiments use distinct seeds.',
            'R10000 gives worst-case population MCSE at most 0.005 for individual probabilities, 0.01 for paired probability differences, log(2)^2/200 for individual capped loss and log(2)^2/100 for paired capped loss. Empirical ddof1 MCSE may differ slightly. These are pointwise precision bounds, not simultaneous confidence intervals.',
            'Report signed log error, log-error variance, tolerance curves and invalid counts as diagnostics. Retain established population moment warnings for unbounded losses; their finite simulation averages do not establish finite risk.',
            'Relative prediction scores are inapplicable at zero true force. Bounded mass losses and fixed-truth interval coverage remain applicable there.'
        ]
    },
    'provenance': deepcopy(pilot['provenance']),
    'interpretation': [
        'This is a fixed-size descriptive assessment chosen after a pilot, not a confirmatory global-winner claim. No stopping, extension or tuning in response to observed main-study scores.',
        'Pilot used for design: .tools/mass_estimator_studies/assessment-pilot-20260920, study_checkpoint.json SHA256=' + pilot_hash,
        'All numerical output remains condition-specific. Conditions and null-truth diagnostic points share random draws, so they cannot be pooled as independent observations.',
        'At zero excitation the measurements contain no information about mass. High coverage around one noise-defined mass scale cannot demonstrate calibration for every positive mass.',
        'Nominal conditional probability is not a guarantee of fixed-truth frequentist coverage. Compare both 50% and95% interval coverage with width and tail asymmetry.',
        'The core grid covers known independent Gaussian noise and stable readings only. Calibration sensitivity is assessed in a separate exploratory protocol; correlation, changing excitation and non-Gaussian noise remain outside this study.',
        'Preserve all attempted trials under the declared failure policy and retain failed numerical experiments explicitly. Resume only completed experiments after integrity checks; incomplete experiments restart in a new attempt.'
    ],
    'fixed_replicates_per_cell': 10000,
    'primary_metrics': metrics,
    'experiments': []
}
for index, n in enumerate((1, 4)):
    main['experiments'].append({'id': f'stable_n{n}', 'config': {
        'name': f'Main mass-excitation assessment: {n} stable reading(s)',
        'samples': 10000, 'order': 96, 'batch_size': 128,
        'training_seed': 2026092011 + index*2, 'heldout_seed': 2026092012 + index*2,
        'save_trials': True, 'noise_model': 'isotropic_unit_3d',
        'scenario': deepcopy(scenario), 'measurement': {'repeats': n, 'calibration': {'mode': 'known'}}}})

robust = deepcopy(main)
robust.update(name='Mass estimator: calibration sensitivity pilot', phase='exploratory_pilot',
              fixed_replicates_per_cell=1024, experiments=[])
robust['ademp']['aims'] = [
    'Assess sensitivity to supplied noise calibration while keeping Newton II correct and the physical observation model fixed.',
    'Distinguish point-estimator sensitivity, interval coverage and uncertainty about those simulation results. This fixed pilot informs later robustness studies; it is not comprehensive validation.'
]
robust['ademp']['data_generating_mechanisms'] = [
    'Use the same 3 positive masses and4 per-reading acceleration SNRs as the main assessment; true coordinate errors are independent Gaussian SD1 in both channels.',
    'All experiments have n=4 stable readings of the same latent pair, R=1024 independent replicates per condition and independent training/held-out seeds across calibration regimes. Cross-regime comparisons are independent, not paired.',
    'Calibration regimes: known SDs; both channels supplied at0.75 times true SD; both at1.25 times true SD; force at1.25 and acceleration at0.75 times true SD; and estimated coordinate SDs from30 independent calibration observations with estimated mean.',
    'Estimated calibration is independently redrawn per replicate and per coordinate, using Gaussian sample variance with29 degrees of freedom. This represents repeated deployments with fresh calibration, not repeated measurements conditional on one shared calibration.',
    'Estimation uses these supplied SDs as plug-in inputs. Calibration uncertainty is not integrated into the conditional mass law. Unsupported realized calibration or numerical refinement failure fails the experiment explicitly without dropping trials.'
]
robust['ademp']['estimands'] = main['ademp']['estimands'][:3]
robust['ademp']['performance_measures'] = [item for item in main['ademp']['performance_measures'] if not item.startswith('R10000')]
robust['ademp']['performance_measures'][2] = 'Within each regime use paired flat-minus-comparator differences. Across calibration regimes use independent-simulation uncertainty; no common random numbers or paired regime comparisons are claimed.'
robust['ademp']['performance_measures'].append('R1024 is fixed before execution. Worst-case individual probability MCSE is1.5625 percentage points and paired probability MCSE3.125 points. Treat small differences as unresolved and use the pilot for future size planning; no automatic extension.')
robust['interpretation'] = [
    'Exploratory calibration sensitivity only; the physical law remains correct. All settings and seeds were fixed before this study ran.',
    'Factor-two tolerance is an illustrative scale. Report the existing tighter tolerance curves as secondary context rather than a binary declaration of practical adequacy.',
    'Null cells contain no mass information. Calibration can change the estimator preferred mass scale without making mass identifiable.',
    'Repeated n=4 physical observations and the within-regime truth grid share the same conventions as the main protocol. Different calibration regimes use independent Monte Carlo seeds.',
    'Keep all1024 trials and report invalid outputs, explicit experiment failures and both nominal interval levels. No tuning or global ranking.',
    'Known-covariance anisotropy, correlated channels, changing excitation, shared calibration errors and non-Gaussian data are not covered by this pilot.'
]
robust['primary_metrics'] = [
    {'scope':'point', 'name':'capped_squared_log_factor_2', 'target_mcse':.015},
    {'scope':'point', 'name':'outside_factor_2', 'target_mcse':.032},
    {'scope':'uncertainty', 'name':'coverage_95', 'target_mcse':.016}]
calibrations = [
    ('known', {'mode':'known'}),
    ('both_under', {'mode':'scaled', 'force_scale':[.75]*3, 'acceleration_scale':[.75]*3}),
    ('both_over', {'mode':'scaled', 'force_scale':[1.25]*3, 'acceleration_scale':[1.25]*3}),
    ('differential', {'mode':'scaled', 'force_scale':[1.25]*3, 'acceleration_scale':[.75]*3}),
    ('estimated_k30', {'mode':'estimated', 'samples':30, 'seed':2026092055})]
for index, (label, calibration) in enumerate(calibrations):
    robust['experiments'].append({'id':label, 'config': {
        'name':'Calibration sensitivity: ' + label, 'samples':1024,
        'order':96, 'batch_size':128, 'training_seed':2026092031+index*2,
        'heldout_seed':2026092032+index*2, 'save_trials':True,
        'noise_model':'isotropic_unit_3d', 'scenario':deepcopy(scenario),
        'measurement':{'repeats':4, 'calibration':calibration}}})

for filename, spec in [('mass_excitation_main.json', main), ('calibration_sensitivity_pilot.json', robust)]:
    with (studies / filename).open('x', encoding='utf-8') as handle:
        json.dump(spec, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(filename)

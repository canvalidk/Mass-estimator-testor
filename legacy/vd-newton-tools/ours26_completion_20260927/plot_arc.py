"""Scientific illustration of the proposed reference and readout checks."""
from pathlib import Path
import os
import sys
import json

root = Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR'] = str(root / 'mpl_cache')
sys.path.insert(0, str(root.parent / 'assessment_plot_dependencies'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from arc_followthrough import logshape

data = json.loads((root / 'arc_results.json').read_text())
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.titleweight': 'bold'})
fig, axs = plt.subplots(1, 2, figsize=(12.8, 5.1), constrained_layout=True)
m = np.logspace(-2.5, 2.5, 1501)
for record, label, color in [(data['priors'][0], 'Both reference ratios = 1', '#2563a6'), (data['priors'][3], 'Reference ratios = 0.1 and 10', '#cf6729')]:
    density = np.log(10)*logshape(np.log(m), record['ratios'], record['counts'])/record['arc_length']
    axs[0].plot(m, density, label=label, color=color, lw=2.8)
axs[0].set_xscale('log')
axs[0].set(xlabel='Mass in the chosen units (log scale)', ylabel='Prior density per unit of log10 mass', title='One centre can become two mass bands')
axs[0].legend(frameon=False, fontsize=10, loc='upper right')
axs[0].grid(alpha=.16)
record = data['no_information_padding_existing_readout'][0]
noise = [x['extra_noise_sd'] for x in record['added_independent_zero_push']]
mass = [x['mass'] for x in record['added_independent_zero_push']]
axs[1].axhline(record['original_mass'], label='Before the extra push', color='#2563a6', lw=2)
axs[1].plot(noise, mass, 'o-', color='#8d3d63', lw=2.8, ms=6, label='After a separate zero reading')
axs[1].axhline(record['infinite_noise_limit_mass'], color='#8d3d63', ls=':', lw=1.5, label='Infinite-noise limit')
axs[1].set_xscale('log')
axs[1].set(xlabel='Both noise SDs of the extra push', ylabel='Reported mass', title='Same mass law, different norm readout', ylim=(1.15, 1.37))
axs[1].legend(frameon=False, fontsize=10, loc='center right')
axs[1].grid(alpha=.16)
fig.suptitle('Following the arc proposal through ours26', fontsize=17, weight='bold')
fig.savefig(root / 'arc_followthrough.png', dpi=150, facecolor='white')
print(str(root / 'arc_followthrough.png'))

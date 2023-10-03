import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
font = {'family' : 'arial',
        'size'   : 15}
text = {'color':'black'}
matplotlib.rc('font', **font)
matplotlib.rc('text', **text)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

shapes = ['gaussian', 'constant', 'rectangular']
for shape in shapes:
    sim = pd.read_csv('../../simulations/simulation_data/dfe_successive_mutations/' + shape + '/inference_results/summary.csv', index_col = 0)
    control = pd.read_csv('../../simulations/simulation_data/dfe_successive_mutations_control/' + shape + '/inference_results/summary.csv', index_col = 0)
    
    merged = sim.merge(control, on = 'Epiweek', suffixes = ('_dfe', '_control'))

    ratio = merged['Netau_HMM_median_dfe']/merged['Netau_HMM_median_control']

    cummean = []
    for i in range(len(ratio)):
        cummean.append(np.nanmean(ratio[:i]))
        
    plt.figure()
    plt.plot(merged['Epiweek'],cummean)
    plt.ylim([0,1])
    plt.xlabel('Epiweek (t)')
    plt.ylabel('Cumulative mean\n$N_e^{DFE}(t)/N_e^{ctrl}(t)$')
    plt.title(shape)
    plt.tight_layout()
    plt.savefig('../figure_outputs/simulation_inferred_Netau_dfe_successive_mutations_' + shape + '_ratio_with_control.png', dpi = 300)
    plt.savefig('../figure_outputs/simulation_inferred_Netau_dfe_successive_mutations_' + shape + '_ratio_with_control.pdf')
    plt.show()
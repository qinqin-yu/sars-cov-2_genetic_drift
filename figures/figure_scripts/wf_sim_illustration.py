#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 15:31:31 2022

@author: qinqinyu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 17:05:47 2022

@author: qinqinyu
"""
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas as pd
import glob, os
from epiweeks import Week, Year
from matplotlib.pyplot import cm
#plt.style.use('ggplot')
font = {'family' : 'arial',
        'size'   : 15}
text = {'color':'black'}
matplotlib.rc('font', **font)
matplotlib.rc('text', **text)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

def counts_to_frequency(counts, total_counts):
    # Function for getting frequencies from counts, frac_neutral, and total_counts
    epiweeks = counts.drop(['lineage'], axis = 1).columns
    neutral_frequency = (counts[epiweeks].values/total_counts[epiweeks].values)
    neutral_frequency = pd.DataFrame(data = neutral_frequency, columns = epiweeks)
    neutral_frequency['lineage'] = counts['lineage']
    columns = neutral_frequency.columns.tolist()
    columns = columns[-1:] + columns[:-1]
    neutral_frequency = neutral_frequency[columns]
    return neutral_frequency

path_folder_all = '../../simulations/simulation_data/wf_sim_for_paper_fig1/'

sim_names = ['Net5000_Nseq1000_c5', 'Net500_Nseq1000_c5']
fig, ax = plt.subplots(2, 2, figsize = (10, 6), sharey = True)
i = 0
for sim_name in sim_names:
    path_folder = path_folder_all + sim_name
    
    Net_label = os.path.basename(path_folder)
    idx = Net_label.find('Net')
    idx2 = Net_label.find('_Nseq')
    Net = Net_label[idx+3:idx2]
    
    counts_actual_all = pd.read_csv(path_folder  + '/counts_actual_all.csv', index_col = 0)
    total_actual_counts = pd.read_csv(path_folder  + '/total_actual_counts.csv', index_col = 0)
    
    counts_observed_all = pd.read_csv(path_folder  + '/counts_observed_all.csv', index_col = 0)
    total_observed_counts = pd.read_csv(path_folder  + '/total_observed_counts.csv', index_col = 0)
    
    epiweeks = counts_actual_all.drop(['lineage'], axis = 1).columns
    frequency_actual = counts_to_frequency(counts_actual_all, total_actual_counts)
    ax[0,i].stackplot(epiweeks.astype('float'), frequency_actual[epiweeks], linewidths = 0)
    
    frequency_observed = counts_to_frequency(counts_observed_all, total_observed_counts)
    ax[1,i].stackplot(epiweeks.astype('float'), frequency_observed[epiweeks], linewidths = 0)

    ax[0,i].set_xlim([0, np.max(epiweeks.astype('float'))])
    ax[0,i].set_ylim([0, 1])
    ax[0,i].set_title('$N_e=$ ' + str(Net) + ', $c_t$ = 0')
    
    ax[1,i].set_xlim([0, np.max(epiweeks.astype('float'))])
    ax[1,i].set_ylim([0, 1])
    ax[1,i].set_title('$N_e=$ ' + str(Net) + ', $c_t$ = 5')
    
    i+=1
ax[1,0].set_xlabel('Weeks')
ax[1,1].set_xlabel('Weeks')

ax[0,0].set_ylabel('True\nfrequency')
ax[1,0].set_ylabel('Observed\nfrequency')

plt.tight_layout()
plt.savefig('../figure_outputs/wf_sim_illustration.png', dpi = 300)
plt.savefig('../figure_outputs/wf_sim_illustration.pdf')

plt.show()

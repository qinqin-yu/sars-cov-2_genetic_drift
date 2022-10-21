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

def counts_to_frequency(counts, total_counts, frac_neutral):
    # Function for getting frequencies from counts, frac_neutral, and total_counts
    epiweeks = counts.drop(['lineage'], axis = 1).columns
    neutral_frequency = (counts[epiweeks].values/total_counts[epiweeks].values)/frac_neutral[epiweeks].values
    neutral_frequency = pd.DataFrame(data = neutral_frequency, columns = epiweeks)
    neutral_frequency['lineage'] = counts['lineage']
    columns = neutral_frequency.columns.tolist()
    columns = columns[-1:] + columns[:-1]
    neutral_frequency = neutral_frequency[columns]
    return neutral_frequency
                               
path_folder_all = '../../simulations/simulated_data/negative_binomial_technical_noise_Ne_change_over_time_with_mutation/'
path_folder = path_folder_all + 'gaussian'

Net_label = os.path.basename(path_folder)
numlin_label = 'numlin500_mu0.01_increasing_seq_burnin_multinomial'
summary_all = pd.read_csv(path_folder + '/' + numlin_label + '/trials/inferred_joint_Ne_c_superlineage_combos_threshold_counts_freq/joint_Ne_c_MSD_initial_guess_profile_likelihood_ci_T_9_summary.csv', index_col = 0)
summary = summary_all[summary_all['simulation_trial']==0]
total_neutral_counts = pd.read_csv(path_folder  + '/' + numlin_label + '/trials/total_counts_trial_0.csv', index_col = 0)
neutral_counts = pd.read_csv(path_folder  + '/' + numlin_label + '/trials/neutral_counts_all_trial_0.csv', index_col = 0)
epiweeks = neutral_counts.drop(['lineage'], axis = 1).columns
ones = pd.DataFrame([[1]*len(epiweeks)], columns = epiweeks)
frequency = counts_to_frequency(neutral_counts, total_neutral_counts, ones)
numlineages = np.count_nonzero(neutral_counts[epiweeks], axis=0)

fig, ax = plt.subplots(4,1, figsize = (11, 11), sharex = True)

ax[0].plot(epiweeks.astype('float'), np.transpose(total_neutral_counts[epiweeks]), color = 'k')
ax[0].set_ylabel('Sequences, $M_t$')
#ax[0].set_title(Net_label + ' ' + numlin_label)
ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax[0].set_ylim([0, 3*10**4])

ax[1].stackplot(epiweeks.astype('float'), frequency[epiweeks], linewidths = 0)#, color = colors)
#     ax[1].set_xlabel('Epiweeks')
ax[1].set_ylabel('Observed \n frequency, $f_t^{obs}$')
ax[1].set_xlim([0, 99])
ax[1].set_ylim([0, 1])
ax[1].grid(False)

p = ax[2].plot(summary['Epiweek'], summary['Netau_HMM_median'], label = 'Inferred HMM')
ax[2].fill_between(summary['Epiweek'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
p = ax[2].plot(summary['Epiweek'], summary['Netau_MSD_median'], label = 'Inferred MSD')
ax[2].fill_between(summary['Epiweek'], summary['Netau_MSD_95%_ci_lower'], summary['Netau_MSD_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[2].plot(summary['Epiweek'], summary['Netau_true'], label = 'True')
ax[2].set_yscale('log')
ax[2].set_ylabel('Scaled effective \n population size, $\\tilde{N}_e(t)$')
ax[2].set_ylim([10**2.25, 10**4.75])
ax[2].set_yscale('log')
ax[2].legend(loc = 'upper right')

p = ax[3].plot(summary['Epiweek'], summary['c_MSD_median'], zorder = 10)#, capsize = 5)
ax[3].fill_between(summary['Epiweek'], summary['c_MSD_95%_ci_lower'], summary['c_MSD_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
p = ax[3].plot(summary['Epiweek'], summary['c_median'], zorder = 10)#, capsize = 5)
ax[3].fill_between(summary['Epiweek'], summary['c_95%_ci_lower'], summary['c_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[3].plot(summary['Epiweek'], summary['c_true'])
ax[3].set_ylim([0, 13])
ax[3].set_xlabel('Week')
ax[3].set_ylabel('Measurement noise \n overdispersion, $c_t$')
ax[3].plot()

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figures_for_paper/inference_simulation_v5_compare_MSD.png', dpi = 300)
plt.savefig('../../figures/figures_for_paper/inference_simulation_v5_compare_MSD.pdf')

plt.show()
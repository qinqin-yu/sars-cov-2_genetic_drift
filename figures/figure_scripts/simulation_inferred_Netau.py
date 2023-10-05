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

path_folder = '../../simulations/simulation_data/neutral/gaussian'

Net_label = os.path.basename(path_folder)
summary = pd.read_csv(path_folder + '/inference_results/summary.csv', index_col = 0)
summary_ct1 = pd.read_csv(path_folder + '/inference_results/summary_ct=1.csv', index_col = 0)

true_Netau = pd.read_csv(path_folder + '/true_Netau.csv', index_col = 0)
true_Netau = pd.DataFrame({'Epiweek':true_Netau.columns.astype('float'), 'Netau_true':true_Netau.values[0]})
true_c = pd.read_csv(path_folder + '/true_c.csv', index_col = 0)
true_c = pd.DataFrame({'Epiweek':true_c.columns.astype('float'), 'c_true':true_c.values[0]})
summary = summary.merge(true_Netau, on = 'Epiweek')
summary = summary.merge(true_c, on = 'Epiweek')

total_neutral_counts = pd.read_csv(path_folder  + '/total_counts_lineages.csv', index_col = 0)
neutral_counts = pd.read_csv(path_folder  + '/counts_lineages.csv', index_col = 0)
epiweeks = neutral_counts.drop(['lineage'], axis = 1).columns
ones = pd.DataFrame([[1]*len(epiweeks)], columns = epiweeks)
frequency = counts_to_frequency(neutral_counts, total_neutral_counts, ones)

fig, ax = plt.subplots(4,1, figsize = (10, 10), sharex = True)

ax[0].plot(epiweeks.astype('float'), np.transpose(total_neutral_counts[epiweeks]))
ax[0].set_ylabel('Sequences')
#ax[0].set_title(Net_label + ' ' + numlin_label)
ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

ax[1].stackplot(epiweeks.astype('float'), frequency[epiweeks], linewidths = 0, alpha = 1.0)
#     ax[1].set_xlabel('Epiweeks')
ax[1].set_ylabel('Observed \n frequency, $f_t^{obs}$')
ax[1].set_xlim([0, 99])
ax[1].set_ylim([0, 1])

p = ax[2].plot(summary['Epiweek'], summary['Netau_HMM_median'], label = 'Inferred')
ax[2].fill_between(summary['Epiweek'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[2].plot(summary['Epiweek'], summary['Netau_true'], label = 'True')
p = ax[2].plot(summary_ct1['Epiweek'], summary_ct1['Netau_HMM_median'], label = 'Inferred, $c_t=1$')
ax[2].fill_between(summary_ct1['Epiweek'], summary_ct1['Netau_HMM_95%_ci_lower'], summary_ct1['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[2].set_yscale('log')
ax[2].set_ylabel('Scaled effective\npopulation size,\n$\\tilde{N}_e(t)$')
ax[2].set_ylim([10**2, 10**4.5])
ax[2].set_yscale('log')
ax[2].legend(loc = 'upper right')

p = ax[3].plot(summary['Epiweek'], summary['c_median'], zorder = 10)#, capsize = 5)
ax[3].fill_between(summary['Epiweek'], summary['c_95%_ci_lower'], summary['c_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[3].plot(summary['Epiweek'], summary['c_true'])
ax[3].set_ylim([0, 13])
ax[3].set_xlabel('Epiweek')
ax[3].set_ylabel('Sampling noise, $c_t$')
ax[3].plot()

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../figure_outputs/simulation_inferred_Netau_gaussian.png', dpi = 300)
plt.savefig('../figure_outputs/simulation_inferred_Netau_gaussian.pdf')

plt.show()
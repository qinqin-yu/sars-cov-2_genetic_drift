#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 17:39:48 2022

@author: qinqinyu
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib.lines import Line2D
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
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

path_folder = '../../simulations/simulation_data/stochastic_seir/N1000000_R02_gammaE0.33_gammaI0.18_numlineages100_labeltime75/'
summary = pd.read_csv(path_folder + '/inference_results/summary_infected_plus_exposed.csv', index_col = 0)

infected = pd.read_csv(path_folder + 'counts_infected.csv', index_col = 0)
exposed = pd.read_csv(path_folder + 'counts_exposed.csv', index_col = 0)
neutral_counts = infected.drop(['lineage'], axis = 1) + exposed.drop(['lineage'], axis = 1)
neutral_counts['lineage'] = infected['lineage']
#    counts_all = counts_all[0:10]

infected_all = pd.read_csv(path_folder + 'total_counts_infected.csv', index_col = 0)
exposed_all = pd.read_csv(path_folder + 'total_counts_exposed.csv', index_col = 0)
infected_and_exposed_all = infected_all + exposed_all

epiweeks = neutral_counts.drop(['lineage'], axis = 1).columns
ones = pd.DataFrame([[1]*len(epiweeks)], columns = epiweeks)
frequency = counts_to_frequency(neutral_counts, infected_and_exposed_all, ones)

fig, ax = plt.subplots(2,1, figsize = (10, 6), sharex = True)

ax[0].stackplot(epiweeks.astype('float'), frequency[epiweeks], edgecolor = 'face')
ax[0].set_ylabel('Frequency')
ax[0].set_xlim([1, 39])
ax[0].set_ylim([0, 1])

p = ax[1].plot(summary['Epiweek'], summary['Netau_HMM_median'], label = 'Inferred $\\tilde{N}_e(t)$', linewidth = 2)
ax[1].fill_between(summary['Epiweek'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
ax[1].plot(infected_all.columns.astype('float'), infected_all.values[0], label = '$I(t)$')
ax[1].plot(infected_and_exposed_all.columns.astype('float'), infected_and_exposed_all.values[0], label = '$I(t) + E(t)$')
ax[1].plot(summary['Epiweek'], summary['Netau_SEIR'], label = '$\\tilde{N}_e^{\mathrm{SEIR, eq}}(t)$ analytical')
ax[1].plot(summary['Epiweek'], summary['Netau_SIR'], label = '$\\tilde{N}_e^{\mathrm{SIR}}(t)$ analytical')
ax[1].plot(summary['Epiweek'], summary['Netau_SEIR_numerical'], label = '$\\tilde{N}_e^{\mathrm{SEIR, eq}}(t)$ numerical')
ax[1].plot(summary['Epiweek'], summary['Netau_SIR_numerical'], label = '$\\tilde{N}_e^{\mathrm{SIR}}(t)$ numerical')
ax[1].set_yscale('log')

ax[1].legend(loc = (1.01,-0.1), fontsize = 12)
ax[1].set_xlabel('Week')
fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../figure_outputs/inference_seir_simulation_infected_exposed.png', dpi = 300)
plt.savefig('../figure_outputs/inference_seir_simulation_infected_exposed.pdf')
plt.show()
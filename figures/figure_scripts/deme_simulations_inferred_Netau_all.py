#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 17:16:22 2022

@author: qinqinyu
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob, os
from scipy.stats import median_absolute_deviation, median_test
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

summary_all = pd.read_csv('../../simulations/simulation_data/deme_simulations/inference_results/summary.csv', index_col = 0)
fig, ax = plt.subplots(1, 2, figsize = (11, 4))
for filled_demes_t0 in np.unique(summary_all['Filled demes t0']):
    summary = summary_all[summary_all['Filled demes t0']==filled_demes_t0]
    ax[0].plot(summary['Deme size'], summary['mean_total_counts'], marker = 'o', linestyle = '', label = filled_demes_t0)
    ax[1].plot(summary['Deme size'], summary['Netau_HMM_median'], marker = 'o', linestyle = '', label = filled_demes_t0)

ax[0].set_xlabel('Number of individuals per deme')
ax[0].set_ylabel('Mean # infecteds per week')
ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax[0].set_xlim([-10, 210])

ax[1].set_xlabel('Number of individuals per deme')
ax[1].set_ylabel('Inferred $\\tilde{N}_e(t)$')
ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax[1].set_xlim([-10, 210])

ax[1].legend(loc = (1.02, 0.1), title = 'Initial #\ninfected\ndemes')#, fontsize = 10)
plt.tight_layout()
plt.savefig('../figure_outputs/deme_simulations_inferred_Netau_all.png', dpi = 300)
plt.savefig('../figure_outputs/deme_simulations_inferred_Netau_all.pdf')
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:45:08 2022

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
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.3)

# Profile likelihood
# ADD FILE NAME
summary_all = pd.read_csv('../../simulations/simulation_data/neutral/gaussian/inference_results/ADD FILE NAME.csv', index_col = 0)

Netau_correct_ci = []
c_correct_ci = []
epiweeks = np.unique(summary_all['Epiweek'])
for epiweek in epiweeks:
    summary = summary_all[summary_all['Epiweek']==epiweek]
    cond1 = summary['Netau_HMM_95%_ci_lower']<=summary['Netau_true']
    cond2 = summary['Netau_HMM_95%_ci_upper']>=summary['Netau_true']
    Netau_correct_ci.append(len(summary[cond1&cond2])/np.sum(~np.isnan(summary['Netau_HMM_median'])))
    
    cond3 = summary['c_95%_ci_lower']<=summary['c_true']
    cond4 = summary['c_95%_ci_upper']>=summary['c_true']
    c_correct_ci.append(len(summary[cond3&cond4])/np.sum(~np.isnan(summary['c_median'])))


correct_ci_df = pd.DataFrame()
Netau_correct_ci_df = pd.DataFrame({'epiweek':epiweeks})
Netau_correct_ci_df['parameter'] = '$\\tilde{N}_e(t)$'
Netau_correct_ci_df['correct_ci'] = Netau_correct_ci
correct_ci_df = correct_ci_df.append(Netau_correct_ci_df)

c_correct_ci_df = pd.DataFrame({'epiweek':epiweeks})
c_correct_ci_df['parameter'] = '$c_t$'
c_correct_ci_df['correct_ci'] = c_correct_ci
correct_ci_df = correct_ci_df.append(c_correct_ci_df)
correct_ci_df.reset_index(drop = True, inplace = True)

fig, ax = plt.subplots(1, 2, figsize = (9, 6), sharey = True, gridspec_kw = {'width_ratios':[5,1]})
ax[0].plot(epiweeks, Netau_correct_ci, label = '$\\tilde{N}_e(t)$')
ax[0].plot(epiweeks, c_correct_ci, label = '$c_t$')
ax[0].set_ylim([-0.1, 1.1])
ax[0].set_xlabel('Epiweek')
ax[0].set_ylabel('Fraction of simulations where\n95% CI included true value')
ax[0].legend()

sns.stripplot(x = 'parameter', y = 'correct_ci', data = correct_ci_df, ax = ax[1], marker = '.', color = 'k')
sns.boxplot(x = 'parameter', y = 'correct_ci', data = correct_ci_df, ax = ax[1])
# sns.violinplot(x = 'parameter', y = 'correct_ci', data = correct_ci_df, ax = ax[1])

ax[1].set_xlabel('')
ax[1].set_ylabel('')
plt.text(-0.35, 1.05, str(np.nanmedian(Netau_correct_ci)))
plt.text(0.65, 1.05, str(np.nanmedian(c_correct_ci)))
plt.tight_layout()
plt.savefig('../figure_outputs/simulation_frac_correct_ci.png', dpi = 300)
plt.savefig('../figure_outputs/simulation_frac_correct_ci.pdf')
plt.show()
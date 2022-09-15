#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 17:05:18 2021

@author: qinqinyu
"""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

colors = ['k', '#44AA99', '#88CCEE', '#DDCC77', '#CC6677']#'#66c2a5', '#fc8d62', '#8da0cb', ]
fig, ax = plt.subplots(2, 1, figsize = (10, 6), sharex = True, gridspec_kw = {'height_ratios':[3,1]})

summary_all = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)

labels = []

merged = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
errcond = merged['Netau SIR error criteria met']

p = ax[0].plot(merged[errcond]['date'], 2*merged[errcond]['Netau SIR'], label = '$N_e\\tau$ SIR', color = 'k', marker = 's', markersize = 4)
ax[0].fill_between(merged[errcond]['date'], 2*merged[errcond]['Netau SIR lower'], 2*merged[errcond]['Netau SIR upper'], color=p[0].get_color(), alpha=0.2)
labels.append('All lineages')

i = 1
for name, summary in summary_all.groupby(['variant']):
    errcond = (summary['Netau SIR variant error criteria met'])&(~np.isnan(summary['Netau_HMM_median']))
    errcond.fillna(False, inplace = True)
    p = ax[0].plot(summary['date'], summary['Netau_HMM_median'], zorder = 10, color = colors[i], marker = 'o', markeredgecolor = 'k')#, capsize = 5)
    ax[0].fill_between(summary['date'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.3)
    
    ax[0].plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant'], zorder = 10, color = colors[i], markeredgecolor = 'k', marker = 's')#, markersize = 4)#, capsize = 5)
    ax[0].fill_between(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant lower'], 2*summary[errcond]['Netau SIR variant upper'], color=p[0].get_color(), alpha=0.2)

    ax[1].plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median'], color = p[0].get_color(), marker = '.')#, markeredgecolor = 'k')
    ax[1].fill_between(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant lower']/summary[errcond]['Netau_HMM_95%_ci_upper'], 2*summary[errcond]['Netau SIR upper']/summary[errcond]['Netau_HMM_95%_ci_lower'], color=p[0].get_color(), alpha=0.3)

    i+=1
    labels.append(name)

lines = [Line2D([0], [0], color=c, linestyle='-') for c in colors]
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'o', markeredgecolor = 'k'))
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 's', markeredgecolor = 'k', alpha = 0.2))
labels.append('Inferred $\\tilde{N}_e(t)$')
labels.append('$\\tilde{N}^{\mathrm{SIR}}(t)$')

order = [0,4,2,1,3,5,6]
ax[0].legend([lines[idx] for idx in order],[labels[idx] for idx in order], loc = (1.01, 0.2))#, fontsize = 12)
ax[0].set_yscale('log')
ax[0].set_ylim([10**1.5, 10**6.5])

ax[1].set_ylabel('$\\frac{\\tilde{N}^{\mathrm{SIR}}(t)}{\\tilde{N}_e^{\mathrm{inf}}(t)}$')
ax[1].set_yscale('log')
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_to_sir_model.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_to_sir_model.pdf')
plt.show()
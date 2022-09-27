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
import datetime
from matplotlib.lines import Line2D
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

colors = ['k', '#44AA99', '#88CCEE', '#DDCC77', '#CC6677']#'#66c2a5', '#fc8d62', '#8da0cb', ]
fig, ax = plt.subplots(2, 1, figsize = (10, 8), sharex = True, gridspec_kw = {'height_ratios':[3,1]})

summary_all = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)

labels = []

merged = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)

p = ax[0].plot(merged['date'], merged['I'], color = 'k', marker = 's', alpha = 0.2)#, markersize = 4)
ax[0].fill_between(merged['date'], merged['I_lower'], merged['I_upper'], color=p[0].get_color(), alpha=0.05)

labels.append('All lineages')

i = 1
for name, summary in summary_all.groupby(['variant']):
    p = ax[0].plot(summary['date'], summary['Netau_HMM_median'], zorder = 10, color = colors[i], marker = 'o', markeredgecolor = 'k')#, capsize = 5)
    ax[0].fill_between(summary['date'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.3)
    cases = summary['I']*summary['frac_variant']
    cases_lower = summary['I_lower']*summary['frac_variant']
    cases_upper = summary['I_upper']*summary['frac_variant']
    
    min_cases_to_plot = 10**3
    cases[cases<min_cases_to_plot] = np.nan
    cases_lower[cases_lower<min_cases_to_plot] = np.nan
    cases_upper[cases_upper<min_cases_to_plot] = np.nan

#    print(np.min(cases/summary['Netau_HMM_median']))
#    print(np.max(cases/summary['Netau_HMM_median']))
    
    ax[0].plot(summary['date'], cases, zorder = 10, color = colors[i], markeredgecolor = 'k', marker = 's', alpha = 0.2)#, markersize = 4)#, capsize = 5)
    ax[0].fill_between(summary['date'], cases_lower, cases_upper, color=p[0].get_color(), alpha=0.1)
    
    ax[1].plot(summary['date'], summary['c_median'], color = p[0].get_color(), marker = '.')#, markeredgecolor = 'k')
    ax[1].fill_between(summary['date'], summary['c_95%_ci_lower'], summary['c_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

    i+=1
    labels.append(name)

lines = [Line2D([0], [0], color=c, linestyle='-') for c in colors]
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'o', markeredgecolor = 'k'))
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 's', markeredgecolor = 'k', alpha = 0.2))
labels.append('Inferred $\\tilde{N}_e(t)$')
labels.append('Positives in community')

order = [0,4,2,1,3,5,6]
ax[0].legend([lines[idx] for idx in order],[labels[idx] for idx in order], loc = 'lower right', fontsize = 12, handlelength = 1)
ax[0].set_yscale('log')
ax[0].set_ylabel('Population size')
ax[0].set_ylim([10**1.5, 10**6.5])
ax[0].text(datetime.datetime(2021, 2, 18), 10**2.1, 'Inferred $\\tilde{N}_e(t)$')
ax[0].text(datetime.datetime(2021, 1, 15), 10**6.2, 'Positives in community', alpha = 0.7)

ax[1].set_ylabel('Measurement noise \n overdispersion , $c_t$')
ax[1].set_ylim([0, 25])
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n%Y'))

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_to_positives.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_to_positives.pdf')
plt.show()
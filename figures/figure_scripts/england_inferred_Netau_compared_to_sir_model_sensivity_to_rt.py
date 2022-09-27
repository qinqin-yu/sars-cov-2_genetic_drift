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
fig, ax = plt.subplots(2, 1, figsize = (13, 7), sharex = True, gridspec_kw = {'height_ratios':[3,1]})

summary_all = pd.read_csv('../../data/epidemiological_models/sensitivity_to_rt/epidemiological_models_Netau_inferred_Netau_by_variant_relRtalpha2.7_relRtdelta1.76.csv', parse_dates = ['date'], index_col = 0)
summary_all2 = pd.read_csv('../../data/epidemiological_models/sensitivity_to_rt/epidemiological_models_Netau_inferred_Netau_by_variant_relRtalpha1.1_relRtdelta2.17.csv', parse_dates = ['date'], index_col = 0)

labels = []

merged = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
errcond = merged['Netau SIR error criteria met']

p = ax[0].plot(merged[errcond]['date'], 2*merged[errcond]['Netau SIR'], label = '$N_e\\tau$ SIR', color = 'k', marker = 's', markersize = 4)
ax[0].fill_between(merged[errcond]['date'], 2*merged[errcond]['Netau SIR lower'], 2*merged[errcond]['Netau SIR upper'], color=p[0].get_color(), alpha=0.2)
labels.append('All lineages')

i = 1
for name in np.unique(summary_all['variant']):
    summary = summary_all[summary_all['variant']==name]
    summary2 = summary_all2[summary_all2['variant']==name]

    errcond = (summary['Netau SIR variant error criteria met'])&(~np.isnan(summary['Netau_HMM_median']))
    errcond.fillna(False, inplace = True)
    p = ax[0].plot(summary['date'], summary['Netau_HMM_median'], zorder = 10, color = colors[i], marker = 'o', markeredgecolor = 'k')#, capsize = 5)
    ax[0].fill_between(summary['date'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.3)

    ax[0].plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant'], zorder = 10, color = colors[i], markeredgecolor = 'k', marker = 's')#, markersize = 4)#, capsize = 5)
    ax[0].fill_between(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant lower'], 2*summary[errcond]['Netau SIR variant upper'], color=p[0].get_color(), alpha=0.2)
    
    ax[1].plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median'], color = p[0].get_color(), markeredgecolor = 'k', marker = 's')#, markeredgecolor = 'k')
    ax[1].fill_between(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant lower']/summary[errcond]['Netau_HMM_95%_ci_upper'], 2*summary[errcond]['Netau SIR upper']/summary[errcond]['Netau_HMM_95%_ci_lower'], color=p[0].get_color(), alpha=0.3)

    ax[0].plot(summary2[errcond]['date'], 2*summary2[errcond]['Netau SIR variant'], zorder = 10, color = colors[i], markeredgecolor = 'k', marker = 'v')#, markersize = 4)#, capsize = 5)
    ax[0].fill_between(summary2[errcond]['date'], 2*summary2[errcond]['Netau SIR variant lower'], 2*summary2[errcond]['Netau SIR variant upper'], color=p[0].get_color(), alpha=0.2)
    
    ax[1].plot(summary2[errcond]['date'], 2*summary2[errcond]['Netau SIR variant']/summary2[errcond]['Netau_HMM_median'], color = p[0].get_color(), markeredgecolor = 'k', marker = 'v')#, markeredgecolor = 'k')
    ax[1].fill_between(summary2[errcond]['date'], 2*summary2[errcond]['Netau SIR variant lower']/summary2[errcond]['Netau_HMM_95%_ci_upper'], 2*summary[errcond]['Netau SIR upper']/summary[errcond]['Netau_HMM_95%_ci_lower'], color=p[0].get_color(), alpha=0.3)

    i+=1
    labels.append(name)

lines = [Line2D([0], [0], color=c, linestyle='-') for c in colors]
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'o', markeredgecolor = 'k'))
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 's', markeredgecolor = 'k', alpha = 0.2))
lines.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'v', markeredgecolor = 'k', alpha = 0.2))
labels.append('Inferred $\\tilde{N}_e(t)$')
labels.append('$\\tilde{N}^{\mathrm{SIR}}(t), \\frac{R_0^{\mathrm{Alpha}}}{R_0^{\mathrm{pre-B.1.1.7}}} = 2.7, \\frac{R_0^{\mathrm{Delta}}}{R_0^{\mathrm{pre-B.1.1.7}}} = 1.76$')
labels.append('$\\tilde{N}^{\mathrm{SIR}}(t), \\frac{R_0^{\mathrm{Alpha}}}{R_0^{\mathrm{pre-B.1.1.7}}} = 1.1, \\frac{R_0^{\mathrm{Delta}}}{R_0^{\mathrm{pre-B.1.1.7}}} = 2.17$')
#print
order = [0,4,2,1,3,5,6,7]
ax[0].legend([lines[idx] for idx in order],[labels[idx] for idx in order], loc = (1.01, 0))#, fontsize = 12)
ax[0].set_yscale('log')
ax[0].set_ylim([10**1.5, 10**6.5])

ax[1].set_ylabel('$\\frac{\\tilde{N}^{\mathrm{SIR}}(t)}{\\tilde{N}_e^{\mathrm{inf}}(t)}$')
ax[1].set_yscale('log')
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_sir_model_sensitivity_to_rt.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compared_sir_model_sensitivity_to_rt.pdf')
plt.show()
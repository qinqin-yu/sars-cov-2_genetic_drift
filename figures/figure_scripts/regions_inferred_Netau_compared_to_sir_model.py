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
import datetime
        
colors = ['#44AA99', '#88CCEE', '#DDCC77', '#CC6677']
          
regions = ['East Midlands',
           'East of England',
           'London',
           'North East',
           'North West',
           'South East',
           'South West',
           'West Midlands',
           'Yorkshire and The Humber']
fig, ax = plt.subplots(6, 3, figsize = (18, 15), sharex = True, gridspec_kw = {'height_ratios':[3,2,3,2,3,2]})

row = 0 
col = 0
for region in regions:
    summary_all = pd.read_csv('../../data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)
    
    labels = []
    
    merged = pd.read_csv('../../data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
    
    i = 0
    for name, summary in summary_all.groupby(['variant']):
        errcond = (summary['Netau SIR variant error criteria met'])&(~np.isnan(summary['Netau_HMM_median']))
        errcond.fillna(False, inplace = True)

f        
        ax[row, col].plot(summary[errcond]['date'], summary[errcond]['Netau SIR'], color = p[0].get_color(), marker = 's', markeredgecolor = 'k', alpha = 0.2)
        ax[row, col].fill_between(summary[errcond]['date'], summary[errcond]['Netau SIR lower'], summary[errcond]['Netau SIR upper'], color=p[0].get_color(), alpha=0.2)
        
        ax[row+1, col].plot(summary[errcond]['date'], summary[errcond]['Netau SIR']/summary[errcond]['Netau_HMM_median'], color = p[0].get_color(), marker = '.')
        ax[row+1, col].fill_between(summary[errcond]['date'], summary[errcond]['Netau SIR lower']/summary[errcond]['Netau_HMM_95%_ci_upper'], summary[errcond]['Netau SIR upper']/summary[errcond]['Netau_HMM_95%_ci_lower'], color=p[0].get_color(), alpha=0.2)

        ax[row, col].set_yscale('log')
        ax[row, col].set_ylim([10**1.5, 10**5.5])
#        ax[row, col].set_xlim([datetime.date(2020, 3, 1), datetime.date(2021, 1, 1)])
        ax[row, col].set_title(region)
        ax[row+1, col].set_yscale('log')
        ax[row+1, col].set_ylim([10**0, 10**3.5])
        i+=1
        
    if col == 0:
#        ax[row, col].set_ylabel('$N_e\\tau$')
        ax[row, col].set_ylabel('$\\tilde{N}_e(t)$')
        ax[row+1, col].set_ylabel('$\\frac{\\tilde{N}_e^{\mathrm{SIR}}(t)}{\\tilde{N}_e^{\mathrm{inf}}(t)}$')
    
    if row == 4:
        ax[row, col].xaxis.set_major_formatter(mdates.DateFormatter('%m\n%y'))

    if col == 2:
        row+=2
        col=0
    else:
        col+=1

lines = [Line2D([0], [0], color='darkgray', linestyle = '-', marker = 'o', markeredgecolor = 'k')]
lines.append(Line2D([0], [0], color='lightgray', linestyle = '-', marker = 's', markeredgecolor = 'k'))
labels.append('$\\tilde{N}_e^{\mathrm{inf}}(t)$')
labels.append('$\\tilde{N}_e^{\mathrm{SIR}}(t)$')

#print
order = [0,1]#,3]
ax[0,0].legend([lines[idx] for idx in order],[labels[idx] for idx in order], loc = 'upper left')

#get handles and labels
handles, labels = ax[0,2].get_legend_handles_labels()
#specify order of items in legend
order = [3,1,0,2]
#add legend to plot
ax[0,2].legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize = 12) 

fig.align_ylabels(ax)
plt.tight_layout()
plt.subplots_adjust(hspace = 0.4, wspace = 0.3)
plt.savefig('../figure_outputs/regions_inferred_Netau_compared_to_sir_model.png', dpi = 300)
plt.savefig('../figure_outputs/regions_inferred_Netau_compared_to_sir_model.pdf')
plt.show()
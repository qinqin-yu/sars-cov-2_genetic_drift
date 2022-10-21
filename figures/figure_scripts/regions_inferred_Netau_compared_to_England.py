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
    summary_all_england = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)
    
    labels = []
    
    merged = pd.read_csv('../../data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
    
    merged_england = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
    
    i = 0
    for name, summary in summary_all.groupby(['variant']):
        summary_england = summary_all_england[summary_all_england['variant']==name]
      
        merged = summary.merge(summary_england, on = ['Epiweek', 'date'])
        errcond_merged = (~np.isnan(merged['Netau_HMM_median_x']))&(~np.isnan(merged['Netau_HMM_median_y']))
        errcond_merged.fillna(False, inplace = True)
        
        p = ax[row, col].plot(summary['date'], summary['Netau_HMM_median'], zorder = 10, color = colors[i], marker = '.', label = name)#, capsize = 5)
        ax[row, col].fill_between(summary['date'], summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
        
        ax[row+1, col].plot(merged[errcond_merged]['date'], merged[errcond_merged]['Netau_HMM_median_y']/merged[errcond_merged]['Netau_HMM_median_x'], color = p[0].get_color(), marker = '.')
        ax[row+1, col].fill_between(merged[errcond_merged]['date'], merged[errcond_merged]['Netau_HMM_95%_ci_lower_y']/merged[errcond_merged]['Netau_HMM_95%_ci_upper_x'], merged[errcond_merged]['Netau_HMM_95%_ci_upper_y']/merged[errcond_merged]['Netau_HMM_95%_ci_lower_x'], color=p[0].get_color(), alpha=0.2)

        ax[row+1, col].axhline(y=1, color='k', linestyle=':')
        ax[row+1, col].set_yscale('log')
        ax[row+1, col].set_ylim([10**-1, 10**2])

        p = ax[row, col].plot(summary_england['date'], summary_england['Netau_HMM_median'], zorder = 10, color = colors[i], marker = '.', alpha = 0.1)#, capsize = 5)
        ax[row, col].fill_between(summary_england['date'], summary_england['Netau_HMM_95%_ci_lower'], summary_england['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.05)
   
        ax[row, col].set_yscale('log')
        ax[row, col].set_ylim([10**1.5, 10**4.5])
        ax[row, col].set_title(region)
        i+=1
        
    if col == 0:
        ax[row, col].set_ylabel('$\\tilde{N}_e(t)$')
        ax[row+1, col].set_ylabel('$\\frac{\\tilde{N}_e(t)^{\mathrm{England}}}{\\tilde{N}_e(t)^{\mathrm{Region}}}$')
    
    if row == 4:
        ax[row, col].xaxis.set_major_formatter(mdates.DateFormatter('%m\n%y'))

    if col == 2:
        row+=2
        col=0
    else:
        col+=1   

lines = [Line2D([0], [0], color='k', linestyle = '-', marker = '.', markeredgecolor = 'k')]
lines.append(Line2D([0], [0], color='lightgray', linestyle = '-', marker = '.'))
labels.append('Region')
labels.append('England')

order = [0,1]
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
plt.savefig('../figure_outputs/regions_inferred_Netau_compared_to_England.png', dpi = 300)
plt.savefig('../figure_outputs/regions_inferred_Netau_compared_to_England.pdf')
plt.show()
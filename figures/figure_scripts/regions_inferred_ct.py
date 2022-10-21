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

fig, ax = plt.subplots(9, 1, figsize = (15, 20), sharex = True, sharey = True)#, gridspec_kw = {'height_ratios':[3,2,3,2,3,2]})

row = 0 
col = 0
j = 0
for region in regions:
    summary_all = pd.read_csv('../../data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)
    
    labels = []
    
    merged = pd.read_csv('../../data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
    
    i = 0
    for name, summary in summary_all.groupby(['variant']):
 
        p = ax[j].plot(summary['date'], summary['c_median'], zorder = 10, color = colors[i], marker = 'o', markeredgecolor = 'k', label = name)#, capsize = 5)
        ax[j].fill_between(summary['date'], summary['c_95%_ci_lower'], summary['c_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

        ax[j].set_ylim([0, 25])
        ax[j].set_title(region)
        i+=1
    j+=1

ax[4].set_ylabel('Measurement noise, $c_t$')
ax[8].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))

#get handles and labels
handles, labels = ax[0].get_legend_handles_labels()
#specify order of items in legend
order = [3,1,0,2]
#add legend to plot
ax[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize = 12) 

fig.align_ylabels(ax)
plt.tight_layout()
plt.subplots_adjust(hspace = 0.4, wspace = 0.3)
plt.savefig('../figure_outputs/regions_inferred_ct.png', dpi = 300)
plt.savefig('../figure_outputs/regions_inferred_ct.pdf')
plt.show()
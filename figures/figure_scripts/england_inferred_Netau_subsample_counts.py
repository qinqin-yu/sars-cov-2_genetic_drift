#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 10:29:47 2022

@author: qinqinyu
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from epiweeks import Week, Year
import datetime
from datetime import timedelta
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)
import os

def epiweeks_to_dates(epiweeks):
    dates = []
    for w in epiweeks:
        w = int(np.floor(float(w)))
        if w<=53:
            week = Week(2020, w)
        elif (w>53)&(w<=105):
            week = Week(2021, w-53)
        else:
            week = Week(2022, w-53-52)
        dates.append(week.startdate()+timedelta(days=3))
    # 2019-01-05
    dates = np.array(dates)
    return dates

path_folders = {'../../data/lineages/pre_B-1-177/microreact/':'pre-B.1.177',
                '../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/':'B.1.177',
                '../../data/lineages/alpha/alpha|2021-06-20|61.5/':'Alpha',
                '../../data/lineages/delta/delta|2022-01-25|49.5+58.5/':'Delta'}

fig, ax = plt.subplots(2, 2, figsize = (14, 8))

i=0
j=0
for path_folder in path_folders:
    label = path_folders[path_folder]
    
    if os.path.exists(path_folder + '/is_pillar_2/England/inference_results/summary_corrected.csv'):
        summary = pd.read_csv(path_folder + '/is_pillar_2/England/inference_results/summary_corrected.csv', index_col = 0)
        summary_random_subsample = pd.read_csv(path_folder + '/is_pillar_2/England_random_subsample_one_tenth/inference_results/summary_corrected.csv', index_col = 0)
    else:
        summary = pd.read_csv(path_folder + '/is_pillar_2/England/inference_results/summary.csv', index_col = 0)
        summary_random_subsample = pd.read_csv(path_folder + '/is_pillar_2/England_random_subsample_one_tenth/inference_results/summary.csv', index_col = 0)
        
    p = ax[i,j].plot(epiweeks_to_dates(summary['Epiweek']), summary['Netau_HMM_median'], zorder = 10, label = 'All sequences')#, capsize = 5)
    ax[i,j].fill_between(epiweeks_to_dates(summary['Epiweek']), summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
    
    p = ax[i,j].plot(epiweeks_to_dates(summary_random_subsample['Epiweek']), summary_random_subsample['Netau_HMM_median'], zorder = 10, label = '10% random subsample \nof sequences')#, capsize = 5)
    ax[i,j].fill_between(epiweeks_to_dates(summary_random_subsample['Epiweek']), summary_random_subsample['Netau_HMM_95%_ci_lower'], summary_random_subsample['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
        
    ax[i,j].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n%y'))
    ax[i,j].set_ylabel('Inferred $\\tilde{N}_e(t)$')
    ax[i,j].set_yscale('log')
    ax[i,j].set_title(label)
    if j==0:
        j+=1
    else:
        i+=1
        j=0

ax[0,1].legend(loc = (1.01, 0.1))
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_subsample_counts_logy.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_subsample_counts_logy.pdf')
plt.show()
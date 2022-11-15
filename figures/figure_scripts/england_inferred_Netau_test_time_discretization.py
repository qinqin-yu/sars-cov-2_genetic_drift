#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:58:56 2022

@author: qinqinyu
"""
import sys
bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import math
import pandas as pd
import misc_useful_functions as muf
import os
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

variant_param_folders = ['../../data/lineages/pre_B-1-177/microreact/is_pillar_2/',
                         '../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/is_pillar_2/',
                         '../../data/lineages/alpha/alpha|2021-06-20|61.5/is_pillar_2/',
                         '../../data/lineages/delta/delta|2022-01-25|49.5+58.5/is_pillar_2/']
variant_labels = ['pre-B.1.177', 'B.1.177', 'Alpha', 'Delta']

fig, ax = plt.subplots(2,2, figsize = (12,8))
row = 0
col = 0
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i in range(len(variant_param_folders)):
    variant_param_folder = variant_param_folders[i]
    variant_label = variant_labels[i]
    if os.path.exists(variant_param_folder + 'England/inference_results/summary_corrected.csv'): 
        delt1 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_corrected.csv', index_col = 0)
        delt2 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_corrected_test_time_discretization.csv', index_col = 0)
    else:
        delt1 = pd.read_csv(variant_param_folder + 'England/inference_results/summary.csv', index_col = 0)
        delt2 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_test_time_discretization.csv', index_col = 0)
    
    delt2['Epiweek']=delt2['Epiweek']
    
    dates = muf.epiweeks2dates(delt1['Epiweek'])
    delt1['date'] = dates
    
    dates = muf.epiweeks2dates(delt2['Epiweek'])
    delt2['date'] = dates

    p = ax[row, col].plot(delt1['date'], delt1['Netau_HMM_median'], zorder = 10, color = 'tab:blue', marker = 'o', markeredgecolor = 'k', label = 'every week')#, capsize = 5)
    ax[row, col].fill_between(delt1['date'], delt1['Netau_HMM_95%_ci_lower'], delt1['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

    p = ax[row, col].plot(delt2['date'], delt2['Netau_HMM_median']*2, zorder = 10, color = 'tab:orange', marker = 'o', markeredgecolor = 'k', label = 'every other week')#, capsize = 5)
    ax[row, col].fill_between(delt2['date'], delt2['Netau_HMM_95%_ci_lower']*2, delt2['Netau_HMM_95%_ci_upper']*2, color=p[0].get_color(), alpha=0.2)

    ax[row, col].set_yscale('log')
    ax[row, col].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n\'%y'))
    ax[row, col].set_ylabel('Inferred $\\tilde{N}_e(t)$')
    ax[row, col].set_title(variant_label)
    
    if col == 1:
        row=1
        col=0
    else:
        col+=1
ax[0,1].legend(title = 'Timepoints used', loc = (1, 0))
plt.tight_layout()
plt.savefig('../figure_outputs/england_inferred_Netau_test_time_discretization.pdf')
plt.savefig('../figure_outputs/england_inferred_Netau_test_time_discretization.png', dpi = 300)

    
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
import numpy as np
import math
import pandas as pd
import misc_useful_functions as muf

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
        cts20 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_corrected.csv', index_col = 0)
        cts10 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_corrected_superlineage_size_10_counts.csv', index_col = 0)
        cts40 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_corrected_superlineage_size_40_counts.csv', index_col = 0)
    else:
        cts20 = pd.read_csv(variant_param_folder + 'England/inference_results/summary.csv', index_col = 0)
        cts10 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_superlineage_size_10_counts.csv', index_col = 0)
        cts40 = pd.read_csv(variant_param_folder + 'England/inference_results/summary_superlineage_size_40_counts.csv', index_col = 0)
    
    dates = muf.epiweeks2dates(cts10['Epiweek'])
    cts10['date'] = dates
    
    dates = muf.epiweeks2dates(cts20['Epiweek'])
    cts20['date'] = dates
    
    dates = muf.epiweeks2dates(cts40['Epiweek'])
    cts40['date'] = dates

    p = ax[row, col].plot(cts10['date'], cts10['Netau_HMM_median'], zorder = 10, color = 'tab:blue', marker = 'o', markeredgecolor = 'k', label = '10')#, capsize = 5)
    ax[row, col].fill_between(cts10['date'], cts10['Netau_HMM_95%_ci_lower'], cts10['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

    p = ax[row, col].plot(cts20['date'], cts20['Netau_HMM_median'], zorder = 10, color = 'tab:orange', marker = 'o', markeredgecolor = 'k', label = '20')#, capsize = 5)
    ax[row, col].fill_between(cts20['date'], cts20['Netau_HMM_95%_ci_lower'], cts20['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

    p = ax[row, col].plot(cts40['date'], cts40['Netau_HMM_median'], zorder = 10, color = 'tab:green', marker = 'o', markeredgecolor = 'k', label = '40')#, capsize = 5)
    ax[row, col].fill_between(cts40['date'], cts40['Netau_HMM_95%_ci_lower'], cts40['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
    
    ax[row, col].set_yscale('log')
    ax[row, col].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n\'%y'))
    ax[row, col].set_ylabel('Inferred $\\tilde{N}_e(t)$')
    ax[row, col].set_title(variant_label)
    
    if col == 1:
        row=1
        col=0
    else:
        col+=1
ax[0,1].legend(title = 'Superlineage \nthreshold counts', loc = (1, 0))
plt.tight_layout()
plt.savefig('../figure_outputs/england_inferred_Netau_compare_superlineage_size.pdf')
plt.savefig('../figure_outputs/england_inferred_Netau_compare_superlineage_size.png', dpi = 300)

    
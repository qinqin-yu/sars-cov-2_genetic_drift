#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:09:16 2022

@author: qinqinyu
"""

import sys
bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
    
import pandas as pd
import numpy as np
import glob
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

import misc_useful_functions as muf

path_folder = '../../data/lineages/'
folders = []

variant_folders = glob.glob(path_folder + '*/')
for variant_folder in variant_folders:
    folders = folders + glob.glob(variant_folder + '*/')
folders = np.sort(np.array(folders))

fig, ax = plt.subplots(2,2, figsize = (11,8))

i = 0
j = 0

folders = ['pre_B-1-177/microreact/',
           'B-1-177/B-1-177|2021-02-22|694.5/',
           'alpha/alpha|2021-06-20|61.5/',
           'delta/delta|2022-01-25|49.5+58.5/']

label_dict = {'pre_B-1-177':'pre-B.1.177',
              'B-1-177':'B.1.177',
              'alpha':'Alpha',
              'delta':'Delta'}

for folder in folders:
    variant_param = folder.split('/')[-3:-1]
    folder_inference_results = path_folder + folder + 'is_pillar_2/England/inference_results/'
    if os.path.exists(folder_inference_results + 'summary_corrected.csv'):
        results = pd.read_csv(folder_inference_results + 'summary_corrected.csv', index_col = 0)
        results_exclude_ben_lin = pd.read_csv(folder_inference_results + 'summary_corrected_exclude_ben_lineages_deterministic.csv', index_col = 0)
        results_exclude_ben_lin_50 = pd.read_csv(folder_inference_results + 'summary_corrected_exclude_ben_lineages_deterministic_50.csv', index_col = 0)
        results_exclude_ben_lin_75 = pd.read_csv(folder_inference_results + 'summary_corrected_exclude_ben_lineages_deterministic_75.csv', index_col = 0)
        results_exclude_ben_lin_90 = pd.read_csv(folder_inference_results + 'summary_corrected_exclude_ben_lineages_deterministic_90.csv', index_col = 0)
    else:
        results = pd.read_csv(folder_inference_results + 'summary.csv', index_col = 0)
        results_exclude_ben_lin = pd.read_csv(folder_inference_results + 'summary_exclude_ben_lineages_deterministic.csv', index_col = 0)
        results_exclude_ben_lin_50 = pd.read_csv(folder_inference_results + 'summary_exclude_ben_lineages_deterministic_50.csv', index_col = 0)
        results_exclude_ben_lin_75 = pd.read_csv(folder_inference_results + 'summary_exclude_ben_lineages_deterministic_75.csv', index_col = 0)
        results_exclude_ben_lin_90 = pd.read_csv(folder_inference_results + 'summary_exclude_ben_lineages_deterministic_90.csv', index_col = 0)

    results['date'] = muf.epiweeks2dates(results['Epiweek'])
    results_exclude_ben_lin['date'] = muf.epiweeks2dates(results_exclude_ben_lin['Epiweek'])
    results_exclude_ben_lin_50['date'] = muf.epiweeks2dates(results_exclude_ben_lin_50['Epiweek'])
    results_exclude_ben_lin_75['date'] = muf.epiweeks2dates(results_exclude_ben_lin_75['Epiweek'])
    results_exclude_ben_lin_90['date'] = muf.epiweeks2dates(results_exclude_ben_lin_90['Epiweek'])

    p = ax[i,j].plot(results['date'], results['Netau_HMM_median'], label = 'All lineages')
    ax[i,j].fill_between(results['date'], results['Netau_HMM_95%_ci_lower'], results['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p[0].get_color())
    
    # p2 = ax[i,j].plot(results_exclude_ben_lin['date'], results_exclude_ben_lin['Netau_HMM_median'], label = 'Excluding non-neutral lineages inferred \nfrom conservative model')
    # ax[i,j].fill_between(results_exclude_ben_lin['date'], results_exclude_ben_lin['Netau_HMM_95%_ci_lower'], results_exclude_ben_lin['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p2[0].get_color())

    p3 = ax[i,j].plot(results_exclude_ben_lin_50['date'], results_exclude_ben_lin_50['Netau_HMM_median'], label = 'Excluding 50th percentile non-neutral lineages')
    ax[i,j].fill_between(results_exclude_ben_lin_50['date'], results_exclude_ben_lin_50['Netau_HMM_95%_ci_lower'], results_exclude_ben_lin_50['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p3[0].get_color())
    
    p4 = ax[i,j].plot(results_exclude_ben_lin_75['date'], results_exclude_ben_lin_75['Netau_HMM_median'], label = 'Excluding 75th percentile non-neutral lineages')
    ax[i,j].fill_between(results_exclude_ben_lin_75['date'], results_exclude_ben_lin_75['Netau_HMM_95%_ci_lower'], results_exclude_ben_lin_75['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p4[0].get_color())
    
    p5 = ax[i,j].plot(results_exclude_ben_lin_90['date'], results_exclude_ben_lin_90['Netau_HMM_median'], label = 'Excluding 90th percentile non-neutral lineages')
    ax[i,j].fill_between(results_exclude_ben_lin_90['date'], results_exclude_ben_lin_90['Netau_HMM_95%_ci_lower'], results_exclude_ben_lin_90['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p5[0].get_color())
    
    ax[i,j].set_yscale('log')
    ax[i,j].set_ylabel('Inferred $N_e\\tau$')
    ax[i,j].set_title(label_dict[variant_param[0]])
    ax[i,j].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n\'%y'))
#    ax[i,j].xaxis.set_major_locator(plt.MaxNLocator(5))
    
    if (i==0)&(j==0):
        ax[i,j].legend(fontsize = 12)
    
    if j == 1:
        i += 1
        j = 0
    else:
        j+=1

plt.tight_layout()
plt.savefig('../figure_outputs/england_inferred_Netau_exclude_ben_lineages_deterministic_summary_fitness_cutoff.png', dpi = 300)
plt.savefig('../figure_outputs/england_inferred_Netau_exclude_ben_lineages_deterministic_summary_fitness_cutoff.pdf')

plt.show()
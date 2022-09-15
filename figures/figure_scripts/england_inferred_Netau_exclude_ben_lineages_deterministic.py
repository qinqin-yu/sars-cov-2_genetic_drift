#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:09:16 2022

@author: qinqinyu
"""

import pandas as pd
import numpy as np
import glob
import os
import matplotlib
import matplotlib.pyplot as plt
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

path_folder = '/Users/qinqinyu/Dropbox/sars-cov-2_github/data/lineages/'
folders = []

variant_folders = glob.glob(path_folder + '*/')
for variant_folder in variant_folders:
    folders = folders + glob.glob(variant_folder + '*/')
folders = np.sort(np.array(folders))

fig, ax = plt.subplots(4,3, figsize = (15,13))
i = 0
j = 0
for folder in folders:
    variant_param = folder.split('/')[-3:-1]
    folder_inference_results = folder + 'is_pillar_2/England/inference_results/'
    if os.path.exists(folder_inference_results + 'summary_corrected.csv'):
        results = pd.read_csv(folder_inference_results + 'summary_corrected.csv', index_col = 0)
        results_exclude_ben_lin = pd.read_csv(folder_inference_results + 'summary_corrected_exclude_ben_lineages_deterministic.csv', index_col = 0)
    else:
        results = pd.read_csv(folder_inference_results + 'summary.csv', index_col = 0)
        results_exclude_ben_lin = pd.read_csv(folder_inference_results + 'summary_exclude_ben_lineages_deterministic.csv', index_col = 0)

    p = ax[i,j].plot(results['Epiweek'], results['Netau_HMM_median'], label = 'All lineages')
    ax[i,j].fill_between(results['Epiweek'], results['Netau_HMM_95%_ci_lower'], results['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p[0].get_color())
    
    p2 = ax[i,j].plot(results_exclude_ben_lin['Epiweek'], results_exclude_ben_lin['Netau_HMM_median'], label = 'Exclude ben lineages, deterministic model')
    ax[i,j].fill_between(results_exclude_ben_lin['Epiweek'], results_exclude_ben_lin['Netau_HMM_95%_ci_lower'], results_exclude_ben_lin['Netau_HMM_95%_ci_upper'], alpha = 0.2, color = p2[0].get_color())
    
    ax[i,j].set_yscale('log')
    ax[i,j].set_xlabel('Epiweek')
    ax[i,j].set_ylabel('Inferred $N_e\\tau$')
    
    if variant_param[1] == 'microreact':
        ax[i,j].set_title(variant_param[0])
    else:
        ax[i,j].set_title(variant_param[1])
            
    if (i==0)&(j==0):
        ax[i,j].legend(fontsize = 12)
    
    if j == 2:
        i += 1
        j = 0
    else:
        j+=1


plt.tight_layout()
plt.savefig('../figure_outputs/england_inferred_Netau_exclude_ben_lineages_deterministic.png', dpi = 300)
plt.savefig('../figure_outputs/england_inferred_Netau_exclude_ben_lineages_deterministic.pdf')

plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 10:49:54 2022

@author: qinqinyu
"""

import numpy as np
import pandas as pd
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from epiweeks import Week, Year
import datetime
from datetime import timedelta
import re
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

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

def figure_label_tree_cut(label):
    label_sep = label.split('|')
    variant = label_sep[0]
    tree_date = label_sep[1]
    dcut_str = label_sep[2]
    dcut_str_sep = dcut_str.split('+')
    variant_label = variant_labels[variant]
    dcut_labels = []
    for dcut_str_sep_i in dcut_str_sep:
        dcut_label = float(dcut_str_sep_i)*(3.34*10**(-5))
        dcut_labels.append("{:.3e}".format(dcut_label))
    if len(dcut_labels) == 1: 
        label_for_figure = '{' + tree_date + ', ' + variant_label + ', $d_{\mathrm{cut}}$=' + str(dcut_labels[0]) + '}'
    elif len(dcut_labels) == 2: 
        label_for_figure = '{' + tree_date + ', ' + variant_label + ', $d_{\mathrm{cut}}^{\mathrm{(1)}}$=' + str(dcut_labels[0]) + ', $d_{\mathrm{cut}}^{\mathrm{(2)}}$=' + str(dcut_labels[1]) + '}'
    return label_for_figure, variant, tree_date, dcut_labels

path_folder_all = '../../data/lineages/'
path_folder_variants = glob.glob(path_folder_all + '*/')
path_folder_variants.remove('../../data/lineages/pre_B-1-177/')

fig, ax = plt.subplots(1, 1, figsize = (15, 4), sharex = True)

for path_folder_variant in path_folder_variants:
    path_folders = glob.glob(path_folder_variant + '*/')    
    path_folders = np.sort(path_folders)
    
    colors = {'B-1-177':['#88CCEE','#8A8AED'], 'alpha':['#017504', '#85A546', '#44AA99', '#147575'], 'delta':['#D4DB79','#D88F79','#DDCC77', '#D15252']}
    variant_labels = {'B-1-177':'B.1.177', 'alpha':'Alpha', 'delta':'Delta'}
    linestyles = ['-', '--', ':', '-.']
    linestyle_idx = {'B-1-177':0, 'alpha':0, 'delta':0}
    dcut_idx = {'B-1-177':0, 'alpha':0, 'delta':0}
    for path_folder in path_folders:
        
        label = path_folder.split('/')[-2]
        label_for_figure, variant, tree_date, dcut_labels = figure_label_tree_cut(label)
        
        summary = pd.read_csv(path_folder + '/is_pillar_2/England/inference_results/summary_corrected.csv', index_col = 0)
        
        color = colors[variant][dcut_idx[variant]]
        
        p = ax.plot(epiweeks_to_dates(summary['Epiweek']), summary['Netau_HMM_median'], zorder = 10, color = color, label = label_for_figure)#, capsize = 5, linestyle = linestyles[linestyle_idx[variant]])
        ax.fill_between(epiweeks_to_dates(summary['Epiweek']), summary['Netau_HMM_95%_ci_lower'], summary['Netau_HMM_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)
        linestyle_idx[variant] = linestyle_idx[variant] + 1
        dcut_idx[variant]+=1
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    ax.set_ylabel('Inferred $\\tilde{N}_e(t)}$')
    ax.set_yscale('log')
    ax.legend(loc = (1.01, -0.15), fontsize = 12, title = '{Tree date, Variant, $d_{\mathrm{cut}}$ (substitutions per site)}')
    
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compare_tree.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_Netau_compare_tree.pdf')
plt.show()
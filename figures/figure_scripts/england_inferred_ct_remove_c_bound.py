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
        'size'   : 15}
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

path_folder_all = '../../data/lineages/'
path_folder_variants = glob.glob(path_folder_all + '*/')
path_folder_variants.remove('../../data/lineages/pre_B-1-177/')

fig, ax = plt.subplots(1, 1, figsize = (13, 4), sharex = True)

path_folders = ['../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5',
                '../../data/lineages/alpha/alpha|2021-06-20|61.5',
                '../../data/lineages/delta/delta|2022-01-25|49.5+58.5']
labels = ['B.1.177', 'Alpha', 'Delta']
colors = ['#88CCEE', '#44AA99', '#DDCC77']

i = 0 

for path_folder in path_folders:

    summary = pd.read_csv(path_folder + '/is_pillar_2/England/inference_results/summary_remove_c_bound.csv', index_col = 0)

    p = ax.plot(epiweeks_to_dates(summary['Epiweek']), summary['c_median'], marker = '.', color = colors[i], zorder = 10, label = labels[i])
    ax.fill_between(epiweeks_to_dates(summary['Epiweek']), summary['c_95%_ci_lower'], summary['c_95%_ci_upper'], color=p[0].get_color(), alpha=0.2)

    i+=1
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d\n%Y'))
ax.set_ylabel('Measurement noise\noverdispersion, $c_t$')
ax.set_ylim([0,20])
ax.legend()
    
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_inferred_ct_remove_c_bound.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/england_inferred_ct_remove_c_bound.pdf')
plt.show()
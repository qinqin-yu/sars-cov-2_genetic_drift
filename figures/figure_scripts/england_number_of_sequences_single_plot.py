#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 2023

@author: qinqinyu
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import math
import pandas as pd
from epiweeks import Week, Year
import datetime
from datetime import timedelta
import seaborn as sns
#plt.style.use('ggplot')
font = {'family' : 'arial',
        'size'   : 15}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

fig, ax = plt.subplots(1, 1, figsize = (15, 4))

path_folders = {'../../data/lineages/pre_B-1-177/microreact/':'pre-B.1.177',
               '../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/':'B.1.177',
               '../../data/lineages/alpha/alpha|2021-06-20|61.5/':'Alpha',
                '../../data/lineages/delta/delta_sublineage50.5+58.5/':'Delta'}
               # '../../data/lineages/delta/delta|2022-01-25|49.5+58.5/':'Delta'}

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
i=0

for path_folder in path_folders:
    total_assigned_counts = pd.read_csv(path_folder + 'is_pillar_2/England/total_counts_lineages.csv', index_col = 0)
    total_variant_counts = pd.read_csv(path_folder + 'is_pillar_2/England/total_counts_variant_metadata.csv', index_col = 0)
    total_metadata_counts = pd.read_csv(path_folder + 'is_pillar_2/England/total_counts_metadata.csv', index_col = 0)
    epiweeks = total_assigned_counts.columns
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
    dates = np.array(dates)
    # ax.plot(dates, total_assigned_counts[epiweeks].values[0], color = colors[i], linestyle = '-', label = 'Used for analysis')
    ax.plot(dates, total_variant_counts[epiweeks].values[0], color = colors[i], linestyle = '-', label = path_folders[path_folder])

    i+=1

ax.plot(dates, total_metadata_counts[epiweeks].values[0], color = 'k', linestyle = '-', label = 'Total')

ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m\n%d\n%y'))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax.set_ylim([0, 4*10**4])
ax.set_yscale('log')
ax.set_ylim([10**0, 10**5])
ax.legend(loc = (1.01, 0))
ax.set_ylabel('Sequences')
# ax.set_title(path_folders[path_folder])
 
fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_number_of_sequences_single_plot.pdf')
plt.savefig('../../figures/figure_outputs/england_number_of_sequences_single_plot.png', dpi = 300)
plt.show()

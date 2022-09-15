#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 22:13:14 2021

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

fig, ax = plt.subplots(2, 2, figsize = (10, 6))

path_folders = {'../../data/lineages/pre_B-1-177/microreact/':'pre-B.1.177',
               '../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/':'B.1.177',
               '../../data/lineages/alpha/alpha|2021-06-20|61.5/':'Alpha',
               '../../data/lineages/delta/delta|2022-01-25|49.5+58.5/':'Delta'}

i = 0
j = 0

for path_folder in path_folders:
    total_assigned_counts = pd.read_csv(path_folder + 'is_pillar_2/England/total_counts_lineages.csv', index_col = 0)
    total_variant_counts = pd.read_csv(path_folder + 'is_pillar_2/England/total_counts_variant_metadata.csv', index_col = 0)
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
    ax[i,j].plot(dates, total_assigned_counts[epiweeks].values[0], color = 'k', linestyle = '-', label = 'Used for analysis')
    ax[i,j].plot(dates, total_variant_counts[epiweeks].values[0], color = 'b', linestyle = '--', label = 'Total surveillance')
    ax[i,j].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%y'))
    ax[i,j].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax[i,j].set_yticklabels([])
    
    if (i==0)&(j==0):
        ax[0,0].set_xlim([dates[0], dates[-1]])
    elif (i==0)&(j==1):
        ax[0,1].set_xlim([datetime.date(2020, 8, 1), datetime.date(2021, 2, 15)])
    elif (i==1)&(j==0):
        ax[1,0].set_xlim([datetime.date(2020, 11, 1), datetime.date(2021, 6, 1)])
    else:
        ax[1,1].set_xlim([datetime.date(2021, 4, 1), dates[-1]])

    ax[i,j].set_ylim([0, 4*10**4])
    ax[i,j].set_yscale('log')
    ax[i,j].set_ylim([10**0, 10**5])
    ax[i,j].legend()
    ax[i,j].set_ylabel('Sequences')
    ax[i,j].set_title(path_folders[path_folder])
    if j==1:
        i+=1
        j=0
    else:
        j+=1
fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_number_of_sequences.pdf')
plt.savefig('../../figures/figure_outputs/england_number_of_sequences.png', dpi = 300)
plt.show()

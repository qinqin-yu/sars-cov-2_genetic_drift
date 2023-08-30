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
import pandas as pd
from epiweeks import Week
import datetime
from datetime import timedelta

font = {'family' : 'arial',
        'size'   : 15}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

def counts_to_frequency(counts, total_counts):
    # Function for getting frequencies from counts and total_counts
    epiweeks = counts.drop(['lineage'], axis = 1).columns
    neutral_frequency = (counts[epiweeks].values/total_counts[epiweeks].values)
    neutral_frequency = pd.DataFrame(data = neutral_frequency, columns = epiweeks)
    neutral_frequency['lineage'] = counts['lineage']
    columns = neutral_frequency.columns.tolist()
    columns = columns[-1:] + columns[:-1]
    neutral_frequency = neutral_frequency[columns]
    return neutral_frequency

fig, ax = plt.subplots(3, 2, figsize = (10, 8))#, sharey = True)

# B-1-177
path_folder = '../../data/lineages/B-1-177/B-1-177|2021-02-22|695.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[0,0].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[0,0].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))
ax[0,0].set_xlim([datetime.date(2020, 8, 10), datetime.date(2021, 2, 15)])
ax[0,0].set_ylim([0, 1])
ax[0,0].set_ylabel('Observed\nfrequency')
ax[0,0].set_title('B.1.177, $d_{\mathrm{cut}}$=2.323e-02')

# B-1-177
path_folder = '../../data/lineages/B-1-177/B-1-177_sublineage692.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[0,1].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[0,1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[0,1].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))
ax[0,1].set_xlim([datetime.date(2020, 8, 10), datetime.date(2021, 2, 15)])
ax[0,1].set_ylim([0, 1])
ax[0,1].set_ylabel('Observed\nfrequency')
ax[0,1].set_title('B.1.177, $d_{\mathrm{cut}}$=2.313e-02')

# Alpha
path_folder = '../../data/lineages/alpha/alpha|2021-06-20|61.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[1,0].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[1,0].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[1,0].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))
ax[1,0].set_xlim([datetime.date(2020, 11, 1), datetime.date(2021, 5, 20)])
ax[1,0].set_ylim([0, 1])
ax[1,0].set_ylabel('Observed\nfrequency')
ax[1,0].set_title('Alpha, $d_{\mathrm{cut}}$=2.054e-03')

# Alpha
path_folder = '../../data/lineages/alpha/alpha_sublineage59.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[1,1].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[1,1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[1,1].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))
ax[1,1].set_xlim([datetime.date(2020, 11, 1), datetime.date(2021, 5, 20)])
ax[1,1].set_ylim([0, 1])
ax[1,1].set_ylabel('Observed\nfrequency')
ax[1,1].set_title('Alpha, $d_{\mathrm{cut}}$=1.987e-03')

# Delta
path_folder = '../../data/lineages/delta/delta|2022-01-25|50.5+58.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[2,0].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[2,0].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[2,0].set_xlim([datetime.date(2021, 4, 1), dates[-1]])
ax[2,0].set_ylim([0, 1])
ax[2,0].set_ylabel('Observed\nfrequency')
ax[2,0].set_title('2022-01-25, Delta, \n$d_{\mathrm{cut}}^{(1)}$=1.687e-03, $d_{\mathrm{cut}}^{(2)}$=1.954e-03')
ax[2,0].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))

# Delta
path_folder = '../../data/lineages/delta/delta_sublineage50.5/is_pillar_2/England/'
neutral_counts = pd.read_csv(path_folder + 'counts_lineages.csv', index_col = 0)
total_assigned_counts = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
neutral_frequency = counts_to_frequency(neutral_counts, total_assigned_counts)
epiweeks = neutral_frequency.drop(['lineage'], axis = 1).columns
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

ax[2,1].stackplot(dates, neutral_frequency[epiweeks], edgecolor = 'face')
ax[2,1].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
ax[2,1].set_xlim([datetime.date(2021, 4, 1), dates[-1]])
ax[2,1].set_ylim([0, 1])
ax[2,1].set_ylabel('Observed\nfrequency')
ax[2,1].set_title('2022-03-25, Delta, \n$d_{\mathrm{cut}}$=1.687e-03')
ax[2,1].xaxis.set_major_locator(mdates.MonthLocator(interval = 2))

fig.align_ylabels(ax)
plt.tight_layout()
plt.savefig('../../figures/figure_outputs/england_lineage_frequencies_tree_cut.pdf')
plt.savefig('../../figures/figure_outputs/england_lineage_frequencies_tree_cut.png', dpi = 300)
plt.show()

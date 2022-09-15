#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 11:29:25 2022

@author: qinqinyu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 17:05:18 2021

@author: qinqinyu
"""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib.lines import Line2D
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

colors = ['#44AA99', '#88CCEE', '#DDCC77', '#CC6677']#'#66c2a5', '#fc8d62', '#8da0cb', ]
fig, ax = plt.subplots(1, 1, figsize = (15, 5))

summary_all = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_inferred_Netau_by_variant.csv', parse_dates = ['date'], index_col = 0)

labels = []

merged = pd.read_csv('../../data/epidemiological_models/epidemiological_models_Netau_overall.csv', parse_dates = ['date'], index_col = 0)
errcond = merged['Netau SIR error criteria met']

i = 0
for name, summary in summary_all.groupby(['variant']):
    errcond = (summary['Netau SIR variant error criteria met'])&(~np.isnan(summary['Netau_HMM_median']))
    errcond.fillna(False, inplace = True)

    ax.plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median'], color = colors[i], label = name)#, markeredgecolor = 'k')    
    p = ax.plot(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median'], color = colors[i], marker = 'o', markeredgecolor = 'k')#, markeredgecolor = 'k')
    ax.fill_between(summary[errcond]['date'], 2*summary[errcond]['Netau SIR variant lower']/summary[errcond]['Netau_HMM_95%_ci_upper'], 2*summary[errcond]['Netau SIR upper']/summary[errcond]['Netau_HMM_95%_ci_lower'], color=p[0].get_color(), alpha=0.3)
    print(2*np.min(summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median']))
    print(2*np.max(summary[errcond]['Netau SIR variant']/summary[errcond]['Netau_HMM_median']))

    i+=1

# Plotting literature values
p = ax.plot([datetime.date(2020, 1, 1), datetime.date(2020, 2, 27)], [65, 65], label = '$\sigma^2$, Worldwide excluding China (Endo et al., \nWellcome Open Research, 2020)', color = 'k', linestyle = ':')
ax.fill_between([datetime.date(2020, 1, 1), datetime.date(2020, 2, 27)], [33.75, 33.75], [127.5, 127.5], color=p[0].get_color(), alpha=0.3)

p = ax.plot([datetime.date(2020, 8, 1), datetime.date(2020, 9, 1)], [7.07, 7.07], label = '$\sigma^2$, UK (Quilty et al., \nReport for SPI-M-O and SAGE, 2021)', color = 'k', linestyle = '--')
ax.fill_between([datetime.date(2020, 8, 1), datetime.date(2020, 9, 1)], [2.65, 2.65], [44.5, 44.5], color=p[0].get_color(), alpha=0.3)

p = ax.plot([datetime.date(2021, 1, 1), datetime.date(2021, 2, 1)], [1.42, 1.42], color = 'k', linestyle = '--')
ax.fill_between([datetime.date(2021, 1, 1), datetime.date(2021, 2, 1)], [0.66, 0.66], [9.19, 9.19], color=p[0].get_color(), alpha=0.3)

handles, labels = plt.gca().get_legend_handles_labels()
handles.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'o', markeredgecolor = 'k'))
labels.append('$\\frac{\\tilde{N}(t)^{\mathrm{SIR}}}{\\tilde{N}_e(t)^{\mathrm{inf}}}$' )
order = [6, 3, 1, 0, 2, 4, 5]
ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = (1.01, 0), fontsize = 18)
ax.set_yscale('log')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))

plt.tight_layout()
plt.savefig('../../figures/figure_outputs/overdispersion_comparison_to_literature_v5.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/overdispersion_comparison_to_literature_v5.pdf')
plt.show()
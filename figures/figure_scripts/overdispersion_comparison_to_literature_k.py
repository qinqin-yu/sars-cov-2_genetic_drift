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

summary_all['sigma2'] = 2*summary_all['Netau SIR variant']/summary_all['Netau_HMM_median']
summary_all['sigma2_lower'] = 2*summary_all['Netau SIR variant lower']/summary_all['Netau_HMM_95%_ci_upper']
summary_all['sigma2_upper'] = 2*summary_all['Netau SIR upper']/summary_all['Netau_HMM_95%_ci_lower']

summary_all['k'] = summary_all['Rt_variant']/(summary_all['sigma2']/summary_all['Rt_variant']-1)
summary_all['k_lower'] = summary_all['Rt_variant_lower']/(summary_all['sigma2_upper']/summary_all['Rt_variant_lower']-1)
summary_all['k_upper'] = summary_all['Rt_variant_upper']/(summary_all['sigma2_lower']/summary_all['Rt_variant_upper']-1)

i = 0
for name, summary in summary_all.groupby(['variant']):
    errcond = (summary['Netau SIR variant error criteria met'])&(~np.isnan(summary['Netau_HMM_median']))
    errcond.fillna(False, inplace = True)

    ax.plot(summary[errcond]['date'], summary[errcond]['k'], color = colors[i], label = name)#, markeredgecolor = 'k')    
    p = ax.plot(summary[errcond]['date'], summary[errcond]['k'], color = colors[i], marker = 'o', markeredgecolor = 'k')#, markeredgecolor = 'k')
    ax.fill_between(summary[errcond]['date'], summary[errcond]['k_lower'], summary[errcond]['k_upper'], color=p[0].get_color(), alpha=0.3)
    print(np.min(summary[errcond]['k']))
    print(np.max(summary[errcond]['k']))

    i+=1

# Plotting literature values
p = ax.plot([datetime.date(2020, 1, 1), datetime.date(2020, 2, 27)], [0.1, 0.1], label = '$k$, Worldwide excluding China (Endo et al., \nWellcome Open Research, 2020)', color = 'k', linestyle = ':')
ax.fill_between([datetime.date(2020, 1, 1), datetime.date(2020, 2, 27)], [0.05, 0.05], [0.2, 0.2], color=p[0].get_color(), alpha=0.3)

p = ax.plot([datetime.date(2020, 8, 1), datetime.date(2020, 9, 1)], [0.25, 0.25], label = '$k$, UK (Quilty et al., \nReport for SPI-M-O and SAGE, 2021)', color = 'k', linestyle = '--')
ax.fill_between([datetime.date(2020, 8, 1), datetime.date(2020, 9, 1)], [0.15, 0.15], [0.39, 0.39], color=p[0].get_color(), alpha=0.3)

p = ax.plot([datetime.date(2021, 1, 1), datetime.date(2021, 2, 1)], [0.33, 0.33], color = 'k', linestyle = '--')
ax.fill_between([datetime.date(2021, 1, 1), datetime.date(2021, 2, 1)], [0.23, 0.23], [0.39, 0.39], color=p[0].get_color(), alpha=0.3)

handles, labels = plt.gca().get_legend_handles_labels()
handles.append(Line2D([0], [0], color='gray', linestyle = '', marker = 'o', markeredgecolor = 'k'))
#labels.append('$\\frac{\\tilde{N}(t)^{\mathrm{SIR}}}{\\tilde{N}_e(t)^{\mathrm{inf}}}$' )
labels.append('Inferred overdispersion parameter, $k$')
order = [6, 3, 1, 0, 2, 4, 5]
ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = (1.01, 0), fontsize = 18)
ax.set_yscale('log')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))

plt.tight_layout()
plt.savefig('../../figures/figure_outputs/overdispersion_comparison_to_literature_k.png', dpi = 300)
plt.savefig('../../figures/figure_outputs/overdispersion_comparison_to_literature_k.pdf')
plt.show()
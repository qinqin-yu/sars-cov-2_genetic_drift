#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:04:29 2022

@author: qinqinyu
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from epiweeks import Week, Year
import datetime
from datetime import timedelta
font = {'family' : 'arial',
        'size'   : 18}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.3)

summary_all = pd.read_csv('../../simulations/simulation_data/deme_simulations/inference_results/summary.csv', index_col = 0)
fig, ax = plt.subplots(1, 1, figsize = (5, 4))
for filled_demes_t0 in np.unique(summary_all['Filled demes t0'])[1:]:
    summary = summary_all[summary_all['Filled demes t0']==filled_demes_t0]
    ax.plot(summary['Deme size'], summary['mean_total_counts']/summary['Netau_HMM_median'], marker = 'o', linestyle = '', label = int(filled_demes_t0))

ax.set_xlabel('Number of individuals per deme')
ax.set_ylabel('Ratio of mean infected to\ninferred $\\tilde{N}_e(t)$, $\\frac{\\langle I\\rangle}{\\tilde{N}_e(t)^{\mathrm{inf}}}$')

ax.legend(title = 'Initial #\ninfected\ndemes', loc = 'lower right')
#plt.xlim([-10, 210])
#plt.ylim([-1, 31])
plt.xlim([0, 210])
plt.ylim([0, 31])
plt.tight_layout()
plt.savefig('../figure_outputs/deme_simulations_inferred_Netau.png', dpi = 300)
plt.savefig('../figure_outputs/deme_simulations_inferred_Netau.pdf')
plt.show()

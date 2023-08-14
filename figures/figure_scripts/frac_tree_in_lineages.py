import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from epiweeks import Week
import datetime
from datetime import timedelta
import matplotlib
import matplotlib.dates as mdates
font = {'family' : 'arial',
        'size'   : 15}
matplotlib.rc('font', **font)
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)

variant_path_folders = ['../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/',
                       '../../data/lineages/alpha/alpha|2021-06-20|61.5/',
                       '../../data/lineages/delta/delta|2022-01-25|50.5+58.5/']
labels = ['B.1.177', 'Alpha', 'Delta']

fig, ax = plt.subplots(1,3, figsize = (20,5))

i = 0

for variant_path_folder in variant_path_folders: 

    counts_lineages = pd.read_csv(variant_path_folder + '/is_pillar_2/England/total_counts_lineages.csv', index_col = 0)
    counts_tree = pd.read_csv(variant_path_folder + '/is_pillar_2/England/total_counts_variant_tree.csv', index_col = 0)
    
    epiweeks = counts_lineages.columns
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
    ax[i].plot(dates, (counts_lineages/counts_tree).values[0])
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%y'))
    ax[i].set_title(labels[i])
    i+=1

ax[0].set_ylabel('Fraction of sequences in tree\nassigned to a lineage')
plt.tight_layout()
plt.savefig('../figure_outputs/frac_tree_in_lineages.png', dpi = 300)
plt.savefig('../figure_outputs/frac_tree_in_lineages.pdf')

plt.show()
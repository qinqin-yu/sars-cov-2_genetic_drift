#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 14:55:00 2022

@author: qinqinyu
"""

import pandas as pd
import numpy as np
import math

# Particular combinations of variants/tree date/tree cut depth, focusing on England
variant_param_folders = ['../data/lineages/pre_B-1-177/microreact/',
                         '../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/',
                         '../data/lineages/alpha/alpha|2021-06-20|61.5/',
                         '../data/lineages/delta/delta|2022-01-25|49.5+58.5/']

total_counts =[]
for variant_param_folder in variant_param_folders:
    path_folder = variant_param_folder + '/is_pillar_2/England/'
    df = pd.read_csv(path_folder + 'total_counts_lineages.csv')
    total_counts_variant = df.values[0].sum()
    total_counts.append(total_counts_variant)
    print(variant_param_folder + ', total counts: ' + str(total_counts_variant))

print('total counts overall: ' + str(np.array(total_counts).sum()))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:58:58 2022

@author: qinqinyu
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import os
sns.set(style = 'whitegrid', font_scale = 1.5)

def wright_fisher(Ne, fj1):
    nj2 = np.random.binomial(Ne, fj1)
    fj2 = nj2/Nei
    return nj2, fj2

def sampling_noise(Nseq, c, fj):    
    if fj == 0:
        nobsj = 0
    elif fj == 1:
        nobsj = Nseq
    else:
        m = fj*Nseq
        v = m*c
        p = m/v
        n = m**2/(v-m)
        nobsj = np.random.negative_binomial(n,p)
    fobsj = nobsj/Nseq
    return nobsj, fobsj

def wright_fisher_all_lineages(Ne, f1):
    n2 = np.empty(f1.shape)
    f2 = np.empty(f1.shape)
    for j in range(len(f1)):
        fj1 = f1[j]
        nj2, fj2 = wright_fisher(Ne, fj1)
        n2[j] = nj2
        f2[j] = fj2
    return n2, f2
        
def sampling_noise_all_lineages(Nseq, c, f):
    nobs = np.empty(f.shape)
    fobs = np.empty(f.shape)
    for j in range(len(f)):
        fj = f[j]
        nobsj, fobsj = sampling_noise(Nseq, c, fj)
        nobs[j] = nobsj
        fobs[j] = fobsj
    return nobs, fobs

# Parameters
Neis = [500, 5000]
Nseqi = 1000#10**3
c = 5

total_epiweeks = 20
t = np.arange(0, total_epiweeks)
t_str = t.astype('float').astype('str')
numlineages = 50
f0 = np.array([1/numlineages]*numlineages)

path_folder = '../simulation_data/wf_sim_for_paper_fig1/'
   
### START
for Nei in Neis:
    n0 = f0*Nei
    output_folder = path_folder + '/Net' + str(Nei) + '_Nseq' + str(Nseqi) + '_c' + str(c) + '/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
            
    counts_sim_all = pd.DataFrame()
    counts_sim_all['lineage'] = np.array(range(len(f0)))
    
    counts_sim_actual_all = pd.DataFrame()
    counts_sim_actual_all['lineage'] = np.array(range(len(f0)))
    
    lineage_counter = len(f0)
    #initialize
    f = f0
    n = n0
    #for each timepoint
    for j in range(0, len(t)):
        epiweek = t_str[j]
        if j == 0:
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
        else:
            n, f = wright_fisher_all_lineages(Nei, f)
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
        #save to df
        counts_sim_all[epiweek] = nobs
        counts_sim_actual_all[epiweek] = n
    
    counts_sim_all.fillna(0, inplace = True)
    counts_sim_actual_all.fillna(0, inplace = True)
    counts_sim_all.reset_index(drop = True, inplace = True)
    counts_sim_actual_all.reset_index(drop = True, inplace = True)
    
    counts_sim_all = counts_sim_all.copy()
    total_neutral_counts_sim = pd.DataFrame(counts_sim_all[t_str].sum(axis = 0))
    total_neutral_counts_sim = total_neutral_counts_sim.transpose()

    total_neutral_counts_sim_actual = pd.DataFrame(counts_sim_actual_all[t_str].sum(axis = 0))
    total_neutral_counts_sim_actual = total_neutral_counts_sim_actual.transpose()
    
    counts_sim_all.to_csv(output_folder + '/counts_observed_all.csv')
    counts_sim_actual_all.to_csv(output_folder + '/counts_actual_all.csv')
    total_neutral_counts_sim.to_csv(output_folder + '/total_observed_counts.csv')
    total_neutral_counts_sim_actual.to_csv(output_folder + '/total_actual_counts.csv')
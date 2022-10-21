#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:21:45 2021

@author: qinqinyu
"""
import sys

bin_path = '../../bin'
if bin_path not in sys.path:
    sys.path.insert(1, '../../bin')
    
import pandas as pd
import scipy as sp 
import numpy as np 
import re
import time
import glob
import os
import matplotlib.pyplot as plt
import format_data_v2 as fd
import analytical_likelihoods_no_emission as lh

# HMM INFERENCE

# Parameters
mincount = 20 # Minimum superlineage size
p_transition=lambda x2, x1, f, beta: lh.normal(x2, x1, f, beta) # Normal 
#p_transition=lambda x2, x1, f, beta: lh.sqrt_normal(x2, x1, f, beta) #Sqrt normal
    
path_folder = '../../simulation_data/deme_simulations/'
output_folder = path_folder + 'inference_results/'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)   

# Importing data
filenames = glob.glob(path_folder + 'counts_t0_42_tl_50/demes_560000.0_lineages_200_filleddemes0_*_demeN_*_counts_t0_42_tl_50.csv')
for filename in filenames:
    #filename = 'counts_t0_42_tl_50/demes_560000.0_lineages_200_filleddemes0_5000_demeN_100_counts_t0_42_tl_50.csv'
    filename_only = os.path.basename(filename)
    filename_no_ending = filename_only[:filename_only.find('.csv')]
    counts=pd.read_csv(filename, index_col = 0)
    counts_original = counts.copy()
    
    epiweeks = counts.drop(['lineage'], axis = 1).columns
    
    d1=epiweeks[0]
    dl=epiweeks[-1]
    
#    counts = counts[0:10]
    total_counts = counts[epiweeks].sum(axis = 0)
    total_counts = total_counts.to_frame().transpose()
    #frac_neutral = pd.read_csv(path_folder + 'frac_neutral.csv', index_col = 0)
    
    numtrials = 100
    Netau_HMM_all = []
    Netau_HMM_lower_all = []
    Netau_HMM_upper_all = []

    for i in range(numtrials):
        print(i)
        # Creating superlineages and formatting the data
        counts = fd.create_superlineages(counts_original, d1, dl, mincount, seed = np.random.randint(10**5))
        data_list = fd.create_data_list(counts, total_counts, epiweeks)
        
        # Inference
        shareddata = lh.getshareddata(total_counts)
        hmm_maxlikelihood = lambda x: lh.hmm_maxlikelihood_Ne(x[0], x[1])
        
        Netau_HMM = hmm_maxlikelihood([data_list, shareddata])
        Netau_HMM_lower, Netau_HMM_upper = confidence_interval_Ne(data_list, shareddata)
        Netau_HMM_all.append(Ne)
        Netau_HMM_lower_all.append(Ne_lower)
        Netau_HMM_upper_all.append(Ne_upper)

    df = pd.DataFrame(data = {'Netau_HMM':Netau_HMM_all, 'Netau_HMM_lower':Netau_HMM_lower_all, 'Netau_HMM_upper':Netau_HMM_upper_all})
    df.to_csv(output_folder + filename_no_ending + '_raw.csv')
    
# SAVE SUMMARY OF RESULTS
path_folder = '../../simulation_data/deme_simulations/inference_results/'
files = glob.glob(path_folder + 'demes_*')
summary = pd.DataFrame()

for file in files:
    filename = os.path.basename(file)
    idx = filename.find('filleddemes0_')
    idx2 = filename.find('demeN')
    filled_demes0 = float(filename[idx+13:idx2-1])
    idx3 = filename.find('_counts_')
    demeN = float(filename[idx2+6:idx3])
    df = pd.read_csv(file, index_col = 0)
    if len(df)>0:
#         Ne = df['Ne']
        cond = df['Ne']<10**5
        Ne_median = np.median(df[cond]['Netau_HMM'])
        Ne_2_5_percentile = np.nanmedian(df[cond]['Netau_HMM_lower'])
        Ne_97_5_percentile = np.nanmedian(df[cond]['Netau_HMM_upper'])
        if Ne_2_5_percentile==Ne_97_5_percentile:
            Ne_median = np.nan
            Ne_2_5_percentile = np.nan
            Ne_97_5_percentile = np.nan
    else:
        Ne_median = np.nan
        Ne_2_5_percentile = np.nan
        Ne_97_5_percentile = np.nan
        

    counts=pd.read_csv('../../simulation_data/deme_simulations/counts_t0_42_tl_50/demes_560000.0_lineages_200_filleddemes0_' + str(int(filled_demes0)) + '_demeN_' + str(int(demeN)) + '_counts_t0_42_tl_50.csv', index_col = 0)
    epiweeks = counts.drop(['lineage'], axis = 1).columns
    total_counts = counts[epiweeks].sum(axis = 0).to_frame().transpose()
    mean_total_counts = np.mean(total_counts[epiweeks].values[0])
    
    summary = summary.append(pd.DataFrame({'Filled demes t0':[filled_demes0], \
                                           'Deme size':[demeN], \
                                           'Netau_HMM_median':[Ne_median], \
                                           'Netau_HMM_95%_ci_lower':[Ne_2_5_percentile], \
                                           'Netau_HMM_95%_ci_upper':[Ne_97_5_percentile], \
                                           'mean_total_counts':[mean_total_counts]}))
    
summary = summary.sort_values(by = ['Filled demes t0'])
summary.reset_index(drop = True, inplace = True)
summary.to_csv(path_folder + '/summary.csv')
     

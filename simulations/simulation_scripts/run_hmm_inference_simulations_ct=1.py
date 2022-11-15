#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 16:40:05 2022

@author: qinqinyu
"""

import sys

bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
    
import os
import glob
from multiprocessing import Pool
import pandas as pd
import numpy as np
import analytical_likelihoods as lh
import format_data as fd
import technical_noise_functions as tnf

import hmm_inference

def infer_Ne_c1(path_folder, \
                      counts_filename, \
                      total_counts_filename, \
                      output_folder, \
                      output_filename, \
                      excluded_lineages = [], lineage_col_name = 'lineage', T=9, mincount = 20, minfreq = 0.01, numtrials = 100, delta_t = 1):
    
    # Currnetly, moving window must be an odd number (can update in future instatiations of code)
    if T%2 == 0:
        print('Moving window (T) must be an odd number') 
        return
    
    # Importing data
    counts_all=pd.read_csv(path_folder + counts_filename, index_col = 0)
    counts_all = counts_all[~counts_all[lineage_col_name].isin(excluded_lineages)]
    counts_all.reset_index(drop = True, inplace = True)
    
    total_sampled_counts_all = pd.read_csv(path_folder + total_counts_filename, index_col = 0)
#    counts_all = counts_all[0:10]

    epiweeks_all = counts_all.drop([lineage_col_name], axis = 1).columns
    
    # Drop times that are not every delta t interval
    epiweeks_all_updated = epiweeks_all[::delta_t]
    epiweeks_to_drop = np.array([epiweek for epiweek in epiweeks_all if epiweek not in epiweeks_all_updated])
    counts_all.drop(labels = epiweeks_to_drop, axis = 'columns', inplace =True)
    total_sampled_counts_all.drop(labels = epiweeks_to_drop, axis = 'columns', inplace =True)
    epiweeks_all = epiweeks_all_updated
    T = int((T-1)/delta_t+1)

    # Renaming epiweek format to floats (easier to work with later on)
    epiweeks_all_rename_dict = {}
    for epiweek in epiweeks_all:
        epiweeks_all_rename_dict[epiweek] = str(float(epiweek))
    counts_all.rename(columns = epiweeks_all_rename_dict, inplace = True)
    total_sampled_counts_all.rename(columns = epiweeks_all_rename_dict, inplace = True)
    epiweeks_all = counts_all.drop([lineage_col_name], axis = 1).columns
#    epiweeks_all = epiweeks_all[70:70+T]

    # Initialize variables to save results to
    df_out = pd.DataFrame()
    
    # Looping through time windows
    for i in range(len(epiweeks_all)-(T-1)):
        epiweeks = epiweeks_all[i:i+T]
        print(str(i+1) + ' out of ' + str(len(epiweeks_all)-(T-1)) + ' times')
        d1=epiweeks[0]
        dl=epiweeks[-1]
        
        # Get data from this time window
        counts = counts_all[epiweeks]
        counts = counts.loc[~(counts==0).all(axis=1)] # Drop any lineages that have 0 counts throughout the time window
        
        counts_original = counts.copy()
        total_sampled_counts = total_sampled_counts_all[epiweeks]
        
        # Set c = 1 for all timepoints
        c_all = pd.DataFrame(data = [[1]*len(epiweeks)], columns = epiweeks)
        shareddata = lh.getshareddata(c_all,total_sampled_counts)
        
        # Looping through trials of bootstrapping lineages
        for j in range(numtrials):
            print('\t superlineage combo ' + str(j+1) + ' out of ' + str(numtrials))
            df_out_j = pd.DataFrame()
            
            # Creating superlineages and formatting the data
            counts, frequency = fd.create_superlineages_counts_and_freq(counts_original, total_sampled_counts, d1, dl, mincount, minfreq, seed = np.random.randint(10**5))
            
            if len(counts)>1: # Check that there is at least one superlineage
                # Inference
                data_list = fd.create_data_list(counts, total_sampled_counts, epiweeks)
                k_df = tnf.get_kappa(counts, epiweeks, total_sampled_counts, mincount = mincount)
                if np.all(k_df['stderr']>0):
                    hmm_maxlikelihood = lambda x: lh.hmm_maxlikelihood_Ne(x[0], x[1])
                    Netau_HMM = hmm_maxlikelihood([data_list, shareddata])
                    
                    # Confidence interval estimation of Ne using profile likelihood
                    Netau_HMM_lower, Netau_HMM_upper = lh.confidence_interval_Ne_fixed_c(data_list, shareddata)
                
                else:
                    Netau_HMM = np.nan
                    Netau_HMM_lower = np.nan
                    Netau_HMM_upper = np.nan
            
            else:
                Netau_HMM = np.nan
                Netau_HMM_lower = np.nan
                Netau_HMM_upper = np.nan
            
            df_out_j['epiweek_start'] = [d1]
            df_out_j['epiweek_middle'] = [float(d1)+(float(dl)-float(d1))/2]
            df_out_j['epiweek_end'] = [dl]
            df_out_j['superlineage_combo'] = [j]
            df_out_j['Netau_HMM'] = [Netau_HMM]
            df_out_j['Netau_HMM_lower'] = [Netau_HMM_lower]
            df_out_j['Netau_HMM_upper'] = [Netau_HMM_upper]
            df_out_j['T'] = (T-1)*delta_t+1
            df_out_j['delta_t'] = delta_t
            df_out = df_out.append(df_out_j)
    df_out.reset_index(inplace = True, drop = True)
    df_out.to_csv(output_folder + output_filename)

def summarize_results_c1(raw_filename, output_folder, output_filename):
    df_all = pd.read_csv(raw_filename, index_col = 0)
    epiweeks_recorded = np.unique(df_all[['epiweek_start', 'epiweek_middle', 'epiweek_end']])
    epiweeks = np.arange(np.min(epiweeks_recorded), np.max(epiweeks_recorded)+1)
    T = df_all['T'].values[0]
    delta_t = df_all['delta_t'].values[0]
    T = int((T-1)/delta_t+1)
    
    summary = pd.DataFrame()
    for epiweek, df in df_all.groupby('epiweek_middle'):
        e1 = df['epiweek_start'].values[0]
        el = df['epiweek_end'].values[0]

        cond = df['Netau_HMM']<10**5
        Ne_median = np.nanmedian(df[cond]['Netau_HMM'])
        Ne_2_5_percentile = np.nanmedian(df[cond]['Netau_HMM_lower'])
        Ne_97_5_percentile = np.nanmedian(df[cond]['Netau_HMM_upper'])

        summary = summary.append(pd.DataFrame({'Epiweek':[epiweek], \
                                               'Netau_HMM_median':[Ne_median], \
                                               'Netau_HMM_95%_ci_lower':[Ne_2_5_percentile], \
                                               'Netau_HMM_95%_ci_upper':[Ne_97_5_percentile]}))

    summary['T'] = (T-1)*delta_t+1
    summary['delta_t'] = delta_t
    summary.sort_values('Epiweek', inplace = True)
    summary.reset_index(inplace = True, drop = True)
    summary.dropna(subset = ['Netau_HMM_median'], how = 'all', inplace = True)
    summary.to_csv(output_folder + output_filename)
    
# Filenames - change to fit needs
output_filename = 'raw_ct=1.csv'
counts_filename = 'counts_lineages.csv'
total_counts_filename = 'total_counts_lineages.csv'

# GET LIST OF PARAMETERS TO RUN ANALYSES ON
# All combinations of variants/tree date/tree cut depth for England
params_all = []
path_folder = '../simulation_data/neutral/gaussian/'
output_folder = path_folder + 'inference_results/'
params = (path_folder, counts_filename, total_counts_filename, output_folder, output_filename)
   
#    params_all = params_all[0:1] # for testing, comment out when doing full run

# MAKE OUTPUT FOLDERS
output_folder = params[3]
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
  
# RUN INFERENCE PIPELINE

infer_Ne_c1(path_folder, \
      counts_filename, \
      total_counts_filename, \
      output_folder, \
      output_filename)
        
# SUMMARIZE RESULTS
    
raw_filename = output_folder + output_filename
output_filename = 'summary_ct=1.csv'
summarize_results_c1(raw_filename, output_folder, output_filename)
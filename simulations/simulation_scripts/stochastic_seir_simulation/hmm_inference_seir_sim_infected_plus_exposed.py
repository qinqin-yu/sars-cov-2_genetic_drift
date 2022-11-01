#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:21:45 2021

@author: qinqinyu

Jointly infers Ne and c using MLE from WF simulations. Simulations have burn-in
time, multinomial transition probabilties, and negative binomial emission probabilities.
Initial guess for Ne and c for optimization fromm MSD method. 
Estimates confidence interval of Ne from profile likelihood and confidence
interval of c from multiple superlineage combinations. Data preprocessing include
creation of superlineages based on count and frequency threshold. 

"""
import sys

bin_path = '../../bin'
if bin_path not in sys.path:
    sys.path.insert(1, '../../bin')
    
import pandas as pd
import numpy as np 
import os
import glob
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import Pool, get_context
import format_data_v3 as fd
import analytical_likelihoods_no_emission_v2 as lh
import technical_noise_functions_v2 as tnf

def run_inference(path_folder, output_folder, T=9, mincount = 20, minfreq = 0.01, numtrials = 20):
#    trial_num_idx = filename_trial.find('_trial_')
#    trial_num_idx2 = filename_trial.find('.csv')
#    trial_num = int(filename_trial[trial_num_idx+7:trial_num_idx2])
                
    # Importing data
    infected = pd.read_csv(path_folder + 'counts_infected.csv', index_col = 0)
    exposed = pd.read_csv(path_folder + 'counts_exposed.csv', index_col = 0)
    counts_all = infected.drop(['lineage'], axis = 1) + exposed.drop(['lineage'], axis = 1)
    counts_all['lineage'] = infected['lineage']
#    counts_all = counts_all[0:10]
    
    infected_all = pd.read_csv(path_folder + 'total_counts_infected.csv', index_col = 0)
    exposed_all = pd.read_csv(path_folder + 'total_counts_exposed.csv', index_col = 0)
    total_sampled_counts_all = infected_all + exposed_all
    
    epiweeks_all = counts_all.drop(['lineage'], axis = 1).columns
#    epiweeks_all = epiweeks_all[60:60+T]
    df_out = pd.DataFrame()
    c_index = np.array(range(T))
    c_index_str = ['c' + str(i) for i in c_index]
    c_MSD_index_str = ['c' + str(i) + '_MSD' for i in c_index]

    for i in range(len(epiweeks_all)-(T-1)):
        epiweeks = epiweeks_all[i:i+T]
#        if i%10 == 0:
        print(str(i+1) + ' out of ' + str(len(epiweeks_all)-(T-1)) + ' times')
        
        d1=epiweeks[0]
        dl=epiweeks[-1]
        
        # Get data from this time window
        counts = counts_all[epiweeks]
        counts_original = counts.copy()
        total_sampled_counts = total_sampled_counts_all[epiweeks]
        
        # Looping through trials
        for j in range(numtrials):
#            if j%10==0:
            print('\t superlineage combo ' + str(j+1) + ' out of ' + str(numtrials))
            df_out_j = pd.DataFrame()
            # Creating superlineages and formatting the data
            counts, frequency = fd.create_superlineages_counts_and_freq(counts_original, total_sampled_counts, d1, dl, mincount, minfreq, seed = np.random.randint(10**5))
            if len(counts)>1: # Check that there is at least one superlineage
                # Inference
                data_list = fd.create_data_list(counts, total_sampled_counts, epiweeks)        
                k_df = tnf.get_kappa(counts, epiweeks, total_sampled_counts, mincount = mincount)

                if k_df['stderr'][0]>0:
                    # first estimate of Ne and c using MSD method
#                    c_all, summary_infer, k_df, df_all_trials = tnf.fit_time_varying_technical_noise_bounded(k_df, epiweeks, numtrials = 0)
#                    Netau_MSD = summary_infer['Netau'].values[0]
                    shareddata = lh.getshareddata(total_sampled_counts)     
                    
                    # use as initial guess for joint inference of Ne and c using MLE
#                    x0 = np.concatenate((np.array([1/Netau_MSD]), c_all.values[0]))
                    hmm_maxlikelihood = lambda x: lh.hmm_maxlikelihood_Ne(x[0], x[1])
                    Netau_HMM = hmm_maxlikelihood([data_list, shareddata])
                    
                    # confidence interval estimation using profile likelihood
                    Netau_HMM_lower, Netau_HMM_upper = lh.hmm_plot_error_find_root(data_list, shareddata)
                
                else:
                    Netau_HMM = np.nan
                    Netau_HMM_lower = np.nan
                    Netau_HMM_upper = np.nan
            else:
                Netau_HMM = np.nan
                Netau_HMM_lower = np.nan
                Netau_HMM_upper = np.nan
                
            df_out_j['epiweek_start'] = [d1]
            df_out_j['epiweek_middle'] = [float(d1)+(T-1)/2]
            df_out_j['epiweek_end'] = [dl]
            df_out_j['superlineage_combo'] = [j]
            df_out_j['Netau_HMM'] = [Netau_HMM]
            df_out_j['Netau_HMM_lower'] = [Netau_HMM_lower]
            df_out_j['Netau_HMM_upper'] = [Netau_HMM_upper]
            df_out = df_out.append(df_out_j)
        df_out.reset_index(inplace = True, drop = True)
        output_filename = 'raw_infected_plus_exposed.csv'
        df_out.to_csv(output_folder + output_filename)

path_folder = '../../simulation_data/stochastic_seir/N1000000_R02_gammaE0.33_gammaI0.18_numlineages100_labeltime75/'
output_folder = path_folder + '/inferred_Ne_superlineage_combos_threshold_counts_freq/'    
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
run_inference(path_folder, output_folder)

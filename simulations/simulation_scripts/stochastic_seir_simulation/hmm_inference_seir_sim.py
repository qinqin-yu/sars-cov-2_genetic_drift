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

bin_path = '../../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../../functions')
    
import pandas as pd
import numpy as np 
import os
import glob
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import Pool, get_context
import format_data as fd
import analytical_likelihoods_no_emission as lh
import technical_noise_functions as tnf
import epidemiological_models_functions as emf

def run_inference(path_folder, output_folder, T=9, mincount = 20, minfreq = 0.01, numtrials = 20):
#    trial_num_idx = filename_trial.find('_trial_')
#    trial_num_idx2 = filename_trial.find('.csv')
#    trial_num = int(filename_trial[trial_num_idx+7:trial_num_idx2])
                
    # Importing data
    counts_all=pd.read_csv(path_folder + 'counts_infected.csv', index_col = 0)
#    counts_all = counts_all[0:10]
    total_sampled_counts_all = pd.read_csv(path_folder + 'total_counts_infected.csv', index_col = 0)
    
    epiweeks_all = counts_all.drop(['lineage'], axis = 1).columns
#    epiweeks_all = epiweeks_all[60:60+T]
    df_out = pd.DataFrame()

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
                    Netau_HMM_lower, Netau_HMM_upper = lh.confidence_interval_Ne(data_list, shareddata)
                
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
            df_out_j['T'] = T
            df_out_j['delta_t']=1
            df_out = df_out.append(df_out_j)
        df_out.reset_index(inplace = True, drop = True)
        output_filename = 'raw.csv'
        df_out.to_csv(output_folder + output_filename)

def summarize_results(path_folder, raw_filename, output_folder, output_filename):
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
#    summary.dropna(subset = ['Netau_HMM_median'], how = 'all', inplace = True)

    epiweeks_str = summary['Epiweek'].values.astype('str')

    # Calculate epideiological model info
    
    infected = pd.read_csv(path_folder  + 'total_counts_infected.csv', index_col = 0)
    exposed = pd.read_csv(path_folder  + 'total_counts_exposed.csv', index_col = 0)
    susceptible = pd.read_csv(path_folder  + 'total_counts_susceptible.csv', index_col = 0)    
    
    idx = path_folder.find('N')
    idx2 = path_folder.find('R0')
    idx3 = path_folder.find('gammaE')
    idx4 = path_folder.find('gammaI')
    idx5 = path_folder.find('numlineages')
    idx6 = path_folder.find('labeltime')
    
    N = int(path_folder[idx+1:idx2-1])
    R0 = float(path_folder[idx2+2:idx3-1])
    gamma_E = float(path_folder[idx3+6:idx4-1])*7
    gamma_I = float(path_folder[idx4+6:idx5-1])*7
    numlineages = int(path_folder[idx5+11:idx6-1])
    labeltime = float(path_folder[idx6+9:-1])

    # Analytical calculation of epidemiological model Ne
    Netau_SEIR = (exposed[epiweeks_str]+infected[epiweeks_str])**2/(2*(R0*susceptible[epiweeks_str]/N)*gamma_I*infected[epiweeks_str])
    Netau_SEIR = Netau_SEIR.values[0]
    
    Netau_SIR = (infected[epiweeks_str])/(2*(R0*susceptible[epiweeks_str]/N)*gamma_I)
    Netau_SIR = Netau_SIR.values[0]

    # Numerical calculation of epidemiological model Ne
    I = infected[epiweeks_str].values[0]
    Rt = R0*susceptible[epiweeks_str].values[0]/N

    merged = pd.DataFrame({'Epiweek':epiweeks_str, 'I':I, 'I_lower':I, 'I_upper':I, 'Rt':Rt, 'Rt_lower':Rt, 'Rt_upper':Rt})
    merged = emf.get_sir_model_Netau(merged, gamma_I)
    merged = emf.get_seir_model_Netau(merged, gamma_E, gamma_I)

    summary['Netau_SEIR_numerical'] = merged['Netau SEIR'].values
    summary['Netau_SIR_numerical'] = merged['Netau SIR'].values
    
    summary['Netau_SEIR'] = Netau_SEIR
    summary['Netau_SIR'] = Netau_SIR

    summary.to_csv(output_folder + output_filename)
    
path_folder = '../../simulation_data/stochastic_seir/N1000000_R02_gammaE0.33_gammaI0.18_numlineages100_labeltime75/'
output_folder = path_folder + '/inference_results/'    
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
run_inference(path_folder, output_folder)
summarize_results(path_folder, output_folder+'raw.csv', output_folder, 'summary.csv')

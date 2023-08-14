#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:43:01 2022

@author: qinqinyu
"""

import sys

bin_path = '../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../functions')
    
import pandas as pd
import numpy as np

import format_data as fd
import analytical_likelihoods as lh
import technical_noise_functions as tnf

def infer_Ne_c(path_folder, \
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
    c_index = np.arange(0, T*delta_t, delta_t) #np.array(range(T))
    c_index_str = ['c' + str(i) for i in c_index]
    c_MSD_index_str = ['c' + str(i) + '_MSD' for i in c_index]
    
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
                    
                    # First, estimate Ne and c using MSD method
                    c_all, summary_infer, k_df, df_all_trials = tnf.fit_time_varying_technical_noise_bounded(k_df, epiweeks, numtrials = 0)
                    Netau_MSD = summary_infer['Netau'].values[0]
                    shareddata = lh.getshareddata(c_all,total_sampled_counts)
                    
                    # Then feed this initial guess into HMM method
                    x0 = np.concatenate((np.array([1/Netau_MSD]), c_all.values[0]))
                    hmm_maxlikelihood = lambda x: lh.hmm_maxlikelihood_Ne_c(x[0], x[1], x[2])
                    output = hmm_maxlikelihood([data_list, shareddata, x0])
                    
                    # Confidence interval estimation of Ne using profile likelihood
                    Netau_HMM = output[0]
                    c_all_HMM = output[1]
                    Netau_HMM_lower, Netau_HMM_upper = lh.confidence_interval_Ne(Netau_HMM, data_list, shareddata, c_all_HMM)
                
                else:
                    c_all = pd.DataFrame([[np.nan]*len(epiweeks)], columns = epiweeks)
                    c_all_HMM = np.array([np.nan]*len(epiweeks))
                    Netau_MSD = np.nan
                    Netau_HMM = np.nan
                    Netau_HMM_lower = np.nan
                    Netau_HMM_upper = np.nan
            
            else:
                c_all = pd.DataFrame([[np.nan]*len(epiweeks)], columns = epiweeks)
                c_all_HMM = np.array([np.nan]*len(epiweeks))
                Netau_MSD = np.nan
                Netau_HMM = np.nan
                Netau_HMM_lower = np.nan
                Netau_HMM_upper = np.nan
            
            df_out_j['epiweek_start'] = [d1]
            df_out_j['epiweek_middle'] = [float(d1)+(float(dl)-float(d1))/2]
            df_out_j['epiweek_end'] = [dl]
            df_out_j['superlineage_combo'] = [j]
            df_out_j[c_MSD_index_str] = [c_all.values[0]]
            df_out_j['Netau_MSD'] = [Netau_MSD]
            df_out_j[c_index_str] = [c_all_HMM]
            df_out_j['Netau_HMM'] = [Netau_HMM]
            df_out_j['Netau_HMM_lower'] = [Netau_HMM_lower]
            df_out_j['Netau_HMM_upper'] = [Netau_HMM_upper]
            df_out_j['T'] = (T-1)*delta_t+1
            df_out_j['delta_t'] = delta_t
            df_out = df_out.append(df_out_j)
    df_out.reset_index(inplace = True, drop = True)
    df_out.to_csv(output_folder + output_filename)
    
def infer_Ne_c_remove_c_bound(path_folder, \
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
    c_index = np.arange(0, T*delta_t, delta_t) #np.array(range(T))
    c_index_str = ['c' + str(i) for i in c_index]
    c_MSD_index_str = ['c' + str(i) + '_MSD' for i in c_index]
    
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
                    
                    # First, estimate Ne and c using MSD method
                    c_all, summary_infer, k_df, df_all_trials = tnf.fit_time_varying_technical_noise_unbounded(k_df, epiweeks, numtrials = 0)
                    Netau_MSD = summary_infer['Netau'].values[0]
                    shareddata = lh.getshareddata(c_all,total_sampled_counts)
                    
                    # Then feed this initial guess into HMM method
                    x0 = np.concatenate((np.array([1/Netau_MSD]), c_all.values[0]))
                    hmm_maxlikelihood = lambda x: lh.hmm_maxlikelihood_Ne_c(x[0], x[1], x[2], bounds_c = (0, 100))
                    output = hmm_maxlikelihood([data_list, shareddata, x0])
                    
                    # Confidence interval estimation of Ne using profile likelihood
                    Netau_HMM = output[0]
                    c_all_HMM = output[1]
                    Netau_HMM_lower, Netau_HMM_upper = lh.confidence_interval_Ne(Netau_HMM, data_list, shareddata, c_all_HMM)
                
                else:
                    c_all = pd.DataFrame([[np.nan]*len(epiweeks)], columns = epiweeks)
                    c_all_HMM = np.array([np.nan]*len(epiweeks))
                    Netau_MSD = np.nan
                    Netau_HMM = np.nan
                    Netau_HMM_lower = np.nan
                    Netau_HMM_upper = np.nan
            
            else:
                c_all = pd.DataFrame([[np.nan]*len(epiweeks)], columns = epiweeks)
                c_all_HMM = np.array([np.nan]*len(epiweeks))
                Netau_MSD = np.nan
                Netau_HMM = np.nan
                Netau_HMM_lower = np.nan
                Netau_HMM_upper = np.nan
            
            df_out_j['epiweek_start'] = [d1]
            df_out_j['epiweek_middle'] = [float(d1)+(float(dl)-float(d1))/2]
            df_out_j['epiweek_end'] = [dl]
            df_out_j['superlineage_combo'] = [j]
            df_out_j[c_MSD_index_str] = [c_all.values[0]]
            df_out_j['Netau_MSD'] = [Netau_MSD]
            df_out_j[c_index_str] = [c_all_HMM]
            df_out_j['Netau_HMM'] = [Netau_HMM]
            df_out_j['Netau_HMM_lower'] = [Netau_HMM_lower]
            df_out_j['Netau_HMM_upper'] = [Netau_HMM_upper]
            df_out_j['T'] = (T-1)*delta_t+1
            df_out_j['delta_t'] = delta_t
            df_out = df_out.append(df_out_j)
    df_out.reset_index(inplace = True, drop = True)
    df_out.to_csv(output_folder + output_filename)
    
def summarize_results(raw_filename, output_folder, output_filename):
    df_all = pd.read_csv(raw_filename, index_col = 0)
    epiweeks_recorded = np.unique(df_all[['epiweek_start', 'epiweek_middle', 'epiweek_end']])
    epiweeks = np.arange(np.min(epiweeks_recorded), np.max(epiweeks_recorded)+1)
    T = df_all['T'].values[0]
    if 'delta_t' in df_all.columns:
        delta_t = df_all['delta_t'].values[0]
    else:
        delta_t = 1
    T = int((T-1)/delta_t+1)
    
    c_index = c_index = np.arange(0, T*delta_t, delta_t) #np.array(range(T))

    c_dict = {}
    c_columns = np.array(['c' + str(i) for i in c_index])

    c_dict_MSD = {}
    c_columns_MSD = np.array(['c' + str(i) + '_MSD' for i in c_index])
 
    summary = pd.DataFrame()
    for epiweek, df in df_all.groupby('epiweek_middle'):
        e1 = df['epiweek_start'].values[0]
        el = df['epiweek_end'].values[0]

        cond = df['Netau_HMM']<10**5
        Ne_median = np.nanmedian(df[cond]['Netau_HMM'])
        Ne_2_5_percentile = np.nanmedian(df[cond]['Netau_HMM_lower'])
        Ne_97_5_percentile = np.nanmedian(df[cond]['Netau_HMM_upper'])

        cond_MSD = df['Netau_MSD']<10**5
        Ne_median_MSD = np.nanmedian(df[cond_MSD]['Netau_MSD'])
        Ne_2_5_percentile_MSD = np.nanpercentile(df[cond_MSD]['Netau_MSD'], 2.5)
        Ne_97_5_percentile_MSD = np.nanpercentile(df[cond_MSD]['Netau_MSD'], 97.5)

        for c_column in c_columns:
            del_idx = int(c_column[-1])
            epiweek_start_idx = np.where(epiweeks==e1)[0][0]
            epiweek_idx = epiweek_start_idx+del_idx
            c_epiweek = epiweeks[epiweek_idx]
            if c_epiweek not in c_dict: 
                c_dict[c_epiweek] = df[c_column].values
            else:
                c_dict[c_epiweek] = np.append(c_dict[c_epiweek], df[c_column].values)

        for c_column_MSD in c_columns_MSD:
            del_idx = int(c_column_MSD[-5])
            epiweek_start_idx = np.where(epiweeks==e1)[0][0]
            epiweek_idx = epiweek_start_idx+del_idx
            c_epiweek = epiweeks[epiweek_idx]
            if c_epiweek not in c_dict_MSD: 
                c_dict_MSD[c_epiweek] = df[c_column_MSD].values
            else:
                c_dict_MSD[c_epiweek] = np.append(c_dict_MSD[c_epiweek], df[c_column_MSD].values)

        summary = summary.append(pd.DataFrame({'Epiweek':[epiweek], \
                                               'Netau_HMM_median':[Ne_median], \
                                               'Netau_HMM_95%_ci_lower':[Ne_2_5_percentile], \
                                               'Netau_HMM_95%_ci_upper':[Ne_97_5_percentile], \
                                               'Netau_MSD_median':[Ne_median_MSD], \
                                               'Netau_MSD_95%_ci_lower':[Ne_2_5_percentile_MSD], \
                                               'Netau_MSD_95%_ci_upper':[Ne_97_5_percentile_MSD]}))

    for c_epiweek in c_dict:
        c_median = np.nanmedian(np.array(c_dict[c_epiweek]))
        c_2_5_percentile = np.nanpercentile(np.array(c_dict[c_epiweek]), 2.5)
        c_97_5_percentile = np.nanpercentile(np.array(c_dict[c_epiweek]), 97.5)
        if np.isin(c_epiweek, summary['Epiweek'].values):
            summary.loc[summary['Epiweek'] == c_epiweek, 'c_median'] = c_median
        else:
            summary = summary.append(pd.DataFrame({'Epiweek':[c_epiweek], 'c_median':c_median}))
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_95%_ci_lower'] = c_2_5_percentile
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_95%_ci_upper'] = c_97_5_percentile
    for c_epiweek in c_dict_MSD:
        c_median = np.nanmedian(np.array(c_dict_MSD[c_epiweek]))
        c_2_5_percentile = np.nanpercentile(np.array(c_dict_MSD[c_epiweek]), 2.5)
        c_97_5_percentile = np.nanpercentile(np.array(c_dict_MSD[c_epiweek]), 97.5)
        if np.isin(c_epiweek, summary['Epiweek'].values):
            summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_median'] = c_median
        else:
            summary = summary.append(pd.DataFrame({'Epiweek':[c_epiweek], 'c_MSD_median':c_median}))
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_95%_ci_lower'] = c_2_5_percentile
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_95%_ci_upper'] = c_97_5_percentile
    summary['T'] = (T-1)*delta_t+1
    summary['delta_t'] = delta_t
    summary.sort_values('Epiweek', inplace = True)
    summary.reset_index(inplace = True, drop = True)
    summary.dropna(subset = ['Netau_HMM_median', 'c_median', 'Netau_MSD_median', 'c_MSD_median'], how = 'all', inplace = True)
    summary.to_csv(output_folder + output_filename)

def summarize_results_keep_high_Ne_runs(raw_filename, output_folder, output_filename):
    df_all = pd.read_csv(raw_filename, index_col = 0)
    epiweeks_recorded = np.unique(df_all[['epiweek_start', 'epiweek_middle', 'epiweek_end']])
    epiweeks = np.arange(np.min(epiweeks_recorded), np.max(epiweeks_recorded)+1)
    T = df_all['T'].values[0]
    if 'delta_t' in df_all.columns:
        delta_t = df_all['delta_t'].values[0]
    else:
        delta_t = 1
    T = int((T-1)/delta_t+1)
    
    c_index = c_index = np.arange(0, T*delta_t, delta_t) #np.array(range(T))

    c_dict = {}
    c_columns = np.array(['c' + str(i) for i in c_index])

    c_dict_MSD = {}
    c_columns_MSD = np.array(['c' + str(i) + '_MSD' for i in c_index])
 
    summary = pd.DataFrame()
    for epiweek, df in df_all.groupby('epiweek_middle'):
        e1 = df['epiweek_start'].values[0]
        el = df['epiweek_end'].values[0]

        cond = df['Netau_HMM']<np.inf#10**5
        Ne_median = np.nanmedian(df[cond]['Netau_HMM'])
        Ne_2_5_percentile = np.nanmedian(df[cond]['Netau_HMM_lower'])
        Ne_97_5_percentile = np.nanmedian(df[cond]['Netau_HMM_upper'])

        cond_MSD = df['Netau_MSD']<np.inf#10**5
        Ne_median_MSD = np.nanmedian(df[cond_MSD]['Netau_MSD'])
        Ne_2_5_percentile_MSD = np.nanpercentile(df[cond_MSD]['Netau_MSD'], 2.5)
        Ne_97_5_percentile_MSD = np.nanpercentile(df[cond_MSD]['Netau_MSD'], 97.5)

        for c_column in c_columns:
            del_idx = int(c_column[-1])
            epiweek_start_idx = np.where(epiweeks==e1)[0][0]
            epiweek_idx = epiweek_start_idx+del_idx
            c_epiweek = epiweeks[epiweek_idx]
            if c_epiweek not in c_dict: 
                c_dict[c_epiweek] = df[c_column].values
            else:
                c_dict[c_epiweek] = np.append(c_dict[c_epiweek], df[c_column].values)

        for c_column_MSD in c_columns_MSD:
            del_idx = int(c_column_MSD[-5])
            epiweek_start_idx = np.where(epiweeks==e1)[0][0]
            epiweek_idx = epiweek_start_idx+del_idx
            c_epiweek = epiweeks[epiweek_idx]
            if c_epiweek not in c_dict_MSD: 
                c_dict_MSD[c_epiweek] = df[c_column_MSD].values
            else:
                c_dict_MSD[c_epiweek] = np.append(c_dict_MSD[c_epiweek], df[c_column_MSD].values)

        summary = summary.append(pd.DataFrame({'Epiweek':[epiweek], \
                                               'Netau_HMM_median':[Ne_median], \
                                               'Netau_HMM_95%_ci_lower':[Ne_2_5_percentile], \
                                               'Netau_HMM_95%_ci_upper':[Ne_97_5_percentile], \
                                               'Netau_MSD_median':[Ne_median_MSD], \
                                               'Netau_MSD_95%_ci_lower':[Ne_2_5_percentile_MSD], \
                                               'Netau_MSD_95%_ci_upper':[Ne_97_5_percentile_MSD]}))

    for c_epiweek in c_dict:
        c_median = np.nanmedian(np.array(c_dict[c_epiweek]))
        c_2_5_percentile = np.nanpercentile(np.array(c_dict[c_epiweek]), 2.5)
        c_97_5_percentile = np.nanpercentile(np.array(c_dict[c_epiweek]), 97.5)
        if np.isin(c_epiweek, summary['Epiweek'].values):
            summary.loc[summary['Epiweek'] == c_epiweek, 'c_median'] = c_median
        else:
            summary = summary.append(pd.DataFrame({'Epiweek':[c_epiweek], 'c_median':c_median}))
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_95%_ci_lower'] = c_2_5_percentile
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_95%_ci_upper'] = c_97_5_percentile
    for c_epiweek in c_dict_MSD:
        c_median = np.nanmedian(np.array(c_dict_MSD[c_epiweek]))
        c_2_5_percentile = np.nanpercentile(np.array(c_dict_MSD[c_epiweek]), 2.5)
        c_97_5_percentile = np.nanpercentile(np.array(c_dict_MSD[c_epiweek]), 97.5)
        if np.isin(c_epiweek, summary['Epiweek'].values):
            summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_median'] = c_median
        else:
            summary = summary.append(pd.DataFrame({'Epiweek':[c_epiweek], 'c_MSD_median':c_median}))
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_95%_ci_lower'] = c_2_5_percentile
        summary.loc[summary['Epiweek'] == c_epiweek, 'c_MSD_95%_ci_upper'] = c_97_5_percentile
    summary['T'] = (T-1)*delta_t+1
    summary['delta_t'] = delta_t
    summary.sort_values('Epiweek', inplace = True)
    summary.reset_index(inplace = True, drop = True)
    summary.dropna(subset = ['Netau_HMM_median', 'c_median', 'Netau_MSD_median', 'c_MSD_median'], how = 'all', inplace = True)
    summary.to_csv(output_folder + output_filename)
    
def correct_for_assigned_sequences(total_counts_lineages_full_filename, total_counts_variant_tree_full_filename, summary_full_filename, summary_corrected_full_filename):
    counts_assigned = pd.read_csv(total_counts_lineages_full_filename, index_col = 0)
    counts_tree = pd.read_csv(total_counts_variant_tree_full_filename, index_col = 0)
    frac_assigned = counts_assigned/counts_tree
    frac_assigned_df = pd.DataFrame(data = {'Epiweek':frac_assigned.columns.astype('float'), 'frac_assigned':frac_assigned.values[0]})
    summary = pd.read_csv(summary_full_filename, index_col = 0)
    summary = summary.merge(frac_assigned_df, on = 'Epiweek', how = 'left')

    summary = summary.rename(columns = {'Netau_HMM_median':'Netau_HMM_assigned_median',
                                        'Netau_HMM_95%_ci_lower':'Netau_HMM_assigned_95%_ci_lower',
                                        'Netau_HMM_95%_ci_upper':'Netau_HMM_assigned_95%_ci_upper',
                                        'Netau_MSD_median':'Netau_MSD_assigned_median',
                                        'Netau_MSD_95%_ci_lower':'Netau_MSD_assigned_95%_ci_lower',
                                        'Netau_MSD_95%_ci_upper':'Netau_MSD_assigned_95%_ci_upper'})

    summary['Netau_HMM_median'] = summary['Netau_HMM_assigned_median']/summary['frac_assigned']
    summary['Netau_HMM_95%_ci_lower'] = summary['Netau_HMM_assigned_95%_ci_lower']/summary['frac_assigned']
    summary['Netau_HMM_95%_ci_upper'] = summary['Netau_HMM_assigned_95%_ci_upper']/summary['frac_assigned']
    summary['Netau_MSD_median'] = summary['Netau_MSD_assigned_median']/summary['frac_assigned']
    summary['Netau_MSD_95%_ci_lower'] = summary['Netau_MSD_assigned_95%_ci_lower']/summary['frac_assigned']
    summary['Netau_MSD_95%_ci_upper'] = summary['Netau_MSD_assigned_95%_ci_upper']/summary['frac_assigned']

    summary.to_csv(summary_corrected_full_filename)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:56:16 2021

Up to date cleaned-up functions for inferring effective population size and 
technical noise using a Hidden Markov Model and a mean squared displacement 
model.

@author: qinqinyu
"""
import pandas as pd
import numpy as np

# FUNCTIONS FOR FORMATTING DATA

def counts_to_frequency(counts, total_counts):
    # Function for getting frequencies from counts and total_counts
    epiweeks = counts.columns
    neutral_frequency = (counts[epiweeks].values/total_counts[epiweeks].values)
    neutral_frequency = pd.DataFrame(data = neutral_frequency, columns = epiweeks)
    return neutral_frequency

def create_superlineages_counts_and_freq(counts, total_counts, d1, dl, mincount, minfreq, seed = 0):
    '''
    Removes trajectories that start at 0 (likely to be newly emerging) and 
    sums lineages until their starting and ending count AND frequency is above a threshold 
    value.
    
    INPUTS:
        
        counts: Dataframe with lineage counts over time. Each row is a 
        lineage and each column is a timepoint. One additional column gives the 
        lineage name.
        
        d1: The column names for the first timepoint.
        
        dl: The column name for the last timepoint.
        
        mincount: The threshold number of counts required at the beginning and 
        end of a superlineage.
        
        minfreq: The threshold frequency required at the beginning and 
        end of a superlineage.
    
    OUTPUTS:
        
        counts: Data frame with superlineage counts over time. Same format as
        input counts variable. Superlineages are named with asecending integer numbers.
    
    '''
#    # Removing trajectories that start at 0:
#    counts = counts.loc[(counts[d1]!=0)]
#    counts.reset_index(drop = True, inplace = True)
    
    # Shuffling rows (so that it's not ordered by abundance)
    counts = counts.sample(frac = 1, random_state = seed)
    counts.reset_index(drop = True, inplace = True)
    
    # Get frequencies
    frequency = counts_to_frequency(counts, total_counts)

    # Merge low count lineages into a super lineage
    # Note that poisson technical noise adds (sum of two poissons is a poisson)
    droprows=[]
    droprows_running = []
    
    # If want to plot original data
    #plt.figure()
    #plt.plot((np.transpose((counts[epiweeks].values/total_counts.values)/frac_neutral.values)))
    #plt.yscale('log')
    
    new_rows = pd.DataFrame()
    new_rows_freq = pd.DataFrame()
    for i,row in frequency.iterrows():
        if row[d1]<minfreq or row[dl]<minfreq or counts.loc[i,d1]<mincount or counts.loc[i,dl]<mincount:
            droprows.append(i)
            droprows_running.append(i)
            summed_rows = np.sum(counts.loc[droprows_running])
            summed_rows_freq = np.sum(frequency.loc[droprows_running])
            if summed_rows_freq[d1]>=minfreq and summed_rows_freq[dl]>=minfreq and summed_rows[d1]>=mincount and summed_rows[dl]>=mincount:
                new_rows = new_rows.append(summed_rows, ignore_index = True)
                new_rows_freq = new_rows_freq.append(summed_rows_freq, ignore_index = True)
                droprows_running = []
    counts.drop(index=droprows,inplace=True)
    counts.reset_index(drop = True, inplace = True)
    counts = counts.append(new_rows, ignore_index = True, sort = True)
    
    frequency.drop(index=droprows,inplace=True)
    frequency.reset_index(drop = True, inplace = True)
    frequency = frequency.append(new_rows_freq, ignore_index = True, sort = True)
    
    counts['lineage'] = np.arange(0, len(counts), 1)
    frequency['lineage'] = np.arange(0, len(frequency), 1)

    return counts, frequency

def create_data_list(counts, total_sampled_counts, epiweeks):
    '''
    Formats data for input into HMM inference method.
    
    INPUTS:
        counts: Dataframe with lineage counts over time. Each row is a 
        lineage and each column is a timepoint. One additional column gives the 
        lineage name.
        
        total_sampled_counts: Dataframe with the total number of sampled counts across all 
        lineages as a function of time. Columns are each timepoint.

        epiweeks: Array of timepoints (same as timepoint column names in other
        dataframes).
        
    OUTPUTS:
        data_list: a list of dataframes, one for each lineage, containing the
        Day, number of counts, total number of counts across all lineages, and 
        fraction of counts from neutral lineages.
        
    '''
    data_list = []
    for g,gene in counts.groupby(by=['lineage']):	
        data_df = pd.DataFrame()
        for week in epiweeks:
            data={
                    'Day':np.array(week),
                    'nobs':np.array(gene[week]),
                    'N_samp':np.array(total_sampled_counts[week])}
            d00 = pd.DataFrame.from_dict(data)
            
            data_df = data_df.append(d00)
            data_df.reset_index(drop = True, inplace = True)    
        data_list.append(data_df)
    return data_list
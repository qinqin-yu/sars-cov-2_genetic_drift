#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:43:01 2022

@author: qinqinyu
"""

import sys

bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
    
import os
import glob
import pandas as pd
import numpy as np
from multiprocessing import Pool

import hmm_inference

if __name__ == '__main__':

    # SPECIFY GLOBAL PARAMETERS
    # Parameters for processing - change to fit needs
    parallel = True # Whether to use parallel processing (best for a large number of data files, on cluster)
    processes = 12 # If using parallel processing, the number of parallel processes to run

    # Filenames - change to fit needs
    counts_filename = 'counts_lineages.csv'
    total_counts_filename = 'total_counts_lineages.csv'
    total_counts_variant_tree_filename = 'total_counts_variant_tree.csv'

    fitness_cutoffs = [0.09199183983679672, 0.1645132902658053, 0.26873537470749415]
    percentiles = [50, 75, 90]
    #50th, 75th, and 90th percentile
    
    # GET EXCLUDED BENEFICIAL LINEAGES
    fitness = pd.read_csv('../../data/fitness/concat_fitness.csv', index_col = 0)
    
    # GET LIST OF PARAMETERS TO RUN ANALYSES ON
    # All combinations of variants/tree date/tree cut depth for England
    params_all = []
    
    variant_param_folders = ['../../data/lineages/pre_B-1-177/microreact/',
                         '../../data/lineages/B-1-177/B-1-177|2021-02-22|694.5/',
                         '../../data/lineages/alpha/alpha|2021-06-20|61.5/',
                         '../../data/lineages/delta/delta|2022-01-25|49.5+58.5/']

    for variant_param_folder in variant_param_folders:
        variant_param = os.path.basename(variant_param_folder[:-1])
        path_folder = variant_param_folder + '/is_pillar_2/England/'
        output_folder = path_folder + 'inference_results/'

        i = 0
        for fitness_cutoff in fitness_cutoffs:
            fitness_sig = fitness[np.abs(fitness['s']) > fitness_cutoff]
            excluded_lineages = fitness_sig[fitness_sig['subset'] == variant_param]['lineage'].values
            
            print(variant_param)
            print(fitness_cutoff)
            print(excluded_lineages)
            
            output_filename = 'raw_exclude_ben_lineages_deterministic_' + str(percentiles[i]) + '.csv'

            params = (path_folder, counts_filename, total_counts_filename, output_folder, output_filename, excluded_lineages)
            params_all.append(params)
            i+=1
#    params_all = params_all[0:5] # comment out for testing
    
    # MAKE OUTPUT FOLDERS
    for params in params_all:
        output_folder = params[3]
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
      
    # RUN INFERENCE PIPELINE
    if parallel:
        with Pool(processes = processes) as pool:
            result = pool.starmap(hmm_inference.infer_Ne_c, params_all)
    else:
        for params in params_all:
            path_folder = params[0]
            print(path_folder)
            counts_filename = params[1]
            total_counts_filename = params[2]
            output_folder = params[3]
            output_filename = params[4]
            excluded_lineages = params[5]

            hmm_inference.infer_Ne_c(path_folder, \
                  counts_filename, \
                  total_counts_filename, \
                  output_folder, \
                  output_filename, \
                  excluded_lineages)
            
    # SUMMARIZE RESULTS
    for params in params_all:
        path_folder = params[0]
        total_counts_filename = params[2]
        output_folder = params[3]
        raw_filename = params[4]
        
        percentile = raw_filename.split('_')[-1][:-4]
        
        raw_filename = output_folder + raw_filename
        output_filename = 'summary_exclude_ben_lineages_deterministic_' + percentile + '.csv'
        hmm_inference.summarize_results(raw_filename, output_folder, output_filename)

        total_counts_lineages_full_filename = path_folder + total_counts_filename
        total_counts_variant_tree_full_filename = path_folder + total_counts_variant_tree_filename
        summary_full_filename = output_folder + output_filename
        summary_corrected_full_filename = output_folder + 'summary_corrected_exclude_ben_lineages_deterministic_' + percentile + '.csv'

        if os.path.exists(total_counts_variant_tree_full_filename):
            hmm_inference.correct_for_assigned_sequences(total_counts_lineages_full_filename, total_counts_variant_tree_full_filename, summary_full_filename, summary_corrected_full_filename)

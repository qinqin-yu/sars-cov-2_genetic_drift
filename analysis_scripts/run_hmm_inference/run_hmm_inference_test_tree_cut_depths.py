#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: qinqinyu
"""

import sys

bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
    
import os
import glob
import numpy as np
from multiprocessing import Pool

import hmm_inference

if __name__ == '__main__':

    # SPECIFY GLOBAL PARAMETERS
    # Parameters for processing - change to fit needs
    parallel = False # Whether to use parallel processing (best for a large number of data files, on cluster)
    processes = 16 # If using parallel processing, the number of parallel processes to run

    # Filenames - change to fit needs
    output_filename = 'raw.csv'
    counts_filename = 'counts_lineages.csv'
    total_counts_filename = 'total_counts_lineages.csv'
    # total_counts_variant_tree_filename = 'total_counts_variant_tree.csv'

    # GET LIST OF PARAMETERS TO RUN ANALYSES ON
    params_all = []
    path_folders = ['../../data/lineages/B-1-177/B-1-177_sublineage692.5/is_pillar_2/England/',
                    '../../data/lineages/alpha/alpha_sublineage59.5/is_pillar_2/England/',
                    '../../data/lineages/delta/delta_sublineage50.5/is_pillar_2/England/',
                   ]
    for path_folder in path_folders:
        output_folder = path_folder + 'inference_results/'
        params = (path_folder, counts_filename, total_counts_filename, output_folder, output_filename)
        params_all.append(params)
 
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

            hmm_inference.infer_Ne_c(path_folder, \
                  counts_filename, \
                  total_counts_filename, \
                  output_folder, \
                  output_filename)
            
    # SUMMARIZE RESULTS
    for params in params_all:
        path_folder = params[0]
        total_counts_filename = params[2]
        output_folder = params[3]
        raw_filename = params[4]
        
        raw_filename = output_folder + raw_filename
        output_filename = 'summary.csv'
        hmm_inference.summarize_results(raw_filename, output_folder, output_filename)

        # total_counts_lineages_full_filename = path_folder + total_counts_filename
        # total_counts_variant_tree_full_filename = path_folder + total_counts_variant_tree_filename
        # summary_full_filename = output_folder + output_filename
        # summary_corrected_full_filename = output_folder + 'summary_corrected.csv'
        
        # if os.path.exists(total_counts_variant_tree_full_filename):
        #     hmm_inference.correct_for_assigned_sequences(total_counts_lineages_full_filename, total_counts_variant_tree_full_filename, summary_full_filename, summary_corrected_full_filename)

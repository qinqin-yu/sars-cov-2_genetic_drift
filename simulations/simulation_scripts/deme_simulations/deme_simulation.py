#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:11:01 2022

@author: qinqinyu
"""
import sys

bin_path = '../../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../../functions')
    
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
import os
import deme_simulation_functions as dsf

output_path = '../../simulation_data/deme_simulations/counts/'

# SEIR model
R0 = 10
gamma_E = 1/2.5
gamma_I = 1/6.5
T = 1000 # total number of timesteps for deterministic SEIR
delta_t = 0.1 # timestep in daysfor deterministic SEIR

# Deme model
num_demes = 5.6*10**5 # England population size ### SHOULD CHANGE TO REFLECT DEME SIZE
num_lineages = 200

Ns = [10, 50, 100, 200] #  100 # individual per deme
num_demes_filled_initially_all = [1000, 2000, 5000]#[100, 1000, 2000, 5000]

if not os.path.exists(output_path):
    os.makedirs(output_path)
    
for N in Ns:
    initial_conditions = [N-1, 1, 0, 0] # Number of S, E, I, R
    time, I_cdf, Id = dsf.deterministic_seir(initial_conditions, R0=R0, gamma_E = gamma_E, gamma_I = gamma_I, N = N, T = T, delta_t = delta_t)
    for num_demes_filled_initially in num_demes_filled_initially_all:
        df_lineages, df, t_transmits = dsf.deme_simulation(num_demes_filled_initially, num_demes, num_lineages, time, I_cdf)
        df_abundances, time_long = dsf.get_abundances(df_lineages, df, T, delta_t, Id)
        dsf.save_results(df_abundances, num_demes, num_lineages, num_demes_filled_initially, N, time_long, output_path)
        
# Resaving data to look at only epiweeks 42-50

epiweeks = np.arange(42, 51)
columns = np.append(np.array('lineage'), epiweeks)
path_folder = '../../simulation_data/deme_simulations/'
output_path = path_folder + 'counts_t0_42_tl_50/'
if not os.path.exists(output_path):
    os.makedirs(output_path)
filenames = glob.glob(path_folder + 'counts/demes_' + str(num_demes) + '_lineages_200_filleddemes0_*_demeN_*')
for filename in filenames:
    filename_only = os.path.basename(filename)
    filename_prefix = filename_only[:filename_only.find('.csv')]
    df = pd.read_csv(filename, index_col = 0)
    df_42_50 = df[columns]
    df_42_50.to_csv(output_path + filename_prefix + '_counts_t0_42_tl_50.csv')
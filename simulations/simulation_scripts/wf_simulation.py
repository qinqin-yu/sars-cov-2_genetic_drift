#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:33:43 2022

@author: qinqinyu
"""
import sys
bin_path = '../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../functions')
import pandas as pd
import numpy as np
import os
import wf_simulation_functions as wsf

# PARAMETERS
# General
total_epiweeks = 100 # weeks
total_burnin_time = 10 # weeks
numlineages = 500
numtrials = 1

# Fitness
s_mean = 0
s_std = 0

# Mutation rate
mu = 0.01

# Number of sampled sequences
Nseqmin = 500
Nseqmax = 2*10**4

# Netau
Netau_mu = 50
Netau_sig = 10
Netau_min = 500
Netau_max = 10**4

# Sampling noise
m = 0.5  # mean minus 1
v = 2 # approximate variance

# PREP INPUTS FOR SIMULATIONS
t = np.arange(0, total_epiweeks)
Net_gaussian = wsf.gaussian(t, Netau_mu, Netau_sig, Netau_min, Netau_max)
Net_rectangular = wsf.rectangular(t, Netau_mu-Netau_sig, Netau_mu+Netau_sig, Netau_min, Netau_max)
Net_constant = np.array([Netau_min + (Netau_max-Netau_min)/2]*total_epiweeks)

Net_all = [Net_gaussian, Net_rectangular, Net_constant]
Net_all_labels = ['gaussian', 'rectangular', 'constant']

theta = v/m
k = m**2/v
seed = 2
rng = np.random.default_rng(seed)
ct = rng.gamma(k, scale = theta, size = total_epiweeks)+1

Nseq_mean = ((Nseqmax-Nseqmin)/total_epiweeks)*t + Nseqmin
Nseq_noise = np.random.normal(scale = Nseq_mean*0.2)
Nseq = Nseq_mean + Nseq_noise
Nseq[Nseq<=Nseqmin] = Nseqmin
Nseq = Nseq.astype('int')

# NEUTRAL SIMULATIONS

path_folder = '../simulation_data/'

# Different shapes of Netau

for i in range(len(Net_all)):
    Net = Net_all[i]
    Net_label = Net_all_labels[i]
    
    output_folder = path_folder + 'neutral/' + Net_label + '/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    ct_df = pd.DataFrame([ct], columns = t)
    Net_df = pd.DataFrame([Net], columns = t)
    
    ct_df.to_csv(output_folder + 'true_c.csv')
    Net_df.to_csv(output_folder + 'true_Netau.csv')

    counts_output_filename = 'counts_lineages.csv'
    total_counts_output_filename = 'total_counts_lineages.csv'
    fitness_filename = 'fitnesses.csv'
    wsf.run_simulation(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, s_mean, s_std, mu, numlineages)

# Multipel trials of gaussian-shaped Netau
Net = Net_gaussian
Net_label = 'gaussian'
numtrials = 20

output_folder = path_folder + 'neutral/' + Net_label + '/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
ct_df = pd.DataFrame([ct], columns = t)
Net_df = pd.DataFrame([Net], columns = t)

ct_df.to_csv(output_folder + 'true_c.csv')
Net_df.to_csv(output_folder + 'true_Netau.csv')

for trial in range(numtrials):
    counts_output_filename = 'counts_lineages_trial_' + str(trial) + '.csv'
    total_counts_output_filename = 'total_counts_lineages_trial_' + str(trial) + '.csv'
    fitness_filename = 'fitnesses_trial_' + str(trial) + '.csv'
    wsf.run_simulation(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, s_mean, s_std, mu, numlineages)

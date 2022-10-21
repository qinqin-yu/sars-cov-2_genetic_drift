#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:56:38 2022

@author: qinqinyu
"""
import sys

bin_path = '../../../functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../../../functions')

import os
import numpy as np
import stochastic_seir_simulation_functions as ssf

R0 = 2#10
gamma_E = 1/3
gamma_I = 1/5.5
N = 10**6 # total number of individuals
# timesteps = 1000
beta = R0*gamma_I

initial_conditions = [N-1,1,0,0]#[N-1, 1, 0, 0]
num_lineages = 100
label_time = 75

num_trials = 1

output_path = '../simulation_data/stochastic_seir/'
folder = 'N' + str(N) + '_R0' + str(R0)+ '_gammaE' + str(round(gamma_E,2)) + '_gammaI' + str(round(gamma_I,2)) + '_numlineages' + str(num_lineages) + '_labeltime' + str(label_time) + '/'

if not os.path.exists(output_path + folder):
    os.makedirs(output_path + folder)
    
num_established = 0
for i in range(num_trials):
    result, lineages_state_initial = ssf.stochastic_seir(initial_conditions, num_lineages, label_time, R0=R0, gamma_E = gamma_E, gamma_I = gamma_I, N = N)
    df_S_total, df_E, df_I, df_R = ssf.get_counts_from_sim_output(result, lineages_state_initial, num_lineages, label_time)
    df_E_total, df_I_total, df_R_total = ssf.save_results(output_path + folder, df_S_total, df_E, df_I, df_R)
    
#    print(df_S_total + df_E_total + df_I_total + df_R_total)
    
    if np.nanmin(result['S'])==0:
        num_established += 1
#    plt.plot(result['time'], result['S'], color = 'b', label = 'S', alpha = 0.3)
#    plt.plot(result['time'], result['E'], color = 'orange', label = 'E', alpha = 0.3)
#    plt.plot(result['time'], result['I'], color = 'g', label = 'I', alpha = 0.3)
#    plt.plot(result['time'], result['R'], color = 'r', label = 'R', alpha = 0.3)
#    plt.plot(result['time'], result['S'] + result['E'] + result['I'] + result['R'], color = 'k', label = 'N', alpha = 0.3)
#
#plt.legend()
#plt.show()
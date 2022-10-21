#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:45:04 2022

@author: qinqinyu
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
import os

def deterministic_seir(initial_conditions, R0=10, gamma_E = 1/2.5, gamma_I = 1/6.5, N = 100, T = 1000, delta_t = 0.1):
    
    '''
    Deterministic SEIR solution to feed into metapopulation model
    Time measured in days

    Default parameters:
    R0 = 10
    gamma_E = 1/2.5
    gamma_I = 1/6.5
    N = 100 # number of individuals per deme
    T = 1000 # total number of days simulated
    delta_t = 0.1 # timestep

    '''
    print('Running deterministic SEIR model')
    
    St0 = initial_conditions[0]#99
    Et0 = initial_conditions[1]#1
    It0 = initial_conditions[2]#0
    Rt0 = initial_conditions[3]#0
    
    beta = R0*gamma_I

    St = St0
    Et = Et0
    It = It0
    Rt = Rt0

    Sd = np.empty(T)
    Sd.fill(np.nan)

    Ed = np.empty(T)
    Ed.fill(np.nan)

    Id = np.empty(T)
    Id.fill(np.nan)

    Rd = np.empty(T)
    Rd.fill(np.nan)

    Sd[0] = St
    Ed[0] = Et
    Id[0] = It
    Rd[0] = Rt

    for t in range(1, T):
        dS = (-beta*It*St/N)*delta_t
        dE = (beta*It*St/N - gamma_E*Et)*delta_t
        dI = (gamma_E*Et - gamma_I*It)*delta_t
        dR = (gamma_I*It)*delta_t

        St = St+dS
        Et = Et+dE
        It = It+dI
        Rt = Rt+dR

        Sd[t] = St
        Ed[t] = Et
        Id[t] = It
        Rd[t] = Rt

    I_cdf = np.cumsum(Id)/np.max(np.cumsum(Id))
    time = np.array(range(T))*delta_t
    
    return time, I_cdf, Id

def deme_simulation(num_demes_filled_initially, num_demes, num_lineages, time, I_cdf):
    # Simulating deme model

    print('Running deme simulation for num_demes_filled_initial = ' + str(num_demes_filled_initially))
    
    num_demes = int(num_demes)
    df_lineages = pd.DataFrame({'lineage':np.array(range(num_lineages)), 'probability':np.ones(num_lineages)/num_lineages})

    df = pd.DataFrame(columns=['deme', 'filled', 'lineage', 't0'])
    df['deme'] = range(num_demes)
    df['filled'] = False
    active_demes = np.array(range(num_demes_filled_initially))

    # Fill initially populated demes each with a lineage
    for i in range(num_demes_filled_initially):
        u = np.random.uniform()
        lineage_prob = np.cumsum(df_lineages['probability'].values)
        lineage = df_lineages['lineage'][np.where(u<lineage_prob)[0][0]]
        df.loc[df['deme']==i, ['lineage', 'filled', 't0']] = [lineage, True, 0]

    # Simulate migration between lineages
    total_num_demes_affected = num_demes_filled_initially
    t_transmits = []
    while (len(active_demes)>0)and(total_num_demes_affected<num_demes-100):
        if total_num_demes_affected%500==0:
            print('\t Total number of demes affected = ' + str(total_num_demes_affected) + '/' + str(num_demes))
    #     print('active demes = ' + str(active_demes))
        host_deme = active_demes[0]
    #     print('host_deme = ' + str(host_deme))
        # Determine number of stochastic migration events that happen
        num_migrants = np.random.poisson(lam = 1)
    #     print('num_migrants = ' + str(num_migrants))
        t0_host = df[df['deme']==host_deme]['t0'].values[0]
        lineage = df[df['deme']==host_deme]['lineage'].values[0]
    #     print('t0_host = ' + str(t0_host))

        for k in range(num_migrants):
            # Determine when the stochastic migration events happen
            v = np.random.uniform()
            t_transmit = time[np.where(v<I_cdf)[0][0]]
            t_transmits.append(t_transmit)
    #         print('t_transmit = ' + str(t_transmit))
            t0_recipient = t_transmit + t0_host
            recipient_deme = df[~df['filled']]['deme'].values[0]
            df.loc[df['deme']==recipient_deme, ['lineage', 'filled', 't0']] = [lineage, True, t0_recipient]
    #         print('recipient_deme = ' + str(recipient_deme))
            active_demes = np.append(active_demes, recipient_deme)
            total_num_demes_affected = total_num_demes_affected + 1
        host_deme_idx = np.where(active_demes == host_deme)[0][0]
        active_demes = np.delete(active_demes, host_deme_idx)
    #     print(active_demes)
    #     print('total_num_demes_affected = ' + str(total_num_demes_affected))
    return df_lineages, df, t_transmits

def get_abundances(df_lineages, df, T, delta_t, Id):
    print('Getting abundaces from deme simulation')
    total_time_long = T*200*delta_t
    total_timesteps = int(total_time_long/delta_t)
    time_long = np.array(range(total_timesteps))*delta_t
    df_abundances = pd.DataFrame({'time':time_long})
    counter = 0
    for lineage in df_lineages['lineage']:
        if counter%20 == 0:
            print('\t Processing lineage number = ' + str(counter) + '/' + str(len(df_lineages)))
        abundance = np.zeros(len(time_long))
        demes = df[df['lineage'] == lineage]['deme'].values
        for deme in demes:
#             if deme%10000 == 0:
#                 print(deme)
#                 print(t0)
#                 print(np.where(t0<time_long))
            t0 = df[df['deme']==deme]['t0'].values[0]
            t0_idx = np.where(t0<time_long)[0][0]
            if (t0_idx+len(Id))<len(abundance):
                abundance[t0_idx:t0_idx+len(Id)] += Id
        df_abundances[lineage] = abundance
        counter+=1
#         plt.plot(time_long, abundance)
    # plt.xlim([0, 150])
#     plt.xlabel('Days')
#     plt.ylabel('Abundance')
#     plt.show()
    return df_abundances, time_long

def save_results(df_abundances, num_demes, num_lineages, num_demes_filled_initially, N, time_long, output_path = '/Users/qinqinyu/Documents/hallatschek_lab/scripts/sars-cov-2/simulations/simulated_data/seir_metapopulation_model/'):
    print('Saving results')
    
    # Remove times with no counts
    df_abundances = df_abundances.loc[~(df_abundances[np.array(range(num_lineages))]==0).all(axis=1)]

    # Convert to weeks
    total_time_long = np.max(df_abundances['time'])
    days_to_weeks = np.arange(0, total_time_long, 7)
    time_long_idxs = []
    for week in days_to_weeks:
        time_long_idxs.append(np.where(week<=time_long)[0][0])
    df_abundances_week = df_abundances.iloc[time_long_idxs]
    df_abundances_week['week'] = np.array(range(len(days_to_weeks)))

    # Save in .csv in same format as other data and plot
    df_abundances_week_transposed = pd.DataFrame(columns = df_abundances_week['week'].values, data = np.transpose(df_abundances_week[np.array(range(num_lineages))]).values)
    df_abundances_week_transposed['lineage'] = np.array(range(num_lineages))
    cols = df_abundances_week_transposed.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_abundances_week_transposed = df_abundances_week_transposed[cols]
    df_abundances_week_transposed.to_csv(output_path + '/demes_' + str(num_demes) + '_lineages_' + str(num_lineages) + '_filleddemes0_' + str(num_demes_filled_initially) + '_demeN_' + str(N) + '.csv')

    plt.plot(np.array(range(len(days_to_weeks))), np.transpose(df_abundances_week_transposed[np.array(range(len(days_to_weeks)))].values))
    plt.xlabel('Week')
    plt.ylabel('Lineage abundance')
    plt.show()


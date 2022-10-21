#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:55:25 2022

@author: qinqinyu
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os

def probabilities(state, params):
    beta = params[0]
    gamma_E = params[1]
    gamma_I = params[2]
    
    S = state[0]
    E = state[1]
    I = state[2]
    R = state[3]
    N = S+E+I+R
    probs = [beta*I*S/N, gamma_E*E, gamma_I*I]
#     'Sm1_Ep1', 'Em1_Ip1', 'Im1_Rp1'
    return probs

def process(probs):
    r = np.random.rand()*np.sum(probs)
    s = 0
    for i in range(len(probs)):
        s = s+probs[i]
        if r<s:
            which = i
            return which
        
def stochastic_seir(initial_conditions, num_lineages, label_time, R0=10, gamma_E = 1/2.5, gamma_I = 1/6.5, N = 100):
    state = initial_conditions.copy()
    timesteps = (N-state[3])  + (N-state[2]) + (N-state[1]) + 1
    states = np.empty((timesteps, 4))
    states.fill(np.nan)
    states[0,:] = np.array(state) # Number of S, E, I, R
    
    lineages_state = np.zeros((num_lineages, 4))
    lineages_state_initial = lineages_state.copy()

    beta = R0*gamma_I
    params = np.array([beta, gamma_E, gamma_I])
    time = np.empty(timesteps)
    time.fill(np.nan)
    time[0] = 0

    lineages = np.empty(timesteps)
    lineages.fill(np.nan)

    labeled = 0
    for i in range(timesteps-1):
        probs = probabilities(state, params)
        total_prob = np.sum(probs)
        if total_prob == 0:
            break
        tau = np.random.exponential(1/total_prob)
        if ((time[i] + tau)>label_time)&(not labeled):
            lineages_state[:,1] = np.random.multinomial(state[1], num_lineages*[1/num_lineages])
            lineages_state[:,2] = np.random.multinomial(state[2], num_lineages*[1/num_lineages])
            lineages_state_initial = lineages_state.copy()
            labeled = 1
        
        which = process(probs)
        if which==0:
            state[0] = state[0]-1
            state[1] = state[1]+1
            if np.sum(lineages_state[:,2])>0:
                lineage = np.random.choice(num_lineages, p=lineages_state[:,2]/np.sum(lineages_state[:,2]))
#                 lineages_state[lineage,0] = lineages_state[lineage,0]-1
                lineages_state[lineage,1] = lineages_state[lineage,1]+1
            else:
                lineage = np.nan
        elif which == 1:
            state[1] = state[1]-1
            state[2] = state[2]+1
            if np.sum(lineages_state[:,1])>0:
                lineage = np.random.choice(num_lineages, p=lineages_state[:,1]/np.sum(lineages_state[:,1]))
                lineages_state[lineage,1] = lineages_state[lineage,1]-1
                lineages_state[lineage,2] = lineages_state[lineage,2]+1
            else:
                lineage = np.nan
        elif which == 2:
            state[2] = state[2]-1
            state[3] = state[3]+1       
            if np.sum(lineages_state[:,2])>0:
                lineage = np.random.choice(num_lineages, p=lineages_state[:,2]/np.sum(lineages_state[:,2]))
                lineages_state[lineage,2] = lineages_state[lineage,2]-1
                lineages_state[lineage,3] = lineages_state[lineage,3]+1
            else:
                lineage = np.nan       
        else:
            print('Process ' + str(which) + ' is not an option')
        time[i+1] = time[i] + tau
        states[i+1,:] = state
        lineages[i+1] = lineage
    result = pd.DataFrame(states, columns = ['S', 'E', 'I', 'R'])
    result['lineage'] = lineages
    result['time'] = time
    return result, lineages_state_initial

def get_counts_from_sim_output(result, lineages_state_initial, num_lineages, label_time):
    result_diff = result.iloc[1:].copy()
    result_diff[['S', 'E', 'I', 'R']] = result[['S', 'E', 'I', 'R']].diff(periods=1, axis=0)

    weeks = np.arange(0, np.max(result['time'])+7, 7)
    I0 = lineages_state_initial[:,2]#result.loc[0, 'I']
    E0 = lineages_state_initial[:,1]#result.loc[0, 'E']
    S0 = lineages_state_initial[:,0]
    R0 = lineages_state_initial[:,3]

    I0_unlabeled = result.loc[0, 'I']
    E0_unlabeled = result.loc[0, 'E']
    S0_unlabeled = result.loc[0, 'S']
    R0_unlabeled = result.loc[0, 'R']

    delI = np.zeros((num_lineages, len(weeks)))
    delE = np.zeros((num_lineages, len(weeks)))
    delR = np.zeros((num_lineages, len(weeks)))

    delI_unlabeled = np.zeros(len(weeks))
    delE_unlabeled = np.zeros(len(weeks))
    delS_unlabeled = np.zeros(len(weeks))
    delR_unlabeled = np.zeros(len(weeks))
    for time, result_diff_i in result_diff.groupby('time'):
        idx = np.where(time>=weeks)[0][-1]+1
        lineage = result_diff_i['lineage'].values[0]
        if np.isnan(lineage):
            delI_unlabeled[idx]+=result_diff_i['I'].values[0]
            delE_unlabeled[idx]+=result_diff_i['E'].values[0]
            delS_unlabeled[idx]+=result_diff_i['S'].values[0]
            delR_unlabeled[idx]+=result_diff_i['R'].values[0]
        else:
            lineage = int(lineage)
            delI[lineage, idx]+=result_diff_i['I'].values[0]
            delE[lineage, idx]+=result_diff_i['E'].values[0]
            delR[lineage, idx]+=result_diff_i['R'].values[0]

    I = np.cumsum(delI, axis = 1)
    E = np.cumsum(delE, axis = 1)#+np.transpose(np.tile(E0, (len(weeks),1)))
    R = np.cumsum(delR, axis = 1)#+np.transpose(np.tile(E0, (len(weeks),1)))

    for i in range(int(np.floor(label_time/7)+1), len(weeks)):
        I[:,i]+=I0#np.transpose(np.tile(I0, (len(weeks),1)))
        E[:,i]+=E0#np.transpose(np.tile(I0, (len(weeks),1)))
        R[:,i]+=R0#np.transpose(np.tile(I0, (len(weeks),1)))

    I_unlabeled = np.cumsum(delI_unlabeled)
    E_unlabeled = np.cumsum(delE_unlabeled)
    S_unlabeled = np.cumsum(delS_unlabeled)
    R_unlabeled = np.cumsum(delR_unlabeled)

    I_unlabeled+=I0_unlabeled
    E_unlabeled+=E0_unlabeled
    S_unlabeled+=S0_unlabeled
    R_unlabeled+=R0_unlabeled

    I_unlabeled[int(np.floor(label_time/7)+1):]-=I0.sum()
    E_unlabeled[int(np.floor(label_time/7)+1):]-=E0.sum()
    S_unlabeled[int(np.floor(label_time/7)+1):]-=S0.sum()
    R_unlabeled[int(np.floor(label_time/7)+1):]-=R0.sum()
    ###
    df_I = pd.DataFrame(I, columns = weeks/7)
    df_I['lineage'] = np.arange(num_lineages)

    df_I_unlabeled = pd.DataFrame([I_unlabeled], columns = weeks/7)
    df_I_unlabeled['lineage'] = np.nan

    df_I = df_I.append(df_I_unlabeled)
    df_I.reset_index(inplace = True, drop = True)
    ###
    df_E = pd.DataFrame(E, columns = weeks/7)
    df_E['lineage'] = np.arange(num_lineages)

    df_E_unlabeled = pd.DataFrame([E_unlabeled], columns = weeks/7)
    df_E_unlabeled['lineage'] = np.nan

    df_E = df_E.append(df_E_unlabeled)
    df_E.reset_index(inplace = True, drop = True)
    ###
    df_S_unlabeled = pd.DataFrame([S_unlabeled], columns = weeks/7)
    df_S_unlabeled.reset_index(inplace = True, drop = True)
    ###
    df_R = pd.DataFrame(R, columns = weeks/7)
    df_R['lineage'] = np.arange(num_lineages)

    df_R_unlabeled = pd.DataFrame([R_unlabeled], columns = weeks/7)
    df_R_unlabeled['lineage'] = np.nan

    df_R = df_R.append(df_R_unlabeled)
    df_R.reset_index(inplace = True, drop = True)
    return df_S_unlabeled, df_E, df_I, df_R

def save_results(output_folder, df_S_total, df_E, df_I, df_R):
    df_E_total = pd.DataFrame(df_E.drop(['lineage'], axis = 1).sum()).transpose()
    df_I_total = pd.DataFrame(df_I.drop(['lineage'], axis = 1).sum()).transpose()
    df_R_total = pd.DataFrame(df_R.drop(['lineage'], axis = 1).sum()).transpose()
    
    df_E.to_csv(output_folder + 'counts_exposed.csv')
    df_I.to_csv(output_folder + 'counts_infected.csv')
    df_R.to_csv(output_folder + 'counts_recovered.csv')

    df_S_total.to_csv(output_folder + 'total_counts_susceptible.csv')
    df_E_total.to_csv(output_folder + 'total_counts_exposed.csv')
    df_I_total.to_csv(output_folder + 'total_counts_infected.csv')
    df_R_total.to_csv(output_folder + 'total_counts_recovered.csv')
    
    return df_E_total, df_I_total, df_R_total
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:37:34 2022

@author: qinqinyu
"""

import pandas as pd
import numpy as np
import os

def wright_fisher_all_lineages(Ne, f1):
    n2 = np.random.multinomial(Ne, f1)
    f2 = n2/Ne
    return n2, f2

def sampling_noise(Nseq, c, fj):    
    if fj == 0:
        nobsj = 0
    elif fj == 1:
        nobsj = Nseq
    else:
        m = fj*Nseq
        v = m*c
        p = m/v
        n = m**2/(v-m)
        nobsj = np.random.negative_binomial(n,p)
    fobsj = nobsj/Nseq
    return nobsj, fobsj
        
def sampling_noise_all_lineages(Nseq, c, f):
    nobs = np.empty(f.shape)
    fobs = np.empty(f.shape)
    for j in range(len(f)):
        fj = f[j]
        nobsj, fobsj = sampling_noise(Nseq, c, fj)
        nobs[j] = nobsj
        fobs[j] = fobsj
    return nobs, fobs

def gaussian(x, mu, sig, Nmin, Nmax):
    A = Nmax*(sig*np.sqrt(2*np.pi))-Nmin
    return A*(1/(sig*np.sqrt(2*np.pi)))*np.exp(-(((x-mu)/sig)**2)/2)+Nmin

def rectangular(x, x1, x2, Nmin, Nmax):
    return (Nmax-Nmin)*(np.heaviside(x-x1, 0.5)-np.heaviside(x-x2, 0.5))+Nmin

def fitness(numlineages, s_mean, s_std):
    return np.random.normal(loc=s_mean, scale=s_std, size=numlineages)

def run_simulation(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, s_mean, s_std, mu, numlineages):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    t_burnin = np.arange(0, total_burnin_time).astype('float')
    t = np.arange(0, total_epiweeks).astype('float')
    f0 = np.array([1/numlineages]*numlineages)
    n0 = np.array([1]*numlineages)
    s0 = fitness(numlineages, s_mean, s_std)

    #initialize
    f = f0
    n = n0
    s = s0
    
    #burn-in time
    for k in range(0, len(t_burnin)):
        Nei = Net[0]
        f_selection_unnorm = f*np.exp(s)
        f_selection_norm = f_selection_unnorm/np.sum(f_selection_unnorm)
        n, f = wright_fisher_all_lineages(Nei, f_selection_norm)
        # Mutations
        n_mut = np.random.binomial(n.astype('int'), mu)
        n_mut_total = np.sum(n_mut)
        n = n-n_mut
        n = np.append(n, np.ones(n_mut_total))
        f = n/Nei
        s = np.append(s, fitness(n_mut_total, s_mean, s_std))
        
    # Drop lineages that have already gone extinct
    s = s[f>0]
    f = f[f>0]
    n = n[n>0]
    
    # Initialize dataframes for saving results
    counts_sim_all = pd.DataFrame()
    counts_sim_all['lineage'] = np.array(range(len(f)))
    
    counts_sim_actual_all = pd.DataFrame()
    counts_sim_actual_all['lineage'] = np.array(range(len(f)))
    
    lineage_counter = len(f)
    
    #for each timepoint
    for j in range(0, len(t)):
        epiweek = t[j]
        Nei = Net[j]
        Nseqi = Nseq[j]
        c = ct[j]
        if j == 0:
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
        else:
            f_selection_unnorm = f*np.exp(s)
            f_selection_norm = f_selection_unnorm/np.sum(f_selection_unnorm)
            n, f = wright_fisher_all_lineages(Nei, f_selection_norm)
            
            # Mutations
            n_mut = np.random.binomial(n.astype('int'), mu)
            n_mut_total = np.sum(n_mut)
            n = n-n_mut
            n = np.append(n, np.ones(n_mut_total))
            f = n/Nei
            s = np.append(s, fitness(n_mut_total, s_mean, s_std))
            counts_sim_all = counts_sim_all.append(pd.DataFrame({'lineage':np.arange(lineage_counter, lineage_counter + n_mut_total)}))
            counts_sim_actual_all = counts_sim_actual_all.append(pd.DataFrame({'lineage':np.arange(lineage_counter, lineage_counter + n_mut_total)}))

            lineage_counter += n_mut_total
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
            
        #save to df
        counts_sim_all[epiweek] = nobs
        counts_sim_actual_all[epiweek] = n
    
    
    counts_sim_all.fillna(0, inplace = True)
    counts_sim_actual_all.fillna(0, inplace = True)
    
    counts_sim_all.reset_index(drop = True, inplace = True)
    counts_sim_actual_all.reset_index(drop = True, inplace = True)
    
    counts_sim_all = counts_sim_all.copy()
    total_neutral_counts_sim = pd.DataFrame(counts_sim_all[t].sum(axis = 0))
    total_neutral_counts_sim = total_neutral_counts_sim.transpose()
    
    # Save results
    counts_sim_all.to_csv(output_folder + counts_output_filename)
#    counts_sim_actual_all.to_csv(output_folder + counts_actual_output_filename)
    total_neutral_counts_sim.to_csv(output_folder + total_counts_output_filename)
    
    if (s_mean!=0)&(s_std!=0):
        lineage_fitnesses = pd.DataFrame({'lineage':counts_sim_all['lineage'],  's':s})
        lineage_fitnesses.to_csv(output_folder + fitness_filename)
    
    return counts_sim_all, total_neutral_counts_sim
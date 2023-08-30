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

def run_simulation(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, mu, numlineages, dfe=[0]):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    t_burnin = np.arange(0, total_burnin_time).astype('float')
    t = np.arange(0, total_epiweeks).astype('float')
    f0 = np.array([1/numlineages]*numlineages)
    n0 = np.array([1]*numlineages)
    s0 = np.random.choice(dfe, size=numlineages)

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
        s = np.append(s, np.random.choice(dfe, size=n_mut_total))
        
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
            s = np.append(s, np.random.choice(dfe, size=n_mut_total))
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
    
    if np.array(dfe).any()!=0:
        lineage_fitnesses = pd.DataFrame({'lineage':counts_sim_all['lineage'],  's':s})
        lineage_fitnesses.to_csv(output_folder + fitness_filename)
    
    return counts_sim_all, total_neutral_counts_sim

def run_simulation_nb_offspring_dist(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, mu, numlineages, offspring_var):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Get parameters of negative binomial distribution
    offspring_mean = 1 # Assume mean of 1
    nb_p = offspring_mean/offspring_var
    nb_n = offspring_mean**2/(offspring_var-offspring_mean)
    
    t_burnin = np.arange(0, total_burnin_time).astype('float')
    t = np.arange(0, total_epiweeks).astype('float')
    f0 = np.array([1/numlineages]*numlineages)
    n0 = np.array([1]*numlineages)

    #initialize
    f = f0
    n = n0
    
    #burn-in time
    for k in range(0, len(t_burnin)):
        Nei = Net[0]
        n_offspring = np.empty(len(n))
        for m in range(len(n)):
            n_m = n[m]
            n_m_offspring = np.sum(np.random.negative_binomial(nb_n, nb_p, size = (int(n_m))))
            n_offspring[m] = n_m_offspring
        f_offspring = n_offspring/np.sum(n_offspring)
        n, f = wright_fisher_all_lineages(Nei, f_offspring)
        # Mutations
        n_mut = np.random.binomial(n.astype('int'), mu)
        n_mut_total = np.sum(n_mut)
        n = n-n_mut
        n = np.append(n, np.ones(n_mut_total))
        f = n/Nei
        
    # Drop lineages that have already gone extinct
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
            n_offspring = np.empty(len(n))
            for m in range(len(n)):
                n_m = n[m]
                n_m_offspring = np.sum(np.random.negative_binomial(nb_n, nb_p, size = (int(n_m))))
                n_offspring[m] = n_m_offspring
            f_offspring = n_offspring/np.sum(n_offspring)
            n, f = wright_fisher_all_lineages(Nei, f_offspring)
            
            # Mutations
            n_mut = np.random.binomial(n.astype('int'), mu)
            n_mut_total = np.sum(n_mut)
            n = n-n_mut
            n = np.append(n, np.ones(n_mut_total))
            f = n/Nei
            counts_sim_all = counts_sim_all.append(pd.DataFrame({'lineage':np.arange(lineage_counter, lineage_counter + n_mut_total)}))
            counts_sim_actual_all = counts_sim_actual_all.append(pd.DataFrame({'lineage':np.arange(lineage_counter, lineage_counter + n_mut_total)}))

            lineage_counter += n_mut_total
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
            
        #save to df
        counts_sim_all[epiweek] = nobs
        counts_sim_actual_all[epiweek] = n
    
    
    counts_sim_all.fillna(0, inplace = True)
    counts_sim_actual_all.fillna(0, inplace = True)
    
    # Drop lineage rows with all 0 entries
    counts_sim_all = counts_sim_all.loc[(counts_sim_all[t]!=0).any(axis=1)]
    counts_sim_actual_all = counts_sim_actual_all.loc[(counts_sim_actual_all[t]!=0).any(axis=1)]
    
    counts_sim_all.reset_index(drop = True, inplace = True)
    counts_sim_actual_all.reset_index(drop = True, inplace = True)
    
    counts_sim_all = counts_sim_all.copy()
    total_neutral_counts_sim = pd.DataFrame(counts_sim_all[t].sum(axis = 0))
    total_neutral_counts_sim = total_neutral_counts_sim.transpose()
    
    # Save results
    counts_sim_all.to_csv(output_folder + counts_output_filename)
#    counts_sim_actual_all.to_csv(output_folder + counts_actual_output_filename)
    total_neutral_counts_sim.to_csv(output_folder + total_counts_output_filename)
    
    return counts_sim_all, total_neutral_counts_sim

def run_simulation_dfe(path_folder, output_folder, counts_output_filename, total_counts_output_filename, fitness_filename, total_epiweeks, total_burnin_time, Net, ct, Nseq, mu, numlineages, dfe=[0]):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    t = np.arange(0, total_epiweeks).astype('float')

    # Initialize

    n0 = np.ones(int(Net[0]))
    s0 = np.random.choice(dfe, size=int(Net[0]))
    f0 = n0/Net[0]
    l0 = np.arange(0, Net[0])

    n = n0
    s = s0
    f = f0
    l = l0

    # Burn-in time until the number of lineages reaches a threshold
    while len(np.unique(l))>numlineages:

        Nei = int(Net[0]) # Use the initial population size for the burn-in time

        weights = f*np.exp(s)
        weights = weights/np.sum(weights)

        n = np.random.multinomial(Nei, weights)

        # Mutations
        n_mut = np.random.binomial(n.astype('int'), mu)
        n_mut_total = np.sum(n_mut)
        n = n-n_mut
        n = np.append(n, np.ones(n_mut_total))

        s_mut = np.array([])
        l_mut = np.array([])
        for i in range(len(n_mut)):
            s_mut = np.concatenate((s_mut, s[i] + np.random.choice(dfe, size=n_mut[i])))
            l_mut = np.concatenate((l_mut, np.array([l[i]]*n_mut[i])))
        s = np.concatenate((s, s_mut))
        l = np.concatenate((l, l_mut))
        f = n/Nei

        # Drop zeros
        s = s[n>0]
        f = f[n>0]
        l = l[n>0]
        n = n[n>0]

    # Initialize dataframes for saving results
    counts_sim_all = pd.DataFrame()
    counts_sim_all['lineage'] = np.array(np.unique(l))

    counts_sim_actual_all = pd.DataFrame()
    counts_sim_actual_all['lineage'] = np.array(np.unique(l))

    # For each timepoint
    for j in range(0, len(t)):
        epiweek = t[j]
        Nei = int(Net[j])
        Nseqi = Nseq[j]
        c = ct[j]
        if j == 0:
            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)
        else:
            weights = f*np.exp(s)
            weights = weights/np.sum(weights)

            n = np.random.multinomial(Nei, weights)

            # Mutations
            n_mut = np.random.binomial(n.astype('int'), mu)
            n_mut_total = np.sum(n_mut)
            n = n-n_mut
            n = np.append(n, np.ones(n_mut_total))

            s_mut = np.array([])
            l_mut = np.array([])
            for i in range(len(n_mut)):
                s_mut = np.concatenate((s_mut, s[i] + np.random.choice(dfe, size=n_mut[i])))
                l_mut = np.concatenate((l_mut, np.array([l[i]]*n_mut[i])))
            s = np.concatenate((s, s_mut))
            l = np.concatenate((l, l_mut))
            f = n/Nei

            # Drop zeros
            s = s[n>0]
            f = f[n>0]
            l = l[n>0]
            n = n[n>0]

            nobs, fobs = sampling_noise_all_lineages(Nseqi, c, f)

        # Save to df
        nobs_by_lineage = pd.DataFrame({'l':l, 'nobs':nobs}).groupby('l').sum()
        nobs_lineages = nobs_by_lineage.values.flatten()
        lineages = nobs_by_lineage.index.values
        for i in range(len(lineages)):
            lineage = lineages[i]
            nobs_lineage = nobs_lineages[i]
            counts_sim_all.loc[counts_sim_all['lineage']==lineage,epiweek] = nobs_lineage

        n_by_lineage = pd.DataFrame({'l':l, 'n':n}).groupby('l').sum()
        n_lineages = n_by_lineage.values.flatten()
        lineages = n_by_lineage.index.values
        for i in range(len(lineages)):
            lineage = lineages[i]
            n_lineage = n_lineages[i]
            counts_sim_actual_all.loc[counts_sim_all['lineage']==lineage,epiweek] = n_lineage

    counts_sim_all['lineage'] = np.arange(len(counts_sim_all))
    counts_sim_actual_all['lineage'] = np.arange(len(counts_sim_actual_all))

    counts_sim_all.fillna(0, inplace = True)
    counts_sim_actual_all.fillna(0, inplace = True)

    # Drop lineage rows with all 0 entries
    counts_sim_all = counts_sim_all.loc[(counts_sim_all[t]!=0).any(axis=1)]
    counts_sim_actual_all = counts_sim_actual_all.loc[(counts_sim_actual_all[t]!=0).any(axis=1)]

    counts_sim_all.reset_index(drop = True, inplace = True)
    counts_sim_actual_all.reset_index(drop = True, inplace = True)

    counts_sim_all = counts_sim_all.copy()
    total_counts_sim = pd.DataFrame(counts_sim_all[t].sum(axis = 0))
    total_counts_sim = total_counts_sim.transpose()
       
    # Save results
    counts_sim_all.to_csv(output_folder + counts_output_filename)
#    counts_sim_actual_all.to_csv(output_folder + counts_actual_output_filename)
    total_counts_sim.to_csv(output_folder + total_counts_output_filename)
    
    return counts_sim_all, total_counts_sim

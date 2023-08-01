#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:28:11 2021

@author: qinqinyu
"""

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style = 'whitegrid', font_scale = 1.5)
import format_data as fd
sns.set(font_scale = 1.5)

def get_kappa(counts, epiweeks, total_sampled_counts, mincount = 20, min_freq = 0, max_freq = np.inf):
    '''
    Calculates kappa for each epiweek pair, where kappa is the variance in the 
    square root frequency changes calculated through the mean squared displacement.
    Uses the neutral frequencies of each day pair (calculated by inputting
    neutral counts, total counts, and fraction of neutral counts).
    
    INPUTS:
        
        counts: Dataframe with NEUTRAL lineage counts over time. Each row is a 
        lineage and each column is a timepoint. One additional column gives the 
        lineage name.
        
        epiweeks: Array of timepoints (same as timepoint column names in other
        dataframes to be used).
        
        total_sampled_counts: Dataframe with the total number of sampled counts across all 
        lineages as a function of time. Columns are each timepoint.
    
        min_freq: Minimum frequency to be included in calculating kappa for 
        first epiweek of each epiweek pair. Default is 0.
        
        max_freq: Maximum frequency to be included in calculating kappa for 
        first epiweek of each epiweek pair. Default is infinity.
        
    OUTPUTS:
        
        k_df: Dataframe of calculated kappas for each epiweek pair. kappa is
        the variance in the square root frequency changes calculated through 
        the mean square displacement. 
        
    '''
    #Format data
    d1=epiweeks[0]
    dl=epiweeks[-1]
    
#    counts = fd.create_superlineages(counts, d1, dl, mincount)

    frequency = pd.DataFrame(counts[epiweeks].values/total_sampled_counts[epiweeks].values, columns = epiweeks)
    frequency['lineage'] = counts['lineage']
    
    frequency_int = frequency.copy()
    Nt = total_sampled_counts.copy()
    days = epiweeks
    
    # Get timepairs
    daypairs = []
    for i in range(len(days)):
        day = days[i]
        for j in range(i+1, len(days)):
            daypairs.append((day, days[j]))
    
    # Calculate kappa for each time pair
    kdata={
        'Time 1':[],
        'Time 2':[],
        'kappa':[],
        'N1':[],
        'N2':[],
        'stderr':[]
    }

    for dp in daypairs:
        d1,d2=dp
        N1=int(Nt[d1])
        N2=int(Nt[d2])
    
        cd = frequency_int[[d1,d2]]
    
        cond1=cd[d1]>min_freq
        cond2=cd[d1]<max_freq       
    
        cd_cond=cd[cond1 & cond2].reset_index(drop=True)
        cd_cond_sqrt = np.sqrt(cd_cond)
        kv_cond=cd_cond_sqrt[d2] - cd_cond_sqrt[d1]
        kappa_t = np.mean(kv_cond**2)
        kappa_t_err = np.std(kv_cond**2)/np.sqrt(len(cd_cond))
        
        #If calculating median instead of mean
            
        # kappa_t = 2*np.median(kv_cond**2)
        # mad = median_absolute_deviation(kv_cond**2)
        # kappa_t_err = mad/(0.67449*np.sqrt(len(cd_cond)))
    
        kdata['Time 1'].append(d1)
        kdata['Time 2'].append(d2)
        kdata['kappa'].append(kappa_t)
        kdata['N1'].append(N1)
        kdata['N2'].append(N2)
        kdata['stderr'].append(kappa_t_err)
    
    k_df=pd.DataFrame(kdata)
    return k_df

def fit_time_varying_technical_noise_bounded(k_df, epiweeks, plot = 0, savefig = 0, figname = None, numtrials = 10**3):
    ''' 
    Fits time varying technical/sampling/measurement noise using kappas
    '''
    days = epiweeks
        
    Ne0=1e3 # initial guess for Ne

    def obj_single(kappa,kappa_err,c1,c2,Ne,c1coeff,c2coeff,Necoeff):
        kappa_est = c1*c1coeff + c2*c2coeff + (1/Ne)*Necoeff
        return ((kappa_est - float(kappa))**2/float(kappa_err))
    
    def add_obj(theta,k_df,tv):
        mse=0
        Ne=np.exp(-theta[-1])
        for i in range(len(k_df)):
            k_df_i = k_df.iloc[i]
            d1 = float(k_df_i['Time 1'])
            d2 = float(k_df_i['Time 2'])
            
            d1_idx = np.where(days == str(k_df_i['Time 1']))[0][0]
            d2_idx = np.where(days == str(k_df_i['Time 2']))[0][0]
        
            idxs = np.arange(d1+1, d2+1, 1).astype('str')
            sum_N = len(idxs)
            c1=np.exp(theta[d1_idx])
            c2=np.exp(theta[d2_idx])
            c1coeff = 1/(4*k_df_i['N1'])
            c2coeff = 1/(4*k_df_i['N2'])
            Necoeff = sum_N/4
            mse+=obj_single(k_df_i['kappa'],k_df_i['stderr'],c1,c2,Ne,c1coeff,c2coeff,Necoeff)
        
        return mse
    
    tv=list(k_df['Time 1']) + list(k_df['Time 2'])
    tv=list(set([str(t) for t in tv]))
    bounds=[]
    for i,dl in enumerate(tv):
        bounds.append((np.log(1),np.log(100)))
    bounds.append((np.log(1/1e7), np.log(1)))
    
    theta0=list(np.zeros(len(tv)))
    theta0.append(np.log(1/Ne0))
    
    res=sp.optimize.minimize(lambda theta: add_obj(theta,k_df,tv),theta0, method='L-BFGS-B',bounds=bounds)
    theta = np.exp(res.x)
    c = theta[:-1]
    Netau = 1/theta[-1]
        
    # Calculate the inferred kappa, use difference with true kappa as error 
    # bars to get errors on betas
    
    coef = np.zeros((len(k_df), len(days)+1))
        
    for i in range(len(k_df)):
        k_df_i = k_df.iloc[i]
        d1 = float(k_df_i['Time 1'])
        d2 = float(k_df_i['Time 2'])
        
        d1_idx = np.where(days == str(k_df_i['Time 1']))[0][0]
        d2_idx = np.where(days == str(k_df_i['Time 2']))[0][0]
    
        idxs = np.arange(d1+1, d2+1, 1).astype('str')
        sum_N = len(idxs)
        
        coef[i,d1_idx] = 1/(4*k_df_i['N1'])
        coef[i,d2_idx] = 1/(4*k_df_i['N2'])
        coef[i,-1] = sum_N/4
    
    kappa = k_df['kappa'].values
    kappa_err = k_df['stderr'].values
    
    kappa_infer = np.matmul(coef,theta)
    kappa_err_tot = np.sqrt((kappa-kappa_infer)**2+kappa_err**2)
    
    if numtrials==0:
        summary_infer = pd.DataFrame(data = {'Time':days.astype('float'), 'c':c, 'c_mean':[np.nan]*len(c), 'c_median':[np.nan]*len(c), 'c_err':[np.nan]*len(c), 'c_95_lower':[np.nan]*len(c), 'c_95_upper':[np.nan]*len(c), 'Netau':[Netau]*len(c), 'Netau_mean':[np.nan]*len(c), 'Netau_err':[np.nan]*len(c), 'Ne_95_lower':[np.nan]*len(c), 'Ne_95_upper':[np.nan]*len(c)})
        
        k_df['kappa_inferred'] = kappa_infer
        
        c_all = pd.DataFrame([summary_infer['c'].values], columns = epiweeks)
        
        df_all_trials = pd.DataFrame()
        
    elif numtrials>0:

        # Get errors
        for i in range(numtrials):
            print(i)
            kappa_i = np.empty(len(kappa))
            kappa_i.fill(np.nan)
            for j in range(len(kappa)):
                kappa_i[j] = np.random.normal(loc = kappa[j], scale = kappa_err_tot[j], size = None)
            k_df_i = k_df.copy()
            k_df_i['kappa'] = kappa_i
            res=sp.optimize.minimize(lambda theta: add_obj(theta,k_df_i,tv),theta0, method='L-BFGS-B',bounds=bounds)

            theta = np.exp(res.x)
            if i == 0:
                x_all = theta
            else:
                x_all = np.vstack([x_all, theta])
        x_all_c = x_all[:,:-1]
        x_all_Ne = 1/x_all[:,-1]
        c_95_lower = np.quantile(x_all_c, 0.025, axis = 0)
        c_95_upper = np.quantile(x_all_c, 0.975, axis = 0)
        Ne_95_lower = np.quantile(x_all_Ne, 0.025, axis = 0)
        Ne_95_upper = np.quantile(x_all_Ne, 0.975, axis = 0)
        
        x_mean = np.mean(x_all, axis = 0)
        x_median = np.median(x_all, axis = 0)
        x_err = np.std(x_all, axis = 0)
        
        c_mean = x_mean[:-1]
        c_err = x_err[:-1]
        c_median = x_median[:-1]
        Netau_mean = 1/(x_mean[-1])
        Netau_err = (Netau_mean*x_err[-1])/x_mean[-1]
        
        summary_infer = pd.DataFrame(data = {'Time':days.astype('float'), 'c':c, 'c_mean':c_mean, 'c_median':c_median, 'c_err':c_err, 'c_95_lower':c_95_lower, 'c_95_upper':c_95_upper, 'Netau':[Netau]*len(c), 'Netau_mean':[Netau_mean]*len(c), 'Netau_err':[Netau_err]*len(c), 'Ne_95_lower':[Ne_95_lower]*len(c), 'Ne_95_upper':[Ne_95_upper]*len(c)})
        
        k_df['kappa_inferred'] = kappa_infer
        
        c_all = pd.DataFrame([summary_infer['c'].values], columns = epiweeks)
        
        df_all_trials = pd.DataFrame()
        df_all_trials[epiweeks] = x_all_c
        df_all_trials['Netau'] = x_all_Ne
        df_all_trials['trial'] = range(len(df_all_trials))
    
    if plot:
        plt.errorbar(summary_infer['Time'], summary_infer['c'], summary_infer['c_err'], marker = 'o')
        plt.xlabel('Epiweek')
        plt.ylabel('c')
        plt.tight_layout()
        if savefig:
            if figname is None:
                print('Figure not saved. Specify output path name')
            else:
                plt.savefig(figname, dpi = 300)
        plt.show()

    return c_all, summary_infer, k_df, df_all_trials
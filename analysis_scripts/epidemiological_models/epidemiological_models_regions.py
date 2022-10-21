#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:05:19 2022

@author: qinqinyu
"""
import sys
bin_path = '../..functions'
if bin_path not in sys.path:
    sys.path.insert(1, '../..functions')
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import median_absolute_deviation
import scipy.stats as st
from epiweeks import Week, Year
from datetime import timedelta
import matplotlib.pyplot as plt

import epidemiological_models_functions as emf

######## Define parameters ########

gamma_E = 1/(3/7) # in units of inverse weeks
gamma_I = 1/(5.5/7)

######## Import and format data on the effective reproduction number and the estimated community cases ########

regions = ['East Midlands',
           'East of England',
           'London',
           'North East',
           'North West',
           'South East',
           'South West',
           'West Midlands',
           'Yorkshire and The Humber']

for region in regions:
    # Get the number of positives and effective reproduction number data
    filename = '../..data/effective_reproduction_number/211210_R_and_growth_rate_time_series_for_publication_v1.0_regions_cleaned_up.csv'
    Rt = pd.read_csv(filename)
    Rt.at[Rt['Upper bound']=='.','Upper bound'] = np.nan
    Rt.at[Rt['Lower bound']=='.','Lower bound'] = np.nan
    Rt = Rt.astype({'Upper bound':'float', 'Lower bound':'float'})
    Rt = Rt[Rt['Region']==region]

    cond1 = Rt['Certainty criteria met'].astype('bool')
    cond2 = ~Rt['Estimate based on fewer days or lower quality data'].astype('bool')
    Rt_cond = Rt[cond1&cond2]
    Rt_cond.reset_index(inplace = True, drop = True)

    filename = '../..data/cases/covid19infectionsurveydatasets20211210england1_regions_estimated_number_of_people_testing_positive_cleaned_up_community_pop_size.csv'
    positives = pd.read_csv(filename)
    positives = positives.astype({'Estimated number of people testing positive for COVID-19':'float', 
               '95% Upper confidence interval number of people':'float',
               '95% Lower confidence interval number of people':'float'})
    positives = positives[positives['Region']==region]
    
    # Interpolate number of cases
    epiweeks = positives['Epiweek'][:-1]+1
    date = np.full(len(epiweeks), 'Interpolated')
    I_region = positives['Region'][:-1]
    method = np.full(len(epiweeks), 'Interpolated')
    population_size = positives['Population size'][:-1]
    positives_interpolated = positives.append(pd.DataFrame({'Epiweek':epiweeks, 'Time period':date, 'Region':I_region, 'Method':method, 'Population size':population_size}))
    positives_interpolated = positives_interpolated.sort_values(['Epiweek'])
    positives_interpolated = positives_interpolated.interpolate()
    positives = positives_interpolated.copy()
    positives.reset_index(inplace = True, drop = True)

    merged = Rt_cond.merge(positives, on = ['Epiweek', 'Region'])
    merged.drop(columns = ['Estimated % testing positive for COVID-19',
                 '95% Lower confidence interval',
                 '95% Upper confidence interval',
                 'Population size',
                 'Region'], inplace = True)

    merged.rename(columns={'Date':'Rt_date',
                           'Lower bound':'Rt_lower',
                           'Upper bound':'Rt_upper',
                           'Certainty criteria met': 'Rt_estimate_certainty_criteria_met',
                           'Estimate based on fewer days or lower quality data': 'Rt_estimate_based_on_fewer_days_or_lower_quality_data',
                           'Time period': 'I_date',
                           'Estimated number of people testing positive for COVID-19': 'I',
                           '95% Upper confidence interval number of people': 'I_upper',
                           '95% Lower confidence interval number of people': 'I_lower',
                           'Method':'I_method'}, inplace=True)
    Rt = (merged['Rt_upper'] + merged['Rt_lower'])/2
    merged['Rt'] = Rt
    
    dates = []
    for w in merged['Epiweek']:
        w = int(np.floor(w))
        if w<=53:
            week = Week(2020, w)
        elif w<=105:
            week = Week(2021, w-53)
        else:
            week = Week(2022, w-105)
        dates.append(week.startdate()+timedelta(days=3))
        
    dates = np.array(dates)
    merged['date'] = dates
    ######## GET SIR AND SEIR MODEL NETAU FOR ALL VARIANTS OR LINEAGES ########
    
    merged = emf.get_sir_model_Netau(merged, gamma_I)
    merged = emf.get_seir_model_Netau(merged, gamma_E, gamma_I)
    merged.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_overall.csv')

######## GET SIR AND SEIR MODEL NETAU FOR BY VARIANTS OR LINEAGE ########

for region in regions:
        
    merged = pd.read_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_overall.csv', index_col = 0)
    merged_original = merged.copy()

    # Get the fraction of each variant
    frac_neutral = pd.read_csv('../..data/lineages/pre_B-1-177/microreact/is_pillar_2/' + region + '/frac_variant_with_b_1_177.csv', index_col = 0)
    frac_neutral_df = pd.DataFrame(data = {'Epiweek':frac_neutral.columns.astype('float'), 'frac_nonvariant_with_b_1_177':frac_neutral.values[0]})
    merged = merged.merge(frac_neutral_df, on = 'Epiweek', how = 'left')
    merged['frac_nonvariant_with_b_1_177'] = merged['frac_nonvariant_with_b_1_177'].fillna(0)

    frac_neutral = pd.read_csv('../..data/lineages/pre_B-1-177/microreact/is_pillar_2/' + region + '/frac_variant.csv', index_col = 0)
    frac_neutral_df = pd.DataFrame(data = {'Epiweek':frac_neutral.columns.astype('float'), 'frac_nonvariant':frac_neutral.values[0]})
    merged = merged.merge(frac_neutral_df, on = 'Epiweek', how = 'left')
    merged['frac_nonvariant'] = merged['frac_nonvariant'].fillna(0)

    frac_neutral = pd.read_csv('../..data/lineages/B-1-177/B-1-177|2021-02-22|694.5/is_pillar_2/'+ region + '/frac_variant.csv', index_col = 0)
    frac_neutral_df = pd.DataFrame(data = {'Epiweek':frac_neutral.columns.astype('float'), 'frac_B_1_177':frac_neutral.values[0]})
    merged = merged.merge(frac_neutral_df, on = 'Epiweek', how = 'left')
    merged['frac_B_1_177'] = merged['frac_B_1_177'].fillna(0)

    frac_neutral = pd.read_csv('../..data/lineages/alpha/alpha|2021-06-20|61.5/is_pillar_2/' + region + '/frac_variant.csv', index_col = 0)
    frac_neutral_df = pd.DataFrame(data = {'Epiweek':frac_neutral.columns.astype('float'), 'frac_alpha':frac_neutral.values[0]})
    merged = merged.merge(frac_neutral_df, on = 'Epiweek', how = 'left')
    merged['frac_alpha'] = merged['frac_alpha'].fillna(0)

    frac_neutral = pd.read_csv('../..data/lineages/delta/delta|2022-01-25|50.5+58.5/is_pillar_2/' + region + '/frac_variant.csv', index_col = 0)
    frac_neutral_df = pd.DataFrame(data = {'Epiweek':frac_neutral.columns.astype('float'), 'frac_delta':frac_neutral.values[0]})
    merged = merged.merge(frac_neutral_df, on = 'Epiweek', how = 'left')
    merged['frac_delta'] = merged['frac_delta'].fillna(0)

    # Calculate the Rt of each variant
    Rt_rel_nonvariant = 1 
    Rt_rel_alpha = 1.7 # Volz et al., Nature, 2021
    Rt_rel_delta = 1.97 # Campbell et al., Euro Surveill, 2021
    
    denom = (merged['frac_nonvariant_with_b_1_177']*Rt_rel_nonvariant + 
                                  merged['frac_alpha']*Rt_rel_alpha + 
                                  merged['frac_delta']*Rt_rel_delta)
    
    
    # Nonvariant
    merged_variant = merged_original.copy()
    merged_variant = emf.get_sir_model_Netau_by_variant(merged_variant, merged['frac_nonvariant'], Rt_rel_nonvariant, denom, gamma_I)
    merged_variant = emf.get_seir_model_Netau_by_variant(merged_variant, merged['frac_nonvariant'], Rt_rel_nonvariant, denom, gamma_E, gamma_I)
    merged_variant.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_pre_B-1-177.csv')
    
    # B.1.177
    merged_variant = merged_original.copy()
    merged_variant = emf.get_sir_model_Netau_by_variant(merged_variant, merged['frac_B_1_177'], Rt_rel_nonvariant, denom, gamma_I)
    merged_variant = emf.get_seir_model_Netau_by_variant(merged_variant, merged['frac_B_1_177'], Rt_rel_nonvariant, denom, gamma_E, gamma_I)
    merged_variant.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_B-1-177.csv')
    
    # Alpha
    merged_variant = merged_original.copy()
    merged_variant = emf.get_sir_model_Netau_by_variant(merged_variant, merged['frac_alpha'], Rt_rel_alpha, denom, gamma_I)
    merged_variant = emf.get_seir_model_Netau_by_variant(merged_variant, merged['frac_alpha'], Rt_rel_nonvariant, denom, gamma_E, gamma_I)
    merged_variant.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_alpha.csv')
    
    # Delta
    merged_variant = merged_original.copy()
    merged_variant = emf.get_sir_model_Netau_by_variant(merged_variant, merged['frac_delta'], Rt_rel_delta, denom, gamma_I)
    merged_variant = emf.get_seir_model_Netau_by_variant(merged_variant, merged['frac_delta'], Rt_rel_nonvariant, denom, gamma_E, gamma_I)
    merged_variant.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_delta.csv')

######## MERGE WITH INFERENCE RESULTS ########

for region in regions:
    summary_all = pd.DataFrame()
    
    df_Netau_variant = pd.read_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_pre_B-1-177.csv', index_col = 0)
    path_folder = '../..data/lineages/pre_B-1-177/microreact/is_pillar_2/' + region + '/'
    summary = pd.read_csv(path_folder + '/inference_results/summary.csv', index_col = 0)
    summary['variant'] = 'pre-B.1.177'
    total_counts_lineages = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
    total_counts_metadata = pd.read_csv(path_folder + 'total_counts_metadata.csv', index_col = 0)
    num_sequences = pd.DataFrame({'Epiweek':total_counts_lineages.columns.astype('float'), 'sequences_analyzed':total_counts_lineages.values[0], 'sequences_total':total_counts_metadata.values[0]})
    summary = summary.merge(num_sequences, how = 'left', on = 'Epiweek')
    summary = summary.merge(df_Netau_variant, on = 'Epiweek', how = 'left')
    summary_all = summary_all.append(summary)
    
    df_Netau_variant = pd.read_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_B-1-177.csv', index_col = 0)
    path_folder = '../..data/lineages/B-1-177/B-1-177|2021-02-22|694.5/is_pillar_2/' + region + '/'
    summary = pd.read_csv(path_folder + 'inference_results/summary_corrected.csv', index_col = 0)
    summary['variant'] = 'B.1.177'
    total_counts_lineages = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
    total_counts_metadata = pd.read_csv(path_folder + 'total_counts_metadata.csv', index_col = 0)
    num_sequences = pd.DataFrame({'Epiweek':total_counts_lineages.columns.astype('float'), 'sequences_analyzed':total_counts_lineages.values[0], 'sequences_total':total_counts_metadata.values[0]})
    summary = summary.merge(num_sequences, how = 'left', on = 'Epiweek')
    summary = summary.merge(df_Netau_variant, on = 'Epiweek', how = 'left')
    summary_all = summary_all.append(summary)
    
    df_Netau_variant = pd.read_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_alpha.csv', index_col = 0)
    path_folder = '../..data/lineages/alpha/alpha|2021-06-20|61.5/is_pillar_2/'+ region + '/'
    summary = pd.read_csv(path_folder + 'inference_results/summary_corrected.csv', index_col = 0)
    summary['variant'] = 'Alpha'
    total_counts_lineages = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
    total_counts_metadata = pd.read_csv(path_folder + 'total_counts_metadata.csv', index_col = 0)
    num_sequences = pd.DataFrame({'Epiweek':total_counts_lineages.columns.astype('float'), 'sequences_analyzed':total_counts_lineages.values[0], 'sequences_total':total_counts_metadata.values[0]})
    summary = summary.merge(num_sequences, how = 'left', on = 'Epiweek')
    summary = summary.merge(df_Netau_variant, on = 'Epiweek', how = 'left')
    summary_all = summary_all.append(summary)
    
    df_Netau_variant = pd.read_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_delta.csv', index_col = 0)
    path_folder = '../..data/lineages/delta/delta|2022-01-25|49.5+58.5/is_pillar_2/' + region + '/'
    summary = pd.read_csv(path_folder + 'inference_results/summary_corrected.csv', index_col = 0)
    summary['variant'] = 'Delta'
    total_counts_lineages = pd.read_csv(path_folder + 'total_counts_lineages.csv', index_col = 0)
    total_counts_metadata = pd.read_csv(path_folder + 'total_counts_metadata.csv', index_col = 0)
    num_sequences = pd.DataFrame({'Epiweek':total_counts_lineages.columns.astype('float'), 'sequences_analyzed':total_counts_lineages.values[0], 'sequences_total':total_counts_metadata.values[0]})
    summary = summary.merge(num_sequences, how = 'left', on = 'Epiweek')
    summary = summary.merge(df_Netau_variant, on = 'Epiweek', how = 'left')
    summary_all = summary_all.append(summary)
    summary_all = summary_all[summary_all['Epiweek']>0]
    
    dates = []
    for w in summary_all['Epiweek']:
        w = int(np.floor(w))
        if w<=53:
            week = Week(2020, w)
        elif w<=105:
            week = Week(2021, w-53)
        else:
            week = Week(2022, w-105)
        dates.append(week.startdate()+timedelta(days=3))
        
    dates = np.array(dates)
    summary_all['date'] = dates
    
    summary_all.reset_index(inplace = True, drop = True)
    summary_all.to_csv('../..data/epidemiological_models/regions/' + region + '_epidemiological_models_Netau_inferred_Netau_by_variant.csv')
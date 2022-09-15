#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:05:15 2022

@author: qinqinyu
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import median_absolute_deviation
import scipy.stats as st
from epiweeks import Week, Year
from datetime import timedelta

def calculate_num_exposed(I, I_upper, I_lower, gamma_E = 1/(3/7), gamma_I = 1/(5.5/7)):
    dIdt_manual = np.empty(len(I))
    dIdt_manual.fill(np.nan)
    dIdt_upper_manual = np.empty(len(I))
    dIdt_upper_manual.fill(np.nan)
    dIdt_lower_manual = np.empty(len(I))
    dIdt_lower_manual.fill(np.nan)

    for i in range(1,len(I)-1):
        Ii = I[i]
        dIdt_manual[i] = (I[i+1]-I[i-1])/2
        dIdt_upper_manual[i] = (I_upper[i+1]-I_lower[i-1])/2
        dIdt_lower_manual[i] = (I_lower[i+1]-I_upper[i-1])/2

    E = (dIdt_manual+gamma_I*I)/gamma_E
    E_upper = (dIdt_upper_manual+gamma_I*I_upper)/gamma_E
    E_lower = (dIdt_lower_manual+gamma_I*I_lower)/gamma_E
    return E, E_upper, E_lower

def get_sir_model_Netau(merged, gamma_I):
    merged['Netau SIR'] = merged['I']/(2*merged['Rt']*gamma_I)
    merged['Netau SIR upper'] = merged['I_upper']/(2*merged['Rt_lower']*gamma_I) ### Improve the estimation of these error bars
    merged['Netau SIR lower'] = merged['I_lower']/(2*merged['Rt_upper']*gamma_I)

    merged['Netau SIR error criteria met'] = merged['Netau SIR upper']-merged['Netau SIR lower']<(3*merged['Netau SIR'])
    return merged

def get_seir_model_Netau(merged, gamma_E, gamma_I):
    I = merged['I']
    I_upper = merged['I_upper']
    I_lower = merged['I_lower']

    E, E_upper, E_lower = calculate_num_exposed(I, I_upper, I_lower, gamma_E = gamma_E, gamma_I = gamma_I)
    
    merged['E'] = E
    merged['E_upper'] = E_upper
    merged['E_lower'] = E_lower

    Netau_SEIR = (merged['E']+merged['I'])**2/(2*merged['Rt']*gamma_I*merged['I'])
    Netau_SEIR_upper = (merged['E_upper']+merged['I_upper'])**2/(2*merged['Rt_lower']*gamma_I*merged['I_lower'])
    Netau_SEIR_lower = (merged['E_lower']+merged['I_lower'])**2/(2*merged['Rt_upper']*gamma_I*merged['I_upper'])

    merged['Netau SEIR'] = Netau_SEIR
    merged['Netau SEIR upper'] = Netau_SEIR_upper
    merged['Netau SEIR lower'] = Netau_SEIR_lower

    merged['Netau SEIR error criteria met'] = (merged['Netau SEIR upper']-merged['Netau SEIR lower'])<(3*merged['Netau SEIR'])
    
    return merged

def get_sir_model_Netau_by_variant(merged_variant, frac, Rt_rel, denom, gamma_I):
    if 'frac_variant' not in merged_variant.columns:
        merged_variant['frac_variant'] = frac
    
    if 'Rt_variant' not in merged_variant.columns:
        merged_variant['Rt_variant'] = merged_variant['Rt']*Rt_rel/denom
        merged_variant['Rt_variant_lower'] = merged_variant['Rt_lower']*Rt_rel/denom
        merged_variant['Rt_variant_upper'] = merged_variant['Rt_upper']*Rt_rel/denom

    merged_variant['Rt_variant'] = merged_variant['Rt']*Rt_rel/denom
    merged_variant['Rt_variant_lower'] = merged_variant['Rt_lower']*Rt_rel/denom
    merged_variant['Rt_variant_upper'] = merged_variant['Rt_upper']*Rt_rel/denom

    merged_variant['Netau SIR variant'] = merged_variant['I']*merged_variant['frac_variant']/(2*merged_variant['Rt_variant']*gamma_I)
    merged_variant['Netau SIR variant upper'] = merged_variant['I_upper']*merged_variant['frac_variant']/(2*merged_variant['Rt_variant_lower']*gamma_I)
    merged_variant['Netau SIR variant lower'] = merged_variant['I_lower']*merged_variant['frac_variant']/(2*merged_variant['Rt_variant_upper']*gamma_I)

    merged_variant['Netau SIR variant error criteria met'] = merged_variant['Netau SIR variant upper']-merged_variant['Netau SIR variant lower']<(3*merged_variant['Netau SIR variant'])
    return merged_variant

def get_seir_model_Netau_by_variant(merged_variant, frac, Rt_rel, denom, gamma_E, gamma_I):    
    if 'frac_variant' not in merged_variant.columns:
        merged_variant['frac_variant'] = frac
    
    if 'Rt_variant' not in merged_variant.columns:
        merged_variant['Rt_variant'] = merged_variant['Rt']*Rt_rel/denom
        merged_variant['Rt_variant_lower'] = merged_variant['Rt_lower']*Rt_rel/denom
        merged_variant['Rt_variant_upper'] = merged_variant['Rt_upper']*Rt_rel/denom

    I = merged_variant['I']*merged_variant['frac_variant']
    I_upper = merged_variant['I_upper']*merged_variant['frac_variant']
    I_lower = merged_variant['I_lower']*merged_variant['frac_variant']

    E, E_upper, E_lower = calculate_num_exposed(I, I_upper, I_lower, gamma_E = gamma_E, gamma_I = gamma_I)

    merged_variant['Netau SEIR variant'] = (E+I)**2/(2*merged_variant['Rt_variant']*gamma_I*I)
    merged_variant['Netau SEIR variant upper'] = (E_upper + I_upper)**2/(2*merged_variant['Rt_variant_lower']*gamma_I*I_lower)
    merged_variant['Netau SEIR variant lower'] = (E_lower+I_lower)**2/(2*merged_variant['Rt_variant_upper']*gamma_I*I_upper)

    merged_variant['Netau SEIR variant error criteria met'] = (merged_variant['Netau SEIR variant upper']-merged_variant['Netau SEIR variant lower'])<(3*merged_variant['Netau SEIR variant'])
    
    return merged_variant
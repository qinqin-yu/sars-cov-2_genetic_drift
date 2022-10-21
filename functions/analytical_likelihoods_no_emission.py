'''
9.13.20
Joao Ascensao
3.25.21
QinQin Yu edit

Functions for priors, likelihoods, posteriors

'''
import pandas as pd
import scipy as sp 
import numpy as np 
from scipy import optimize
from scipy.optimize import minimize
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import format_data as fd
import itertools
sns.set(font_scale = 1.5)

def getshareddata(N):
    days = N.columns
    N = N.values[0]
    dd={
        'day':days,
        'N':N}
    return dd

def log_normal(x,mu,v):
    return  -((x - mu)**2)/(2*v) - 0.5*np.log(2*np.pi*v)

def hmm_forward_flat(data):
    nu,phi = data

    ll = 0
#    print(phi)
#    print(nu)
    for i in range(1,len(phi)):
#        print(i)
#        print(phi[i])
#        print(phi[i-1])
#        print(nu[0:i])
        ll += log_normal(phi[i], phi[i-1], nu)
    return ll

def hmm_add_likelihoods(theta, data_list, shareddata):
    ll=0

    N = shareddata['N']
    
    nu = theta/(4)
            
    for _,data in enumerate(data_list):
        
        nobs = data['nobs'].values
        phi = np.sqrt(nobs/(N))
        
        dd = (nu,phi)
        ll += hmm_forward_flat(dd)
    
    if pd.isnull(ll):
        ll=-1e7
    return ll*-1
    
def hmm_maxlikelihood_Ne(data_list, shareddata, bounds=(1e-11,1)):
    objective = lambda theta: hmm_add_likelihoods(theta, data_list, shareddata)
    res=optimize.minimize_scalar(objective,method='Bounded',bounds=bounds)
    if res.success:
        return 1/res.x
    else:
        return np.nan
    
### CALCULATING CIs WITH LIKELIHOOD FUNCTION ###
        
def p_value_array(Nes, likelihood_norm, norm, ll_unnorm_max, dNes, data_list, shareddata):
    thetas = 1/Nes
    objective = lambda theta: hmm_add_likelihoods(theta, data_list, shareddata)
    p_value = np.empty(len(thetas))
    for i in range(len(thetas)):
        theta = thetas[i]
        ll_unnorm_scalar = -1*objective(theta)
        ll_unnorm_scalar = ll_unnorm_scalar - ll_unnorm_max
        likelihood_unnorm_scalar = np.exp(ll_unnorm_scalar)
        likelihood_norm_scalar = likelihood_unnorm_scalar/norm
        p_value[i] = np.sum((likelihood_norm_scalar>likelihood_norm)*likelihood_norm*dNes)
    return p_value
    
def confidence_interval_Ne(data_list, shareddata, p_thresh = 0.05):
    
    objective = lambda theta: hmm_add_likelihoods(theta, data_list, shareddata)
    Nes = np.logspace(0, 11, 10**3)
    dNes = np.diff(Nes)
    Nes = Nes[:-1]
    thetas = 1/Nes

    # get the normalization for the likelihood function
    ll_unnorm = np.empty(len(thetas))
    for i in range(len(thetas)):
        ll_unnorm[i] = -1*objective(thetas[i])
    ll_unnorm_max = np.max(ll_unnorm)
    ll_unnorm = ll_unnorm - ll_unnorm_max
    likelihood_unnorm = np.exp(ll_unnorm)
    norm = np.trapz(likelihood_unnorm,x=Nes)
    likelihood_norm = likelihood_unnorm/norm

    error_objective_array = lambda Nes: (p_value_array(Nes, likelihood_norm, norm, ll_unnorm_max, dNes, data_list, shareddata)-p_thresh)
    
    Nemax = Nes[np.argmax(likelihood_norm)]
    Ne_lower0 = 0.9*Nemax
    Ne_upper0 = 1.1*Nemax 
    Ne_lower = optimize.fsolve(error_objective_array, Ne_lower0, factor = 0.1)[0]
    Ne_upper = optimize.fsolve(error_objective_array, Ne_upper0, factor = 0.1)[0]
    
#    plt.plot(Nes, error_objective_array(Nes))
#    plt.xscale('log')
#    plt.xlim([2*10**3, 10**4])
#    plt.show()
    
    return Ne_lower, Ne_upper
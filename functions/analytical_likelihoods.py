'''
Functions for priors, likelihoods, posteriors
Updates from previous version: 
    -Joint inference of Ne and c using MLE
    -Confidence interval estimation for Ne using profile likelihood
    
Updated 3/10/2022
'''
import pandas as pd
import numpy as np 
from scipy import optimize
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale = 1.5)

def getshareddata(c,N_samp):
    if np.isscalar(c):
        c = np.array([c]*len(N_samp))
    if isinstance(c, pd.DataFrame):
        c = c.values[0]
    c[c<=1] = 1
    days = N_samp.columns
    N_samp = N_samp.values[0]
    dd={
        'c':c,
        'day':days,
        'N_samp':N_samp}
    return dd

def log_normal(x,mu,v):
    return  -((x - mu)**2)/(2*v) - 0.5*np.log(2*np.pi*v)

def combo_mu(mu1,mu2,v1,v2):
    return (mu1*v2 + mu2*v1)/(v1 + v2)

def combo_v(v1,v2):
    return v1*v2/(v1 + v2)

def hmm_forward_flat(data):
    nu,beta,phi = data

    ll = 0
    mu = phi[0]
    v = beta[0]
    for i in range(1,len(phi)):
        v = v + nu
        ll += log_normal(phi[i], mu, v + beta[i])
        if i!=len(phi):
            mu = combo_mu(phi[i], mu, beta[i], v)
            v = combo_v(v, beta[i])
    return ll

### INFERENCE OF ONLY NE ###
    
def hmm_add_likelihoods(theta, data_list, shareddata):
    ll=0

    c = shareddata['c']
    N_samp = shareddata['N_samp']
    
    beta = c/(4*N_samp)
    nu = theta/(4)
            
    for _,data in enumerate(data_list):
        
        nobs = data['nobs'].values
        phi = np.sqrt(nobs/(N_samp))
        
        dd = (nu,beta,phi)
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

### JOINT INFERENCE OF NE AND C ###
# Can use variable T, where T is the amount time when Ne is assumed to be constant.
# However, T>=3 and must be an odd number.

def hmm_add_likelihoods_Ne_c(theta, c0, data_list, shareddata, *argv):
    ll=0
    
    N_samp = shareddata['N_samp']
    c = np.array([c0])
    for arg in argv:
        c = np.append(c,arg)
    beta = c/(4*N_samp)
    nu = theta/(4)
            
    for _,data in enumerate(data_list):
        
        nobs = data['nobs'].values
        phi = np.sqrt(nobs/(N_samp))
        
        dd = (nu,beta,phi)
        ll += hmm_forward_flat(dd)
    
    if pd.isnull(ll):
        ll=-1e7
    return ll*-1

def hmm_maxlikelihood_Ne_c(data_list, shareddata, x0, bounds_theta=(1e-11,1), bounds_c = (1.0, 100)):
    T = len(x0)-1

    def objective(theta):
        if T>1:
            args = theta[2:]
            output = hmm_add_likelihoods_Ne_c(theta[0], theta[1], data_list, shareddata, *args)
        else:
            output = hmm_add_likelihoods_Ne_c(theta[0], theta[1], data_list, shareddata)
        return output
        
    bounds = [bounds_theta, bounds_c]
    if T>1:
        for i in range(T-1):
            bounds.append(bounds_c)
    res=optimize.minimize(objective,x0,bounds=bounds)
    if res.success:
        Netau = 1/res.x[0]
        c = res.x[1:]
        return Netau, c  
    else:
        print(x0)
        print(T)
        return np.nan, np.array([np.nan]*T)

### PROFILE LIKELIHOOD FOR GETTING CONFIDENCE INTERVAL OF Ne###
        
def hmm_maxlikelihood_c(Ne, data_list, shareddata, x0, bounds_theta=(1e-11,1), bounds_c = (1.0, 100)):
    T = len(x0)
    def objective(theta):
        if T>1:
            args = theta[1:]
            output = hmm_add_likelihoods_Ne_c(1/Ne, theta[0], data_list, shareddata, *args)
        else:
            output = hmm_add_likelihoods_Ne_c(1/Ne, theta[0], data_list, shareddata)
        return output
        
    bounds = [bounds_c]
    if T>1:
        for i in range(T-1):
            bounds.append(bounds_c)
    res=optimize.minimize(objective,x0,bounds=bounds)
    if res.success:
        return res.x  
    else:
        return np.array([np.nan]*T)
        
def profile_likelihood(Ne, data_list, shareddata, x0):
    # At a given value of Ne, get the max likelihood across c's
    T = len(x0)
    c_all = hmm_maxlikelihood_c(Ne, data_list, shareddata, x0)
    if T>1:
        args = c_all[1:]
        pll = hmm_add_likelihoods_Ne_c(1/Ne, c_all[0], data_list, shareddata, *args)
    else:
        pll = hmm_add_likelihoods_Ne_c(1/Ne, c_all[0], data_list, shareddata)
    return pll
    
def p_value(Nes, likelihood_norm, norm, ll_unnorm_max, dNes, data_list, shareddata, x0):
    # Gets p value from profile likelihood
    thetas = 1/Nes
    objective = lambda theta: profile_likelihood(1/theta, data_list, shareddata, x0)
    p_value = np.empty(len(thetas))
    for i in range(len(thetas)):
        theta = thetas[i]
        ll_unnorm_scalar = -1*objective(theta)
        ll_unnorm_scalar = ll_unnorm_scalar - ll_unnorm_max
        likelihood_unnorm_scalar = np.exp(ll_unnorm_scalar)
        likelihood_norm_scalar = likelihood_unnorm_scalar/norm
        p_value[i] = np.sum((likelihood_norm_scalar>likelihood_norm)*likelihood_norm*dNes)
    return p_value
    
def confidence_interval_Ne(Netau_HMM, data_list, shareddata, x0, p_thresh = 0.05, plot = False):
        
    objective = lambda theta: profile_likelihood(1/theta, data_list, shareddata, x0)
        
    Nes_middle = np.linspace(Netau_HMM/1.5, Netau_HMM*1.5, 40)
    Nes_before = np.linspace(Netau_HMM/10, Netau_HMM/1.49, 5)
    Nes_after = np.linspace(Netau_HMM*1.51, Netau_HMM*10, 5)
    Nes = np.append(Nes_before, Nes_middle)
    Nes = np.append(Nes, Nes_after)
    
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
    error_objective_array = lambda Nes: (p_value(Nes, likelihood_norm, norm, ll_unnorm_max, dNes, data_list, shareddata, x0)-p_thresh)
    
    Nemax = Nes[np.argmax(likelihood_norm)]
    if plot:
        plt.plot(Nes, likelihood_norm, marker = '.')
        plt.show()
    Ne_lower0 = 0.9*Nemax
    Ne_upper0 = 1.1*Nemax 
    Ne_lower = optimize.fsolve(error_objective_array, Ne_lower0, factor = 0.1)[0]
    Ne_upper = optimize.fsolve(error_objective_array, Ne_upper0, factor = 0.1)[0]
    
    return Ne_lower, Ne_upper
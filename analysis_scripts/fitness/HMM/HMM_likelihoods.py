'''
9.13.20
Joao Ascensao

Functions for priors, likelihoods, posteriors


'''

import pandas as pd
import scipy as sp 
import numpy as np 
from scipy import optimize
from scipy.optimize import minimize
from scipy import interpolate
import scipy.special
import copy
import re
import random
from statsmodels.stats import multitest
import time


largenum = 1e9

# clean up all data for a gene
def getgenedata(countdf,betas):
	gene_data=betas.copy()
	r1=[]

	for _,row in gene_data.iterrows():
		d1=str(row['Day label'])
		if re.match(r'^-?\d+(?:\.\d+)?$', d1) is not None:
			d10=np.int(np.float(d1))
			d1=str(d10)
		r1.append(np.float(countdf[d1]))
	r2=copy.deepcopy(r1)
	
	#get rid of timepoints after the bc goes extinct
	if r1[-1]==0:
		x=np.array(r1)
		xa=np.flip(np.short(x==0))
		xb=np.flip(np.cumprod(xa))
		wh = np.where(xb==1)
		xb[wh[0][0]]=0
		r1 = x[xb==0][:-1]

	
	gene_data_dic = gene_data[['beta + Ne R','R']].to_dict(orient='series')
	gene_data_dic['r']=r1
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))

	r = gene_data_dic2['r']
	k = gene_data_dic2['beta + Ne R'][0:len(r)]
	R = gene_data_dic2['R'][0:len(r)]

	fls=[]
	fus=[]
	for i,ri in enumerate(r):
		if ri==0:
			continue
		rv=np.linspace(1e-2,ri*100,10**4)
		icdf = sqrt_desai_cdf(rv,ri,k[i])
		fls.append(icdf(0.01)/R[i])
		fus.append(icdf(0.99)/R[i])

	fl = np.min(fls)
	fu = np.max(fus)

	eps=1e-6
	if fl<=eps:
		fl=eps
	if fu>1:
		fu=1
	elif fu<=eps:
		fu=50/R[0]
	
	gene_data_dic2['f0range']=(fl,fu)

	return {'f0range':gene_data_dic2['f0range'], 'r':gene_data_dic2['r']}

# clean up data that is shared between genes, i.e. variance parameters and mean fitness changes
def getshareddata(betas):
	gene_data=betas.copy()
	
	gene_data['xbar']=np.zeros(len(gene_data))
	gene_data_dic = gene_data[['beta','Ne','R','gens','xbar']].to_dict(orient='series')
	
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))
	return gene_data_dic2

# stirlings approx
def gamma_stirlings(x):
	return 0.5*np.log(2*np.pi) + (x - 0.5)*np.log(x) - x + 1/(12*x)

def sqrt_transition(f,a,Ne):
	return np.sqrt( np.sqrt(a)/(2 * np.pi * (f**(3/2)) / Ne) )*np.exp( - (2*Ne)*(( np.sqrt(f) - np.sqrt(a))**2 ) )

def sqrt_desai_cdf(r,a,beta):
	b = 0.5*beta
	pdf = np.sqrt( np.sqrt(a)/(4 * b * np.pi * (r**(3/2))) )*np.exp( - (1/b)*(( np.sqrt(r) - np.sqrt(a))**2 ) )
	pdf = pdf/np.trapz(pdf,x=r)
	return interpolate.interp1d(np.cumsum(pdf)*(r[1]-r[0]), r, fill_value='extrapolate')

def NB_transition(x,mu,c):
	if c==1:
		c=1 + 1e-6

	g1=x + mu/(c-1)
	g2=mu/(c-1)
	g3=x+1
	
	b1=np.all(g1>1)
	b2=np.all(g2>1)

	if b1: #stirlings approx
		dg1 = gamma_stirlings(g1)
	else:
		dg1 = np.longdouble(sp.special.loggamma(np.float64(g1)))
	
	if b2: #stirlings approx
		dg2 = gamma_stirlings(g2)
	else:
		dg2 = np.longdouble(sp.special.loggamma(np.float64(g2)))
	dg3 = gamma_stirlings(g3)
	
	ll = dg1 - dg2 - np.log(c)*g2 + x*np.log(c - 1) - x*np.log(c) - dg3
	return np.exp(ll)

def hmm_init(s,r,shareddata,f,t,p_emission):
	R = shareddata['R'][t]
	beta = shareddata['beta'][t]

	alpha = p_emission(r,f*R,beta)
	return alpha

def hmm_recursion(s,r,shareddata,f,alpha,t,p_transition,p_emission,Ne0):

	R = shareddata['R'][t]
	beta = shareddata['beta'][t]
	Ne = Ne0[t]

	a=np.exp(s)

	alpha_new = np.zeros_like(alpha)

	for j,fj in enumerate(f):
		ai=alpha*p_transition(fj,a*f,Ne)*p_emission(r,fj*R,beta)
		alpha_new[j]=np.trapz(ai,x=f)

	return alpha_new


# returns log-likelihood
def hmm_forward_algorithm(s, data, shareddata, Ne):
	p_transition = lambda f,a,Ne: sqrt_transition(f,a,Ne)
	p_emission = lambda x,mu,c: NB_transition(x,mu,c)


	s=np.longdouble(s)
	r=data['r']

	lower,upper = data['f0range']

	


	#we need to start at a later time if the reads are initially 0
	mR=[]
	for t,ri in enumerate(r):
		mR.append(shareddata['R'][t])
		if ri!=0:
			t0=t
			break
	delta_f = 1/(2*np.max(mR))
	f=np.arange(lower,upper,delta_f)

	for t in range(t0,len(r)):
		if t==t0:
			alpha = hmm_init(s,r[t],shareddata,f,t,p_emission)
		else:
			alpha = hmm_recursion(s,r[t],shareddata,f,alpha,t,p_transition,p_emission,Ne)

	return np.log(np.trapz(alpha,x=f))


def s_maxlikelihood(data, shareddata, Ne_factor):
	Ne=shareddata['Ne']*Ne_factor

	objective = lambda s: hmm_forward_algorithm(s, data, shareddata, Ne)*-1
	res=optimize.minimize_scalar(objective)

	ll = np.float(res.fun*-1) # log likelihood
	ll0 = np.float(objective(0)*-1)

	# get pvalue and Delta AIC (positive values mean selective model is preferred)
	D_AIC = 2*(ll - ll0 - 1)

	LR=2*(ll - ll0)

	pval = 1-sp.stats.chi2.cdf(LR,1)
	if res.success:
		s=np.float(res.x)
		stderr=std_from_likelihood(s, np.float(res.fun), objective)
		return s, np.float(D_AIC), pval, stderr, True
	else:
		return np.nan, np.nan, np.nan, np.nan, False

# BH FDR correaction
def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr


def std_from_likelihood(optim_s, optim_ll, get_ll, h=1e-5):
	deriv2 = (-get_ll(optim_s + 2*h) + 16*get_ll(optim_s + h) - 30*optim_ll + 16*get_ll(optim_s - h) - get_ll(optim_s - 2*h))/(12*(h**2))
	std = 1/np.sqrt(deriv2)
	return std


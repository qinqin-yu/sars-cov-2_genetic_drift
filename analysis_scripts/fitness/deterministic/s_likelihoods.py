'''
9.13.20
Joao Ascensao

Functions for knockout fitness likelihoods & helper functions

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

	k = gene_data_dic2['beta + Ne R']
	r = gene_data_dic2['r']
	R = gene_data_dic2['R']


	f0guess = r[0]/R[0]
	fl = (r[0] - 2.75*r[0]*np.sqrt(k[0]))/R[0]
	fu = (r[0] + 2.75*r[0]*np.sqrt(k[0]))/R[0]
	if fl<=0:
		fl=1e-10
	if fu>1:
		fu=1
	elif fu<=0:
		fu=50/R[0]
	
	gene_data_dic2['f0range']=(f0guess,fl,fu)

	return {'f0range':gene_data_dic2['f0range'], 'r':gene_data_dic2['r']}

# clean up data that is shared between genes, i.e. variance parameters and mean fitness changes
def getshareddata(betas):
	gene_data=betas.copy()
	xbar=[]
	for _,row in gene_data.iterrows():
		d1=str(row['Day label'])
		if re.match(r'^-?\d+(?:\.\d+)?$', d1) is not None:
			d10=np.int(np.float(d1))
			d1=str(d10)

		

		xbar.append(0)
	gene_data['xbar']=np.array(xbar)
	gene_data_dic = gene_data[['beta + Ne R','R','gens','xbar']].to_dict(orient='series')
	
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))
	return gene_data_dic2

# stirlings approx
def gamma_stirlings(x):
	return 0.5*np.log(2*np.pi) + (x - 0.5)*np.log(x) - x + 1/(12*x)

# likelihood for a single barcode, single timepoint
def og_likelihood0(gs):
	g1,g2,g3,b1,b2,kappa,r=gs
	if b1: #stirlings approx
		dg1 = gamma_stirlings(g1)
	else:
		dg1 = np.longdouble(sp.special.loggamma(np.float64(g1)))
	
	if b2: #stirlings approx
		dg2 = gamma_stirlings(g2)
	else:
		dg2 = np.longdouble(sp.special.loggamma(np.float64(g2)))
	dg3 = gamma_stirlings(g3)
	return dg1 - dg2 - np.log(kappa)*g2 + r*np.log(kappa - 1) - r*np.log(kappa) - dg3

# add all log likelihoods across time for a single barcode, then integrate over nuisance intercept parameter
def int_likelihood(data):
	f0,R,s,xbar_v,t00,kappa,r=data

	bb=0
	c=0
	t=-1
	xbar=0
	for i,ri in enumerate(r):
		if ri==0 and bb==0:
			continue
		bb=1
		t+=1
		if ri==0:
			continue

		m=f0*R[i]*np.exp((s-xbar)*t)

		g1=ri + m/(kappa[i]-1)
		g2=m/(kappa[i]-1)
		g3=ri+1
		
		b1=np.all(g1>1)
		b2=np.all(g2>1)

		gs=(g1,g2,g3,b1,b2,kappa[i],ri)
		if c==0:
			ll=og_likelihood0(gs)
			c+=1
		else:
			ll+=og_likelihood0(gs)
		

	
	return ll

# add likelihoods over all barcodes within a gene
def add_likelihoods(s,logf0,data_list,shareddata):

	s=np.longdouble(s)
	f0 = np.longdouble(10**(logf0))
	kappa = shareddata['beta + Ne R']
	R = shareddata['R']
	xbar = shareddata['xbar']
	
	#print(shareddata,data_list)
	data = data_list[0]

	r = data['r']
	fg = np.where(np.array(r)!=0)[0][0]
	t = shareddata['gens']# - shareddata['gens'][fg]


	dd = (f0,R,s,xbar,t,kappa,r)
	lls = int_likelihood(dd)
		
	return lls


# p-value calculation
def posterior_pval(svec,ll):
	
	log_post_int= interpolate.interp1d(svec,ll,kind='cubic')
	svec1 = np.linspace(min(svec),max(svec),5*10**4)
	ll0 = log_post_int(0)
	log_post2 = log_post_int(svec1)

	lam = np.nan_to_num(ll0-log_post2)
	mm=np.max(log_post2)
	post2=np.nan_to_num(np.exp(log_post2-mm))
	post3=post2/np.sum(post2)

	pval=np.sum(post3[lam>0])


	amx=np.argmax(log_post2)
	s_hat = svec1[amx]
	if s_hat==0.3 or s_hat==-0.3:
		return np.nan, np.nan,np.nan, np.nan, False
	else:
		success=True

	pvals=[]
	for si in svec:
		ll0 = log_post_int(si)
		lam = np.nan_to_num(ll0-log_post2)
		mm=np.max(log_post2)
		post2=np.nan_to_num(np.exp(log_post2-mm))
		post3=post2/np.sum(post2)
		pvals.append(np.sum(post3[lam>0]) - 0.3174)
	

	w=np.where(np.diff(np.sign(pvals))!=0)[0]
	if len(w)!=2:
		sl=np.nan
		su=np.nan
	else:
		sl=svec[w[0]]
		su=svec[w[1]]

	return s_hat, pval, sl, su, success


# get maximum likelihood estimate of fitness
def s_maxlikelihood(data_list, shareddata, boot=False, bounds=(-1.5,0.8), boot_bracket=False):
	s=np.linspace(-0.5,0.5,10**3)
	f0=np.linspace(-6,-0.5,10**3)
	s_v, f0_v = np.meshgrid(s, f0)
	ll = add_likelihoods(s_v,f0_v,data_list,shareddata)


	prof_ll=np.max(ll,axis=0)
	#print(prof_ll)

	s_hat, pval,sl,su,success=posterior_pval(s,prof_ll)


	return s_hat, su-sl, pval, success
	

# BH FDR correaction
def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr

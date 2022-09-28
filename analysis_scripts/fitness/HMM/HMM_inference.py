'''
3.31.21
Joao Ascensao

MLE to infer strength of selection on knockouts, per biological replicate

python s_inference.py dir/count dir/out dir/technical_noise dir/outliers dir/meanfitness experiment,replicate mincount

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import time
import HMM_likelihoods as lh
from multiprocessing import Pool, TimeoutError
from functools import partial
from itertools import compress
import os
pd.set_option('mode.chained_assignment', None)


dirs = {
'A':'alpha|2021-06-01|36.5',
'B':'alpha|2021-06-01|37.5',
'C':'alpha|2021-06-20|61.5',
'D':'alpha|2021-06-20|62.5',

'E':'B-1-177|2021-02-22|694.5',
'F':'B-1-177|2021-02-22|695.5',

'G':'delta|2022-01-25|49.5+58.5',
'H':'delta|2022-01-25|49.5+59.5',
'I':'delta|2022-01-25|50.5+58.5',
'J':'delta|2022-01-25|50.5+59.5',
'K':'microreact',
}


for experiment in dirs:
	print(experiment,'HMM')
	countdir='../'+dirs[experiment]+'/data/bc_counts/'
	outdir='../'+dirs[experiment]+'/data/fitness/'
	ktdir='../'+dirs[experiment]+'/data/kt_errors/'

	replicate=0
	counts=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate))
	total_counts=pd.read_csv(countdir + '/{}{}_Rtot.csv'.format(experiment,replicate)).drop(columns='Unnamed: 0')
	betas=pd.read_csv(ktdir+ '/{}{}_betas.csv'.format(experiment,replicate))
	tv = list(total_counts.columns)


	#counts = counts[counts['gene_symbol']=='lacY']

	shareddata=lh.getshareddata(betas)

	start=time.time()

	def get_s_est(ggene):
		g,gene = ggene

		gene = gene[tv]

		if len(gene)<1:
			return None
		
		data_list=[]
		for _,row in gene.iterrows():
			d00=lh.getgenedata(row,betas)
			if len(d00['r'])>2:
				data_list.append(d00)
		if len(data_list)==0 or np.sum(data_list[0]['r'])<10:
			return None

		s_v=[]
		D_AIC_v=[]
		s_pval_v=[]
		std_v=[]
		success=True

		for Ne_factor in [1]:
			s, D_AIC, s_pval, std, success0 = lh.s_maxlikelihood(data_list[0], shareddata, Ne_factor)
			s_v.append(s)
			D_AIC_v.append(D_AIC)
			s_pval_v.append(s_pval)
			std_v.append(std)
			if not success0:
				success=False

		if success:
			return pd.DataFrame({
					'lineage':[g],
					's':s_v,
					'pval':s_pval_v,
					'D_AIC':D_AIC_v,
					'stderr':std_v,
					'Ne factor':[1]
				})
		else:
			return None
	


	batches = list(counts.groupby(by=['lineage']))
	with Pool(processes=4) as pool:
		s_data = pool.map(get_s_est, batches)

	s_data.append(None)
	cond = pd.notnull(s_data)
	s_data = list(compress(s_data, list(cond)))

	df_save=pd.concat(s_data,ignore_index=True)
	s_signi,s_pval_corr=lh.FDR_correction(df_save['pval'])
	df_save['s significant']=s_signi
	df_save['s corrected pval']=s_pval_corr
	df_save.to_csv(outdir+'/{}{}_fitness_HMM.csv'.format(experiment,replicate))
	print('Wrote '+outdir+'/{}{}_fitness.csv'.format(experiment,replicate))
	print((time.time()-start)/3600)

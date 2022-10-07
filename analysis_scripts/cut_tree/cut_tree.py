#!/usr/bin/env python
# coding: utf-8

# This code identifies which lineage each metadata sample belongs to by cutting a phylogenetic tree at a certain depth, d.
# 
# The input files are the tree data (in a newick format) and the metadata file.
# While the same metadata file is used for all the variants, the tree data is different for each variant.
# 
# By executing all the cells, the file 'create_sublineages/sublin'+metaname+'with_sublineages.csv' is generated,
# wherein each column shows which lineage each metadata sample belongs to.

# # Modules

import os
import random
from statistics import harmonic_mean
import shutil
from io import StringIO
from Bio import Phylo
import numpy as np
#import time
import copy
import csv
import pandas as pd
import math
#import seaborn as sns
import matplotlib.pyplot as plt 
from collections import Counter

from datetime import datetime
from datetime import timedelta


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    if len(node_path)>=2:
        res=node_path[-2]
    else: 
        res="root"
    return res

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = "%d_%s" % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names


def remove_nodes(tree, listofnames):
    newtree = copy.deepcopy(tree)
    for name in listofnames:
        newtree.collapse(name)
    return newtree

def print_child(tree):
    for clade in tree.find_clades():
        print("parent = ",clade.name)
        for child in clade:
            print("         child = ",child)
            
def rm(strings):
    return strings[1:-1]

def to_int(lst):
    return [int(i) for i in lst]

def flatten(t):
    return [item for sublist in t for item in sublist]


'''
Returns the set of descendants of an internal node.
(eg. treestring='((L1L,L2L,L3L)I2I,L0L)I0I' and int_node='I2I')
'''
def get_descendants(treestring, int_node):
    pos_rangle=treestring.find(int_node)-1 # position of ) next to IxI
    counter=1
    shift=1
    # Move leftward until the number of ( becomes equal to the number of). 
    while counter>0:
        if treestring[pos_rangle-shift]==')':
            counter+=1
        elif treestring[pos_rangle-shift]=='(':
            counter-=1
        shift+=1
    pos_langle=pos_rangle-shift+1

    letters_inside=treestring[pos_langle:pos_rangle+1] 
    # Remark:the 2nd argument needed to be shifted by +1 in slicing
    Lpos=[pos for pos, char in enumerate(letters_inside) if char == 'L']
    numLL = int(len(Lpos)/2)
    externals_inside=[letters_inside[Lpos[2*i]:Lpos[2*i+1]+1] for i in range(numLL)]
    return externals_inside

def create_hist(data):
    hist=Counter(data)
    list_elements=list(hist.keys())
    list_elements.sort()
    list_counts=[ hist[i] for i in list_elements]
    return [list_elements,list_counts]

def epiweek_date(epiweek):
    return (datetime.strptime('2019-12-29', "%Y-%m-%d") + timedelta(days=7*(epiweek-1))).strftime("%Y/%m/%d")

def date_epiweek(d):
    d1= datetime.strptime(d, "%Y-%m-%d")
    d2 = datetime.strptime('2019-12-29', "%Y-%m-%d")

    return (d1 -d2).days//7+1


def rename_pango(pango):
   
    pango = str(pango)
    
    if pango =='nan':
        aux ='nan'
        
        
    elif pango=='B.1.1.7'or pango[0:2]=='Q.':
        aux='alpha'
        
    elif pango=='AY.4.2' or pango[0:7]=='AY.4.2.': # remove so-called delta plus
        aux='delta_plus'
    elif pango=='B.1.617.2' or pango[0:3]=='AY.':
        aux='delta'
        
    elif pango=='B.1.1.529' or pango=='BA.1':
        aux='omicron'
        
    elif pango[0:3]=='BA.':
        aux='omicron_other'
        
    elif pango=='B.1.177' or pango[0:8]=='B.1.177.':
        aux='B-1-177'
        
    elif pango=='B.1.351' or pango[0:8]=='B.1.351.':
        aux='beta'
        
    elif pango=='P.1' or pango[0:4]=='P.1.':
        aux='gamma'
        
    elif pango=='B.1.427' or pango=='B.1.429':
        aux='epsilon'    
        
    elif pango=='B.1.525':
        aux='eta'
    
    elif pango=='B.1.526':
        aux='iota'
        
    elif pango=='B.1.617.1':
        aux='kappa'
        
    elif pango=='B.1.617.1':
        aux='kappa'
    
    elif pango=='B.1.621' or pango==' B.1.621.1':
        aux='mu'
        
    else:
        aux='others'
    return aux


from pathlib import Path


import sys
sys.setrecursionlimit(10000)


Path('create_sublineages').mkdir(parents=True, exist_ok=True)
Path('output').mkdir(parents=True, exist_ok=True)


################################################################
# # The newick files of trees should be placed in tree_dir
# # The metadata file should be placed in meta_dir
tree_dir='../../data/trees/'
meta_dir='../../data/metadata/'
metaname='cog_england_2022-01-16_'
meta_df = pd.read_csv(meta_dir+metaname+'metadata.csv',low_memory=False, index_col=0)
meta_seq=list(meta_df.sequence_name)

################################################################
# # For a given tree, extract depths of internal and leaf nodes


for treedate in ['2021-02-22','2021-06-01','2021-06-20','2022-01-25']:

    treename='cog_global_'+treedate+'_'

    print("TREE NAME: "+treename)

    tree = Phylo.read(tree_dir+treename+'tree.newick', 'newick')
    tree.rooted=True 
    # if len(tree.get_terminals( order='preorder'))<100:
    #     Phylo.draw(tree)

    #count nodes
    num_ext=len(tree.get_terminals())
    num_int=len(tree.get_nonterminals())

    print('num_ext: ',num_ext)
    print('num_int: ',num_int)

    #record sample names
    samplename=[ clade.name for clade in tree.get_terminals( order='preorder')]

    #rename external nodes & name internal nodes
    for idx, clade in enumerate(tree.get_terminals( order='preorder')):
                clade.name="L"+str(idx)+'L'
    for idx, clade in enumerate(tree.get_nonterminals( order='preorder')):
                clade.name="I"+str(idx)+'I'

    # Store labels
    ext_name=[i.name for i in tree.get_terminals( order='preorder')]
    int_name=[i.name for i in tree.get_nonterminals( order='preorder')]

    #Store branch lengths
    ext_brln=[i.branch_length for i in tree.get_terminals( order='preorder')]
    int_brln=[i.branch_length for i in tree.get_nonterminals( order='preorder')]

    #Compute&store depths
    tree_dep=tree.depths() # compute depths of all nodes
    dict_name_dep=dict() # dict: node name -> depth
    for dep in tree_dep:
        dict_name_dep[dep.name]=round(tree_dep[dep],7)
    ext_dep=[dict_name_dep[ 'L' + str(i)+'L'] for i in range(num_ext)]
    int_dep=[dict_name_dep[ 'I' + str(i)+'I'] for i in range(num_int)]

    Phylo.write(tree, 'create_sublineages/'+treename+"renamed.txt","newick")

    # Extract topological information to save the memory
    # 1.Set all branch lengths to be 0.00000
    for clade in tree.find_clades():
        clade.branch_length = None
    Phylo.write(tree, 'create_sublineages/'+treename+"aux.txt","newick")

    #2. Read the tree_aux.txt as string and remove :0.0000
    file = open('create_sublineages/'+treename+"aux.txt","r")
    tree_str= file.read().replace(':0.00000',"").replace('\n',"").replace(';',"")
    file.close()
    #3. Write the tree
    file = open('create_sublineages/'+treename+'tree_topology.newick.txt', 'w')
    file.write(tree_str)
    file.close()
    os.remove('create_sublineages/'+treename+"aux.txt")

    ext_df =  pd.DataFrame(
        {'samplename':samplename,
         'label': ext_name,
         'depth': ext_dep,
         'brlen':ext_brln
        })
    int_df =  pd.DataFrame(
        {
         'label': int_name,
         'depth': int_dep,
         'brlen':int_brln
        })

    ext_df.to_csv('create_sublineages/'+treename+'ext_df.csv')
    int_df.to_csv('create_sublineages/'+treename+'int_df.csv')

    # # Extract tree data for metadata sammples
    dict_name_label=dict()
    for idx, i in enumerate(samplename):
        dict_name_label[i]=idx

    print('metaname',metaname)
    df=pd.DataFrame()
    label_in_tree=['L'+str(dict_name_label.get(i))+'L' if dict_name_label.get(i)!=None else np.nan for i in meta_seq] 
    depth_in_tree=[ext_dep[dict_name_label.get(i)] if dict_name_label.get(i)!=None else  np.nan for i in meta_seq] 
    df['depth_in_tree'] = depth_in_tree
    df['leaf_label']=label_in_tree
    df.to_csv('create_sublineages/metadata_'+metaname+treename+'with_depth.csv',index=False)
    del df

    # #Attach metadata to ext_df
    dict_seqname_row=dict()
    for idx, i in enumerate(meta_df['sequence_name']):
        dict_seqname_row[i]=idx
    ext_df=pd.read_csv('create_sublineages/'+treename+'ext_df.csv',low_memory=False, index_col=0)

    #Attach the following lists to ext_df
    region_ext_df=[]
    lineage_ext_df=[]
    ew_ext_df=[]
    date_ext_df=[]
    for i in ext_df['samplename']:
        if i[0:8]=='England/':
            if i in dict_seqname_row:
                region_ext_df.append(meta_df['region'].iloc[dict_seqname_row[i]])  
                lineage_ext_df.append(meta_df['lineage'].iloc[dict_seqname_row[i]])  
                ew_ext_df.append( date_epiweek(meta_df['sample_date'].iloc[dict_seqname_row[i]] ) )
                date_ext_df.append( meta_df['sample_date'].iloc[dict_seqname_row[i]] ) 
            else:
                region_ext_df.append('nan') 
                lineage_ext_df.append('nan') 
                ew_ext_df.append('nan') 
                date_ext_df.append('nan') 
        else:
            region_ext_df.append('nan') 
            lineage_ext_df.append('nan') 
            ew_ext_df.append('nan')
            date_ext_df.append('nan') 

    #ext_df['region']= region_ext_df
    ext_df['lineage']= lineage_ext_df
    #ext_df['epi_week']=ew_ext_df
    #ext_df['sample_date']=date_ext_df
    ext_df =ext_df.drop(columns=['samplename'])
    ext_df.to_csv('create_sublineages/'+treename+'ext_df_with_meta.csv')

    os.remove('create_sublineages/'+treename+'ext_df.csv')
################################################################
# # Display the distribution of deppths for leaf nodes
#meta_df = pd.read_csv(meta_dir+metaname+'metadata.csv',low_memory=False, index_col=0)

for focal_variant in ['B-1-177','alpha', 'delta']:
    
    if focal_variant=='delta':
        tree_date = '2022-01-25'
    if focal_variant=='B-1-177':
        tree_date = '2021-02-22'
    if focal_variant=='alpha':
        tree_date ='2021-06-20'
  
    treename='cog_global_'+tree_date+'_'
    meta=pd.read_csv('create_sublineages/metadata_'+metaname+treename+'with_depth.csv')
    meta['epi_week']= meta_df['epi_week']
    meta['variant']= meta_df['variant']
  
    dunit=3.34e-05
    plt.figure(figsize=(6,5))
    aux=meta[(meta['depth_in_tree']>0) & (meta['variant']==focal_variant) ]


    if focal_variant=='B-1-177':
        dmin=690
        dmax=700
        ewmin=30
        ewmax=70
    if focal_variant=='alpha':
        dmin=40
        dmax=70
        ewmin=40
        ewmax=80
    if focal_variant=='delta':
        dmin=40
        dmax=80
        ewmin=60
        ewmax=110

    fs=16
    x_edges=np.linspace(ewmin,ewmax,ewmax-ewmin+1)
    y_edges=np.linspace(dmin,dmax,dmax-dmin+1)
    plt.hist2d(aux['epi_week'],aux['depth_in_tree']/dunit, bins=[x_edges,y_edges])
    plt.xlabel('Epiweek',fontsize=fs)
    plt.ylabel('Depth '+r'$/ 3.34\times 10^{-5}$',fontsize=fs)
    plt.title(focal_variant,fontsize=20)

    
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    plt.colorbar().set_label(label='Number of sequences',size=14)
    plt.savefig('output/depthdistribution_'+focal_variant+'.pdf',bbox_inches='tight')
    #plt.show()



################################################################
# # Define sublineages by cutting the tree at a specified depth

variantlist=['alpha','delta','B-1-177']
sunlin_df=pd.DataFrame()
sunlin_df['sequence_name']=meta_df['sequence_name']
del meta_df
for focal_variant in variantlist:
    print('focal_variant',focal_variant)
    
    
    if focal_variant=='delta':
        tree_date = '2022-01-25'
        dCL_rescaled_list=[50.5,58.5]
    if focal_variant=='B-1-177':
        tree_date = '2021-02-22'
        dCL_rescaled_list=[694.5]
    if focal_variant=='alpha':
        # tree_date ='2021-06-01'
        # dCL_rescaled_list=[36.5,37.5]
        tree_date ='2021-06-20'
        dCL_rescaled_list=[61.5]   
        
        
    treename='cog_global_'+tree_date+'_'
    print('create_sublineages/metadata_'+metaname+treename+'with_depth.csv', dCL_rescaled_list)
    meta = pd.read_csv('create_sublineages/metadata_'+metaname+treename+'with_depth.csv',low_memory=False)
    
    
    print('treename',treename)
    ext_df=pd.read_csv('create_sublineages/'+treename+'ext_df_with_meta.csv',low_memory=False,index_col=0 )
    int_df=pd.read_csv('create_sublineages/'+treename+'int_df.csv',low_memory=False,index_col=0)

    ext_name=list(ext_df['label'])
    int_name=list(int_df['label'])

    num_ext=len(ext_name)
    num_int=len(int_name)

    ext_depth=list(ext_df['depth'])
    int_depth=list(int_df['depth'])

    ext_brln=list(ext_df['brlen'])
    int_brln=list(int_df['brlen'])

    ext_lin=[rename_pango(i) for i in ext_df.lineage]


    dunit=3.34e-05

    for dCL_rescaled in dCL_rescaled_list:
        dCL=dunit*dCL_rescaled # Depth at which the tree is cut 
        print('dCL',dCL, ', dCL_rescaled',dCL_rescaled)

        #Edges cut by dCL. Each of such edges defines a lineage
        branch_CL=[]
        for i in range(len(int_name)):
            if int_depth[i]-int_brln[i]<dCL and dCL<int_depth[i]:
                branch_CL.append(int_name[i])

        file = open('create_sublineages/'+treename+'tree_topology.newick.txt','r')
        treestring=file.read()
        file.close()

        #Write out the descendants for each internal node at dCL 
        counts_intsubtrees=0
        counts_extsubtrees=0
        counts_intnodes=0
        counts_extnodes=0

        descendants_CL=[] # The i-th element represents the list of leaf nodes belonging to the i-th lineage
        subtree_root=[] #node 
        for idx,i in enumerate(branch_CL):

            des_anylin=get_descendants(treestring,i)
            if focal_variant != 'ALL':
                   des_focal=[des for des in des_anylin if ext_lin[int(rm(des))] == focal_variant]
            else:
                des_focal=des_anylin

            if len(des_focal)>0:
                subtree_root.append(branch_CL[idx])
                counts_intsubtrees+=1
                counts_intnodes+=len(des_focal)
                descendants_CL.append(des_focal)

        #UK-external branches cut by dCL
        for i in range(len(ext_name)):
            if ext_depth[i]-ext_brln[i]<dCL and dCL<ext_depth[i] and ext_lin[i]==focal_variant:
                # Remark: The condition of being sampled in UK is unnecessary to be imposed because if the leaf node has a lineage name, it is automatically a UK leaf. 
                counts_extsubtrees+=1
                counts_extnodes+=1
                subtree_root.append(ext_name[i])
                descendants_CL.append([ext_name[i]])
        file.close()

        print("#subtrees", counts_intsubtrees+counts_extsubtrees, 'int/ext', counts_intsubtrees, counts_extsubtrees)
        print("#nodes", counts_intnodes+counts_extnodes,'int/ext',counts_intnodes,counts_extnodes)

        dict_leaflabel_CL=dict()
        for CLidx, members in enumerate(descendants_CL):
            for i in members:
                dict_leaflabel_CL[i]=CLidx

        sublineages=[ 'SL.'+str(dict_leaflabel_CL[i]) if i in dict_leaflabel_CL else np.nan for i in meta['leaf_label']]  
        sunlin_df[focal_variant+'_sublineage'+str(dCL_rescaled)]=sublineages

sunlin_df.to_csv('create_sublineages/sublin'+metaname+'with_sublineages.csv')

################################################################
# ## For Delta, cutting the tree at the two depths
# - After the first cut, one of the lienages has a large frequency ("dominant lineage").
# - To obtain many lineages, divide the dominant lineage further by doing the 2nd cut. The 2nd cut is applied to the dominant lineage only.

focal_variant ='delta'
treename='cog_global_2022-01-25_'

sublin_df= pd.read_csv('create_sublineages/sublin'+metaname+'with_sublineages.csv',low_memory=False, index_col=0)

first_cuts=['delta_sublineage50.5']
second_cuts=['delta_sublineage58.5']

print(first_cuts)
print(second_cuts)

for c0 in first_cuts:
    for c1 in second_cuts:
        slcounter=Counter([i for i in sublin_df[c0] if str(i)[:3]=='SL.'])
        slcounter_df=pd.DataFrame()
        slcounter_df['sl']=[ i for i in slcounter]
        slcounter_df['cts']=[ slcounter[i] for i in slcounter]
        slcounter_df=slcounter_df.sort_values(by=['cts'])
        dominantSL=slcounter_df['sl'].iloc[-1] # Find the dominant lineage after the 1st cutting

        
        # For the dominant lineage, decompose further by using the 2nd cut.
        combined=[]
        aux=0
        for i in range(len(sublin_df)):
            if str(sublin_df[c0].iloc[i])!=dominantSL:
                combined.append(sublin_df[c0].iloc[i])
            elif str(sublin_df[c1].iloc[i]) !='nan':
                combined.append(dominantSL+'-'+sublin_df[c1].iloc[i]) #  Give new names to the dominant-lineage samples
                aux+=1
            else:
                combined.append(np.nan)

        #Simplify lineage names
        combined_remove=set([i for i in combined if str(i)!='nan'])
        dict_combined_renamed=dict()
        for idx,i in enumerate(combined_remove):
            dict_combined_renamed[i] = 'SL.'+str(idx)
        dict_combined_renamed[np.nan]=np.nan
        combined_renamed=[ dict_combined_renamed[i] for i in combined]
        sublin_df['delta_sublineage'+c0[16:]+'+'+c1[16:]]=combined_renamed
        
sunlin_df.to_csv('output/sublin'+metaname+'with_sublineages.csv')

################################################################
shutil.rmtree('create_sublineages')



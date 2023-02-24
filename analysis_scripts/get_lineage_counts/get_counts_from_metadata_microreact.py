#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:06:14 2022

@author: qinqinyu
"""

import pandas as pd
import numpy as np
import os
from zipfile import ZipFile

def save_counts_from_metadata_microreact(metadata_pillar2, variant_column, output_folder):
    # INPUTS:
    # metadata_pillar2: all metadata of the geographical region to be considered
    # variant_column: the heading in the metadata corresponding to assigned lineages the variant/tree date/tree cut depth of interest
    # tree column: the corresponding heading in the metadata indicating whether the sequence was in a tree from a particular date
    # output_folder: the folder to save output counts files to
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    epiweeks = np.unique(metadata_pillar2['epi_week'])
       
    # Drop sequences that were not assigned a lineage for this variant/date/depth combo
    metadata = metadata_pillar2.dropna(axis=0, subset=[variant_column])
    lineages = np.unique(metadata[variant_column])
    
    # Get and save the number of sequences in each lineage
    counts = pd.DataFrame(columns = epiweeks, index = lineages)
    for (epiweek, lineage), df in metadata.groupby(['epi_week', variant_column]):
        counts.at[lineage, epiweek] = len(df)
    counts.fillna(0, inplace = True)
    counts.reset_index(level=0, inplace=True)
    counts.rename(columns = {'index':'lineage'}, inplace = True)
    counts.to_csv(output_folder + '/counts_all_lineages.csv')
    
    # Exclude B.1.1.7 and B.1.177, and all sublineages
    to_exclude = ['B.1.1.7', 'B.1.177']
    to_exclude_with_sublineages = []
    droprows = []
    for excluded in to_exclude:
        for i, row in counts.iterrows():
            lineage = row['lineage']
            if lineage == excluded:
                to_exclude_with_sublineages.append(lineage)
                droprows.append(i)
            elif lineage.startswith(excluded+'.'):
                to_exclude_with_sublineages.append(lineage)
                droprows.append(i)
    counts_neutral = counts.drop(index = droprows)
    counts_neutral.to_csv(output_folder + '/counts_lineages.csv')

    # Get total number of sequences assigned to lineages
    total_counts_lineages = counts_neutral[epiweeks].sum().to_frame().transpose()
    total_counts_lineages.to_csv(output_folder + '/total_counts_lineages.csv')
    
    # Gotal total number of sequences of a given variant in the metadata
    # Since all sequences are assigned a lineage, it's the same as the total number of sequences assigned to lineages
    total_counts_variant_metadata = total_counts_lineages.copy()
    total_counts_variant_metadata.to_csv(output_folder + '/total_counts_variant_metadata.csv')
    
    # Get total number of sequences in the metadata
    total_counts_metadata = pd.DataFrame(columns = epiweeks, index = [0])
    total_counts_metadata[epiweeks] = 0
    
    for epiweek, df in metadata_pillar2.groupby(['epi_week']):
        total_counts_metadata.at[0, epiweek] = len(df)
    total_counts_metadata.to_csv(output_folder + '/total_counts_metadata.csv')
    
    # Get fraction of sequences in the metadata that are of this variant
    frac_variant = total_counts_variant_metadata/total_counts_metadata
    frac_variant.to_csv(output_folder + '/frac_variant.csv')

    # Get fraction of sequences in the metadata that only exclude B.1.1.7 (includes B.1.177)
    to_exclude = ['B.1.1.7']
    to_exclude_with_sublineages = []
    droprows = []
    for excluded in to_exclude:
        for i, row in counts.iterrows():
            lineage = row['lineage']
            if lineage == excluded:
                to_exclude_with_sublineages.append(lineage)
                droprows.append(i)
            elif lineage.startswith(excluded+'.'):
                to_exclude_with_sublineages.append(lineage)
                droprows.append(i)
    counts_neutral_with_b_1_177 = counts.drop(index = droprows)
    
    total_counts_lineages_with_b_1_177 = counts_neutral_with_b_1_177[epiweeks].sum().to_frame().transpose()
    total_counts_lineages_with_b_1_177.to_csv(output_folder + '/total_counts_lineages_with_b_1_177.csv')

    frac_variant_with_b_1_177 = total_counts_lineages_with_b_1_177/total_counts_metadata
    frac_variant_with_b_1_177.to_csv(output_folder + '/frac_variant_with_b_1_177.csv')
    
path_folder = '../../data/lineages/'

# Unzip metadata file
if not os.path.exists(path_folder + 'metadata_cog_england_2022-01-16_with_sublineages_with_istree.csv'):
    with ZipFile(path_folder + 'metadata_cog_england_2022-01-16_with_sublineages_with_istree.csv.zip', 'r') as zip_ref:
        zip_ref.extractall(path_folder)

# Load metadata file
metadata_original = pd.read_csv(path_folder + 'metadata_cog_england_2021_microreact_isocode_region.csv', index_col = 0)

# Filter for pillar2 (community cases)
metadata_pillar2 = metadata_original[metadata_original['pillar_2']]

variant = 'pre_B-1-177'
variant_column = 'lineage'

# GET COUNTS FOR ALL OF ENGLAND

# Create folder for output counts files
output_folder = path_folder + variant + '/microreact/is_pillar_2/England/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

save_counts_from_metadata_microreact(metadata_pillar2, variant_column, output_folder)

# GET COUNTS FOR REGIONS
for region, metadata_pillar2_region in metadata_pillar2.groupby('region'):
    variant_column_params = variant_column.split('|')
    
    # Create folder for output counts files
    output_folder = path_folder + variant + '/microreact/is_pillar_2/' + region + '/'
    
    save_counts_from_metadata_microreact(metadata_pillar2_region, variant_column, output_folder)





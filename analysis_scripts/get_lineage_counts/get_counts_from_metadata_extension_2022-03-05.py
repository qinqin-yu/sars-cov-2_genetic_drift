#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: qinqinyu
"""

import pandas as pd
import numpy as np
import os
from zipfile import ZipFile

def save_counts_from_metadata(metadata_pillar2, variant_column, output_folder):
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
    counts.to_csv(output_folder + '/counts_lineages.csv')

    # Get total number of sequences assigned to lineages
    total_counts_lineages = counts[epiweeks].sum().to_frame().transpose()
    total_counts_lineages.to_csv(output_folder + '/total_counts_lineages.csv')
    
    # Gotal total number of sequences of a given variant in the metadata
    variant_pango_lineage_names = np.unique(metadata['pango_lineage'])
    metadata_variant = metadata_pillar2[metadata_pillar2['pango_lineage'].isin(variant_pango_lineage_names)]
    
    total_counts_variant_metadata = pd.DataFrame(columns = epiweeks, index = [0])
    total_counts_variant_metadata[epiweeks] = 0
    for epiweek, df in metadata_variant.groupby(['epi_week']):
        total_counts_variant_metadata.at[0, epiweek] = len(df)
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

#     # Get the number of sequences of the variant in the tree
#     metadata_variant_tree = metadata_variant[metadata_variant[tree_column]]
    
#     total_counts_variant_tree = pd.DataFrame(columns = epiweeks, index = [0])
#     total_counts_variant_tree[epiweeks] = 0
    
#     for epiweek, df in metadata_variant_tree.groupby(['epi_week']):
#         total_counts_variant_tree.at[0, epiweek] = len(df)
#     total_counts_variant_tree.to_csv(output_folder + '/total_counts_variant_tree.csv')
   
########
     
path_folder = '../../data/lineages/'
metadata_path_folder = '../../data/metadata/'

# Unzip metadata files
if not os.path.exists(path_folder + 'sublincog_england_2022-03-25_with_sublineages_created2023Jul30.csv'):
    with ZipFile(path_folder + 'sublincog_england_2022-03-25_with_sublineages_created2023Jul30.csv.zip', 'r') as zip_ref:
        zip_ref.extractall(path_folder)

if not os.path.exists(metadata_path_folder + 'cog_england_2022-03-25_metadata.csv'):
    with ZipFile(metadata_path_folder + 'cog_england_2022-03-25_metadata.csv.zip', 'r') as zip_ref:
        zip_ref.extractall(metadata_path_folder)
    
# Load metadata files
sublineages = pd.read_csv(path_folder + 'sublincog_england_2022-03-25_with_sublineages_created2023Jul30.csv', index_col = 0)
metadata = pd.read_csv(path_folder + 'cog_england_2022-03-25_metadata.csv', index_col = 0)

metadata_original = metadata.merge(sublineages, on = 'sequence_name')

if 'lineage' in metadata_original.columns:
    metadata_original.rename(columns={'lineage':'pango_lineage'}, inplace = True)

variant_columns = ['alpha_sublineage59.5', 'B-1-177_sublineage692.5', 'delta_sublineage50.5', 'delta_sublineage50.5+58.5']

# Filter for pillar2 (community cases)
metadata_pillar2 = metadata_original[metadata_original['is_pillar_2']=='Y']

# GET COUNTS FOR ALL OF ENGLAND
# Loop through variant/tree date/tree cut depth combinations
for variant_column in variant_columns:
    variant_column_params = variant_column.split('_sublineage')
    variant = variant_column_params[0]
    # tree_date = variant_column_params[1]
    # tree_cut_depth = variant_column_params[1][11:]
    # tree_column = 'is_tree_' + tree_date
    # Create folder for output counts files
    output_folder = path_folder + variant + '/' + variant_column + '/is_pillar_2/England/'

    save_counts_from_metadata(metadata_pillar2, variant_column, output_folder)

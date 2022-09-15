#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 12:06:41 2022

@author: qinqinyu
"""
import pandas as pd

metadata_isocode = pd.read_csv('../data/lineages/metadata_cog_england_2021_microreact_isocode.csv')

iso_code_to_utla = pd.read_csv('../data/geographical/iso_code_legend_better_format_hand_annotated.csv', index_col = 0)
utla_to_region = pd.read_csv('../data/geographical/utla_to_region.csv', index_col = 0)

merged = utla_to_region.merge(iso_code_to_utla, left_on = 'UTLANM', right_on = 'Subdivision name', how = 'outer')
merged_utla_region = merged[['3166-2 code', 'RGNNM','UTLANM']]
merged_utla_region = merged_utla_region.rename(columns={'3166-2 code': 'iso_3166_code', 'RGNNM': 'region', 'UTLANM':'name'})

merged_utla_region = merged_utla_region[~pd.isnull(merged_utla_region['region'])]
merged_utla_region = merged_utla_region.sort_values('region')

metadata_isocode_region = metadata_isocode.merge(merged_utla_region, on = 'iso_3166_code', how = 'left')
metadata_isocode_region.to_csv('../data/lineages/metadata_cog_england_2021_microreact_isocode_region.csv')

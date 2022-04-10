#!/usr/bin/python

# Writen in 20 of April 2020 (Finished on 25 of April) by L. Felipe Benites
# Inspired by https://github.com/peterthorpe5/public_scripts/tree/master/Lateral_gene_transfer_prediction_tool and https://github.com/reubwn/hgt
# CHECK IF THERE IS NO STRANGE LINES IN THE FILE (either way in the calculation will show "None" instead of "NaN")
# This should be the csv expected headers "Query,Gene_Identifier,Scientific_name,Phylum,Kingdom,Superkingdom,Evalue,Identity,Bitscore"
# This script performs analyses in the Diamond output and customize to contain species names axa (eukarya,bacteria and etc as kingdoms) and categories (metazoa nonmetazoa) to calculate indexes to predict candidate sequences with horizontal gene transfer (HGT) signals #changed from "metazoan" to "target" and "nonmetazoa" to nontarget
# Usage: python Get_HGT_indexes_V5.py taxon_blast_file.csv

import pandas as pd
import numpy
import numpy as np
import sys
import csv

# import and load csv file
f = open(sys.argv[1])

data = pd.read_csv((f), delimiter=',', index_col = 0)

# insert new column with category (nontarget and target) according to data in Superkingdom column and according with desired target group (Eukaryota or Bacteria etc.) and nontarget (Viruses or Bacteria etc.)
data.loc[data.Superkingdom.str.contains("Viruses",na=False), 'category'] = 'nontarget'

data.loc[data.Superkingdom.str.contains("Eukaryota|Bacteria|Archaea",na=False), 'category'] = 'target'

# remove Unknown taxid and Symbiodinium
data = data[data.Scientific_name.str.contains("Unknown taxid|Symbiodinium") == False]

# get best Evalues
BEST_EV= pd.pivot_table(data, index= 'Query', columns='category', values='Evalue', aggfunc='min')

# replace NaN by 1 (in case there is no target or nontarget
BEST_EV_REPLACE= BEST_EV.replace(np.nan, 1, regex=True)

# get best bitscores
BEST_BIT= pd.pivot_table(data, index= 'Query', columns='category', values='Bitscore', aggfunc='max')

# replace NaN by 1 (in case there is no target or nontarget
BEST_BIT_REPLACE= BEST_BIT.replace(np.nan, 0, regex=True)

# Calculate Alien Index (AI) - AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-Metazoa) + e-200).
"""Based on AI, genes are classified as foreign (AI>30), indeterminate (0<AI<30), or target (AI<0) ... Entries with incomplete taxonomy,
    such as X-ray structures, or hits belonging to the same phylum as the query sequence (
    i.e Rotifera, filter_out_tax_id or Nematoda), were excluded from the analysis"""

BEST_EV_REPLACE["Alien index (AI)"] = (numpy.log((BEST_EV_REPLACE["target"])+1e-200) - (numpy.log((BEST_EV_REPLACE["nontarget"]+ 1e-200))))

# Calculate HGT (hU) index = (best_bitscore_nontarget - best_bitscore_target)
"""hU is a measure of how well a given sequence matches to one set of taxa (eg. Metazoa) relative to another, mutually exclusive set of taxa (eg. non-Metazoa). It uses best-hit bitscores and is defined as: (best-hit bitscore for OUTGROUP) - (best-hit bitscore for INGROUP). See Boschetti et al. 2012 for more details. hU >= 30"""

BEST_BIT_REPLACE["HGT index (hU)"] = BEST_BIT_REPLACE["nontarget"] - BEST_BIT_REPLACE["target"]

# concat two objects from AI and Hu (BEST_EV and BEST_BIT)
concat_BEST_EV_BIT = pd.concat([BEST_EV_REPLACE,BEST_BIT_REPLACE], axis=1)

# delete target and nontarget columns
drop_met_nonmet = concat_BEST_EV_BIT.drop(['target', 'nontarget'], axis=1) # more than one value needs brackets

# get best Evalue for indexed sequences (Query=query name)
Best_Evalues1= pd.pivot_table(data, index= 'Query', values='Evalue', aggfunc='min')

# get species name for indexed sequences (Query=query name)
Best_sp_name1= pd.pivot_table(data, index=['Query','Evalue'], values='Scientific_name', aggfunc='min')

# group by species name get species name (first)
Best_sp_name2= Best_sp_name1.groupby('Query').first()

# concat tables of alien index, HGT index, e-value and species name
cols_final = pd.concat( [drop_met_nonmet, Best_Evalues1,Best_sp_name2], axis=1)

# change col names 'species name'
cols_final.rename(columns={'Scientific_name':'Top species Evalue'}, inplace=True)

# print final table final - AI = TRUE if foreign (AI>30) - hU = TRUE if >=30 #
cols_final.loc[cols_final['Alien index (AI)'] < 30, 'AI > 30'] = 'False'
cols_final.loc[cols_final['Alien index (AI)'] > 30, 'AI > 30'] = 'True'

cols_final.loc[cols_final['HGT index (hU)'] < 30, 'hU >= 30'] = 'False'
cols_final.loc[cols_final['HGT index (hU)'] >= 30, 'hU >= 30'] = 'True'

# reorder cols
cols_final = cols_final[['Alien index (AI)','AI > 30','HGT index (hU)','hU >= 30','Evalue','Top species Evalue']]

# print final table raw files
cols_final.to_csv(sys.argv[1] + '_HGT_index_raw.csv')

# drop cols
drop_cols = cols_final.drop(['Alien index (AI)', 'HGT index (hU)','Evalue'], axis=1) # more than one value needs brackets

# rename cols
drop_cols.rename(columns={'AI > 30':'Alien index (AI)','hU >= 30':'HGT index (hU)'}, inplace=True)

# print final table good files
drop_cols.to_csv(sys.argv[1] + '_HGT_index_final.csv')

#!/usr/bin/python

# This script takes a diamond blast output table taxonomically annotated with a perl script, and compares the sum of the bitscores related to a target group (ie. Symbiodiniaceae+bacteria and etc.) with a nontarget group (viruses).
# Writen in the last week of May/2020 by L. Felipe Benites.
# Usage: python Get_HGT_bitscore_sum.py blast_taxon_table.csv

import pandas as pd
import numpy
import numpy as np
import sys
import csv

# import and load csv file
f = open(sys.argv[1])

data = pd.read_csv((f), delimiter=',', index_col = 0)

# remove "Unknown tax and Symbiodinium"
data = data[data.Scientific_name.str.contains("Unknown taxid|Symbiodinium") == False]

# String to be searched
search1 ="Eukaryota|Bacteria|Archaea"

data["count_tax_target"]= data["Superkingdom"].str.count(search1)

# String to be searched
search2 ="Viruses"

# count occurrence distribution of string (word) and create new column (count_tax, search=virus) based on species name
data["count_tax_nontarget"]= data["Superkingdom"].str.count(search2)

COUNT_NONTARGET= pd.pivot_table(data, index= 'Query', values='count_tax_nontarget', aggfunc='sum')

SUM_BIT_NONTARGET= pd.pivot_table(data, index= 'Query', columns='count_tax_nontarget', values='Bitscore', aggfunc='sum')

COUNT_TARGET= pd.pivot_table(data, index= 'Query', values='count_tax_target', aggfunc='sum')

concat_TARGET_NONTARGET_BITSUM = pd.concat([COUNT_NONTARGET,COUNT_TARGET,SUM_BIT_NONTARGET], axis=1)

concat_TARGET_NONTARGET_BITSUM.rename(columns={1:'Bitscore_sum_tax_nontarget',0:'Bitscore_sum_tax_target'}, inplace=True)

# reorder cols
concat_TARGET_NONTARGET_BITSUM = concat_TARGET_NONTARGET_BITSUM[['count_tax_nontarget','Bitscore_sum_tax_nontarget','count_tax_target','Bitscore_sum_tax_target']]

# Replace NaN with zeroes (0) (when thereis no nontarget or targe)
concat_TARGET_NONTARGET_BITSUM.fillna(0, inplace=True)

# Here if 'Result_nontarget' (viruses) has higher sum and distribution than 'Result_target' (symbio genome) is 'True', other is 'False'
concat_TARGET_NONTARGET_BITSUM.loc[concat_TARGET_NONTARGET_BITSUM['Bitscore_sum_tax_nontarget'] > concat_TARGET_NONTARGET_BITSUM['Bitscore_sum_tax_target'], 'Bitscore_sum_tax_nontarget > Bitscore_sum_tax_target?'] = 'True'
concat_TARGET_NONTARGET_BITSUM.loc[concat_TARGET_NONTARGET_BITSUM['Bitscore_sum_tax_nontarget'] < concat_TARGET_NONTARGET_BITSUM['Bitscore_sum_tax_target'], 'Bitscore_sum_tax_nontarget > Bitscore_sum_tax_target?'] = 'False'

# print raw table
concat_TARGET_NONTARGET_BITSUM.to_csv(sys.argv[1] + '_raw_bitscore_sum.csv')

# drop columns from row table
drop_concat_TARGET_NONTARGET_BITSUM = concat_TARGET_NONTARGET_BITSUM.drop(['count_tax_nontarget', 'count_tax_target', 'Bitscore_sum_tax_nontarget', 'Bitscore_sum_tax_target'], axis=1) # more than one value needs brackets

# rename cols - same as "Bitscore_sum_tax_nontarget > Bitscore_sum_tax_target?"
drop_concat_TARGET_NONTARGET_BITSUM.rename(columns={'Bitscore_sum_tax_nontarget > Bitscore_sum_tax_target?':'Higher_nontarget_bit_sum'}, inplace=True)

# print final files
drop_concat_TARGET_NONTARGET_BITSUM.to_csv(sys.argv[1] + '_final_bitscore_sum.csv')






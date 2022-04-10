#!/usr/bin/python

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
search ="Viruses"

# remove duplicated scientific names ocurrence
data = data.drop_duplicates('Scientific_name', keep='first')

# count occurrence distribution of string (word) and create new column (count_tax, search=virus) based on species name
data["count_tax"]= data["Superkingdom"].str.count(search)

# create pivot dataframe (df) in this case "pivot" with 3 cols Query, count (number of lines (hits) of Query), and sum (number of times of searched word - string was found in this case "virus")
pivot = data.pivot_table(index= 'Query', values= "count_tax", aggfunc={'sum','count'})

# create new col in pivot df with name "next" and get results from number of lines of seq_tag divided by 2 (half of total of lines in Scientific_name #
pivot['next'] = pivot['count']/2

# create new col in pivot df with name "Taxon dist HGT" that compares sum(number of times that word string was found) and get true and false to keep and discard (keep if TRUE and discard if FALSE)
"""KEEP: n of viruses - count tax is > (greater) than n of lines (of hits)
   DISCARD: n of viruses - count tax is > (smaller) than n of lines (of hits)"""

pivot['Taxon dist HGT (Viruses)'] = pivot['sum'] > pivot['next']

# change col names from 'count, sum, next, keep/discard' to 'n of hits', 'n of viruses','half n of hits', 'Taxon dist HGT'. This should be raw table#
pivot.columns = ['n of hits', 'n of viruses','half n of hits', 'Taxon dist HGT (Viruses)']

# print this raw table
pivot.to_csv(sys.argv[1] + '_full_HGT_dist_raw.csv')

# drop columns from row table
drop_cols = pivot.drop(['n of hits', 'n of viruses','half n of hits'], axis=1) # more than one value needs brackets

# print final files
drop_cols.to_csv(sys.argv[1] + '_full_HGT_dist_final.csv')

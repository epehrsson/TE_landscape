#! /usr/bin/env python

# state_sharing_inter_meth.py
# For TEs ever in each state, calculates the number of samples the TE is in each state
# Performs transformation directly on TE x sample average methylation matrix using pandas
# Erica Pehrsson, 2019

# Load required packages
import sys
import os
import numpy as np
import pandas as pd

# List of methylation states
states = ['Hypomethylated','Intermediate','Hypermethylated','Missing']

# Load matrix of TE average methylation per sample
print("Loading matrix")
meth = pd.read_table(sys.argv[1])

# Count the number of samples each TE is in each state
print("Counting states")
meth['Hypomethylated']=meth.iloc[:,7:44].apply(lambda x: (x < 0.3).sum(),axis=1)
meth['Intermediate']=meth.iloc[:,7:44].apply(lambda x: ((x >= 0.3) & (x <= 0.7)).sum(),axis=1)
meth['Hypermethylated']=meth.iloc[:,7:44].apply(lambda x: (x > 0.7).sum(),axis=1)
meth['Missing']=meth.iloc[:,7:44].isnull().sum(axis=1)

# Remove average methylation columns
meth.drop(meth.columns[7:44],inplace=True,axis=1)

# Load matrix of number of CpGs per TE
print("Removing TEs without CpGs")
cpgs = pd.read_table(sys.argv[2],header=None,names=['chromosome','start','stop','subfamily','class','family','strand','CpGs'])

# Merge matrices by TE to remove TEs without CpGs
meth = meth.merge(cpgs,on=['chromosome','start','stop','subfamily','family','class','strand'],how='inner')

# For TEs ever in each state, sum the number of samples the TEs are in each state 
print("State switching profiles")
inter = meth.loc[:,states].apply(lambda x: meth.loc[x > 0,states].sum()).T

# Count the number of TEs ever in each state
inter['Total'] = meth.loc[:,states].apply(lambda x: (x > 0).sum(),axis=0)

# Output: Write table to file
print("Write table")
inter.to_csv(path_or_buf=sys.argv[3],sep='\t')

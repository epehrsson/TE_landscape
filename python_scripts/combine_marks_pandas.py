#! /usr/bin/env python

# combine_marks_pandas.py 
# For each TE x sample, combine chromHMM, WGBS, DNase, H3K27ac, and RNA state, using pandas
# Takes as input metric files split by sample
# Prints combined state matrix by sample
# Erica Pehrsson, 2018

# Load required packages
import sys
import os
import numpy as np
import pandas as pd

# Load list of chromHMM samples
print("Loading samples")
samples = [line.rstrip() for line in open(sys.argv[1], 'r')]

print("Loading dataframes")

# For each Roadmap sample,
for sample in samples:
    print(sample)

    # chromHMM: Loads matrix of chromHMM state assignments per TE for that sample
    print("chromHMM")
    chromHMM = pd.read_table('chromHMM/'+sample,header=None,names=['chromosome','start','stop','subfamily','class','family','strand','Sample','chromHMM_bp','chromHMM_State','Category'])
    print(chromHMM.shape)

    # H3K27ac: Loads matrix of H3K27ac peaks overlapping the TE for that sample, if applicable
    if os.path.isfile('H3K27ac/'+sample):
        print("H3K27ac")
        H3K27ac = pd.read_table('H3K27ac/'+sample,header=None,names=['chromosome','start','stop','subfamily','class','family','strand','Sample','H3K27ac_peaks'])
        print(H3K27ac.shape)

        # Combine chromHMM and H3K27ac matrices
        print("Merging dataframes")
        chromHMM = pd.merge(chromHMM, H3K27ac, how='outer', on=['chromosome','start','stop','subfamily','class','family','strand','Sample'], sort=False)

        # Add a column for whether a TE overlaps an H3K27ac peak (True/False)
        print("Adding state columns")
        chromHMM = chromHMM.assign(H3K27ac_State=(chromHMM.H3K27ac_peaks > 0))

    # DHS: Loads matrix of DHS peaks overlapping the TE for that sample, if applicable
    if os.path.isfile('DNase/'+sample):
        print("DNase")
        DNase = pd.read_table('DNase/'+sample,header=None,names=['chromosome','start','stop','subfamily','class','family','strand','Sample','DNase_peaks'])
        print(DNase.shape)

        # Combine master matrix with DHS matrix
        print("Merging dataframes")
        chromHMM = pd.merge(chromHMM, DNase, how='outer', on=['chromosome','start','stop','subfamily','class','family','strand','Sample'], sort=False)

        # Add a column for whether a TE overlaps a DHS peak (True/False)
        print("Adding state columns")
        chromHMM = chromHMM.assign(DNase_State=(chromHMM.DNase_peaks > 0))

    # RNA: Loads matrix of average RPKM per TE and sample, if applicable
    if os.path.isfile('RNA/'+sample):
        print("RNA")
        RNA = pd.read_table('RNA/'+sample,header=None,names=['chromosome','start','stop','subfamily','family','class','strand','Sample','Expression'])
        print(RNA.shape)

        # Add a column for whether a TE is expressed RPKM > 1
        print("Adding state columns")
        RNA = RNA.assign(RNA_State=(RNA.Expression > 1))

        # Merge master matrix with RNA-seq matrix
        print("Merging dataframes")
        chromHMM = pd.merge(chromHMM, RNA, how='outer', on=['chromosome','start','stop','subfamily','class','family','strand','Sample'], sort=False)

    # WGBS: Loads matrix of average methylation and methylation state per TE and sample, if applicable
    if os.path.isfile('WGBS/'+sample):
        print("WGBS")
        WGBS = pd.read_table('WGBS/'+sample,header=None,names=['chromosome','start','stop','subfamily','class','family','strand','Sample','Methylation','WGBS_State'])
        print(WGBS.shape)

        # Merge master matrix with methylation matrix
        print("Merging dataframes")
        chromHMM = pd.merge(chromHMM, WGBS, how='outer', on=['chromosome','start','stop','subfamily','class','family','strand','Sample'], sort=False)

    # Print merged dataframe
    print(chromHMM.shape)
    chromHMM.to_csv(path_or_buf=sample,sep='\t',index=False,na_rep="NA")

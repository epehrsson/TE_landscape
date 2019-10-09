#! /usr/bin/env python

# aggregate_marks_pandas.py
# Count TE x sample in each chromHMM, WGBS, DNase, H3K27ac, and RNA state combination, using pandas
# Also counts unique TE x sample for samples with non-chromHMM annotations
# Erica Pehrsson, 2018

# Load required packages
import sys
import os
import numpy as np
import pandas as pd

# Load list of samples with chromHMM
print("Loading samples")
samples = [line.rstrip() for line in open(sys.argv[1], 'r')]

# List of states across five metrics
categories = ['chromHMM_State','H3K27ac_State','DNase_State','RNA_State','WGBS_State']
summaries = {} # Number of TEs per state combination, all metrics
unique = {} # Number of TEs per state combination, excluding chromHMM (unique TEs)

# Input: Matrix of TE x sample with combined states, by sample
for sample in samples:
    print(sample)

    # Load combined state matrix
    combined = pd.read_table(sample,header=0)
    
    # Identifies the metrics present in that sample
    keys = list(set(combined.columns.tolist()).intersection(set(categories)))

    # Group by all available metrics, counting the number of TEs with that state combination in that sample
    # Allows NA as a level
    summaries[sample] = combined[keys+['chromosome']].fillna('missing').groupby(keys, as_index=False).count()

    # Group by non-chromHMM metrics, for samples with additional metrics
    # Counts unique TEs (removes multiple chromHMM state entries per TE)
    if len(keys) > 1:
        keys.remove('chromHMM_State')
        unique[sample] = combined[['chromosome','start','stop','subfamily','class','family','strand']+keys].drop_duplicates()[keys+['chromosome']].fillna('missing').groupby(keys, as_index=False).count()

# Sums counts of state combinations across samples and writes to file
pd.concat(summaries.values(),ignore_index=True).fillna('missing').groupby(categories, as_index=False).sum().to_csv(path_or_buf="combine_marks_counts.txt",sep='\t',na_rep="NA")

# Sums counts of state combinations, unique TEs only, excluding chromHMM
categories.remove('chromHMM_State')
pd.concat(unique.values(),ignore_index=True).fillna('missing').groupby(categories, as_index=False).sum().to_csv(path_or_buf="combine_marks_counts_unique.txt",sep='\t',na_rep="NA")

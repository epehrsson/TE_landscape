#! /usr/bin/env python

#shared_meth_CpG.py
#Summarize number of CpG states per TE x sample
#Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

start = time.time()

# Number of methylation states per TE x sample
# Initialized with 0 instances per count
counts = {k: 0 for k in range(1,5)}

# Input: Output of count_CpG_state.py, table of TE x sample x number of CpGs per methylation state
inFile = open(sys.argv[1],'r')
inFile.readline() # Header
count = 1
# For each TE x sample entry, count the number of methylation states represented by any CpG
for line in inFile:
    line = line.rstrip().split('\t')  
    states = 0
    for column in range(8,12):
        if int(line[column]) > 0:
            states += 1 
    # Then add a count to the matrix of number of states
    counts[states] += 1
    count += 1
    if count % 10000000 == 0:
        print time.time() - start
inFile.close()
print time.time() - start

# Output: Number of TE x sample instances with each number of states represented
outFile = open(sys.argv[2],'w')
for k in range(1,5):
    print >> outFile, str(k)+'\t'+str(counts[k])
outFile.close()

#! /usr/bin/env python

#state_sharing_intra_lite.py
#Count the number of states per TE x sample
#Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

start = time.time()

# Input: Bed file of TEs x sample x state, total bp of overlap. Must be sorted by position and sample.
inFile = open(sys.argv[1],'r')

# Output: Table of number of states per TE x sample
outFile = open(sys.argv[2],'w')

# Dictionary with number of states per TE x sample
counts = {k: 0 for k in range(1,16)}
sample_column = int(sys.argv[3]) # Column with Roadmap sample
state_column = int(sys.argv[4]) # Column with chromHMM state

# For each TE x sample, creates a list of chromHMM states annotating the TE
line = inFile.readline().rstrip().split('\t')
TE = '\t'.join(line[0:7])
sample = line[sample_column]
state_list = [line[state_column]]
count = 1

for line in inFile:
    line = line.rstrip().split('\t')
    # Skips TEs that do not overlap the center of a 200bp chromHMM annotation window
    # (Can only be annotated with one state)
    if line[10] == "majority":
        next
    # Add additional states to the list for that TE x sample
    if ('\t'.join(line[0:7]) == TE) & (line[sample_column] == sample):
        state_list.append(line[state_column])
    # When a new TE x sample is reached,
    else:
        # Count the total number of states for that TE x instance
        counts[len(state_list)] +=1
        # Reset TE x sample
        TE = '\t'.join(line[0:7])
        sample = line[sample_column]
        state_list = [line[state_column]]
    count += 1
    if count % 10000000 == 0:
        print time.time() - start
# Last entry
counts[len(state_list)] +=1

inFile.close()
print time.time() - start

# Print the number of states per TE x sample
for k in range(1,16):
    print >> outFile, str(k)+'\t'+str(counts[k])
outFile.close()

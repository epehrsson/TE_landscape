#! /usr/bin/env python

#state_sharing_intra_max.py
#Calculates the maximum number of states each TE is annotated with in the same sample
#Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

start = time.time()

# Input: Bed file of TEs x sample x state, total bp overlap. Must be sorted by position and sample
inFile = open(sys.argv[1],'r')

# Output: File listing each TE (row) and the maximum number of chromHMM states per sample
outFile = open(sys.argv[2],'w')

sample_column = int(sys.argv[3]) # Column listing Roadmap sample

# For each TE, create a dictionary with the number of chromHMM states the TE is annotated with by sample
line = inFile.readline().rstrip().split('\t')
TE = '\t'.join(line[0:7])
sample = line[sample_column]
state_count = {sample:1}
count = 1

for line in inFile:
    line = line.rstrip().split('\t')
    # For the same TE and sample, increase the number of states 
    if ('\t'.join(line[0:7]) == TE) & (line[sample_column] == sample):
        state_count[sample] +=1 
    # For the same TE, new sample, add an entry to the sample-state dictionary for that TE
    elif ('\t'.join(line[0:7]) == TE):
        sample = line[sample_column]
        state_count[sample] = 1
    # When a new TE is reached, print out the previous TE with the maximum value of states per sample
    else:
        outFile.write(TE+'\t'+str(max(state_count.values()))+'\n')
        # Reset the TE and sample dictionary
        TE = '\t'.join(line[0:7])
        sample = line[sample_column]
        state_count = {sample:1}
    count += 1
    if count % 10000000 == 0:
        print time.time() - start

# Write last entry
outFile.write(TE+'\t'+str(max(state_count.values()))+'\n')

inFile.close()
outFile.close()

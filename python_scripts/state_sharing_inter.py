#! /usr/bin/env python

#state_sharing_inter.py
#Calculate the frequency with which different states annotate the same TE in different samples
#Erica Pehrsson, 2016; updated 2019

# Load required packages
import os
import sys
import itertools
import time

start = time.time()

# Input: Bedfile of TEs x sample x state, total bp of overlap. Must be sorted by position and sample. 
inFile = open(sys.argv[1],'r')

# List of chromHMM states
states = [line.rstrip() for line in open(sys.argv[2],'r')] 

# Initialize a matrix of chromHMM states (row and column)
freqs = {stateA: {stateB: 0 for stateB in states} for stateA in states}

# Output: Matrix of chromHMM states
outFile = open(sys.argv[3],'w')

sample_column = int(sys.argv[4]) # Column with Roadmap sample
state_column = int(sys.argv[5]) # Column with chromHMM state
width = int(sys.argv[6]) # Width of the feature

# For each TE, create a list of states with which the TE is annotated in each sample
line = inFile.readline().rstrip().split('\t')
TE = '\t'.join(line[0:width])
sample = line[sample_column]
state_list_TE = [] # List for TE, each sample
state_list_sample = [line[state_column]] # List by sample
count = 1

for line in inFile:
    line = line.rstrip().split('\t')
    # For same TE and sample, append the list of states for that sample
    if ('\t'.join(line[0:width]) == TE) & (line[sample_column] == sample):
        state_list_sample.append(line[state_column])
    # When a new sample is reached, add the sample list to the TE list and begin a new sample list
    elif ('\t'.join(line[0:width]) == TE) & (line[sample_column] != sample):
        state_list_TE.append(state_list_sample)
        sample = line[sample_column]
        state_list_sample = [line[state_column]]
    # When a new TE is reached,
    else:
        state_list_TE.append(state_list_sample)

        # Find each state the TE is ever annotated with
        TE_states = list(set(itertools.chain(*state_list_TE)))

        # For each state the TE is ever annotated with (row)
        # Count the number of samples the TE is in each other state
        # And add a count for those states (column) for each row state
        for state in TE_states:
            for m in range(0,len(state_list_TE)):
                for n in range(0,len(state_list_TE[m])):
                    freqs[state][state_list_TE[m][n]] +=1

        # Reset TE and sample
        TE = '\t'.join(line[0:width])
        sample = line[sample_column]
        state_list_sample = [line[state_column]]
        state_list_TE = []
    count += 1
    if count % 10000000 == 0:
        print(time.time() - start)

# Last entry
TE_states = list(set(itertools.chain(*state_list_TE)))
for state in TE_states:
    for m in range(0,len(state_list_TE)):
        for n in range(0,len(state_list_TE[m])):
            freqs[state][state_list_TE[m][n]] +=1

inFile.close()
print(time.time() - start)

# Write out a matrix of chromHMM states (row: TEs ever in each state, column: number of samples in state)
outFile.write('\t'+'\t'.join(states)+'\n')
for i in states:
    outFile.write(i)
    for j in states:
        outFile.write('\t'+str(freqs[i][j]))
    outFile.write('\n')
outFile.close()

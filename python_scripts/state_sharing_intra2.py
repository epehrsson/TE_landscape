#! /usr/bin/env python

#state_sharing_intra2.py
#Count the number of TE instances annotated with one chromHMM state also annotated with a second chromHMM state in the same sample
#Erica Pehrsson, 2016

# Load required packages
import os
import sys
import itertools
import time

start = time.time()

# Input: Bed file of TEs x sample x state, total bp of overlap. Must be sorted by position and sample. 
inFile = open(sys.argv[1],'r') 

states = [line.rstrip() for line in open(sys.argv[2],'r')] #Load list of chromHMM states

# Initialize a matrix of chromHMM state combinations (row and column)
freqs = {stateA: {stateB: 0 for stateB in states} for stateA in states}

# Count number of TE x sample instances annotated with one one state, by state
single_counts = {state: 0 for state in states}

# Count the number of states per TE instance
counts = {k: 0 for k in range(1,16)}

sample_column = int(sys.argv[4]) # Column with Roadmap sample
state_column = int(sys.argv[5]) # Column with chromHMM state

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
        # For TE x sample annotated with a single state, record the state
        if len(state_list) == 1:
            single_counts[state_list[0]] += 1
        # For all pairwise combinations of states, add a count to the state matrix
        # Including an identity count (diagonal)
        # Restricted to the upper right corner
        for i in range(0,len(state_list)):
             for j in range(i,len(state_list)):
                 if states.index(state_list[i]) < states.index(state_list[j]):
                     freqs[state_list[i]][state_list[j]] +=1
                 else:
                     freqs[state_list[j]][state_list[i]] +=1
        # Reset TE x sample
        TE = '\t'.join(line[0:7])
        sample = line[sample_column]
        state_list = [line[state_column]]
    count += 1
    if count % 10000000 == 0:
        print time.time() - start

# Last entry
counts[len(state_list)] +=1
if len(state_list) == 1:
    single_counts[state_list[0]] += 1
for i in range(0,len(state_list)):
    for j in range(i,len(state_list)):
        freqs[state_list[i]][state_list[j]] +=1

inFile.close()
print time.time() - start

# Outfile: Matrix of chromHMM states (row and column) with counts of TE instances with each pairwise state combination
outFile = open(sys.argv[3],'w')
print >> outFile, '\t'+'\t'.join(states)
for i in states:
    outFile.write(i)
    for j in states:
        outFile.write('\t'+str(freqs[i][j])) #Print number of TEs x samples where pair is observed 
    outFile.write('\n')
outFile.close()

# Print the number of states per TE x sample
print "Number of states per TE x sample:"
for k in range(1,16):
    print str(k)+'\t'+str(counts[k])

# Print the number of TEs x sample annotated with a single state
print "Number of instances with single state:"
for state in states:
    print state+'\t'+str(single_counts[state])

#! /usr/bin/env python

#pick_majority.py
#Adapted from identify_no_summit.py. For TEs that do not overlap the center of a 200bp chromHMM annotation window
#Identify the state that overlaps the majority of the TE
#Erica Pehrsson, 2019

# Load required packages
import os
import sys

# Input: Bedfile of TEs x state, total bp of overlap, for a single sample. Must be sorted by position. 
# Prefiltered to those that do not overlap a 200bp window center
inFile = open(sys.argv[1],'r') 

# Output: For each TE, the state assignment for that sample
outFile = open(sys.argv[2],'w')

state_column = int(sys.argv[3]) # Column with chromHMM state
bp_column = int(sys.argv[4]) # Column with length of overlap

# For each TE, creates a dictionary with the number of bp the TE is in the state
entries = {}
line = inFile.readline().rstrip().split('\t')
TE = '\t'.join(line[0:7])
entries[line[state_column]] = line[bp_column]

for line in inFile:
    line = line.rstrip().split('\t')
    # Add states to the dictionary for the previous TE
    if '\t'.join(line[0:7]) == TE:
        entries[line[state_column]] = line[bp_column]
    # When the file reaches a new TE, identify the state covering the majority of the previous TE (largest bp)
    else:
        max_state = max(entries, key=lambda k: entries[k])
        # Write out the TE, the state, and the bp overlap
        outFile.write(TE+'\t'+max_state+'\t'+str(entries[max_state])+'\n')
        # Reset the TE
        entries = {}
        TE = '\t'.join(line[0:7])
        entries[line[state_column]] = line[bp_column]

# Last entry
max_state = max(entries, key=lambda k: entries[k])
outFile.write(TE+'\t'+max_state+'\t'+str(entries[max_state])+'\n')

inFile.close()
outFile.close()

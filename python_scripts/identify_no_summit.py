#! /usr/bin/env python

#identify_no_summit.py
#For TEs that do not overlap the center of a 200bp chromHMM annotation window,
#Identify the state the covers the majority of the TE
#Erica Pehrsson, 2018

# Load required packages
import os
import sys

# Input: Bed file of TEs x state, total bp of overlap, for a single sample. MUST be sorted by position.
inFile = open(sys.argv[1],'r')

# Load list of TEs overlapping 200bp windows in that sample
TEs = set(['\t'.join(line.rstrip().split('\t')[0:7]) for line in open(sys.argv[2], 'r')])

# Output: For TEs that do not overlap 200bp window centers, the TE and its state assignment for that sample
outFile = open(sys.argv[3],'w')

state_column = int(sys.argv[4]) # Column with chromHMM state
bp_column = int(sys.argv[5]) # Column with length of overlap between TE and chromHMM state
print len(TEs)

# For each TE, creates a dictionary with the number of bp the TE is in the state
entries = {}
line = inFile.readline().rstrip().split('\t')
TE = '\t'.join(line[0:7])
entries[line[state_column]] = line[bp_column]
count = 1

for line in inFile:
    line = line.rstrip().split('\t')
    # Add states to the dictionary for the previous TE
    if '\t'.join(line[0:7]) == TE:
        entries[line[state_column]] = line[bp_column]
    # When the file reaches a new TE,
    else:
        # If the previous TE does not overlap the center of a 200bp annotation window,
        if TE not in TEs:
            # Identify the state covering the majority of the TE (largest bp)
            max_state = max(entries, key=lambda k: entries[k])
            # Write out the TE, the state, and the bp overlap
            outFile.write(TE+'\t'+max_state+'\t'+str(entries[max_state])+'\n')
            # For TEs in 2 states, write out the discarded state
            if len(entries) > 1:
                print TE, entries, max_state
        # Reset the TE
        entries = {}
        TE = '\t'.join(line[0:7])
        entries[line[state_column]] = line[bp_column]
        count += 1
    if count % 100000 == 0:
        print count

# Last entry
if TE in TEs:
    next
else:
    max_state = max(entries, key=lambda k: entries[k])
    outFile.write(TE+'\t'+max_state+'\t'+str(entries[max_state])+'\n')
    if len(entries) > 1:
        print TE, entries, max_state

inFile.close()
outFile.close()

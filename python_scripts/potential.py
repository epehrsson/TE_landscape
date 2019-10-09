#! /usr/bin/env python

#potential.py 
#Finds the number of samples each TE/promoter is annotated with each chromHMM state
#Erica Pehrsson, 2016

# Load required packages
import os
import sys
import itertools
import time

# Input: Matrix of TE/promoter x sample x state, total bp overlap. Must be sorted by position. 
inFile = open(sys.argv[1],'r')

# Load list of features under consideration
TEs = [line.rstrip() for line in open(sys.argv[2], 'r')]
# Width of the feature
width = len(TEs[0].split('\t'))

states = [line.rstrip() for line in open(sys.argv[3],'r')] #Load list of chromHMM states
samples = [line.rstrip() for line in open(sys.argv[4],'r')] #Load list of samples under consideration

# Output: Matrix of feature (row) by chromHMM states (column), listing the number of samples each feature is annotated with each state
outFile = open(sys.argv[5],'w')
threshold = float(sys.argv[6]) # Minimum length of overlap
sample_column = int(sys.argv[7]) # Column listing Roadmap sample
state_column = int(sys.argv[8]) # Column listing chromHMM state
overlap_column = int(sys.argv[9]) # Column listing bp overlap

start = time.time()

#Initialize matrix where all TE-state tuples have 0 count
matrix = dict.fromkeys(list(itertools.product(TEs,states)),0) 
print time.time() - start
count = 0
# For each entry,
for line in inFile:
    line = line.rstrip().split('\t')
    # If sample is in list of samples under consideration and length of overlap is greater than threshold, add a count for that state
    if (line[sample_column] in samples) & (float(line[overlap_column])/(float(line[2])-float(line[1])) > threshold): 
        matrix[('\t'.join(line[0:width]),line[state_column])] +=1 
    count +=1
    if count % 1000000 == 0:
        print time.time() - start
inFile.close()

# Print output column names
print >> outFile, 'chromosome\tstart\tstop\tsubfamily\tclass\tfamily\tstrand\t'+'\t'.join(states)
# For each feature (row), print the number of samples the feature is in each state
for TE in TEs:
    outFile.write(TE)
    for state in states:
        outFile.write('\t'+str(int(matrix[(TE,state)])))
    outFile.write('\n')
outFile.close()

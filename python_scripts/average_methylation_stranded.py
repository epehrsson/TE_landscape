#! /usr/bin/env python

# average_methylation_stranded.py
# From TE-CpG fractional methylation level intersection, find average methylation level across TE
# Uses as input mouseENCODE fractional methylation files, which include read coverage
# Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

# Input: Bedtools intersection of TEs with CpG methylation levels (two columns for each sample)
inFile = open(sys.argv[1],'r')

# Load list of TEs
TEs = [line.rstrip() for line in open(sys.argv[2], 'r')]

# Load list of samples
samples = [line.rstrip() for line in open(sys.argv[3], 'r')]
width = len(samples) # Number of samples

start = time.time()

# For each TE, initiate an array for each sample
matrix = {TE: [[] for i in range(width)] for TE in list(TEs)}
print time.time() - start

# For each entry, append the methylation value for each sample column to its array for that TE
# Only adding the value if read coverage is >3
count = 0
for line in inFile:
    line = line.rstrip().split('\t')
    [matrix['\t'.join(line[0:7])][i].append(float(line[(i*2)+8])) for i in range(0,width) if float(line[(i*2)+7]) > 3]
    count +=1
    if count % 1000000 == 0:
        print time.time() - start
inFile.close()
print matrix

# Output: For each TE (row), write the average methylation level in each sample (column)
outFile = open(sys.argv[4],'w')
print >> outFile, 'chromosome\tstart\tstop\tsubfamily\tclass\tfamily\tstrand\t'+'\t'.join(samples)
for TE in matrix.keys():
    outFile.write(TE)
    for i in range(0,width):	
        # If the TE overlaps any CpG with sufficient read coverage in that sample,
        # Output the average methylation level over all CpGs overlapping the TE for that sample
        if len(matrix[TE][i]) > 0: 
            outFile.write('\t'+str(sum(matrix[TE][i])/float(len(matrix[TE][i]))))
        # Or add a missing value for that sample
        else:
            outFile.write('\tNA')
    outFile.write('\n')
outFile.close()

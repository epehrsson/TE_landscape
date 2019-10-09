#! /usr/bin/env python

#average_methylation_promoter.py
#From promoter-CpG fractional methylation level intersection, find the average methylation level across each promoter
#Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

# Input: Bedtools intersection between RefSeq promoters and CpG methylation levels (column for each sample)
inFile = open(sys.argv[1],'r')

# Load list of promoters
promoters = [line.rstrip() for line in open(sys.argv[2], 'r')] 

# Load list of samples
samples = [line.rstrip() for line in open(sys.argv[3], 'r')]
width = len(samples) # Number of samles

start = time.time()

# For each promoter, initialize an empty array for each sample
matrix = {promoter: [[] for i in range(width)] for promoter in list(promoters)}
print time.time() - start

# For each entry, append the methylation value for each sample column to its array for that promoter
# Only adding the value if it is not -1 (read coverage 3 reads or less)
count = 0
for line in inFile:
    line = line.rstrip().split('\t')
    if line[7] != '.': #Skips promoters that do not overlap a CpG (-wao flag)
        [matrix['\t'.join(line[0:4])][i].append(float(line[i+7])) for i in range(0,width) if float(line[i+7]) >= 0]
    count +=1
    if count % 1000000 == 0:
        print time.time() - start
inFile.close()

# Output: For each promoter (row), the average methylation level in each sample (column)
outFile = open(sys.argv[4],'w')
print >> outFile, 'chromosome\tstart\tstop\tstrand\t'+'\t'.join(samples)
for promoter in matrix.keys():
    outFile.write(promoter)
    for i in range(0,width):	
        # If the promoter overlaps any CpG with sufficient read coverage in that sample,
        # Output the average methylation level over all CpGs overlapping the promoter for that sample
        if len(matrix[promoter][i]) > 0: 
            outFile.write('\t'+str(sum(matrix[promoter][i])/float(len(matrix[promoter][i]))))
        else:
        # Or add a missing value for that sample
            outFile.write('\tNA')
    outFile.write('\n')
outFile.close()

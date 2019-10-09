#! /usr/bin/env python

# calculate_bin_average_meth.py
# Given a bed file of intersection between bins and methylation, will calculated the weighted average for each bin
# Averages methylation level over each bin instance first, to account for different numbers of overlapping CpGs
# Erica Pehrsson, 2016

# Load required packages
import os
import sys
import math

# Input: Intersection between bins and fractional CpG methylation level
inFile = open(sys.argv[1],'r')

# Output: Average methylation per bin
outFile = open(sys.argv[2],'w')

bins = {} # Dictionary of bin (keys) and total methylation
bin_length = {} # Dictionary of bin (keys) and length of intersections

# Reading in the first line
line = inFile.readline().rstrip().split('\t')
bin = line[4] # Bin
position = '\t'.join(line[0:3]) # Location of bin
level = float(line[8]) # Methylation level
CpGs = 1 # Number of CpGs per bin

# For each intersection, 
for line in inFile:
    line = line.rstrip().split('\t')
    # If the bin is the same, add the methylation level to the total and add an additional CpG
    if '\t'.join(line[0:3]) == position:    
        level = level + float(line[8])
        CpGs += 1
    # When a new bin is reached, calculate the average methylation level across all CpGs overlapping that bin
    # Then add the average methylation to the total for that bin and count the number of bins
    else:    
        average = level/float(CpGs)
        if bin not in bins.keys():    
            bins[bin] = average
            bin_length[bin] = 1
        else:
            bins[bin] += average
            bin_length[bin] += 1
        # Reset the bin
        bin = line[4]
        position = '\t'.join(line[0:3])
        level = float(line[8])
        CpGs = 1
inFile.close()
# Finish last instance
average = level/float(CpGs)
if bin not in bins.keys():
    bins[bin] = average
    bin_length[bin] = 1
else:
    bins[bin] += average
    bin_length[bin] += 1

# For each bin, write the bin, average methylation over all instances, and number of contributing instances
for bin in bins.keys():
   outFile.write(bin+'\t'+str(bins[bin]/float(bin_length[bin]))+'\t'+str(bin_length[bin])+'\n')
outFile.close()

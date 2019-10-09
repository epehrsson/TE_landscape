#! /usr/bin/env python

# calculate_bin_average_all.py
# Given a file of average signal value per bin for many samples, will calculated the weighted average for each bin
# Erica Pehrsson, 2016

# Load required packages
import os
import sys
import math

# Input: Average signal per bin by sample
inFile = open(sys.argv[1],'r')

# Output: Averaged signal per bin across many samples
outFile = open(sys.argv[2],'w')

bins = {} # Dictionary of bins (keys) and total signal
bin_length = {} # Dictionary of bins (keys) and length of coverage

# For each bin, calculate the total coverage (integral of signal strenth) and total length of coverage
for line in inFile:
    line = line.rstrip().split('\t')
    bin = line[0] # Bin number
    if bin not in bins.keys():
        bins[bin] = float(line[1])*float(line[2])
        bin_length[bin] = int(line[2])
    else:
        bins[bin] += float(line[1])*float(line[2])
        bin_length[bin] += int(line[2])
inFile.close()

# For each bin, output the bin number, average coverage across samples, and length of coverage
for r in range(1,len(bins.keys())+1):
   outFile.write("Bin_"+str(r)+'\t'+str(bins["Bin_"+str(r)]/float(bin_length["Bin_"+str(r)]))+'\t'+str(bin_length["Bin_"+str(r)])+'\n')
outFile.close()

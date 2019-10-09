#! /usr/bin/env python

# calculate_bin_average.py
# Given a bed file of intersection between bins and signal value, will calculated the weighted average for each bin
# Erica Pehrsson, 2016

# Load required packages
import os
import sys
import math

# Input: Intersection between bins and signal (fold change enrichment)
# Because signal is a bed file, different signal strengths are represented as separate intersections of varying length
inFile = open(sys.argv[1],'r')

# Output: Average signal per bin
outFile = open(sys.argv[2],'w')

bins = {} # Dictionary of bin (keys) and total signal coverage
bin_length = {} # Dictionary of bin (keys) and length of signal coverage

# For each intersection, add the total coverage (integral) and length of coverage to dictionaries for that bin
for line in inFile:
    line = line.rstrip().split('\t')
    bin = line[4]
    if bin not in bins.keys():
        bins[bin] = float(line[8])*float(line[9])
        bin_length[bin] = int(line[9])
    else:
        bins[bin] += float(line[8])*float(line[9])
        bin_length[bin] += int(line[9])
inFile.close()

# For each bin, write out bin number, average signal (total coverage/length of coverage), and length of coverage
for r in range(1,len(bins.keys())+1):
   outFile.write("Bin_"+str(r)+'\t'+str(bins["Bin_"+str(r)]/float(bin_length["Bin_"+str(r)]))+'\t'+str(bin_length["Bin_"+str(r)])+'\n')
outFile.close()

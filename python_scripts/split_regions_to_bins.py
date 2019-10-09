#! /usr/bin/env python

# split_regions_to_bins.py
# Given a bed file, will split each region into x-sized bins
# If the length of the region is not evenly divided by the bin size, the last bin will be shorter
# Erica Pehrsson, 2016

# Load required packages
import os
import sys
import math

# Input: Bedfile of regions
inFile = open(sys.argv[1],'r')

# Output: Regions split into bins
outFile = open(sys.argv[2],'w')

# Bin size
bin_size = int(sys.argv[3])

# For each region,
for line in inFile:
    line = line.rstrip().split('\t')
    # Find the total length of the region and number of bins that will be created
    bins = math.ceil((float(line[2])-float(line[1]))/bin_size)
    start = int(line[1]) # Start position
    b = 0 # Bin number
    # For each bin,
    while (b < bins):
        # For the last bin, output chromosome, bin start, region stop position, and other fields
        if (b == bins-1) and (int(line[2]) < start+bin_size):
            outFile.write(line[0]+'\t'+str(start)+'\t'+line[2]+'\t'+'\t'.join(line[3:])+'\tBin_'+str(b+1)+'\n')
            break
        # For all other bins, output chromosome, bin start, bin stop, and other fields
        else:
            outFile.write(line[0]+'\t'+str(start)+'\t'+str(start+bin_size)+'\t'+'\t'.join(line[3:])+'\tBin_'+str(b+1)+'\n')
            # Update the start position of the next bin and the bin number
            start = start+bin_size
            b = b+1
inFile.close()
outFile.close()

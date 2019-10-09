#! /usr/bin/env python

#filter_summit.py 
#Filter intersections of individual TEs and chromHMM to those overlapping the center of the 200bp window
#Erica Pehrsson, 2018

# Load required packages
import os
import sys
import itertools
import time

# Input: Intersection of TEs and chromHMM annotations, using bedtools intersect with -wo flag
inFile = open(sys.argv[1],'r')

feature_start = int(sys.argv[2]) # Column with start coordinate of TE
peak_start = int(sys.argv[3]) # Column with start coordinate of chromHMM segment

# Output: Same columns as input
outFile = open(sys.argv[4],'w')

metric = sys.argv[5] # Epigenetic technique, chromHMM only

# By line:
for line in inFile:
    line = line.rstrip().split('\t')
    if metric != "chromHMM":
        next
    # For chromHMM, determines whether the center of a 200bp window overlaps the TE
    # Annotation blocks are multiples of 200bp
    else:
        start = int(line[peak_start])
        stop = int(line[peak_start + 1])
        while (start+200 <= stop):
            summit = start + 99.5
            if (summit >= int(line[feature_start])) and (summit < int(line[feature_start + 1])):
                print >> outFile, '\t'.join(line)
                break
            start += 200

inFile.close()
outFile.close()

#! /usr/bin/env python

#subset_file.py
#Subset a file based on line numbers
#Erica Pehrsson, 2016

# Load required packages
import os
import sys

# Input: file to subset
inFile = open(sys.argv[1],'r')

# Load line numbers
numbers = [int(line.rstrip()) for line in open(sys.argv[2], 'r')]

# Output: Subset file
outFile = open(sys.argv[3],'w')

# Print line if line number is in list of line numbers to subset
count = 0
for line in inFile:
    count += 1
    if count in numbers:
        outFile.write(line)
inFile.close()
outFile.close()

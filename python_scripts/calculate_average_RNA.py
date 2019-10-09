#! /usr/bin/env python

#calculate_average_RNA.py
#Given a bed file of intersection between TEs and RNA-seq read counts, calculates the average read coverage for each TE
#Erica Pehrsson, 2016

# Load required packages
import os
import sys
import math

# Input: Bedtools intersection between TEs/exons and RNA-seq read coverage (strand agnostic) 
inFile = open(sys.argv[1],'r')

# Load list of TEs/exons to be profiled
TEs = [line.rstrip() for line in open(sys.argv[2], 'r')]
width = len(TEs[1].split('\t')) # Width of the feature

# For each TE/exon, create a dictionary with the length of the feature, total read coverage, and total length of intersection
matrix = {TE: [int(TE.split('\t')[2])-int(TE.split('\t')[1]),0,0] for TE in list(TEs)}

# For each entry, add total read coverage and length of overlap to dictionary for the feature
for line in inFile:
    line = line.rstrip().split('\t')
    matrix['\t'.join(line[0:width])][1] += float(line[(width+4)])*float(line[(width+5)])
    matrix['\t'.join(line[0:width])][2] += int(line[(width+5)])
inFile.close()

# Output: For each feature (row), total read coverage, total length of overlap, and average read coverage
outFile = open(sys.argv[3],'w')
for TE in TEs:
   outFile.write(TE+'\t'+str(matrix[TE][1])+'\t'+str(matrix[TE][2])+'\t'+str(matrix[TE][1]/float(matrix[TE][0]))+'\n')
outFile.close()

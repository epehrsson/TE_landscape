#! /usr/bin/env python

#count_CpG_state.py 
#From TE-CpG fractional methylation level intersection, find number of CpGs in each state per TE x sample
#Erica Pehrsson, 2017

# Load required packages
import os
import sys
import itertools
import time

# Input: Intersection of TEs with CpG methylation levels
inFile = open(sys.argv[1],'r')

# Load list of TEs
TEs = ['\t'.join(line.rstrip().split('\t')[0:7]) for line in open(sys.argv[2], 'r')]

# Load list of samples
samples = [line.rstrip() for line in open(sys.argv[3], 'r')]

start = time.time()

# For each TE x sample, initialize an array with 0 for each methylation state
matrix = {TE: {sample: [0,0,0,0] for sample in list(samples)} for TE in list(TEs)}
print time.time() - start

# For each TE x CpG intersection, for each sample, add a value to the appropriate methylation state
count = 0
for line in inFile:
    line = line.rstrip().split('\t')
    TE = '\t'.join(line[0:7])
    for i in range(0,len(samples)):
        sample = samples[i]
        sample_CpG = float(line[i+10]) #Skip TE and CpG coordinates
        if sample_CpG == -1: # Missing methylation data
            matrix[TE][sample][3] += 1
        elif sample_CpG < 0.3: # Hypomethylated
            matrix[TE][sample][0] += 1
        elif sample_CpG >= 0.3 and sample_CpG <= 0.7: # Intermediately methylated
            matrix[TE][sample][1] += 1
        elif sample_CpG > 0.7: # Hypermethylated
            matrix[TE][sample][2] += 1
    count +=1
    if count % 1000000 == 0:
        print time.time() - start
inFile.close()

# Output: For each TE x sample, the number of CpGs per methylation state
outFile = open(sys.argv[4],'w')
print >> outFile, 'chromosome\tstart\tstop\tsubfamily\tclass\tfamily\tstrand\tSample\tHypomethylated\tIntermediate\tHypermethylated\tMissing'
TE_samples = [(TE,sample) for TE in TEs for sample in samples]
for TE_sample in TE_samples:
    TE = TE_sample[0]
    sample = TE_sample[1]
    print >> outFile, TE+'\t'+sample+'\t'+str(matrix[TE][sample][0])+'\t'+str(matrix[TE][sample][1])+'\t'+str(matrix[TE][sample][2])+'\t'+str(matrix[TE][sample][3])
outFile.close()

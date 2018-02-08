#!/bin/bash

#****************************************#
#   Master code for TE landscape paper   #
#      Written by Erica Pehrsson         #
#          January 10, 2018              #
#****************************************#

# Sample list locations
chromHMM_samples=sample_lists/mnemonics.txt
WGBS_samples=samples_lists/WGBS_samples.txt
DNase_samples=sample_lists/DNase_samples.txt
H3K27ac_samples=sample_lists/H3K27ac_samples.txt

# Raw data locations 
chromHMM_raw=raw_data/chromHMM
WGBS_raw=
DNase_raw=raw_data/DNase/DNase_narrow_peaks
H3K27ac_raw=raw_data/H3K27ac/H3K27ac_narrow_peaks

# Functions
{
  ## Intersect feature with Roadmap data
  ### 4/19/2016, 4/25/2016, 5/3/2016, 5/5/2016, 5/19/2016, 1/26/2017, 2/2/2017, 2/3/2017, 2/6/2017, 2/8/2017, 3/2/2017, 3/6/2017, 5/8/2017, 5/10/2017, 5/11/2017, 5/22/2017, 5/29/2017, 5/30/2017, 5/31/2017, 6/5/2017, 6/12/2017, 6/15/2017, 6/19/2017, 7/4/2017, 8/2/2017, 8/4/2017, 8/5/2017, 8/7/2017, 8/18/2017, 8/25/2017, 8/28/2017, 8/29/2017

  input="$1"
  subdir="$2"
  output="$3" 
   
  ### chromHMM
  while read line  
  do
    bedtools intersect -wo -a $input -b $chromHMM_raw/$line\_15_coreMarks_mnemonics.bed | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> chromHMM/intersect/$subdir/$output\_chromHMM.bed
  done

  ### DNase
  ###TE_landscape/DNase/intersect_sample_list_DNase_peak.sh
  while read line
  do
    bedtools intersect -wo -a $input -b $DNase_raw/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> DNase/intersect/$subdir/$output\_DNase.bed
  done < $DNase_samples

  ### H3K27ac
  while read line
  do 
    bedtools intersect -wo -a $input -b $H3K27ac_raw/$line\-H3K27ac.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> H3K27ac/intersect/$subdir/$output\_H3K27ac.bed
  done < $H3K27ac_samples

  ### WGBS (CpGs) - UPDATE
  split -l 1000000 ~/TE_landscape/all_CpG_Meth.bed
  for file in xa*
  do
    echo $file
    bedtools intersect -wo -a $input -b $file >> WGBS/$output\_CpG_Meth.bed
  done
}

#### Main ####

# Intersect features with Roadmap metrics
# TEs
# Individual TEs
intersect_roadmap features/TEs/rmsk_TEother.txt TEs TE

# Merged TEs
intersect_roadmap features/TEs/merge/rmsk_TEother_merge.txt TEs TE_merge

# Merged TE classes
while read line
do
  class=$line
  intersect_roadmap TEs/class TE_$class
done < features/TEs/class/TEother_class.txt

# Merged TE subfamilies
while read line
do
  subfam=$line
  intersect_roadmap TEs/subfamily TE_$subfam
done < features/TEs/subfamily/subfamilies.txt


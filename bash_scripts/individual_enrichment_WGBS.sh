#!/bin/bash
# TEs in state when subfamily is enriched, WGBS
# 9/30/2016, 5/4/2017, 5/24/2017, 5/29/2017, 5/30/2017, 5/31/2017

input_matrix=$1 #Three columns: subfamily, sample, state
combined_output=$2

# Get individual TEs in the state from candidate subfamilies, filtered to enriched samples
python ~/bin/TE_landscape/pull_individual_TEs.py $input_matrix WGBS/TE_WGBS_state_sorted.txt $combined_output\_temp 8 10 

# Filter to subfamily-sample-state pairs
while read subfamily sample state
do
  awk -v OFS='\t' -v subfam=$subfamily -v state=$state -v sample=$sample '{if(($4 == subfam) && ($10 == state) && ($8 == sample)) print $0}' $combined_output\_temp >> $combined_output
done < $input_matrix

# Split by subfamily, reformat to bed file, unique TEs in state when enriched
while read subfamily state
do
  awk -v OFS='\t' -v subfam=$subfamily -v state=$state '{if(($4 == subfam) && ($10 == state)) a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i], sep[7];}}' $combined_output > enrichment/$subfamily\_$state\_enriched.bed
done < <(awk -v OFS='\t' '{print $1, $3}' $input_matrix | sort | uniq)

rm $combined_output\_temp

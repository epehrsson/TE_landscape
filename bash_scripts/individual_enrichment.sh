#!/bin/bash
# TEs in state when subfamily is enriched
# 9/30/2016, 5/4/2017, 5/24/2017, 5/29/2017, 5/30/2017, 5/31/2017

input_matrix=$1 #Three columns: subfamily, sample, state
combined_output=$2

# Filter subfamily-state to subfamily-sample-state trios 
while read subfamily sample state
do
  awk -v OFS='\t' -v sample=$sample '{if($10 == sample) print $0}' chromHMM/subfamily/by_state/$subfamily\_$state\.txt >> $combined_output
done < $input_matrix

# Split by subfamily, reformat to bed file, unique TEs in state when enriched
while read subfamily state
do
  awk -v OFS='\t' -v subfam=$subfamily -v state=$state '{if(($4 == subfam) && ($8 == state)) a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i], sep[7];}}' $combined_output > enrichment/$subfamily\_$state\_enriched.bed
done < <(awk -v OFS='\t' '{print $1, $3}' $input_matrix | sort | uniq)

# For 8_ZNF/Rpts
#  while read subfamily state; do awk -v OFS='\t' -v subfam=$subfamily '{if($4 == subfam) a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i], sep[7];}}' enrichment/candidate_8_ZNF.Rpts_chromHMM.txt > enrichment/$subfamily\_8_ZNF.Rpts_enriched.bed; done < <(awk -v OFS='\t' '{print $1, $3}' enrichment/candidate_8_ZNF.Rpts_enrich.txt | sort | uniq)

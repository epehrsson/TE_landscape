#!/bin/bash
# TEs in state when subfamily is enriched
# 9/30/2016, 5/4/2017, 5/24/2017, 5/29/2017, 5/30/2017, 5/31/2017

input_matrix=$1 #Three columns: subfamily, sample, state
metric=$2

if [ $metric = "chromHMM" ]; then
   split_path=chromHMM/subfamily/by_state/
elif [ $metric = "WGBS" ]; then
   split_path=WGBS/subfamily/by_state/
fi

# Filter subfamily-state to subfamily-sample-state trios, reformat to bedfile of unique TEs in state when enriched
while read subfamily sample state
do
  if [ $state = "8_ZNF/Rpts" ]; then
    $state="8_ZNF.Rpts"
  fi
  # Change so WGBS column is correct, all subfamily x sample processed at the same time into bedfile
  awk -v OFS='\t' -v sample=$sample '{if($10 == sample) print $0}' "$split_path"$subfamily\_$state\.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i], sep[7];}}' - > enrichment/bedfiles/$subfamily\_$state\_enriched.bed
done < $input_matrix

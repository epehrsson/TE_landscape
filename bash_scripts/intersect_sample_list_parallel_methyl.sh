#!/bin/bash

mapfile -t samples < $1
mapfile -t states < $2
allCpG=$3
binfile=$4

declare -i count=3
for i in "${samples[@]}"
do 
  sample=$i
  count=$count+1

  #Filter CpG file to that sample
  awk -v OFS="\t" -v column="$count" '{print $1,$2,$3,$column}' $allCpG > $sample\_$allCpG

  for j in "${states[@]}"
  do 
    TEs=$j
    suffix=${TEs%.txt}
    #Pull out TEs in that sample
    grep $sample $TEs > $suffix\_$sample\.txt
    echo "Pulled out TEs for state: "$TEs" and sample: "$sample
  
    #Get 10kb bins for those TEs
    bedtools intersect -wo -a $suffix\_$sample\.txt -b $binfile -f 1 -r > $suffix\_$sample\_temp.bed #Removed -s for new bedtools
    awk -v OFS="\t" '{print $18,$19,$20}' $suffix\_$sample\_temp.bed > $suffix\_$sample\_10kb.bed
    echo "Got 10kb regions"

    #Split 10kb bins into 50bp
    python /bar/epehrsson/bin/TE_landscape/split_regions_to_bins.py $suffix\_$sample\_10kb.bed $suffix\_$sample\_10kb_bins.bed 50
    echo "Split 10kb regions into 50bp bins"

    #Intersect bins with CpGs
    bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\_$allCpG > $suffix\_$sample\_temp\-$allCpG
    awk -v OFS="\t" '{if($8!="-1") print $0}' $suffix\_$sample\_temp\-$allCpG > $suffix\_$sample\-$allCpG 
    echo "Intersected with methylation"

    #Get average over that file
    python /bar/epehrsson/bin/TE_landscape/calculate_bin_average_meth.py $suffix\_$sample\-$allCpG /bar/epehrsson/TE_landscape/ChIP_histone/$suffix\_$sample\-$allCpG\_average.txt

    rm $suffix\_$sample\.txt
    rm $suffix\_$sample\_temp.bed
    rm $suffix\_$sample\_10kb.bed
    rm $suffix\_$sample\_10kb_bins.bed 
    rm $suffix\_$sample\_temp\-$allCpG 
    rm $suffix\_$sample\-$allCpG
  done

  rm $sample\_$allCpG
done

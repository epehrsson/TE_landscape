#!/bin/bash

mapfile -t samples < $1
mapfile -t states < $2
mark=$3
binfile=$4

for i in "${samples[@]}"
do 
  sample=$i
  #Get histone files and convert to bedGraph
  wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-$mark\.fc.signal.bigwig
  echo "Downloaded bigwig file"

  bigWigToBedGraph $sample\-$mark\.fc.signal.bigwig $sample\-$mark\.fc.signal.bedGraph
  echo "Converted to bedGraph"

  rm $sample\-$mark\.fc.signal.bigwig

  for j in "${states[@]}"
  do 
    TEs=$j
    suffix=${TEs%.txt}
    #Pull out TEs in that sample
    grep $sample $TEs > $suffix\_$sample\.txt
    echo "Pulled out TEs for state: "$TEs" and sample: "$sample
  
    #Get 10kb bins for those TEs
    bedtools intersect -wo -a $suffix\_$sample\.txt -b $binfile -f 1 -r > $suffix\_$sample\_temp.bed #Removed -s for new version of bedtools
    awk -v OFS="\t" '{print $18,$19,$20}' $suffix\_$sample\_temp.bed > $suffix\_$sample\_10kb.bed
    echo "Got 10kb regions"

    #Split 10kb bins into 50bp
    python /bar/epehrsson/bin/TE_landscape/split_regions_to_bins.py $suffix\_$sample\_10kb.bed $suffix\_$sample\_10kb_bins.bed 50
    echo "Split 10kb regions into 50bp bins"

    #Intersect bins with all histone modification files, then remove
    bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-$mark\.fc.signal.bedGraph > $suffix\_$sample\-$mark\.fc.signal.bedGraph
    echo "Intersected with histone modifications"

    #Get average over that file
    python /bar/epehrsson/bin/TE_landscape/calculate_bin_average.py $suffix\_$sample\-$mark\.fc.signal.bedGraph /bar/epehrsson/TE_landscape/ChIP_histone/$suffix\_$sample\-$mark\_average.txt

    rm $suffix\_$sample\-$mark\.fc.signal.bedGraph
    rm $suffix\_$sample\.txt
    rm $suffix\_$sample\_temp.bed
    rm $suffix\_$sample\_10kb.bed
    rm $suffix\_$sample\_10kb_bins.bed  
  done

  #Removing sample files
  rm $sample\-$mark\.fc.signal.bedGraph
done

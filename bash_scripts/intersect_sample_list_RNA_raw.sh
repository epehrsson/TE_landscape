#!/bin/bash

mapfile -t samples < $1
mapfile -t rmsk < $2

for i in "${samples[@]}"
do 
  sample=$i
  #Get unnormalized RNA-seq file and convert to bed file
  wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/stranded/$sample\.gz
  echo "Downloaded wig file"

  gunzip $sample\.gz
  wig2bed < $sample > $sample\.bed
  echo "Converted to bed file"

  rm $sample

  for j in "${rmsk[@]}"
  do 
    TEs=$j
    suffix=${TEs%.txt}

    #Intersect TEs with RNA file, then remove
    bedtools intersect -wo -a $TEs -b $sample\.bed > $suffix\_$sample\.bed
    echo "Intersected with RNA"

    #Get average over that file
    python /bar/epehrsson/bin/TE_landscape/calculate_average_RNA.py $suffix\_$sample\.bed $TEs $suffix\_$sample\_average.txt

    rm $suffix\_$sample\.bed
  done

  #Removing sample files
  rm $sample\.bed
done


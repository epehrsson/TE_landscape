#!/bin/bash

mapfile -t samples < $1

for i in "${samples[@]}"
do 
  sample=$i
  #Get narrow DNase peaks
  wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$sample\.gz
  echo "Downloaded peak file"

  gunzip $sample\.gz

  TEs=$2
  suffix=${TEs%.txt}

  #Intersect TEs with peak file, then remove
  bedtools intersect -wo -a $TEs -b $sample > $suffix\_$sample
  echo "Intersected with peaks"

  #Removing sample files
  rm $sample
done


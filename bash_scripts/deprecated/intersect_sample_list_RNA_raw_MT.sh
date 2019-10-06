#!/bin/bash

mapfile -t samples < $1

for i in "${samples[@]}"
do 
  sample=$i
  #Get unnormalized RNA-seq file and convert to bed file
  if [ "$sample" == "E050.wig" ]; then
    wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/strandagnostic/$sample\.gz
  else
    wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/stranded/$sample\.gz
  fi
  echo "Downloaded wig file"

  gunzip $sample\.gz
  sed -n '/chrMT/,/chrGL/ p' $sample | sed '$d' > $sample\_MT
  wig2bed < $sample\_MT > $sample\.bed
  echo "Converted to bed file"

  rm $sample
  #rm $sample\_MT

  #Intersect TEs with RNA file, then remove
  bedtools intersect -wo -a refseq_exons_chrM.txt.sorted -b $sample\.bed > refseq_exons_chrM_$sample\.bed
  echo "Intersected with RNA"

  #Get average over that file
  python /bar/epehrsson/bin/TE_landscape/calculate_average_RNA.py refseq_exons_chrM_$sample\.bed refseq_exons_chrM.txt.sorted refseq_exons_chrM_$sample\_average.txt

  #rm refseq_exons_chrM_$sample\.bed

  #Removing sample files
  #rm $sample\.bed
done


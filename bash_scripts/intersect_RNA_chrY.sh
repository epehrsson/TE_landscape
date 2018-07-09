#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=1-4%4

ml bedops
ml bedtools

mapfile -t features < $1

sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p RNA_samples_agnostic.txt )

#Get unnormalized RNA-seq file and convert to bed file
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/strandagnostic/$sample\.wig.gz
echo "Downloaded wig file"

gunzip $sample\.wig.gz
wig2bed < $sample.wig > $sample\.bed
echo "Converted to bed file"

rm $sample.wig
awk '{if($1 == "chrY") print $0}' $sample\.bed > $sample
mv $sample $sample\.bed

for feature in "${features[@]}"
do
  #Intersect feature with RNA file
  bedtools intersect -wo -a $feature -b $sample\.bed > $feature\_$sample\.bed -sorted
  echo "Intersected with RNA"

  #Get average over that file
  declare -i columns=$( head -n 1 $feature\_$sample\.bed | awk -v FS='\t' '{print NF}' -)
  if [ $columns -eq 10 ]
  then
    awk -v OFS='\t' -v sample=$sample '{total[$4]+=$9*$10;len[$4]+=$10}END{for(i in total){print sample, i, total[i], len[i]}}' $feature\_$sample\.bed >> $feature\_average_chrY.txt
  else
    awk -v OFS='\t' -v sample=$sample '{total+=$8*$9;len+=$9}END{print sample, total, len}' $feature\_$sample\.bed >> TE_average_chrY.txt
  fi
  echo "Got average expression"

  rm $feature\_$sample\.bed
done

# Average over genome
awk -v OFS='\t' -v sample=$sample '{total+=$5;len+=$3-$2}END{print sample, total, len}' $sample\.bed >> Genome_average_chrY.txt
echo "Got average genome expression"

#Removing sample files
rm $sample\.bed

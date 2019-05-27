#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=1-4%4

# For chrY only, calculates total RNA-seq read coverage and length of coverage for genome, TEs, TE classes, and TE subfamilies

# Load required packages
ml bedops
ml bedtools

# Load feature list
mapfile -t features < $1

# Load sample list
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p RNA_samples_agnostic.txt )

# Get unnormalized RNA-seq read coverage file from Roadmap site
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/strandagnostic/$sample\.wig.gz
echo "Downloaded wig file"

# Unzip and convert to bed file
gunzip $sample\.wig.gz
wig2bed < $sample.wig > $sample\.bed
echo "Converted to bed file"

# Remove wig file
rm $sample.wig

# Restrict the bed file to chrY
awk '{if($1 == "chrY") print $0}' $sample\.bed > $sample
mv $sample $sample\.bed

# For each feature:
for feature in "${features[@]}"
do
  #Intersect feature with RNA-seq read coverage file
  bedtools intersect -wo -a $feature -b $sample\.bed > $feature\_$sample\.bed -sorted
  echo "Intersected with RNA"

  # Get total read coverage and length of coverage over the feature on chrY

  # Check the number of columns in the file to determine feature
  declare -i columns=$( head -n 1 $feature\_$sample\.bed | awk -v FS='\t' '{print NF}' -)

  ## TE classes and subfamilies
  if [ $columns -eq 10 ]
  then
    awk -v OFS='\t' -v sample=$sample '{total[$4]+=$9*$10;len[$4]+=$10}END{for(i in total){print sample, i, total[i], len[i]}}' $feature\_$sample\.bed >> $feature\_average_chrY.txt

  ## TEs
  else
    awk -v OFS='\t' -v sample=$sample '{total+=$8*$9;len+=$9}END{print sample, total, len}' $feature\_$sample\.bed >> TE_average_chrY.txt
  fi
  echo "Got average expression"

  # Remove intersection file
  rm $feature\_$sample\.bed
done

# Total read coverage and length of coverage over genome
awk -v OFS='\t' -v sample=$sample '{total+=$5;len+=$3-$2}END{print sample, total, len}' $sample\.bed >> Genome_average_chrY.txt
echo "Got average genome expression"

# Remove sample file
rm $sample\.bed

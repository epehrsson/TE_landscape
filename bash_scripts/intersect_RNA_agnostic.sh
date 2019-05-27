#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=20-56%3

# Calculates average RNA-seq read coverage over individual features

# Load required packages
ml bedops
ml bedtools

# Load sample list
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p RNA_samples_agnostic.txt )

# Get unnormalized RNA-seq read coverage file from Roadmap site
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/strandagnostic/$sample\.wig.gz
echo "Downloaded wig file"

# Unzip
gunzip $sample\.wig.gz

# Find all chromosomes within file
echo $sample >> RNA_chr.txt
grep 'chrom=' $sample.wig | awk -v FS=" " '{print $2}' | uniq >> RNA_chr.txt

# Convert to bed file
wig2bed < $sample.wig > $sample\.bed
echo "Converted to bed file"

# Remove wig file
rm $sample.wig

# For exons and TEs:
for j in refseq_exons_unique.txt.sorted rmsk_TEother.txt.sorted
do 
  TEs=$j
  suffix=${TEs%.txt}

  # Intersect exons/TEs with RNA-seq read coverage file
  bedtools intersect -wo -a $TEs -b $sample\.bed > $suffix\_$sample\.bed -sorted
  echo "Intersected with RNA"

  # Get average read coverage over each feature
  python calculate_average_RNA.py $suffix\_$sample\.bed $TEs $suffix\_$sample\_average.txt

  # Remove intersection file
  rm $suffix\_$sample\.bed
done

# Remove sample file
rm $sample\.bed

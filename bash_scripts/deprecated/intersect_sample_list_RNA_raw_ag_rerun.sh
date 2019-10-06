#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --array=1-3%1

ml bedops
ml bedtools

sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p RNA_samples_rerun.txt )

#Get unnormalized RNA-seq file and convert to bed file
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/strandagnostic/$sample\.gz
echo "Downloaded wig file"

gunzip $sample\.gz
wig2bed < $sample > $sample\.bed
echo "Converted to bed file"

rm $sample

for j in refseq_exons.txt.sorted rmsk_TEother.txt.sorted
do 
  TEs=$j
  suffix=${TEs%.txt}

  #Intersect TEs with RNA file, then remove
  bedtools intersect -wo -a $TEs -b $sample\.bed > $suffix\_$sample\.bed -sorted
  echo "Intersected with RNA"

  #Get average over that file
  python calculate_average_RNA.py $suffix\_$sample\.bed $TEs $suffix\_$sample\_average.txt

  rm $suffix\_$sample\.bed
done

#Removing sample files
rm $sample\.bed

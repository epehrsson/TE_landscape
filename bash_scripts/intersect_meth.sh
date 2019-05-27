#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --array=1-37%10

# For TEs in each chromHMM state, finds the DNA methylation level over a 10kb region centered on the TE, in bins of 50bp

# Load required packages
ml bedtools

# Load sample list
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p WGBS_samples.txt )

# Load chromHMM states
mapfile -t states < $1

# Filter the CpG methylation level bedfile to the sample
column="$((${SLURM_ARRAY_TASK_ID}+3))"
cut -f1-3,$column CpG_TE_Meth.bed > $sample\_meth.bed

# For each chromHMM state:
for state in "${states[@]}"
do
  # Get extended 10kb regions centered on the TE for TEs in that state in that sample
  bedtools intersect -wo -a ever/rmsk_TEother_$sample\_$state\.txt -b rmsk_TEother_10kb.bed -f 1 -r > rmsk_TEother_$sample\_$state\_temp.bed
  awk -v OFS="\t" '{print $19,$20,$21}' rmsk_TEother_$sample\_$state\_temp.bed > rmsk_TEother_$sample\_$state\_10kb.bed
  echo "Got 10kb regions"

  # Split 10kb regions into 50bp bins
  python split_regions_to_bins.py rmsk_TEother_$sample\_$state\_10kb.bed rmsk_TEother_$sample\_$state\_10kb_bins.bed 50
  echo "Split 10kb regions into 50bp bins"

  # Intersect 50bp bins with CpG methylation file
  bedtools intersect -wo -a rmsk_TEother_$sample\_$state\_10kb_bins.bed -b $sample\_meth.bed | awk -v OFS="\t" '{if($8!="-1") print $0}' - > rmsk_TEother_$sample\_$state\_meth.bed
  echo "Intersected with methylation"

  # Calculate the average methylation level over that bin
  python calculate_bin_average_meth.py rmsk_TEother_$sample\_$state\_meth.bed averages/rmsk_TEother_$sample\_$state\_meth_average.txt
  echo "Calculated average coverage"

  # Remove intersection file and extended TE and bin files
  rm rmsk_TEother_$sample\_$state\_temp.bed
  rm rmsk_TEother_$sample\_$state\_10kb.bed
  rm rmsk_TEother_$sample\_$state\_10kb_bins.bed
  rm rmsk_TEother_$sample\_$state\_meth.bed
done

# Remove sample file
rm $sample\_meth.bed

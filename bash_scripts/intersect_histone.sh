#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH --array=1%1

# For TEs in each chromHMM state, finds the fold change enrichment over input for several histone modifications and DHS
# Over a 10kb region centered on the TE, in bins of 50bp

# Load required packages
ml kentsrc
ml bedtools

# Load sample list
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p mnemonics.txt )

# Load chromHMM states
mapfile -t states < $1

# Load epigenetic marks
mapfile -t marks < $2

# For each epigenetic mark:
for mark in "${marks[@]}"
do
  # Download fold change files from Roadmap site
  wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-$mark\.fc.signal.bigwig
  echo "Downloaded "$mark" bigwig file for sample "$sample

  # Convert to bedGraph file
  bigWigToBedGraph $sample\-$mark\.fc.signal.bigwig $sample\-$mark\.fc.signal.bedGraph
  echo "Converted to bedGraph"

  # Remove bigwig file
  rm $sample\-$mark\.fc.signal.bigwig
done

# For each chromHMM state:
for state in "${states[@]}"
do
  # Get extended 10kb regions centered on the TE for TEs in that state in that sample
  bedtools intersect -wo -a ever/rmsk_TEother_$sample\_$state\.txt -b rmsk_TEother_10kb.bed -f 1 -r > rmsk_TEother_$sample\_$state\_temp.bed
  awk -v OFS="\t" '{print $19,$20,$21}' rmsk_TEother_$sample\_$state\_temp.bed > rmsk_TEother_$sample\_$state\_10kb.bed
  echo "Got 10kb regions"

  # Split the 10kb regions into 50bp bins
  python split_regions_to_bins.py rmsk_TEother_$sample\_$state\_10kb.bed rmsk_TEother_$sample\_$state\_10kb_bins.bed 50
  echo "Split 10kb regions into 50bp bins"

  # For each epigenetic mark
  for mark in "${marks[@]}"
  do
    # Intersect 50bp bins with fold change file
    bedtools intersect -wo -a rmsk_TEother_$sample\_$state\_10kb_bins.bed -b $sample\-$mark\.fc.signal.bedGraph > rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph
    echo "Intersected with histone modifications"

    # Calculate the average fold change over that bin
    python calculate_bin_average.py rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph averages/rmsk_TEother_$sample\_$state\_$mark\_average.txt
    echo "Calculated average coverage"

    # Remove intersection file
    rm rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph
  done

  # Remove extneded TE and bin files
  rm rmsk_TEother_$sample\_$state\_temp.bed
  rm rmsk_TEother_$sample\_$state\_10kb.bed
  rm rmsk_TEother_$sample\_$state\_10kb_bins.bed
done

# Remove sample file
rm $sample\-*.fc.signal.bedGraph

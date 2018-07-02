#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH --array=1%1

ml kentsrc
ml bedtools

sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p mnemonics.txt )

mapfile -t states < $1
mapfile -t marks < $2

# Get histone/DNase fold change files and convert to bedGraph
for mark in "${marks[@]}"
do
  wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-$mark\.fc.signal.bigwig
  echo "Downloaded "$mark" bigwig file for sample "$sample

  bigWigToBedGraph $sample\-$mark\.fc.signal.bigwig $sample\-$mark\.fc.signal.bedGraph
  echo "Converted to bedGraph"

  rm $sample\-$mark\.fc.signal.bigwig
done

for state in "${states[@]}"
do
  # Get 10kb bins for TEs in that state in that sample
  bedtools intersect -wo -a ever/rmsk_TEother_$sample\_$state\.txt -b rmsk_TEother_10kb.bed -f 1 -r > rmsk_TEother_$sample\_$state\_temp.bed #New version (2.25.0) fails with -s
  awk -v OFS="\t" '{print $19,$20,$21}' rmsk_TEother_$sample\_$state\_temp.bed > rmsk_TEother_$sample\_$state\_10kb.bed
  echo "Got 10kb regions"

  # Split 10kb bins into 50bp
  python split_regions_to_bins.py rmsk_TEother_$sample\_$state\_10kb.bed rmsk_TEother_$sample\_$state\_10kb_bins.bed 50
  echo "Split 10kb regions into 50bp bins"

  for mark in "${marks[@]}"
  do
    # Intersect bins with all histone modification/DNase files, then remove
    bedtools intersect -wo -a rmsk_TEother_$sample\_$state\_10kb_bins.bed -b $sample\-$mark\.fc.signal.bedGraph > rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph
    echo "Intersected with histone modifications"

    # Get average over that file
    python calculate_bin_average.py rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph averages/rmsk_TEother_$sample\_$state\_$mark\_average.txt
    echo "Calculated average coverage"

    rm rmsk_TEother_$sample\_$state\_$mark\.fc.signal.bedGraph
  done

  rm rmsk_TEother_$sample\_$state\_temp.bed
  rm rmsk_TEother_$sample\_$state\_10kb.bed
  rm rmsk_TEother_$sample\_$state\_10kb_bins.bed
done

# Removing sample files
rm $sample\-*.fc.signal.bedGraph

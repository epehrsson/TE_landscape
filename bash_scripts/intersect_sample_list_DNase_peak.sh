#!/bin/bash

# Intersects DHS peaks identified in each sample with individual TEs

# Read in list of samples with DHS peaks
mapfile -t samples < $1

# For each sample:
for i in "${samples[@]}"
do 
  sample=$i

  # Download narrow DHS peaks bedfile from Roadmap site
  wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$sample\.gz
  echo "Downloaded peak file"

  # Unzip
  gunzip $sample\.gz

  # TE bedfile
  TEs=$2
  suffix=${TEs%.txt}

  # Intersect TEs with peak file
  bedtools intersect -wo -a $TEs -b $sample > $suffix\_$sample
  echo "Intersected with peaks"

  # Remove DHS peaks sample file
  rm $sample
done

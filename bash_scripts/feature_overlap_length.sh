# Creates a matrix listing the length of each feature and its overlap with TEs

# For each feature, where features is a list of 7 RefSeq features
while read feature
do 
  # For intergenic regions:
  if [ $feature = "intergenic" ]; then
    echo "refseq_intergenic"
    # Total length, no contigs
    awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' ~/genic_features/RefSeq/refseq_intergenic.txt
    # Total length, no chrY or contigs
    awk '{if(($1 !~ /_/) && ($1 != "chrY")) sum+=$3-$2}END{print sum}' ~/genic_features/RefSeq/refseq_intergenic.txt
    # Length of intergenic regions overlapping TEs
    bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk '{sum+=$7}END{print sum}' -   
    bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk '{sum+=$7}END{print sum}' -
  # For genic features
  else
    # Uses merged features (redundant bases removed)
    for file in ~/genic_features/RefSeq/refseq_$feature*_merge.txt
    do
      echo $( basename $file _merge.txt )
      # Total length of feature, no contigs
      awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' $file
      # Total kength of feature, no chrY or contigs
      awk '{if(($1 !~ /_/) && ($1 != "chrY")) sum+=$3-$2}END{print sum}' $file
      # Length of overlap between TEs, Refseq features, unstranded  
      bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b $file | awk '{sum+=$7}END{print sum}' -
      # Length of overlap between TEs, Refseq features, stranded
      awk -v OFS='\t' '{print $1, $2, $3, $4, '0', $7}' features/TEs/rmsk_TEother.txt | bedtools intersect -a - -b ~/genic_features/RefSeq/$( basename $file _merge.txt ).txt -s | sort -k1,1V -k2,2n - | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'
    done
  fi
done < features/features.txt > features/intersect_features/feature_overlap.txt

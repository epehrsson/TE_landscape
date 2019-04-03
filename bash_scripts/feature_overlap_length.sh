# Feature overlap
#TE_landscape/features/intersect_features/feature_overlap.txt

while read feature
do 
  if [ $feature = "intergenic" ]; then
    echo "refseq_intergenic"
    awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' ~/genic_features/RefSeq/refseq_intergenic.txt
    awk '{if(($1 !~ /_/) && ($1 != "chrY")) sum+=$3-$2}END{print sum}' ~/genic_features/RefSeq/refseq_intergenic.txt
    bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk '{sum+=$7}END{print sum}' -   
    bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk '{sum+=$7}END{print sum}' -
  else
    for file in ~/genic_features/RefSeq/refseq_$feature*_merge.txt
    do
      echo $( basename $file _merge.txt )
      # Length of feature
      awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' $file
      # Length of feature, no chrY
      awk '{if(($1 !~ /_/) && ($1 != "chrY")) sum+=$3-$2}END{print sum}' $file
      # Length of overlap between TEs, Refseq features, unstranded  
      bedtools intersect -wo -a features/TEs/merge/rmsk_TEother_merge.txt -b $file | awk '{sum+=$7}END{print sum}' -
      # Length of overlap between TEs, Refseq features, stranded
      awk -v OFS='\t' '{print $1, $2, $3, $4, '0', $7}' features/TEs/rmsk_TEother.txt | bedtools intersect -a - -b ~/genic_features/RefSeq/$( basename $file _merge.txt ).txt -s | sort -k1,1V -k2,2n - | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'
    done
  fi
done < features/features.txt > features/intersect_features/feature_overlap.txt

# Where features is a list of 7 genic features

# Feature overlap
# 6/13/2017, 6/14/2017, 6/15/2017

#TE_landscape/features/intersect_features/feature_overlap.txt	
# Length of overlap between TEs, Refseq features	 
for file in ~/genic_features/refseq*merge*txt; do echo $file; awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' $file; done
for file in ~/genic_features/refseq*merge*txt; do bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file | awk '{sum+=$7}END{print sum}' -; done
awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' ~/genic_features/refseq_intergenic.txt
bedtools intersect -wo -a rmsk_TEother_merge.txt -b ~/genic_features/refseq_intergenic.txt | awk '{sum+=$7}END{print sum}' -

# Length of stranded overlap	 
for file in ~/genic_features/refseq_*.txt; do echo $file; awk -v OFS='\t' '{print $1, $2, $3, $4, '0', $7}' rmsk_TEother.txt | bedtools intersect -a - -b $file -s | sort -k1,1V -k2,2n - | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'; done

# Length of coding/non-coding features	 
for file in refseq*pc_merge.txt; do echo $file ; awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' $file ; done
for file in refseq*nc_merge.txt; do echo $file ; awk '{if($1 !~ /_/) sum+=$3-$2}END{print sum}' $file ; done

# Length of overlap between TEs, Refseq features (coding/non-coding)	 
for file in ~/genic_features/refseq*pc_merge.txt; do echo $file; bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file | awk '{sum+=$7}END{print sum}' -; done
for file in ~/genic_features/refseq*nc_merge.txt; do echo $file; bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file | awk '{sum+=$7}END{print sum}' -; done

# Length of stranded overlap between TEs, Refseq features (coding/non-coding)	 
for file in ~/genic_features/refseq_*_pc.txt; do echo $file; awk -v OFS='\t' '{print $1, $2, $3, $4, '0', $7}' rmsk_TEother.txt | bedtools intersect -a - -b $file -s | sort -k1,1V -k2,2n - | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'; done
for file in ~/genic_features/refseq_*_nc.txt; do echo $file; awk -v OFS='\t' '{print $1, $2, $3, $4, '0', $7}' rmsk_TEother.txt | bedtools intersect -a - -b $file -s | sort -k1,1V -k2,2n - | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'; done

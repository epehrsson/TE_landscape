# Generated merged Refseq features
# 1/18/2017, 5/8/2017, 6/13/2017, 6/15/2017, 8/2/2017

# Merged features
#genic_features/refseq_3UTR_merge.txt
#genic_features/refseq_5UTR_merge.txt
#genic_features/refseq_coding_exon_merge.txt
#genic_features/refseq_exons_merge.txt
#genic_features/refseq_introns_merge.txt
#genic_features/refseq_genes_merge.txt
for file in refseq_*.txt ; do sort -k1,1 -k2,2n $file | bedtools merge -i - > $file\_merge; done
rename 's/.txt_merge/_merge.txt/' *.txt_merge

#genic_features/refseq_promoters_merge.txt
bedtools slop -i refseq_up2000.txt -g hg19.genome -l 0 -r 500 -s | sort -k1,1 -k2,2n - | bedtools merge -i - > refseq_promoters_merge.txt

#genic_features/cpgIslandExtUnmasked_merge.txt
sort -k1,1 -k2,2n cpgIslandExtUnmasked.txt | bedtools merge -i - > cpgIslandExtUnmasked_merge.txt

# Merged features, coding/non-coding
#genic_features/refseq_[feature]_pc_merge.txt [7 files]
for file in refseq*.txt; do awk -v OFS='\t' '{if(substr($4,0,2) == "NM") print $1, $2, $3}' $file | sort -k1,1V -k2,2n - | bedtools merge -i - > $(basename "$file" .txt)\_pc_merge.txt; done
#genic_features/refseq_[feature]_nc_merge.txt [6 files]
for file in refseq*.txt; do awk -v OFS='\t' '{if(substr($4,0,2) == "NR") print $1, $2, $3}' $file | sort -k1,1V -k2,2n - | bedtools merge -i - > $(basename "$file" .txt)\_nc_merge.txt; done

# Merged features, no TEs
#genic_features/RefSeq/refseq_[feature]_merge_noTE.txt [7 files]
#genic_features/RefSeq/refseq_[feature]_nc_merge_noTE.txt [6 files]
#genic_features/RefSeq/refseq_[feature]_pc_merge_noTE.txt [7 files]
for file in *merge*; do bedtools subtract -a $file -b ~/TE_landscape/features/TEs/merge/rmsk_TEother_merge.txt > $(basename "$file" .txt)_noTE.txt; done
#genic_features/refseq_intergenic_noTE.txt
bedtools subtract -a refseq_intergenic.txt -b ../TE_landscape/rmsk_TEother_merge.txt > refseq_intergenic_noTE.txt

# Merged features, no TEs or repeats
#genic_features/refseq_intergenic_noTERepeats.txt
bedtools subtract -a refseq_intergenic_noTE.txt -b ../TE_landscape/rmsk_repeats.txt > refseq_intergenic_noTERepeats.txt

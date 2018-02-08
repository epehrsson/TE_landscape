# TE intersection with features
# 5/3/2016, 9/12/2016, 9/20/2016, 10/28/2016, 11/7/2016, 2/2/2017, 2/8/2017, 6/13/2017, 6/14/2017, 8/24/2017, 8/28/2017, 1/7/2018

# Intersection between TEs, RefSeq features	 
#TE_landscape/features/intersect_features/rmsk_TEother_refseq_[feature].txt [8 files]	
for file in ~/genic_features/refseq*merge*txt; do bedtools intersect -wo -a rmsk_TEother.txt -b $file | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TEother_$(basename "$file" .txt)\.txt; done
bedtools intersect -wo -a rmsk_TEother.txt -b ~/genic_features/refseq_intergenic.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TEother_refseq_intergenic.txt

# Intersection between TEs, Refseq features (coding/non-coding)	 
#TE_landscape/features/intersect_features/rmsk_TEother_refseq_[feature]_pc_merge.txt [6 files]	
#TE_landscape/features/intersect_features/rmsk_TEother_refseq_[feature]_nc_merge.txt [5 files]		
for file in ~/genic_features/RefSeq/refseq_*c_merge.txt; do bedtools intersect -wo -a ../TEs/rmsk_TEother.txt -b $file | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TEother_$(basename "$file" .txt)\.txt; done

# Intersection between TEs, CpG islands	 
#TE_landscape/features/intersect_features/rmsk_TEother_cpgIsland.txt	
bedtools intersect -wo -a rmsk_TEother.txt -b ~/genic_features/cpgIslandExtUnmasked_merge.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TEother_cpgIsland.txt

# Intersection between shuffled TEs, Refseq features	 
#TE_landscape/features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_#_cpgIsland.txt [10 files]	
#TE_landscape/features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_#_refseq_intergenic.txt [10 files]		
#TE_landscape/features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_#_refseq_[feature]_merge.txt [70 files]		
#TE_landscape/features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_#_refseq_[feature]_nc_merge.txt [60 files]		
#TE_landscape/features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_#_refseq_[feature]_pc_merge.txt [70 files]		
for i in {1..10}; do for file in ~/genic_features/RefSeq/*_merge.txt; do bedtools intersect -wo -a ../rmsk_TE_shuffle_$i\.txt -b $file | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TE_shuffle_$i\_$(basename "$file" .txt)\.txt; done; bedtools intersect -wo -a ../rmsk_TE_shuffle_$i\.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TE_shuffle_$i\_refseq_intergenic.txt; bedtools intersect -wo -a ../rmsk_TE_shuffle_$i\.txt -b ~/genic_features/cpgIslandExtUnmasked_merge.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > rmsk_TE_shuffle_$i\_cpgIsland.txt; done

# Updated 1/5/18
 for j in {1..10}; do for file in ~/genic_features/RefSeq/*_merge.txt; do bedtools intersect -wo -a features/shuffled_TEs/rmsk_TE_shuffle_$j\.txt -b $file | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_$j\_$(basename "$file" .txt)\.txt; done; bedtools intersect -wo -a features/shuffled_TEs/rmsk_TE_shuffle_$j\.txt -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_$j\_refseq_intergenic.txt; bedtools intersect -wo -a features/shuffled_TEs/rmsk_TE_shuffle_$j\.txt -b ~/genic_features/cpgIslandExtUnmasked_merge.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=$11}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' - > features/shuffled_TEs/intersect_features/rmsk_TE_shuffle_$j\_cpgIsland.txt; done


# TEs that overlap the hg19 blacklist	 
#TE_landscape/features/intersect_features/TE_blacklist.bed	
bedtools intersect -wao -a rmsk_TE.txt -b ../../genomes/hg19/blacklist.bed > TE_blacklist.bed
awk '{if($11 != 0) print $0}' TE_blacklist.bed
#TE_landscape/features/intersect_features/other_blacklist.bed	
bedtools intersect -wao -a rmsk_other.txt -b ../../genomes/hg19/blacklist.bed > other_blacklist.bed

# Vista
# Human, positively validated VISTA enhancers intersected with TEs	 
#TE_landscape/features/intersect_features/vista_enhancers_rmsk_TE.txt	
bedtools intersect -wo -a vista_enhancers_human.bed -b rmsk_TE.txt > vista_enhancers_rmsk_TE.txt
#TE_landscape/features/intersect_features/vista_enhancers_rmsk_other.txt	
bedtools intersect -wo -a vista_enhancers_human.bed -b rmsk_other.txt > vista_enhancers_rmsk_other.txt

# C-Gate
# CGate TE coordinates intersected with all TEs	 
#TE_landscape/features/intersect_features/CGate_TE.txt	
bedtools intersect -wo -a C_gate.txt -b rmsk_TE.txt > CGate_TE.txt
bedtools intersect -wo -a C_gate.txt -b rmsk_other.txt > CGate_other.txt
# Updated 1/7/18 after hg18->hg19 liftover
bedtools intersect -wo -a features/C_gate.txt -b features/TEs/rmsk_TEother.txt > features/intersect_features/CGate_TE.txt

# Intersection between TEs and Segway promoters	 
#TE_landscape/features/intersect_features/TE_promoter_intersect.txt	
bedtools intersect -wo -a rmsk_TE.txt -b promoters.txt > TE_promoter_intersect.txt
#TE_landscape/features/intersect_features/other_promoter_intersect.txt	
bedtools intersect -wo -a rmsk_other.txt -b promoters.txt > other_promoter_intersect.txt

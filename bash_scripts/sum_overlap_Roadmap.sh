# Sum overlap by merged feature
# 4/20/2016, 4/25/2016, 4/26/2016, 5/19/2016, 6/27/2016, 9/5/2016,
# 2/2/2017, 2/3/2017, 2/8/2017, 3/2/2017, 5/8/2017, 5/24/2017, 6/5/2017, 6/12/2017, 7/4/2017, 8/2/2017, 8/4/2017, 8/5/2017, 8/7/2017, 8/18/2017, 9/19/2017

# chromHMM 
# Difference between hg19 chromosome lengths and chromHMM annotation (11/9/17)
## Genome: chromHMM misses 2583 from all chromosomes and 2417 from all but chrY (3/1/18)
awk -v OFS='\t' '{chr[$1]+=$3-$2}END{for(i in chr){print i, chr[i]}}' raw_data/chromHMM/E001_15_coreMarks_mnemonics.bed

# Features: not missing any bases, indicating all unannotated bases are intergenic
for file in chromHMM/Refseq_features/intersect/chromHMM_*.bed; do awk -v OFS='\t' -v feature=$(basename $file .bed) '{a[$9]+=$8}END{for(i in a){print feature, a[i]}}' $file | sort | uniq; done < features/features.txt

# Combines bp counts of state per sample into one file
#TE_landscape/chromHMM/combine_states.sh	

# Whole genome
# Number of bases in each state in each sample	 
#TE_landscape/chromHMM/genome/sample_summaries/E#_15_coreMarks_mnemonics.bed_state [127 files]	
for file in chromHMM_bedfiles/E*_15_coreMarks_mnemonics.bed; do awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' $file > $file\_state; while read line; do grep $line $file | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' - >> $file\_state; done < chromHMM_states.txt ; done

# Table of number of bases in each state in each sample	
#TE_landscape/chromHMM/genome/mnemonics_state.txt	
bash combine_states.sh chromHMM_states.txt mnemonics.txt chromHMM_bedfiles/*state

# By chromosome
 while read line; do awk -v OFS='\t' -v sample=$line '{a[$1,$4]+=$3-$2}END{for(i in a){split(i,sep,SUBSEP); print sample, sep[1], sep[2], a[i];}}' raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed; done < sample_lists/mnemonics.txt > chromHMM/chromosome_states.txt

# Merged TEs
#TE_landscape/chromHMM/TEs/sample_summaries/TE_merge/E#_15_coreMarks_mnemonics.bed_TE_merge_state [127 files]
for file in chromHMM_TE_merge/*; do awk '{SUM+=$8}END{print SUM}' $file > $file\_state; while read line; do grep $line $file | awk '{SUM+=$8}END{print SUM}' - >> $file\_state; done < chromHMM_states.txt ; done
#TE_landscape/chromHMM/TEs/sample_summaries/TEother_merge/E#_15_coreMarks_mnemonics.bed_TEother_merge_state [127 files]
for file in chromHMM_other/*; do awk '{SUM+=$8}END{print SUM}' $file > $file\_state; while read line; do grep $line $file | awk '{SUM+=$8}END{print SUM}' - >> $file\_state; done < chromHMM_states.txt ; done

# Table of number of bases in each state in merged TEs in each sample	
#TE_landscape/chromHMM/mnemonics_TEmerge_states.txt	
bash combine_states.sh chromHMM_states.txt mnemonics.txt chromHMM_TE_merge/*state
#TE_landscape/chromHMM/mnemonics_TEother_merge_states.txt	
bash combine_states.sh chromHMM_states.txt mnemonics.txt chromHMM_other/*state

# Number of bases in each state in each sample in TEs	(redundant bases!) 
#TE_landscape/chromHMM/TEs/sample_summaries/TE/E#_15_coreMarks_mnemonics.bed_TE_state [127 files]	
for file in chromHMM_bedfiles/E*_15_coreMarks_mnemonics.bed_TE; do awk '{SUM+=$12}END{print SUM}' $file > $file\_TE_state; while read line; do grep $line $file | awk '{SUM+=$12}END{print SUM}' - >> $file\_TE_state; done < chromHMM_states.txt ; done

# Table of number of bases in each state in TEs in each sample	
#TE_landscape/chromHMM/mnemonics_TE_states.txt	
bash combine_states.sh chromHMM_states.txt mnemonics.txt chromHMM_TE/*state

# Merged TE classes
# Number of bases in each state in merged TE classes in each sample
 while read line; do awk -v OFS='\t' -v class=$line '{a[$7,$9]+=$8}END{for(i in a) {split (i, sep, SUBSEP); print class, sep[1], sep[2], a[i];}}' ../TEs/intersect/class/rmsk_$line\.txt_chromHMM.bed >> class_state_sample.txt; done < classes.txt

#TE_landscape/chromHMM/class/rmsk_[class].txt_chromHMM.bed_state [13 files]
for file in TE_classes/rmsk_*.txt_chromHMM.bed; do awk 'BEGIN{SUM=0}{SUM+=$8}END{print SUM}' $file > $file\_state; while read line; do grep $line $file | awk 'BEGIN{SUM=0}{SUM+=$8}END{print SUM}' - >> $file\_state; done < chromHMM_states.txt ; done
awk 'BEGIN{SUM=0}{SUM+=$8}END{print SUM}' rmsk_Unconfident.txt_chromHMM.bed > rmsk_Unconfident.txt_chromHMM.bed_state; while read line; do grep $line rmsk_Unconfident.txt_chromHMM.bed | awk 'BEGIN{SUM=0}{SUM+=$8}END{print SUM}' - >> rmsk_Unconfident.txt_chromHMM.bed_state; done < chromHMM_states.txt
#TE_landscape/chromHMM/class/rmsk_Unconfident_RC.txt_chromHMM.bed_state
awk '{SUM+=$8}END{print SUM}' ../TEs/intersect/class/rmsk_Unconfident_RC.txt_chromHMM.bed > rmsk_Unconfident_RC.txt_chromHMM.bed_state; while read line; do grep $line ../TEs/intersect/class/rmsk_Unconfident_RC.txt_chromHMM.bed | awk '{SUM+=$8}END{print SUM}' - >> rmsk_Unconfident_RC.txt_chromHMM.bed_state; done < ../chromHMM_states.txt 
#TE_landscape/chromHMM/class/chromHMM_class_TEother.txt
paste rmsk*state >chromHMM_class_TEother.txt
paste chromHMM_class_TEother.txt rmsk_Unconfident.txt_chromHMM.bed_state > chromHMM_class_TEother.txt
paste chromHMM_class_TEother.txt rmsk_Unconfident_RC.txt_chromHMM.bed_state > chromHMM_class_TEother.txt

# Merged TE subfamilies
#TE_landscape/chromHMM/subfamily/subfamily_state_sample.txt
awk -v OFS='\t' '{a[$4, $8, $10]+=$9;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' subfamily_state_sample.bed > subfamily_state_sample.txt

# Unique peaks overlapping each subfamily
while read line; do awk -v OFS='\t' -v sample=$line '{if($11 != "8_ZNF/Rpts") print > ""$4"_"$11"_"sample".txt"}' ~/TE_landscape/chromHMM/summit/TEs/rmsk_TEother_$line\_chromHMM.bed; while read subfamily; do while read state; do awk -v sample=$line -v subfam=$subfamily -v state=$state -v OFS='\t' '{a[$8, $9, $10]+=1}END{print state, sample, subfam, length(a)}' $subfamily\_$state\_$line\.txt >> subfamily_chromHMM_sample_summit.txt; done < chromHMM_states.txt; done < ~/TE_landscape/features/TEs/subfamily/subfamilies.txt; rm *_$line\.txt; done < mnemonics.txt
while read line; do awk -v OFS='\t' -v sample=$line '{if($11 == "8_ZNF/Rpts") print > ""$4"_8_ZNF.Rpts_"sample".txt"}' ~/TE_landscape/chromHMM/summit/TEs/rmsk_TEother_$line\_chromHMM.bed; while read subfamily; do awk -v sample=$line -v subfam=$subfamily -v OFS='\t' '{a[$8, $9, $10]+=1}END{print "8_ZNF/Rpts", sample, subfam, length(a)}' $subfamily\_8_ZNF.Rpts_$line\.txt >> subfamily_chromHMM_sample_summit_8.txt; done < ~/TE_landscape/features/TEs/subfamily/subfamilies.txt; rm *_$line\.txt; done < ~/TE_landscape/sample_lists/mnemonics.txt
cat subfamily_chromHMM_sample_summit.txt subfamily_chromHMM_sample_summit_8.txt > ~/TE_landscape/chromHMM/subfamily/subfamily_chromHMM_sample_summit.txt

# Merged features
#TE_landscape/chromHMM/Refseq_features/chromHMM_CDS_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_3UTR_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_5UTR_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_exon_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_intron_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_promoter_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_refseq_intergenic_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_repeats_states.txt
awk -v OFS='\t' '{a[$1, $2, $9]+=$10}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], a[i]}}' chromHMM/Refseq_features/chromHMM_features.txt > chromHMM/chromHMM_refseq_features.txt

# Table of number of bases in each state in merged genic features in each sample	 
#TE_landscape/chromHMM/Refseq_features/chromHMM_feature_states.txt	
while read line; do awk -v OFS='\t' -v feature=$line '{print $0, feature}' chromHMM_$line\_states.txt >> chromHMM_feature_states.txt; done < features.txt
awk -v OFS='\t' '{print $0, "genome_noTE"}' chromHMM_genome_noTE_states.txt >> chromHMM_feature_states.txt

# Merged features, no TEs
#TE_landscape/chromHMM/Refseq_features/chromHMM_refseq_[feature]_merge_noTE_states.txt [6 files]
#TE_landscape/chromHMM/Refseq_features/chromHMM_refseq_[feature]_nc_merge_noTE_states.txt [5 files]
#TE_landscape/chromHMM/Refseq_features/chromHMM_refseq_[feature]_pc_merge_noTE_states.txt [6 files]
for file in chromHMM_refseq*; do name=$(basename "$file" .txt); awk -v OFS='\t' '{a[$7,$9]+=$8}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i];}}' $file > ../$name\_states.txt; done
#TE_landscape/chromHMM/Refseq_features/chromHMM_refseq_intergenic_noTE_states.txt
#TE_landscape/chromHMM/Refseq_features/chromHMM_genome_noTE_states.txt
awk -v OFS='\t' '{a[$5,$4]+=$3-$2}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i];}}' intersect/chromHMM_genome_noTE.bed > chromHMM_genome_noTE_states.txt

# Table of number of bases in each state in merged genic features in each sample, no TEs	 
#TE_landscape/chromHMM/Refseq_features/chromHMM_feature_noTE_states.txt	
while read line; do awk -v OFS='\t' -v feature=$line '{print $0, feature}' chromHMM_refseq_$line\_merge_noTE_states.txt >> chromHMM_feature_noTE_states.txt; done < features.txt
awk -v OFS='\t' '{print $0, "intergenic"}' chromHMM_refseq_intergenic_noTE_states.txt >> chromHMM_feature_noTE_states.txt
awk -v OFS='\t' '{print $2, $1, $3, "genome_noTE"}' chromHMM_genome_noTE_states.txt >> chromHMM_feature_noTE_states.txt #Updated 8/7/17

# Merged features, no TEs or repeats
#TE_landscape/chromHMM/Refseq_features/chromHMM_intergenic_states.txt
awk -v OFS='\t' '{a[$7,$9]+=$8}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i];}}' chromHMM_intergenic_noTERepeats.bed > chromHMM_intergenic_states.txt

# Segwey promoters (incorrect)
#TE_landscape/chromHMM/Segway_promoters/chromHMM_promoter.txt
awk -v OFS='\t' '{a[$4, $6]+=$5;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' all_chromHMM_promoter.txt > chromHMM_promoter.txt

# DNase

# Merged TEs
#TE_landscape/DNase/rmsk_TEother_merge_DNase_contribution.txt
awk -v OFS='\t' '{a[$15]+=$14}END{for(i in a){print i, a[i]}}' rmsk_TEother_merge_DNase.txt | sort > rmsk_TEother_merge_DNase_contribution.txt

# Number and width of DNase peaks overall, overlapping TEs per sample; including summit rules	 
#TE_landscape/DNase/DNase_stats.txt	
while read line; do wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$line\.gz; gunzip $line\.gz; echo $line >> DNase_stats.txt; wc -l $line>> DNase_stats.txt; awk '{sum+=$3-$2}END{print sum}' $line >> DNase_stats.txt; done < DNase_peaks.txt
while read line; do awk '{print $8, $9, $10}' rmsk_TEother_$line\-DNase.macs2.narrowPeak | sort | uniq | wc -l >> test; done < ../DNase_samples.txt
while read line; do awk -v OFS='\t' '{print $8, $9, $10}' DNase/summit/TEs/rmsk_TEother_$line\_DNase_summit.txt | sort | uniq | wc -l >> DNase/rmsk_TEother_DNase_summit_stats.txt; done < sample_lists/DNase_samples.txt

# Merged TE classes
#TE_landscape/DNase/class_DNase_sample.txt
awk -v OFS='\t' '{a[$4, $16]+=$15}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i];}}' rmsk_TEother_class_Dnase.txt > class_DNase_sample.txt

# Merged TE subfamilies (unique peaks overlapping subfamily)
while read line; do while read line2; do awk -v sample=$line -v subfam=$line2 -v OFS='\t' '{if($4 == subfam) a[$8, $9, $10]+=1}END{print sample, subfam, length(a)}' DNase/summit/TEs/rmsk_TEother_$line\_DNase_summit.txt >> DNase/subfamily_DNase_sample_summit.txt; done < features/TEs/subfamily/subfamilies.txt; done < sample_lists/DNase_samples.txt

# Merged features
awk -v OFS='\t' '{a[$1, $2]+=$16}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i]}}' DNase/refseq_features_DNase_intersect.txt > DNase/refseq_features_DNase.txt

# Merged features, no TEs
#TE_landscape/DNase/Refseq_features/refseq_features_DNase.txt
for file in refseq_*DNase.txt; do awk -v OFS='\t' -v feature=$(basename "$file" .txt) '{a[$15]+=$14}END{for(i in a){print i, a[i], feature}}' $file | sort >> refseq_features_DNase.txt; done
awk -v OFS='\t' '{a[$11]+=$3-$2}END{for(i in a){print i, a[i], "genome_noTE"}}' genome_noTE_DNase.txt | sort >> refseq_features_DNase.txt #Updated 8/7/17

# H3K27ac

# Merged TEs
#TE_landscape/H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt
awk -v OFS='\t' '{a[$15]+=$14}END{for(i in a){print i, a[i]}}' rmsk_TEother_merge_H3K27ac.txt | sort > rmsk_TEother_merge_H3K27ac_contribution.txt

# Number and width of H3K27ac peaks overall, overlapping TEs per sample; including summit rules
#TE_landscape/H3K27ac/H3K27ac_stats.txt	
while read line; do echo $line >> H3K27ac_stats.txt; wc -l H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak >> H3K27ac_stats.txt; awk '{sum+=$3-$2}END{print sum}' H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak >> H3K27ac_stats.txt; done < H3K27ac_samples.txt
while read line; do awk '{print $8, $9, $10}' H3K27ac_TEs/rmsk_TEother_$line\-H3K27ac.narrowPeak | sort | uniq | wc -l >> test; done < ../H3K27ac_samples.txt #Combine in Excel
while read line; do awk -v OFS='\t' '{print $8, $9, $10}' H3K27ac/summit/TEs/rmsk_TEother_$line\_H3K27ac_summit.txt | sort | uniq | wc -l >> H3K27ac/rmsk_TEother_H3K27ac_summit_stats.txt; done < sample_lists/H3K27ac_samples.txt

# Merged TE classes
#TE_landscape/H3K27ac/class_H3K27ac_sample.txt
awk -v OFS='\t' '{a[$4, $16]+=$15}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i];}}' rmsk_TEother_class_H3K27ac.txt > class_H3K27ac_sample.txt

# Merged TE subfamilies (unique peaks overlapping subfamily)
while read line; do while read line2; do awk -v sample=$line -v subfam=$line2 -v OFS='\t' '{if($4 == subfam) a[$8, $9, $10]+=1}END{print sample, subfam, length(a)}' H3K27ac/summit/TEs/rmsk_TEother_$line\_H3K27ac_summit.txt >> H3K27ac/subfamily_H3K27ac_sample_summit.txt; done < features/TEs/subfamily/subfamilies.txt; done < sample_lists/H3K27ac_samples.txt

# Merged features
awk -v OFS='\t' '{a[$1, $2]+=$16}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], a[i]}}' H3K27ac/refseq_features_H3K27ac_intersect.txt > H3K27ac/refseq_features_H3K27ac.txt

# Merged features, no TEs
#TE_landscape/H3K27ac/Refseq_features/refseq_features_H3K27ac.txt
for file in refseq_*H3K27ac.txt; do awk -v OFS='\t' -v feature=$(basename "$file" .txt) '{a[$15]+=$14}END{for(i in a){print i, a[i], feature}}' $file | sort >> refseq_features_H3K27ac.txt; done
awk -v OFS='\t' '{a[$11]+=$3-$2}END{for(i in a){print i, a[i], "genome_noTE"}}' genome_noTE_H3K27ac.txt | sort >> refseq_features_H3K27ac.txt #Update 8/7/17

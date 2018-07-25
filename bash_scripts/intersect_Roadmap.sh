# Intersect with Roadmap data
# 4/19/2016, 4/25/2016, 5/3/2016, 5/5/2016, 5/19/2016,
# 1/26/2017, 2/2/2017, 2/3/2017, 2/6/2017, 2/8/2017, 3/2/2017, 3/6/2017, 5/8/2017, 5/10/2017, 5/11/2017, 5/22/2017, 5/29/2017, 5/30/2017, 5/31/2017, 6/5/2017, 6/12/2017, 6/15/2017, 6/19/2017, 7/4/2017
# 8/2/2017, 8/4/2017, 8/5/2017, 8/7/2017, 8/18/2017, 8/25/2017, 8/28/2017, 8/29/2017

# chromHMM

# Individual TEs
#TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TE [127 files]
for file in chromHMM_bedfiles/*.bed ; do bedtools intersect -wo -a rmsk_TE.txt -b $file > $file\_TE; done
#TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_other [127 files]
for file in chromHMM_bedfiles/*.bed; do bedtools intersect -wo -a rmsk_other.txt -b $file > $file\_other; done

# Merged TEs
#TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TE_merge [127 files]
for file in chromHMM_bedfiles/E*_15_coreMarks_mnemonics.bed; do output=$(basename $file); bedtools intersect -wo -a rmsk_TE_merge.txt -b $file > chromHMM_TE_merge/$output\_TE_merge; done
#TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TEother_merge [127 files]
for file in chromHMM_bedfiles/E*_15_coreMarks_mnemonics.bed; do output=$(basename $file); bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file > chromHMM_other/$output\_TEother_merge; done

# TE merged classes
#TE_landscape/chromHMM/TEs/intersect/class/rmsk_[class].txt_chromHMM.bed [13 files]
for class in TE_classes/rmsk_*.txt; do for file in chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a $class -b $file | awk -v OFS='\t' -v tag=$suffix '{print $0, tag}' - >> $class\_chromHMM.bed; done; done
#Unconfident class
for file in ../chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a rmsk_Unconfident.txt -b $file | awk -v OFS='\t' -v tag=$suffix '{print $0, tag}' - >> rmsk_Unconfident.txt_chromHMM.bed; done
#TE_landscape/chromHMM/TEs/intersect/class/rmsk_Unconfident_RC.txt_chromHMM.bed
for file in ../../raw_data/chromHMM/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a ../../features/TEs/class/rmsk_Unconfident_RC.txt -b $file | awk -v OFS='\t' -v tag=$suffix '{print $0, tag}' - >> rmsk_Unconfident_RC.txt_chromHMM.bed; done

# TE merged subfamilies
#TE_landscape/chromHMM/subfamily/subfamily_state_sample.bed
for file in ../chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a TEother_subfamily_merge.txt -b $file | awk -v OFS='\t' -v tag=$suffix '{print $0, tag}' - >> subfamily_state_sample.bed; done

# Features
# Merged genic features
#TE_landscape/chromHMM/chromHMM_refseq_features.txt
for file in ~/genic_features/RefSeq/*_merge.txt; do feature=$(basename "$file" _merge.txt); while read line; do bedtools intersect -wo -a $file -b raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed | awk -v OFS='\t' -v sample=$line -v feature=$feature '{print sample, feature, $0}' - >> chromHMM/Refseq_features/chromHMM_features.txt; done < sample_lists/mnemonics.txt; done
while read line; do bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_intergenic.txt -b raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed | awk -v OFS='\t' -v sample=$line '{print sample, "intergenic", $0}' - >> chromHMM/Refseq_features/chromHMM_features.txt; done < sample_lists/mnemonics.txt

#TE_landscape/chromHMM/simple_repeats/chromHMM_repeats.bed
for file in chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a rmsk_repeats_merge.txt -b $file | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM_repeats.bed; done

# Merged genic features, no TEs
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_merge_noTE.txt [6 files]
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_nc_merge_noTE.txt [5 files]
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_pc_merge_noTE.txt [6 files]
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt); for file2 in raw_data/chromHMM/E*.bed; do suffix=$(basename $file2 | cut -d '_' -f1); bedtools intersect -wo -a $file -b $file2 | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM/Refseq_features/intersect/chromHMM_$feature\.txt; done; done
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_intergenic_noTE.txt
for file in raw_data/chromHMM/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_intergenic_noTE.txt -b $file | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM/Refseq_features/intersect/chromHMM_refseq_intergenic_noTE.txt; done
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_genome_noTE.bed
for file in raw_data/chromHMM/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools subtract -a $file -b features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM/Refseq_features/intersect/chromHMM_genome_noTE.bed; done

# Merged genic features, no TEs or repeats
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_intergenic_noTERepeats.bed
for file in chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a ../genic_features/refseq_intergenic_noTERepeats.txt -b $file | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM_features/chromHMM_intergenic_noTERepeats.bed; done

# Shuffled TEs
#TE_landscape/chromHMM/shuffled_TEs/chromHMM_rmsk_TEother_shuffle_*.txt [10 files]
for i in {1..10}; do for file2 in ~/TE_landscape/raw_data/chromHMM/E*.bed; do suffix=$(basename $file2 | cut -d '_' -f1); bedtools intersect -wo -a ../rmsk_TE_shuffle_$i\.txt -b $file2 | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM_rmsk_TE_shuffle_$i.txt; done; done

# Refseq promoters
#TE_landscape/chromHMM/Refseq_promoters/chromHMM_refseq_promoters_unique.txt
for file in ../../raw_data/chromHMM/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a refseq_promoters_unique.txt -b $file | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM_refseq_promoters_unique.txt; done

# Segwey promoters
#TE_landscape/chromHMM/Segway_promoters/intersect/E#_15_coreMarks_mnemonics_promoter.bed [127 files]
 for file in chromHMM_bedfiles/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a $file -b promoters.txt > chromHMM_promoter/$suffix\_15_coreMarks_mnemonics_promoter.bed; done

# WGBS (CpGs)

# TEs
# Individual TEs
#TE_landscape/WGBS/TE_CpG_Meth_new.bed
split -l 1000000 ~/TE_landscape/all_CpG_Meth.bed
for file in xa*; do echo $file; bedtools intersect -wo -a ~/TE_landscape/rmsk_TEother.txt -b $file >>  TE_CpG_Meth_new.bed ; done
#bedtools intersect -wao -a rmsk_TE.txt -b all_CpG_Meth.bed > TE_CpG_Meth.bed #Prints missing value for TEs without CpGs
#bedtools intersect -wao -a rmsk_other.txt -b all_CpG_Meth.bed > other_CpG_Meth.bed

# Merged TEs
#TE_landscape/WGBS/CpG_TE_Meth.bed
bedtools intersect -a all_CpG_Meth.bed -b rmsk_TEother_merge.txt > CpG_TE_Meth.bed &
# Reverse - no TEs
#TE_landscape/WGBS/CpG_noTE_Meth.bed
bedtools intersect -a all_CpG_Meth.bed -b rmsk_TEother_merge.txt -v > CpG_noTE_Meth.bed &

# Merged features
for file in ~/genic_features/RefSeq/*_merge.txt; do echo $file; feature=$(basename "$file" _merge.txt); bedtools intersect -a WGBS/all_CpG_Meth.bed -b $file | awk -v OFS='\t' -v feature=$feature '{print feature, $0}' >> WGBS/feature_CpG_Meth.bed; done
bedtools intersect -a WGBS/all_CpG_Meth.bed -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk -v OFS='\t' '{print "intergenic", $0}' >> WGBS/feature_CpG_Meth.bed

# Merged features, no TEs
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_merge_noTE_Meth.bed [6 files]
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_nc_merge_noTE_Meth.bed [5 files]
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_pc_merge_noTE_Meth.bed [6 files]
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt); bedtools intersect -a ../all_CpG_Meth.bed -b $file > CpG_$feature\_Meth.bed; done
#TE_landscape/WGBS/Refseq_features/CpG_refseq_intergenic_noTE_Meth.bed
bedtools intersect -a ../all_CpG_Meth.bed -b ~/genic_features/RefSeq/refseq_intergenic_noTE.txt > CpG_refseq_intergenic_noTE_Meth.bed

# Refseq promoters
#TE_landscape/WGBS/Refseq_promoters/refseq_promoter_unique_CpG_Meth.bed
bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_promoters_unique.txt -b ../all_CpG_Meth.bed > refseq_promoter_unique_CpG_Meth.bed

# Refseq exons
 split -l 1000000 ~/TE_landscape/WGBS/all_CpG_Meth.bed
 for file in x*; do echo $file; bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_exons_unique.txt -b $file >> ~/TE_landscape/WGBS/refseq_exons_unique_CpG_Meth.bed ; done

# Shuffled TEs
split -l 1000000 ~/TE_landscape/WGBS/all_CpG_Meth.bed
for i in {1..10}; do for file in x*; do bedtools intersect -wo -a ../rmsk_TE_shuffle_$i\.txt -b $file >> rmsk_TE_shuffle_$i\_Meth.bed; done; done

# DNase

# Individual TEs
#TE_landscape/DNase/intersect_sample_list_DNase_peak.sh
#TE_landscape/DNase/intersect/TEs/rmsk_TEother_E#-DNase.macs2.narrowPeak [53 files]

# Merged TEs
#TE_landscape/DNase/rmsk_TEother_merge_DNase.txt
while read line; do bedtools intersect -wo -a DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b ../rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_merge_DNase.txt; done < DNase_samples.txt

# Merged TE classes
#TE_landscape/DNase/rmsk_TEother_class_Dnase.txt
while read line; do bedtools intersect -wo -a ../TE_classes/TEother_class_merge.txt -b DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_class_Dnase.txt ; done < DNase_samples.txt
while read line; do grep "Unconfident_RC" ../features/TEs/class/TEother_class_merge.txt | bedtools intersect -wo -a - -b ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_class_Dnase.txt ; done < ../sample_lists/DNase_samples.txt

# Merged TE subfamilies
#TE_landscape/DNase/rmsk_TEother_subfamily_DNase.txt
while read line; do bedtools intersect -wo -a ../TE_subfamilies/TEother_subfamily_merge.txt -b DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_subfamily_DNase.txt ; done < DNase_samples.txt

# Features
for file in ~/genic_features/RefSeq/*_merge.txt; do feature=$(basename "$file" _merge.txt); echo $feature; while read line; do echo $line; bedtools intersect -wo -a raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b $file | awk -v OFS='\t' -v sample=$line -v feature=$feature '{print sample, feature, $0}' - >> DNase/refseq_features_DNase_intersect.txt; done < sample_lists/DNase_samples.txt; done
while read line; do echo $line; bedtools intersect -wo -a raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk -v OFS='\t' -v sample=$line '{print sample, "intergenic", $0}' - >> DNase/refseq_features_DNase_intersect.txt; done < sample_lists/DNase_samples.txt

# Features, no TEs
#TE_landscape/DNase/Refseq_features/refseq_[feature]_merge_noTE_DNase.txt [6 files]
#TE_landscape/DNase/Refseq_features/refseq_[feature]_nc_merge_noTE_DNase.txt [5 files]
#TE_landscape/DNase/Refseq_features/refseq_[feature]_pc_merge_noTE_DNase.txt [6 files]
#TE_landscape/DNase/Refseq_features/refseq_intergenic_noTE_DNase.txt
#TE_landscape/DNase/Refseq_features/genome_noTE_DNase.txt
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt); echo $feature; while read line; do echo $line; bedtools intersect -wo -a ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b $file | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> $feature\_DNase.txt; done < ../sample_lists/DNase_samples.txt; done
while read line; do bedtools intersect -wo -a ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b ~/genic_features/RefSeq/refseq_intergenic_noTE.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> refseq_intergenic_noTE_DNase.txt; done < ../sample_lists/DNase_samples.txt
while read line; do bedtools subtract -a ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b ../features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> genome_noTE_DNase.txt; done < ../sample_lists/DNase_samples.txt

# Refseq promoters
#TE_landscape/DNase/Refseq_promoters/refseq_promoter_unique_E*-DNase.macs2.narrowPeak [53 files]
 while read line; do bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_promoters_unique.txt -b ../../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak > refseq_promoter_unique_$line\-DNase.macs2.narrowPeak; done < ../../sample_lists/DNase_samples.txt

# Shuffled TEs
#TE_landscape/features/shuffled_TEs/run_intersect.sh

# H3K27ac

# Individual TEs
#TE_landscape/H3K27ac/intersect/TEs/rmsk_TEother_E#-H3K27ac.narrowPeak [98 files]
while read line; do bedtools intersect -wo -a ../rmsk_TEother.txt -b H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak > H3K27ac_TEs/rmsk_TEother_$line\-H3K27ac.narrowPeak; done < H3K27ac_samples.txt

# Merged TEs
#TE_landscape/H3K27ac/rmsk_TEother_merge_H3K27ac.txt
while read line; do bedtools intersect -wo -a H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b ../rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_merge_H3K27ac.txt; done < H3K27ac_samples.txt

# Merged TE classes
#TE_landscape/H3K27ac/rmsk_TEother_class_H3K27ac.txt
 while read line; do bedtools intersect -wo -a ../features/TEs/class/TEother_class_merge.txt -b ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_class_H3K27ac.txt ; done < ../sample_lists/H3K27ac_samples.txt

# Merged TE subfamilies
#TE_landscape/H3K27ac/rmsk_TEother_subfamily_H3K27ac.txt
while read line; do bedtools intersect -wo -a ../TE_subfamilies/TEother_subfamily_merge.txt -b H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_subfamily_H3K27ac.txt ; done < H3K27ac_samples.txt

# Features 
for file in ~/genic_features/RefSeq/*_merge.txt; do feature=$(basename "$file" _merge.txt); echo $feature; while read line; do echo $line; bedtools intersect -wo -a raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b $file | awk -v OFS='\t' -v sample=$line -v feature=$feature '{print sample, feature, $0}' - >> H3K27ac/refseq_features_H3K27ac_intersect.txt; done < sample_lists/H3K27ac_samples.txt; done
while read line; do echo $line; bedtools intersect -wo -a raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b ~/genic_features/RefSeq/refseq_intergenic.txt | awk -v OFS='\t' -v sample=$line '{print sample, "intergenic", $0}' - >> H3K27ac/refseq_features_H3K27ac_intersect.txt; done < sample_lists/H3K27ac_samples.txt

# Features, no TEs
#TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_merge_noTE_H3K27ac.txt [6 files]
#TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_nc_merge_noTE_H3K27ac.txt [5 files]
#TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_pc_merge_noTE_H3K27ac.txt [6 files]
#TE_landscape/H3K27ac/Refseq_features/refseq_intergenic_noTE_H3K27ac.txt
#TE_landscape/H3K27ac/Refseq_features/genome_noTE_H3K27ac.txt
 for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt); echo $feature; while read line; do echo $line; bedtools intersect -wo -a ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b $file | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> $feature\_H3K27ac.txt; done < ../sample_lists/H3K27ac_samples.txt; done
 while read line; do bedtools intersect -wo -a ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b ~/genic_features/RefSeq/refseq_intergenic_noTE.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> refseq_intergenic_noTE_H3K27ac.txt; done < ../sample_lists/H3K27ac_samples.txt
 while read line; do bedtools subtract -a ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b ../features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> genome_noTE_H3K27ac.txt; done < ../sample_lists/H3K27ac_samples.txt

# Refseq promoters
#TE_landscape/H3K27ac/Refseq_promoters/refseq_promoter_unique_E*-H3K27ac.narrowPeak [98 files]
 while read line; do bedtools intersect -wo -a ~/genic_features/RefSeq/refseq_promoters_unique.txt -b ../../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak > refseq_promoter_unique_$line\-H3K27ac.narrowPeak; done < ../../sample_lists/H3K27ac_samples.txt

# Shuffled TEs
#TE_landscape/features/shuffled_TEs/run_intersect.sh

# RNA-seq
# Intersect TEs with RNA-seq bedfiles (run by Daofeng on server)
#TE_landscape/RNAseq/intersect_sample_list_RNA_raw.sh
#TE_landscape/RNAseq/intersect_sample_list_RNA_raw_ag.sh
#TE_landscape/RNAseq/htcf_ecp/refseq_exons.txt.sorted
#TE_landscape/RNAseq/htcf_ecp/rmsk_TEother.txt.sorted
#TE_landscape/RNAseq/htcf_ecp/run.sh
#TE_landscape/RNAseq/htcf_ecp/samples.txt
#TE_landscape/RNAseq/htcf_ecp/samples2.txt

# Mouse

# chromHMM
# mm9 TEs
#TE_landscape/Mouse/chromHMM/intersect/TEs/ENCFF#.bed_TE [15 files]
for file in mouse_chromHMM/*.bed ; do bedtools intersect -wo -a $file -b mm9_rmsk_TE.txt > $file\_TE; done

# WGBS (mm10)
#TE_landscape/Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_ENCFF#.bed [9 files]
#TE_landscape/Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_WGBS.bed
awk -v OFS='\t' '{print $1, $2, $3, $10, $11}' $file | bedtools intersect -wo -a mm10_rmsk_TE.txt -b - > mm10_rmsk_TE_$file;

# DNase
#TE_landscape/Mouse/DNase_mm10/mm10_orthologs_DNase.txt
for file in ENCFF*.bed; do bedtools intersect -wo -a mm10_hg19_TE_intersect_same.bed -b $file | awk -v OFS='\t' -v sample=$file '{print $0, sample}' - >> mm10_orthologs_DNase.txt; done

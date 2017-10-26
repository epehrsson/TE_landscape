# Intersect with Roadmap data
# 4/19/2016, 4/25/2016, 5/3/2016, 5/5/2016, 5/19/2016,
# 1/26/2017, 2/2/2017, 2/3/2017, 2/6/2017, 2/8/2017, 3/2/2017, 3/6/2017, 5/8/2017, 5/10/2017, 5/11/2017, 5/22/2017, 5/29/2017, 5/30/2017, 5/31/2017, 6/5/2017, 6/12/2017, 6/15/2017, 6/19/2017, 7/4/2017
# 8/2/2017, 8/4/2017, 8/5/2017, 8/7/2017, 8/18/2017, 8/25/2017, 8/28/2017, 8/29/2017

# chromHMM 
for file in raw_data/chromHMM/E*_15_coreMarks_mnemonics.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools intersect -wo -a $input -b $file | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> $output; done
# DNase
#TE_landscape/DNase/intersect_sample_list_DNase_peak.sh
while read line; do bedtools intersect -wo -a $input -b ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> $output ; done < ../sample_lists/DNase_samples.txt
# H3K27ac
while read line; do bedtools intersect -wo -a $input -b ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> $output ; done < ../sample_lists/H3K27ac_samples.txt
# WGBS (CpGs)
split -l 1000000 ~/TE_landscape/all_CpG_Meth.bed
for file in xa*; do echo $file; bedtools intersect -wo -a $input -b $file >> $output ; done

# TEs
# Individual TEs
rmsk_TE.txt	TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TE [127 files]
rmsk_other.txt	TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_other [127 files]
	TE_landscape/DNase/intersect/TEs/rmsk_TEother_E#-DNase.macs2.narrowPeak [53 files]
../rmsk_TEother.txt	TE_landscape/H3K27ac/intersect/TEs/rmsk_TEother_E#-H3K27ac.narrowPeak [98 files]
~/TE_landscape/rmsk_TEother.txt	TE_landscape/WGBS/TE_CpG_Meth_new.bed

# Merged TEs
rmsk_TE_merge.txt	TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TE_merge [127 files]
rmsk_TEother_merge.txt	TE_landscape/chromHMM/TEs/intersect/E#_15_coreMarks_mnemonics.bed_TEother_merge [127 files]
../rmsk_TEother_merge.txt	TE_landscape/DNase/rmsk_TEother_merge_DNase.txt
../rmsk_TEother_merge.txt	TE_landscape/H3K27ac/rmsk_TEother_merge_H3K27ac.txt
rmsk_TEother_merge.txt	TE_landscape/WGBS/CpG_TE_Meth.bed #Not -wo

# Merged TE classes
for class in TE_classes/rmsk_*.txt	TE_landscape/chromHMM/TEs/intersect/class/rmsk_[class].txt_chromHMM.bed [13 files]
../../features/TEs/class/rmsk_Unconfident_RC.txt	TE_landscape/chromHMM/TEs/intersect/class/rmsk_Unconfident_RC.txt_chromHMM.bed
../TE_classes/TEother_class_merge.txt	TE_landscape/DNase/rmsk_TEother_class_Dnase.txt
grep "Unconfident_RC" ../features/TEs/class/TEother_class_merge.txt	
../features/TEs/class/TEother_class_merge.txt	TE_landscape/H3K27ac/rmsk_TEother_class_H3K27ac.txt

# Merged TE subfamilies
TEother_subfamily_merge.txt	TE_landscape/chromHMM/subfamily/subfamily_state_sample.bed
../TE_subfamilies/TEother_subfamily_merge.txt	TE_landscape/DNase/rmsk_TEother_subfamily_DNase.txt
../TE_subfamilies/TEother_subfamily_merge.txt	TE_landscape/H3K27ac/rmsk_TEother_subfamily_H3K27ac.txt

# Features
# Merged genic features
../genic_features/refseq_coding_exon_merge.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_CDS.bed
../genic_features/refseq_3UTR_merge.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_3UTR.bed
../genic_features/refseq_5UTR_merge.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_5UTR.bed
../genic_features/refseq_exons_merge.txt	
../genic_features/refseq_introns_merge.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_intron.bed
../genic_features/refseq_promoters_merge.txt	
	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_intergenic.txt
rmsk_repeats_merge.txt	TE_landscape/chromHMM/simple_repeats/chromHMM_repeats.bed

# Merged genic features, no TEs
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt)	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_merge_noTE.txt [6 files]
	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_nc_merge_noTE.txt [5 files]
	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_[feature]_pc_merge_noTE.txt [6 files]
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; feature=$(basename "$file" .txt)	TE_landscape/DNase/Refseq_features/refseq_[feature]_merge_noTE_DNase.txt [6 files]
	TE_landscape/DNase/Refseq_features/refseq_[feature]_nc_merge_noTE_DNase.txt [5 files]
	TE_landscape/DNase/Refseq_features/refseq_[feature]_pc_merge_noTE_DNase.txt [6 files]
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt)	TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_merge_noTE_H3K27ac.txt [6 files]
	TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_nc_merge_noTE_H3K27ac.txt [5 files]
	TE_landscape/H3K27ac/Refseq_features/refseq_[feature]_pc_merge_noTE_H3K27ac.txt [6 files]
for file in ~/genic_features/RefSeq/*_merge_noTE.txt; do feature=$(basename "$file" .txt)	TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_merge_noTE_Meth.bed [6 files] #Not -wo
	TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_nc_merge_noTE_Meth.bed [5 files]
	TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_pc_merge_noTE_Meth.bed [6 files]
~/genic_features/RefSeq/refseq_intergenic_noTE.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_refseq_intergenic_noTE.txt
~/genic_features/RefSeq/refseq_intergenic_noTE.txt	TE_landscape/DNase/Refseq_features/refseq_intergenic_noTE_DNase.txt
~/genic_features/RefSeq/refseq_intergenic_noTE.txt	TE_landscape/H3K27ac/Refseq_features/refseq_intergenic_noTE_H3K27ac.txt
~/genic_features/RefSeq/refseq_intergenic_noTE.txt	TE_landscape/WGBS/Refseq_features/CpG_refseq_intergenic_noTE_Meth.bed #Not -wo

# Merged genic features, no TEs or repeats
../genic_features/refseq_intergenic_noTERepeats.txt	TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_intergenic_noTERepeats.bed

# Shuffled TEs
for i in {1..10}; ../rmsk_TE_shuffle_$i\.txt	TE_landscape/chromHMM/shuffled_TEs/chromHMM_rmsk_TEother_shuffle_*.txt [10 files]
for i in {1..10}; rmsk_TE_shuffle_$i\.txt
for i in {1..10}; rmsk_TE_shuffle_$i\.txt	H3K27ac/rmsk_TE_shuffle_$i\_$line\-H3K27ac.narrowPeak
for i in {1..10}; ../rmsk_TE_shuffle_$i\.txt	rmsk_TE_shuffle_$i\_Meth.bed

# Refseq promoters
refseq_promoters_unique.txt	TE_landscape/chromHMM/Refseq_promoters/chromHMM_refseq_promoters_unique.txt
~/genic_features/RefSeq/refseq_promoters_unique.txt	TE_landscape/DNase/Refseq_promoters/refseq_promoter_unique_E*-DNase.macs2.narrowPeak [53 files]
~/genic_features/RefSeq/refseq_promoters_unique.txt	TE_landscape/H3K27ac/Refseq_promoters/refseq_promoter_unique_E*-H3K27ac.narrowPeak [98 files]
~/genic_features/RefSeq/refseq_promoters_unique.txt	TE_landscape/WGBS/Refseq_promoters/refseq_promoter_unique_CpG_Meth.bed

# Segwey promoters
promoters.txt	TE_landscape/chromHMM/Segway_promoters/intersect/E#_15_coreMarks_mnemonics_promoter.bed [127 files]

# Genome, no TEs
#TE_landscape/chromHMM/Refseq_features/intersect/chromHMM_genome_noTE.bed
for file in raw_data/chromHMM/E*.bed; do suffix=$(basename $file | cut -d '_' -f1); bedtools subtract -a $file -b features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$suffix '{print $0, sample}' - >> chromHMM/Refseq_features/intersect/chromHMM_genome_noTE.bed; done
#TE_landscape/DNase/Refseq_features/genome_noTE_DNase.txt
while read line; do bedtools subtract -a ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak -b ../features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> genome_noTE_DNase.txt; done < ../sample_lists/DNase_samples.txt
#TE_landscape/H3K27ac/Refseq_features/genome_noTE_H3K27ac.txt
 while read line; do bedtools subtract -a ../raw_data/H3K27ac/H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak -b ../features/TEs/merge/rmsk_TEother_merge.txt | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> genome_noTE_H3K27ac.txt; done < ../sample_lists/H3K27ac_samples.txt
  bedtools intersect -a all_CpG_Meth.bed -b rmsk_TEother_merge.txt -v > CpG_noTE_Meth.bed &


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

# mm9 merged TEs
#TE_landscape/Mouse/chromHMM/intersect/TEs/ENCFF#.bed_TEmerge [15 files]
for file in mouse_chromHMM/*.bed ; do bedtools intersect -wo -a $file -b mm9_rmsk_TEmerge.txt > $file\_TEmerge; done
#TE_landscape/Mouse/chromHMM/intersect/TEs/ENCFF#.bed_TEother_merge [15 files]
for file in mouse_chromHMM/*.bed ; do output=$(basename $file); bedtools intersect -wo -a mm9_rmsk_TEother_merge.txt -b $file > mouse_chromHMM_other/$output\_TEother_merge; done

# WGBS (mm10)
#TE_landscape/Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_ENCFF#.bed [9 files]
#TE_landscape/Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_WGBS.bed
awk -v OFS='\t' '{print $1, $2, $3, $10, $11}' $file | bedtools intersect -wo -a mm10_rmsk_TE.txt -b - > mm10_rmsk_TE_$file;

# DNase
#TE_landscape/Mouse/DNase_mm10/mm10_orthologs_DNase.txt
for file in ENCFF*.bed; do bedtools intersect -wo -a mm10_hg19_TE_intersect_same.bed -b $file | awk -v OFS='\t' -v sample=$file '{print $0, sample}' - >> mm10_orthologs_DNase.txt; done

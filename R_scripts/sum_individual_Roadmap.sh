# Sum overlap by individual feature
# 4/20/2016, 4/26/2016, 5/3/2016, 7/6/2016, 1/26/2017, 2/2/2017, 5/24/2017, 7/4/2017, 8/7/2017, 8/23/2017, 8/24/2017, 9/12/2017

# chromHMM

# TEs
#TE_landscape/chromHMM/all_chromHMM_TE.txt
for file in chromHMM_bedfiles/*.bed_TE; do suffix=$(basename $file | cut -d '_' -f1); awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7, $11]+=$12;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], a[i];}}' $file | awk -v x=$suffix 'BEGIN{OFS="\t";}{print $0, x}' - ; done >> all_chromHMM_TE.txt
#TE_landscape/chromHMM/all_chromHMM_TE_sorted.txt
sort -k1,1V -k2,2n -k3,3n -k6,6 -k10,10 all_chromHMM_TE.txt > all_chromHMM_TE_sorted.txt
#TE_landscape/chromHMM/all_chromHMM_other.txt
for file in chromHMM_bedfiles/*.bed_other; do suffix=$(basename $file | cut -d '_' -f1); awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7, $11]+=$12;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], a[i];}}' $file | awk -v x=$suffix 'BEGIN{OFS="\t";}{print $0, x}' - ; done >> all_chromHMM_other.txt
#TE_landscape/chromHMM/all_chromHMM_other_sorted.txt
sort -k1,1V -k2,2n -k3,3n -k6,6 -k10,10 all_chromHMM_other.txt > all_chromHMM_other_sorted.txt

# Number of bases in each TE block in each state in each sample (old!)
#TE_landscape/chromHMM/all_chromHMM_TE_merge.txt	
for file in chromHMM_TE_merge/* ; do suffix=$(basename $file | cut -d '_' -f1); awk -v OFS='\t' '{a[$1, $2, $3, $7]+=$8;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i];}}' $file | awk -v x=$suffix 'BEGIN{OFS="\t";}{print $0, x}' - >> all_chromHMM_TE_merge.txt; done 

# TEs by class
#TE_landscape/chromHMM/all_chromHMM_[class]_sorted.txt [6 files]
while read line ; do awk -v OFS='\t' -v class=$line '{if($5 == class) print $0}' all_chromHMM_TE_sorted.txt > all_chromHMM_$line\_sorted.txt; done < ../features/TEs/class/TE_class.txt
awk -v OFS='\t' '{if($5 == "Other") print $0}' all_chromHMM_other_sorted.txt > all_chromHMM_SVA_sorted.txt
awk -v OFS='\t' '{if($5 != "Other") print $0}' all_chromHMM_other_sorted.txt > all_chromHMM_Other_sorted.txt

# Refseq promoters
#TE_landscape/chromHMM/Refseq_promoters/chromHMM_refseq_promoters_unique_sorted.txt
awk -v OFS='\t' '{a[$1, $2, $3, $4, $8, $10]+=$9;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], a[i];}}' chromHMM_refseq_promoters_unique.txt | sort -k1,1V -k2,2n -k3,3n -k6,6 > chromHMM_refseq_promoters_unique_sorted.txt

# Shuffled TEs
#/scratch/ecp/shuffled/rmsk_TEother_shuffle_#_sorted.txt [10 files]
for i in {1..10}; do awk '{print>$13}' chromHMM_rmsk_TE_shuffle_$i\.txt; for file in E*; do awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7, $11, $13]+=$12;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[9], sep[8], a[i];}}' $file >> rmsk_TE_shuffle_$i\_sum.txt; done; rm E*; sort -k1,1V -k2,2n -k3,3n -k6,6 -k4,4 -k8,8 rmsk_TE_shuffle_$i\_sum.txt > rmsk_TE_shuffle_$i\_sorted.txt; rm rmsk_TE_shuffle_$i\_sum.txt; done

# Segwey promoters
#TE_landscape/chromHMM/Segway_promoters/all_chromHMM_promoter.txt
for file in chromHMM_promoter/*.bed; do suffix=$(basename $file | cut -d '_' -f1); awk -v OFS='\t' '{a[$4, $5, $6, $7]+=$8;}END{for(i in a) {split (i, sep, SUBSEP); print sep[2], sep[3], sep[4], sep[1], a[i];}}' $file | awk -v x=$suffix 'BEGIN{OFS="\t";}{print $0, x}' - ; done >> all_chromHMM_promoter.txt

# DNase

# TEs
#TE_landscape/DNase/rmsk_TEother_DNase_peaks.txt
#TE_landscape/DNase/rmsk_TEother_DNase_peaks_filter.txt

# Refseq promoters
#TE_landscape/DNase/Refseq_promoters/refseq_promoters_unique_DNase_peaks.txt
python ~/bin/TE_landscape/DNase_peaks_promoter.py ~/genic_features/RefSeq/refseq_promoters_unique_std.txt ../../sample_lists/DNase_samples.txt refseq_promoters_unique_DNase_peaks.txt

# Shuffled TEs
#TE_landscape/DNase/shuffled/rmsk_TE_shuffle_#_DNase_peaks.txt [10 files]
for i in {1..10}; do python ~/bin/TE_landscape/DNase_peaks.py rmsk_TE_shuffle_$i\.txt ~/TE_landscape/sample_lists/DNase_samples.txt DNase/rmsk_TE_shuffle_$i\_ DNase/rmsk_TE_shuffle_$i\_DNase_peaks.txt; done

# H3K27ac

# TEs
#TE_landscape/H3K27ac/rmsk_TEother_H3K27ac_peaks.txt
#TE_landscape/H3K27ac/rmsk_TEother_H3K27ac_peaks_filter.txt
python ~/bin/TE_landscape/H3K27ac_peaks.py ../rmsk_TEother.txt H3K27ac_samples.txt rmsk_TEother_H3K27ac_peaks.txt
awk -v OFS='\t' '{if($9 > 0) print $0}' rmsk_TEother_H3K27ac_peaks.txt > rmsk_TEother_H3K27ac_peaks_filter.txt

# Refseq promoters
#TE_landscape/H3K27ac/Refseq_promoters/refseq_promoters_unique_H3K27ac_peaks.txt
python ~/bin/TE_landscape/H3K27ac_peaks_promoter.py ~/genic_features/RefSeq/refseq_promoters_unique_std.txt ../../sample_lists/H3K27ac_samples.txt refseq_promoters_unique_H3K27ac_peaks.txt

# Shuffled TEs
#TE_landscape/H3K27ac/shuffled/rmsk_TE_shuffle_#_H3K27ac_peaks.txt [10 files]
for i in {1..10}; do cat rmsk_TE_shuffle_$i\_E*-H3K27ac.narrowPeak | cut -f1-7 | sort | uniq > test; python ~/bin/TE_landscape/H3K27ac_peaks.py test ~/TE_landscape/sample_lists/H3K27ac_samples.txt rmsk_TE_shuffle_$i\_ rmsk_TE_shuffle_$i\_H3K27ac_peaks.txt; done

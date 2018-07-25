# Sum overlap by individual feature
# 4/20/2016, 4/26/2016, 5/3/2016, 7/6/2016, 1/26/2017, 2/2/2017, 5/24/2017, 7/4/2017, 8/7/2017, 8/23/2017, 8/24/2017, 9/12/2017, 9/14/2017
# Updated 4/24/18 with summit rules

# chromHMM

# TEs
# Filtering to summit
while read line; do python ~/bin/TE_landscape/filter_summit.py chromHMM/TEs/intersect/$line\_15_coreMarks_mnemonics.bed_TE 1 8 chromHMM/summit/TEs/rmsk_TE_$line\_chromHMM.bed chromHMM; python ~/bin/TE_landscape/filter_summit.py chromHMM/TEs/intersect/$line\_15_coreMarks_mnemonics.bed_other 1 8 chromHMM/summit/TEs/rmsk_other_$line\_chromHMM.bed chromHMM; cat chromHMM/summit/TEs/rmsk_TE_$line\_chromHMM.bed chromHMM/summit/TEs/rmsk_other_$line\_chromHMM.bed > chromHMM/summit/TEs/rmsk_TEother_$line\_chromHMM.bed; rm chromHMM/summit/TEs/rmsk_TE_$line\_chromHMM.bed; rm chromHMM/summit/TEs/rmsk_other_$line\_chromHMM.bed; done < sample_lists/mnemonics.txt

# Including those with no summit
while read line; do python identify_no_summit.py $line rmsk_TEother_$line\_chromHMM.bed rmsk_TEother_$line\_chromHMM_noSummit.txt 7 8; done < mnemonics.txt

# Matrix of TE x sample x state
while read line; do awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4, $5, $6, $7, $11]+=$12}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sample, a[i], sep[8], "summit"}}' rmsk_TEother_$line\_chromHMM.bed >> rmsk_TEother_chromHMM.txt; awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4, $5, $6, $7, $8]+=$9}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sample, a[i], sep[8], "majority"}}' rmsk_TEother_$line\_chromHMM_noSummit.txt >> rmsk_TEother_chromHMM.txt; done < mnemonics.txt
sort -k1,1V -k2,2n -k3,3n -k4,4 -k8,8 rmsk_TEother_chromHMM.txt > ~/TE_landscape/chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt

# QC
## Number of states per TE overlapping summit
for file in chromHMM/summit/TEs/rmsk_TEother_E*_chromHMM.bed; do cut -f1-7,11 $file | sort | uniq | awk '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){print a[i]}}' - | awk '{b[$1]+=1}END{for(i in b){print i, b[i]}}' -; done

## TEs not completely covered with summit rules
python ~/bin/TE_landscape/check_missing_bp.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt rmsk_TEother_total_chromHMM.txt

## Missing bp per sample
awk '{if($8 == "summit") a[$9]+=($3-$2-$10)}END{for(i in a){print i, a[i]}}' rmsk_TEother_total_chromHMM.txt
awk '{if($8 == "majority") a[$9]+=($3-$2-$10)}END{for(i in a){print i, a[i]}}' rmsk_TEother_total_chromHMM.txt

## Number of TEs missing bp per sample
awk '{if($8 == "summit") a[$9]+=1}END{for(i in a){print i, a[i]}}' rmsk_TEother_total_chromHMM.txt
awk '{if($8 == "majority") a[$9]+=1}END{for(i in a){print i, a[i]}}' rmsk_TEother_total_chromHMM.txt

# TEs by class
while read line ; do awk -v OFS='\t' -v class=$line '{if($5 == class) print $0}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt > chromHMM/chromHMM_summit_$line\.txt; done < features/TEs/class/TE_class.txt
awk -v OFS='\t' '{if($5 == "Other") print $0}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt > chromHMM/chromHMM_summit_SVA.txt
awk -v OFS='\t' '{if(($5 != "Other") && ($5 != "DNA") && ($5 != "LINE") && ($5 != "LTR") && ($5 != "SINE")) print $0}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt > chromHMM/chromHMM_summit_Other.txt

# TEs by subfamily and state
awk '{if($10 != "8_ZNF/Rpts") print > "chromHMM/subfamily/by_state/"$4"_"$10".txt"; else print > "chromHMM/subfamily/by_state/"$4"_8_ZNF.Rpts.txt"}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt

# Refseq promoters
# Filtering to summit
python ~/bin/TE_landscape/filter_summit.py chromHMM/Refseq_promoters/chromHMM_refseq_promoters_unique.txt 1 5 chromHMM/summit/promoters/chromHMM_refseq_promoters_unique_summit.txt chromHMM
awk -v OFS='\t' '{a[$1, $2, $3, $4, $8, $10]+=$9}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[6], a[i], sep[5]}}' chromHMM/summit/promoters/chromHMM_refseq_promoters_unique_summit.txt | sort -k1,1V -k2,2n -k3,3n -k5,5 - > chromHMM/refseq_promoters_unique_chromHMM_summit_sorted.txt

# Shuffled TEs
#/scratch/ecp/shuffled/rmsk_TEother_shuffle_#_sorted.txt [10 files]
for i in {1..10}; do awk '{print>$13}' chromHMM_rmsk_TE_shuffle_$i\.txt; for file in E*; do awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7, $11, $13]+=$12;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[9], sep[8], a[i];}}' $file >> rmsk_TE_shuffle_$i\_sum.txt; done; rm E*; sort -k1,1V -k2,2n -k3,3n -k6,6 -k4,4 -k8,8 rmsk_TE_shuffle_$i\_sum.txt > rmsk_TE_shuffle_$i\_sorted.txt; rm rmsk_TE_shuffle_$i\_sum.txt; done

# Segwey promoters
#TE_landscape/chromHMM/Segway_promoters/all_chromHMM_promoter.txt
for file in chromHMM_promoter/*.bed; do suffix=$(basename $file | cut -d '_' -f1); awk -v OFS='\t' '{a[$4, $5, $6, $7]+=$8;}END{for(i in a) {split (i, sep, SUBSEP); print sep[2], sep[3], sep[4], sep[1], a[i];}}' $file | awk -v x=$suffix 'BEGIN{OFS="\t";}{print $0, x}' - ; done >> all_chromHMM_promoter.txt

# DNase

# TEs
# Filtering to summit
for file in DNase/intersect/TEs/rmsk_TEother_E*-DNase.macs2.narrowPeak; do python ~/bin/TE_landscape/filter_summit.py $file 1 8 DNase/summit/TEs/$( basename $file -DNase.macs2.narrowPeak)_DNase_summit.txt; done

# Matrix of TE x sample x state
while read line; do awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sample, a[i]}}' DNase/summit/TEs/rmsk_TEother_$line\_DNase_summit.txt >> DNase/rmsk_TEother_DNase_summit.txt; done < sample_lists/DNase_samples.txt

# TEs by subfamily and state
awk '{print $0 > "DNase/subfamily/"$4"_DNase.txt"}' DNase/rmsk_TEother_DNase_summit.txt

# Refseq promoters
for file in DNase/Refseq_promoters/refseq_promoter_unique_E*-DNase.macs2.narrowPeak; do python ~/bin/TE_landscape/filter_summit.py $file 1 5 DNase/summit/promoters/$( basename $file -DNase.macs2.narrowPeak)_DNase_summit.txt; done

while read line; do awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sample, a[i]}}' DNase/summit/promoters/refseq_promoter_unique_$line\_DNase_summit.txt >> DNase/refseq_promoter_unique_DNase_summit.txt; done < sample_lists/DNase_samples.txt

# Shuffled TEs
#TE_landscape/DNase/shuffled/rmsk_TE_shuffle_#_DNase_peaks.txt [10 files]
for i in {1..10}; do python ~/bin/TE_landscape/DNase_peaks.py rmsk_TE_shuffle_$i\.txt ~/TE_landscape/sample_lists/DNase_samples.txt DNase/rmsk_TE_shuffle_$i\_ DNase/rmsk_TE_shuffle_$i\_DNase_peaks.txt; done

# H3K27ac

# TEs
# Filtering to summit
for file in H3K27ac/intersect/TEs/rmsk_TEother_E*-H3K27ac.narrowPeak; do python ~/bin/TE_landscape/filter_summit.py $file 1 8 H3K27ac/summit/TEs/$( basename $file -H3K27ac.narrowPeak)_H3K27ac_summit.txt; done

# Matrix of TE x sample x state
while read line; do awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sample, a[i]}}' H3K27ac/summit/TEs/rmsk_TEother_$line\_H3K27ac_summit.txt >> H3K27ac/rmsk_TEother_H3K27ac_summit.txt; done < sample_lists/H3K27ac_samples.txt

# TEs by subfamily and state
awk '{print $0 > "H3K27ac/subfamily/"$4"_H3K27ac.txt"}' H3K27ac/rmsk_TEother_H3K27ac_summit.txt

# Refseq promoters
for file in H3K27ac/Refseq_promoters/refseq_promoter_unique_E*-H3K27ac.narrowPeak; do python ~/bin/TE_landscape/filter_summit.py $file 1 5 H3K27ac/summit/promoters/$( basename $file -H3K27ac.narrowPeak)_H3K27ac_summit.txt; done

while read line; do awk -v OFS='\t' -v sample=$line '{a[$1, $2, $3, $4]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sample, a[i]}}' H3K27ac/summit/promoters/refseq_promoter_unique_$line\_H3K27ac_summit.txt >> H3K27ac/refseq_promoter_unique_H3K27ac_summit.txt; done < sample_lists/H3K27ac_samples.txt

# Shuffled TEs
#TE_landscape/H3K27ac/shuffled/rmsk_TE_shuffle_#_H3K27ac_peaks.txt [10 files]
for i in {1..10}; do cat rmsk_TE_shuffle_$i\_E*-H3K27ac.narrowPeak | cut -f1-7 | sort | uniq > test; python ~/bin/TE_landscape/H3K27ac_peaks.py test ~/TE_landscape/sample_lists/H3K27ac_samples.txt rmsk_TE_shuffle_$i\_ rmsk_TE_shuffle_$i\_H3K27ac_peaks.txt; done

# Mouse
# chromHMM state for each mm9 TE
#TE_landscape/Mouse/chromHMM/mouse_mm9_chromHMM_TE.txt
while read line ; do awk -v OFS='\t' -v sample=$line '{a[$10, $11, $12, $13, $14, $15, $16, $4, sample]+=$17;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], sep[9], a[i];}}'  mouse_chromHMM_TE/$line\.bed_TE >> mouse_mm9_chromHMM_TE.txt; done < mouse_samples.txt

# DNase
# Total overlap with Dnase per TE
#TE_landscape/Mouse/DNase_mm10/mm10_orthologs_DNase_sum.txt
awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7, $19]+=$18; b[$1, $2, $3, $4, $5, $6, $7, $19]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], a[i], b[i];}}' mm10_orthologs_DNase.txt > mm10_orthologs_DNase_sum.txt

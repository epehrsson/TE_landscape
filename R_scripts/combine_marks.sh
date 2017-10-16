# Combine multiple epigenetic marks
# 1/30/2017, 1/31/2017, 2/9/2017, 5/15/2017, 5/16/2017, 6/6/2017, 6/7/2017, 7/4/2017, 7/20/2017, 7/21/2017, 7/27/2017, 7/28/2017

# Combined file with each TE x sample and chromHMM, WGBS, Dnase, and H3K27ac state	
#TE_landscape/compare_marks/TE_combine_marks.txt	
ln -s /bar/epehrsson/TE_landscape/chromHMM/all_chromHMM_TE_sorted.txt
ln -s ~/TE_landscape/chromHMM/all_chromHMM_other_sorted.txt chromHMM/.
ln -s /bar/epehrsson/TE_landscape/DNase/rmsk_TEother_DNase_peaks.txt DNase/.
ln -s /bar/epehrsson/TE_landscape/WGBS/TE_WGBS_state.txt WGBS/.
ln -s ~/TE_landscape/H3K27ac/rmsk_TEother_H3K27ac_peaks.txt

awk '{print>$10}' all_chromHMM_TE_sorted.txt &
awk '{print>$8}' TE_WGBS_state.txt &
awk '{print>$8}' rmsk_TEother_H3K27ac_peaks.txt &
awk '{print>$8}' rmsk_TEother_DNase_peaks.txt &
awk '{print>>$10}' all_chromHMM_other_sorted.txt &

python ~/bin/TE_landscape/combine_marks.py mnemonics.txt WGBS_samples.txt DNase_samples.txt H3K27ac_samples.txt TE_combine_marks.txt

# Counts of unique TEs x sample in each WGBS/Dnase/H3K27ac state
#TE_landscape/compare_marks/combine_marks_table_counts.txt	
awk -v OFS='\t' '{a[$9,$11,$13,$15]+=1}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i];}}' TE_combine_marks.txt > combine_marks_table_counts.txt

# Counts of TEs x sample in each chromHMM/WGBS/Dnase/H3K27ac state	 
#TE_landscape/compare_marks/combine_marks_counts.txt	
awk '{print>$8}' TE_combine_marks.txt
for file in E*; do cut -f1-8,11,13,15 $file | sort | uniq | awk -v OFS='\t' '{a[$8,$9,$10,$11]+=1}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i];}}' - >> combine_marks_counts.txt; done &
awk -v OFS='\t' '{a[$2,$3,$4]+=$5}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' combine_marks_counts.txt

# ChIP-seq line plots
# 9/8/2016, 9/14/2016, 9/15/2016, 9/19/2016, 11/3/2016, 11/11/2016, 11/17/2016, 11/18/16, 12/13/2016
# 2/8/2017, 2/27/2017, 2/28/2017, 3/1/2017, 3/6/2017, 3/7/2017, 3/8/2017, 7/31/2017, 8/1/2017, 8/4/2017

# Finds level of five histone modifications over 10kb region centered on TEs in bins	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample.sh	
# Reads in list of samples, not single sample	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample_list.sh	
# Processes multiple samples for multiple input TE sets in parallel	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample_list_parallel.sh	
# Processes multiple samples for multiple input TE sets in parallel, histone acetylation (added option for bin file)	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample_list_parallel_acetyl.sh	
# Processes multiple samples for multiple input TE sets in parallel, Dnase (added option for bin file)	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample_list_parallel_DNase.sh	
# Processes multiple samples for multiple input TE sets in parallel, methylation (added option for bin file)	
#TE_landscape/compare_marks/ChIP_histone/intersect_sample_list_parallel_methyl.sh	

# Each TE extended from center to 10kb region, excluding those that would extend past the ends of the chromosomes.	 
#TE_landscape/features/TEs/rmsk_TE_10kb.txt	
while read chr size; do awk -v OFS="\t" -v chr=$chr -v size=$size '{if($1 == chr && $2 > 4999 && ($3+5000 < size)) print $1,int($2+(($3-$2)/2))-5000,int($2+(($3-$2)/2))+5000,$1,$2,$3,$4,$5,$6,$7}' rmsk_TE.txt >> rmsk_TE_10kb.txt; done < hg19_standard.genome
#TE_landscape/features/TEs/rmsk_other_10kb.txt	
while read chr size; do awk -v OFS="\t" -v chr=$chr -v size=$size '{if($1 == chr && $2 > 4999 && ($3+5000 < size)) print $1,int($2+(($3-$2)/2))-5000,int($2+(($3-$2)/2))+5000,$1,$2,$3,$4,$5,$6,$7}' rmsk_other.txt >> rmsk_other_10kb.txt; done < hg19_standard.genome

# 10kb-extended TE bed file, rearranged so original coordinates are first	 
#TE_landscape/features/TEs/rmsk_TE_10kb.bed	
awk -v OFS="\t" '{print $4,$5,$6,$7,$8,$9,$10,$1,$2,$3}' /bar/epehrsson/TE_landscape/rmsk_TE_10kb.txt > rmsk_TE_10kb.bed
#TE_landscape/features/TEs/rmsk_other_10kb.bed	
awk -v OFS="\t" '{print $4,$5,$6,$7,$8,$9,$10,$1,$2,$3}' rmsk_other_10kb.txt > /scratch/ecp/rmsk_other_10kb.bed

# All TEs in all samples in that state	 
#TE_landscape/chromHMM/chromHMM_ever/TE_[state].txt [12 files]	
grep '1_TssA' all_chromHMM_TE_sorted.txt > TE_1TssA.txt
#TE_landscape/chromHMM/chromHMM_ever/other_[state].txt [15 files]	
while read line; do grep $line all_chromHMM_other_sorted.txt > /scratch/ecp/other_$line\.txt ; done < chromHMM_states.txt

# Epigenetic measure level in bins	 
#TE_landscape/compare_marks/ChIP_histone/TE_[state]-[mod]_average.txt [108 files]	
#TE_landscape/compare_marks/ChIP_histone/other_[state]_all-[mod]_average.txt [108 files]		
cat TE_1TssA_*-H3K27me3_average.txt > TE_1TssA-H3K27me3_average.txt

# Methylation level in bins	 
#TE_landscape/compare_marks/ChIP_histone/other_[state]_all-Meth_average.txt [12 files]	
bash intersect_sample_list_parallel_methyl.sh WGBS_samples.txt other_files.txt all_CpG_Meth.bed rmsk_other_10kb.bed #other_files is list of files of TEs in state by sample (chromHMM_ever)
while read state; do cat other_$state\_E*-all_CpG_Meth.bed_average.txt > other_$state\_all-Meth_average.txt; done < chromHMM_states.txt
#TE_landscape/compare_marks/ChIP_histone/TE_[state]_all-Meth_average.txt [12 files]		 
bash intersect_sample_list_parallel_methyl.sh WGBS_samples.txt TE_files.txt all_CpG_Meth.bed rmsk_TE_10kb.bed
while read state; do cat TE_$state\_E*-all_CpG_Meth.bed_average.txt > TE_$state\_all-Meth_average.txt; done < chromHMM_states.txt

# Averages for all TEs x state x modification	 
#TE_landscape/compare_marks/ChIP_histone/TE_[state]_all-[mod]_average.txt [108 files]	
python ../../bin/TE_landscape/calculate_bin_average_all.py TE_1TssA-H3K27me3_average.txt TE_1TssA_all-H3K27me3_average.txt
#TE_landscape/compare_marks/ChIP_histone/TEother_[state]-[mod]_average.txt [108 files]		 
while read mod; do while read a b; do cat other_$a\_all-$mod\_average.txt TE_$b\_all-$mod\_average.txt > temp.txt; python ../../bin/TE_landscape/calculate_bin_average_all.py temp.txt TEother_$a\-$mod\_average.txt; done <chromHMM_states.txt; done < modifications.txt

# Averages for all TEs x state for methylation	 
#TE_landscape/compare_marks/ChIP_histone/TEother_[state]-Meth_average.txt [12 files]	
while read a b; do cat other_$a\_all-Meth_average.txt TE_$b\_all-Meth_average.txt > temp.txt; python ~/bin/TE_landscape/calculate_bin_average_all.py temp.txt TEother_$a\-Meth_average.txt; done <chromHMM_states.txt

# Epigenetic marks for all TEs in that state	 
#TE_landscape/compare_marks/ChIP_histone/TE_[state]_average.txt [12 files]	
while read state; do while read modification; do awk -v OFS='\t' -v mod=$modification '{print $0, mod}' TE_$state\-$modification\_average.txt >> TE_$state\_average.txt; done < modifications.txt ; done < chromHMM_states.txt

# Combined into one file, all TEs	 
#TE_landscape/compare_marks/ChIP_histone/TEother_average.txt	
while read a b; do while read modification; do awk -v OFS='\t' -v mod=$modification -v state=$a '{print $0, mod, state}' TEother_$a\-$modification\_average.txt >> TEother_average.txt; done < modifications.txt ; done < chromHMM_states.txt

# ChIP-seq line plots
# 9/8/2016, 9/14/2016, 9/15/2016, 9/19/2016, 11/3/2016, 11/11/2016, 11/17/2016, 11/18/16, 12/13/2016
# 2/8/2017, 2/27/2017, 2/28/2017, 3/1/2017, 3/6/2017, 3/7/2017, 3/8/2017, 7/31/2017, 8/1/2017, 8/4/2017
# Updated with summit TEs 5/29/18 to 6/19/18 (should be run on htcf)

# Finds level of histone modifications and DNase over 10kb region centered on TEs in bins	
#intersect_histone.sh
# Finds level of DNA methylation over 10kb region centered on TEs, in bins
#intersect_meth.sh

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
# Combined into one file
cat features/TEs/rmsk_*_10kb.bed > features/TEs/rmsk_TEother_10kb.bed

# TEs in state in sample
awk -v OFS='\t' '{print $0 > "ever/rmsk_TEother_"$8".txt"}' rmsk_TEother_chromHMM_summit_sorted.txt
while read line; do echo $line; awk -v OFS='\t' '{if($10 != "8_ZNF/Rpts") print $0 > "ever/rmsk_TEother_"$8"_"$10".txt"}' ever/rmsk_TEother_$line\.txt; done < mnemonics.txt
while read line; do echo $line; awk -v OFS='\t' '{if($10 == "8_ZNF/Rpts") print $0 > "ever/rmsk_TEother_"$8"_8_ZNF.Rpts.txt"}' ever/rmsk_TEother_$line\.txt; done < mnemonics.txt
 
# Epigenetic profiles by state and mark
#TE_landscape/compare_marks/profile_histone/rmsk_TEother_[state]_[mark]_average.txt [108 files]
while read state; do while read mark; do echo $state $mark; cat averages/rmsk_TEother_E*_$state\_$mark\_average.txt > rmsk_TEother_$state\_$mark\_average.tmp; python calculate_bin_average_all.py rmsk_TEother_$state\_$mark\_average.tmp rmsk_TEother_$state\_$mark\_average.txt; rm rmsk_TEother_$state\_$mark\_average.tmp; done < marks.txt; done < chromHMM_states.txt

# Combined histone profiles
while read state; do while read mark; do awk -v OFS='\t' -v mark=$mark -v state=$state '{print $0, mark, state}' rmsk_TEother_$state\_$mark\_average.txt >> TEother_average.txt; done < marks.txt; done < chromHMM_states.txt

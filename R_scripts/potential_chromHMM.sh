# chromHMM potential 
# 4/25/2016, 4/26/2016, 5/3/2016, 5/12/2016, 6/27/2016, 8/30/2016, 1/26/2017, 2/2/2017, 5/10/2017, 5/29/2017, 8/7/2017, 8/23/2017, 9/14/2017, 9/18/2017

# TEs
#TE_landscape/chromHMM/potential/all_chromHMM_TE_potential_0.txt
#TE_landscape/chromHMM/potential/all_chromHMM_TE_potential_0.75.txt
#TE_landscape/chromHMM/potential/all_chromHMM_TE_potential_0_nodiv.txt
#TE_landscape/chromHMM/potential/all_chromHMM_TE_noCancer_potential_0.txt
#TE_landscape/chromHMM/potential/all_chromHMM_other_potential.txt
#TE_landscape/chromHMM/potential/all_chromHMM_other_potential_noCancer.txt
python ../bin/TE_landscape/potential_TE.py all_chromHMM_TE.txt rmsk_TE.txt chromHMM_states.txt mnemonics.txt test.txt 0
 python ../bin/TE_landscape/potential_TE.py all_chromHMM_TE.txt rmsk_TE.txt chromHMM_states.txt mnemonics_noCancer.txt all_chromHMM_TE_noCancer_potential_0.txt 0 &
 python ../bin/TE_landscape/potential_TE.py all_chromHMM_other.txt rmsk_other.txt chromHMM_states.txt mnemonics.txt all_chromHMM_other_potential.txt 0
 python ../bin/TE_landscape/potential_TE.py all_chromHMM_other.txt rmsk_other.txt chromHMM_states.txt mnemonics_noCancer.txt all_chromHMM_other_potential_noCancer.txt 0

# Number of TEs in state in sample
#TE_landscape/chromHMM/subfamily/state_sample_counts.txt
 awk -v OFS='\t' '{a[$4,$5]+=$6}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' subfamily_state_sample.txt > state_sample_counts.txt

# By class
awk -v OFS='\t' '{a[$5,$8,$10]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' all_chromHMM_TE.txt > class_state_sample.txt
awk -v OFS='\t' '{if($5 == "Other") a[$8,$10]+=1}END{for(i in a) {split (i, sep, SUBSEP); print "SVA", sep[1], sep[2], a[i];}}' all_chromHMM_other.txt >> class_state_sample.txt
awk -v OFS='\t' '{if($5 != "Other") a[$8,$10]+=1}END{for(i in a) {split (i, sep, SUBSEP); print "Other", sep[1], sep[2], a[i];}}' all_chromHMM_other.txt >> class_state_sample.txt

# Refseq promoters
#TE_landscape/chromHMM/Refseq_promoters/potential_refseq_promoters_unique.txt
 python ~/bin/TE_landscape/potential_promoter.py chromHMM_refseq_promoters_unique_sorted.txt ~/genic_features/RefSeq/refseq_promoters_unique_std.txt ../chromHMM_states.txt potential_refseq_promoters_unique.txt 0

# Number of promoters in state in sample
#TE_landscape/chromHMM/Refseq_promoters/promoter_state_sample_count.txt
 awk -v OFS='\t' '{if($1 !~ /_/) a[$5, $6]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}'  chromHMM_refseq_promoters_unique_sorted.txt > promoter_state_sample_count.txt

# Shuffled TEs
 for j in {1..10}; do awk -v OFS='\t' '{a[$8,$9]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' rmsk_TE_shuffle_$j\_sorted.txt > state_sample_count_$j\.txt; done

# Segwey promoters
#TE_landscape/chromHMM/potential/all_chromHMM_promoter_potential_0.75.txt
#TE_landscape/chromHMM/potential/all_chromHMM_promoter_potential_0.txt
 python ../bin/TE_landscape/potential_promoter.py all_chromHMM_promoter.txt promoters.txt chromHMM_states.txt all_chromHMM_promoter_potential_0.txt 0

# Mouse
# chromHMM state for each mm9 TE	 
#TE_landscape/Mouse/chromHMM/mouse_mm9_chromHMM_TE.txt	
while read line ; do awk -v OFS='\t' -v sample=$line '{a[$10, $11, $12, $13, $14, $15, $16, $4, sample]+=$17;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], sep[9], a[i];}}'  mouse_chromHMM_TE/$line\.bed_TE >> mouse_mm9_chromHMM_TE.txt; done < mouse_samples.txt

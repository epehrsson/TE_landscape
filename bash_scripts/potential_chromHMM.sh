# chromHMM potential 
# 4/25/2016, 4/26/2016, 5/3/2016, 5/12/2016, 6/27/2016, 8/30/2016, 1/26/2017, 2/2/2017, 5/10/2017, 5/29/2017, 8/7/2017, 8/23/2017, 9/14/2017, 9/15/2017, 9/18/2017
# Updated 5/8/18 with summit rules

# TEs
python ~/bin/TE_landscape/potential.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt features/TEs/rmsk_TEother.txt chromHMM/chromHMM_states.txt sample_lists/mnemonics.txt chromHMM/potential/rmsk_TEother_chromHMM_summit_potential.txt 0 7 9 8
python ~/bin/TE_landscape/potential.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt features/TEs/rmsk_TEother.txt chromHMM/chromHMM_states.txt sample_lists/mnemonics_noCancer.txt chromHMM/potential/rmsk_TEother_chromHMM_summit_potential_noCancer.txt 0 7 9 8

# Number of TEs in state in sample
awk -v OFS='\t' '{a[$8,$10]+=1}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt > chromHMM/state_sample_counts_summit.txt

# By class
for file in chromHMM/chromHMM_summit_*.txt; do awk -v OFS='\t' -v class=$(basename $file .txt) '{a[$8, $10]+=1}END{for(i in a) {split (i, sep, SUBSEP); print class, sep[1], sep[2], a[i];}}' $file >> chromHMM/class_state_sample_summit.txt; done

# By subfamily
awk -v OFS='\t' '{a[$4, $8, $10]+=1}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt > chromHMM/subfamily/subfamily_state_sample_summit.txt

# Maximum states per TE in any sample
python ~/bin/TE_landscape/state_sharing_intra_max.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt chromHMM/rmsk_TEother_chromHMM_summit_max.txt 7

# Refseq promoters
python ~/bin/TE_landscape/potential.py chromHMM/refseq_promoters_unique_chromHMM_summit_sorted.txt ~/genic_features/RefSeq/refseq_promoters_unique_std.txt chromHMM/chromHMM_states.txt sample_lists/mnemonics.txt chromHMM/potential/refseq_promoters_chromHMM_summit_potential.txt 0 4 6 5

# Number of promoters in state in sample
awk -v OFS='\t' '{if($1 !~ /_/) a[$5, $7]+=1}END{for(i in a){split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' chromHMM/refseq_promoters_unique_chromHMM_summit_sorted.txt > chromHMM/promoters_state_sample_counts_summit.txt

# Shuffled TEs
#TE_landscape/chromHMM/shuffled_TEs/rmsk_TE_shuffle_#_max.txt [10 files]
#TE_landscape/chromHMM/shuffled_TEs/rmsk_TE_shuffle_#_potential.txt [10 files]
#TE_landscape/chromHMM/shuffled_TEs/state_sample_count_#.txt [10 files]
 for j in {1..10}; do awk -v OFS='\t' '{a[$8,$9]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' rmsk_TE_shuffle_$j\_sorted.txt > state_sample_count_$j\.txt; done

# Segwey promoters
#TE_landscape/chromHMM/potential/all_chromHMM_promoter_potential_0.75.txt
#TE_landscape/chromHMM/potential/all_chromHMM_promoter_potential_0.txt
 python ../bin/TE_landscape/potential_promoter.py all_chromHMM_promoter.txt promoters.txt chromHMM_states.txt all_chromHMM_promoter_potential_0.txt 0

# Shuffled TEs, epigenetic
#TE_landscape/chromHMM/shuffled_TEs/state_sample_count_#.txt [10 files]
 for j in {1..10}; do awk -v OFS='\t' '{a[$8,$9]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], a[i];}}' rmsk_TE_shuffle_$j\_sorted.txt > state_sample_count_$j\.txt; done
 
#TE_landscape/chromHMM/shuffled_TEs/rmsk_TE_shuffle_#_max.txt [10 files]
for i in {1..10}; do python ~/bin/TE_landscape/state_sharing_intra_max.py rmsk_TE_shuffle_$i\_sorted.txt rmsk_TE_shuffle_$i\_max.txt 7; done

# Intersection of shuffled TEs and mappability
split -l 5000000 ~/TE_landscape/mappability/wgEncodeCrgMapabilityAlign36mer.bedGraph
#TE_landscape/features/shuffled_TEs/run_intersect.sh

# Shuffled TEs intersect with mappability
for i in {1..10}; do for file in mappability/x*; do bedtools intersect -wo -a rmsk_TE_shuffle_$i\.txt -b $file >> mappability/rmsk_TE_shuffle_$i\_mappabililty.bed; done; done

# Average mappability per shuffled TE
#TE_landscape/mappability/shuffled/rmsk_TE_shuffle_#_mappabililty_avg.txt [10 files]
for i in {1..10}; do python ~/bin/TE_landscape/mappability_TE.py rmsk_TE_shuffle_$i\_mappabililty.bed ../rmsk_TE_shuffle_$i\.txt rmsk_TE_shuffle_$i\_mappabililty_avg.txt; done


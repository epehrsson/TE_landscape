# Shuffled mappability
# 8/25/2017, 8/29/2017, 8/30/2017

# Intersection of shuffled TEs and mappability	 
split -l 5000000 ~/TE_landscape/mappability/wgEncodeCrgMapabilityAlign36mer.bedGraph
#TE_landscape/features/shuffled_TEs/run_intersect.sh		

# Average mappability per shuffled TE	 
#TE_landscape/mappability/shuffled/rmsk_TE_shuffle_#_mappabililty_avg.txt [10 files]	
for i in {1..10}; do python ~/bin/TE_landscape/mappability_TE.py rmsk_TE_shuffle_$i\_mappabililty.bed ../rmsk_TE_shuffle_$i\.txt rmsk_TE_shuffle_$i\_mappabililty_avg.txt; done

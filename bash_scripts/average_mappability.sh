# From split genome mappability bedGraph file
# Average mappability by chromHMM state (1/9/18-1/10/18)
while read line; do echo $line; for file in x*; do bedtools intersect -wo -a $file -b ~/TE_landscape/raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed >> mappabililty_$line\.bed; done; python ~/bin/TE_landscape/calculate_average_mappability.py mappabililty_$line\.bed ~/TE_landscape/chromHMM/chromHMM_states.txt mappability_$line\_average.txt; done < mnemonics.txt

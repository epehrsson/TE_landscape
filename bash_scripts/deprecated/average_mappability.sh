# 36bp mappability, genome-wide	 
#TE_landscape/mappability/wgEncodeCrgMapabilityAlign36mer.bedGraph	
wget hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign36mer.bw

# From raw file
~/bin/bigWigToBedGraph wgEncodeCrgMapabilityAlign36mer.bw.2 wgEncodeCrgMapabilityAlign36mer.bedGraph

# Average mappability, whole genome
awk '{sum+=($4*($3-$2));total+=$3-$2;}END{print sum/total}' wgEncodeCrgMapabilityAlign36mer.bedGraph

# Average mappability, TEs 	 
split -l 5000000 ~/TE_landscape/mappability/wgEncodeCrgMapabilityAlign36mer.bedGraph
for file in x*; do bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file >> rmsk_TEother_merge_map.txt; done
awk '{sum+=($7*$8);total+=$8;}END{print sum/total}' rmsk_TEother_merge_map.txt

# From split genome mappability bedGraph file
# Average mappability by chromHMM state (1/9/18-1/10/18)
while read line; do echo $line; for file in x*; do bedtools intersect -wo -a $file -b ~/TE_landscape/raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed >> mappabililty_$line\.bed; done; python ~/bin/TE_landscape/calculate_average_mappability.py mappabililty_$line\.bed ~/TE_landscape/chromHMM/chromHMM_states.txt mappability_$line\_average.txt; done < mnemonics.txt

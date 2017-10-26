# Mappability
# 11/4/2016, 2/2/2017, 8/27/2017

# TE 36bp mappability file	 
#TE_landscape/mappability/rmsk_TE_mappability_36mer.txt	
while read line; do awk -v subfam=$line '{if($4 == subfam) print $0}' rmsk.txt.mapability.36mer >> rmsk_TE_mappability_36mer.txt; done < subfamilies.txt
#TE_landscape/mappability/rmsk_other_mappability_36mer.txt	
while read line; do awk -v subfam=$line '{if($4 == subfam) print $0}' rmsk.txt.mapability.36mer >> rmsk_other_mappability_36mer.txt; done < other_subfamilies.txt

# From raw file
~/bin/bigWigToBedGraph wgEncodeCrgMapabilityAlign36mer.bw.2 wgEncodeCrgMapabilityAlign36mer.bedGraph

# Average mappability, whole genome	 
awk '{sum+=($4*($3-$2));total+=$3-$2;}END{print sum/total}' wgEncodeCrgMapabilityAlign36mer.bedGraph
# Average mappability, TEs 	 
for file in x*; do bedtools intersect -wo -a rmsk_TEother_merge.txt -b $file >> rmsk_TEother_merge_map.txt; done
awk '{sum+=($7*$8);total+=$8;}END{print sum/total}' rmsk_TEother_merge_map.txt

# Shuffled TEs intersect with mappability
for i in {1..10}; do for file in mappability/x*; do bedtools intersect -wo -a rmsk_TE_shuffle_$i\.txt -b $file >> mappability/rmsk_TE_shuffle_$i\_mappabililty.bed; done; done

# Filter repeat file
# 4/19/2016, 4/25/2016, 4/27/2016, 5/9/2016, 5/19/2016, 6/9/2016, 9/12/2016, 9/14/2016, 1/26/2017, 2/2/2017, 2/6/2017, 2/8/2017, 5/10/2017, 6/14/2017

# RepeatMasker file restricted to chr 1-22, X, Y, M	
#TE_landscape/features/rmsk_standard.txt	
awk -v OFS='\t' '{if($6 !~ /_/) print $0}' rmsk.txt

# RepeatMasker file restricted to all TE classes, chr 1-22, X, Y, M	 
#TE_landscape/features/TEs/rmsk_TE.txt	
awk -v OFS='\t' '{if(($12 == "LTR" || $12 == "DNA" || $12 == "SINE" || $12 == "LINE") && ($6 !~ /_/)) print $6, $7, $8, $11, $12, $13, $10}' rmsk.txt > rmsk_TE.txt
#TE_landscape/features/TEs/rmsk_other.txt	
awk -v OFS='\t' '{if(($12=="Unknown"||$12=="Unknown?"||$12=="DNA?"||$12=="LINE?"||$12=="SINE?"||$12=="LTR?"||$12=="Other"||$12=="RC") && ($6 !~ /_/))print $6, $7, $8, $11, $12, $13, $10}' rmsk.txt > rmsk_other.txt
#TE_landscape/features/TEs/rmsk_TEother.txt	
cat rmsk_TE.txt rmsk_other.txt > rmsk_TEother.txt
#Total length
awk '{sum+=$3-$2}END{print sum}' rmsk_TEother.txt #1,391,359,975

# Merged TE RepeatMasker file	 
#TE_landscape/features/TEs/merge/rmsk_TE_merge.txt	
bedtools merge -i rmsk_TE.txt > rmsk_TE_merge.txt
while read line; do grep $line rmsk_TE.txt | sort -k 1,1 -k 2,2n > temp.bed; bedtools merge -i temp.bed | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' - ; done < TE_class.txt
# Merged all TE RepeatMasker file	 
#TE_landscape/features/TEs/merge/rmsk_TEother_merge.txt	
cat rmsk_TE.txt rmsk_other.txt | sort -k1,1 -k2,2n - | bedtools merge -i - > rmsk_TEother_merge.txt
#Total length
awk '{sum+=$3-$2}END{print sum}' rmsk_TEother_merge.txt #1,389,947,349

# All TE RepeatMasker bedfile with unique number for each	 
#TE_landscape/features/TEs/rmsk_TE_numbered.txt	
for each	 awk -v OFS="\t" '{print $1,$2,$3,NR,".",$7}' rmsk_TE.txt > rmsk_TE_numbered.txt
#TE_landscape/features/TEs/rmsk_other_numbered.txt	
for each	 awk -v OFS="\t" '{print $1,$2,$3,NR,".",$7}' rmsk_other.txt > rmsk_other_numbered.txt
#TE_landscape/features/TEs/rmsk_TEother_numbered.txt	
cat rmsk_TE.txt rmsk_other.txt | awk -v OFS="\t" '{print $1,$2,$3,NR,".",$7}' - > rmsk_TEother_numbered.txt

# Simple repeats
# RepeatMasker file restricted to repeats, chr 1-22, X, Y, M	 
#TE_landscape/features/repeats/rmsk_repeats.txt	
awk -v OFS='\t' '{if($12 == "Low_complexity" || $12 == "Satellite" || $12 == "Simple_repeat") print $6, $7, $8, $11, $12, $13, $10}' rmsk.txt > rmsk_repeats.txt

# Merged repeat RepeatMasker file	 
#TE_landscape/features/repeats/rmsk_repeats_merge.txt	
bedtools merge -i rmsk_repeats.txt > rmsk_repeats_merge.txt

# Mouse
# mm9
# Mouse TEs, chromosomes 1-19, X, Y, M only	
#TE_landscape/features/mouse/mm9_rmsk_standard.txt	
awk -v OFS='\t' '{if($6 !~ /_/) print $0}' mm9_rmsk.txt > mm9_rmsk_standard.txt

# RepeatMasker file restricted to all TE classes, chr 1-19, X, Y, and M	
#TE_landscape/features/mouse/TEs/mm9_rmsk_TE.txt	
awk -v OFS='\t' '{if(($12 == "LTR" || $12 == "DNA" || $12 == "SINE" || $12 == "LINE") && ($6 !~ /_/)) print $6, $7, $8, $11, $12, $13, $10}' mm9_rmsk.txt > mm9_rmsk_TE.txt
#TE_landscape/features/mouse/TEs/mm9_rmsk_other.txt	
awk -v OFS='\t' '{if(($12 == "Other" || $12 == "RC" || $12 == "Unknown") && ($6 !~ /_/)) print $6, $7, $8, $11, $12, $13, $10}' mm9_rmsk.txt > mm9_rmsk_other.txt
#TE_landscape/features/mouse/TEs/mm9_rmsk_TEother.txt	

# Merged TE RepeatMasker file	 
#TE_landscape/features/mouse/TEs/mm9_rmsk_TEmerge.txt	
bedtools merge -i mm9_rmsk_TE.txt > mm9_rmsk_TEmerge.txt
#TE_landscape/features/mouse/TEs/mm9_rmsk_TEother_merge.txt	
cat mm9_rmsk_TE.txt mm9_rmsk_other.txt | sort -k1,1 -k2,2n - | bedtools merge -i - > mm9_rmsk_TEother_merge.txt

# mm10
# mm10 mouse TEs, chromosomes 1-19, X, Y, M only	 
#TE_landscape/features/mouse/TEs/mm10_rmsk_TE.txt	
awk -v OFS='\t' '{if(($12 == "LTR" || $12 == "DNA" || $12 == "SINE" || $12 == "LINE" || 12 == "Other" || $12 == "RC" || $12 == "Unknown" || $12 == "DNA?" || $12 == "LINE?" || $12 == "LTR?" || $12 == "SINE?" || $12 == "RC?") && ($6 !~ /_/)) print $6, $7, $8, $11, $12, $13, $10}' rmsk_mm10.txt > mm10_rmsk_TE.txt
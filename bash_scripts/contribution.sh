# Contribution
# 4/26/2016, 4/27/2016, 2/3/2017, 3/2/2017, 8/18/2017

# Number of bases in each state in TEs (old!)
#TE_landscape/chromHMM/chromHMM_TEmerge_all.txt	

# Number of bases in each TE class in each state across all samples (old!)
#TE_landscape/chromHMM/class/chromHMM_class_all.txt	
while read line; do while read line2; do awk -v OFS='\t' -v class=$line -v state=$line2 '{if(($5 == class)&&($8 == state))SUM+=$9}END{print class, state, SUM}' all_chromHMM_TE.txt; done < chromHMM_states.txt; done < TE_class.txt > chromHMM_class_all.txt

# Number of bp in each state, overall/in TEs/by TE class across all tissues (old!)
#TE_landscape/chromHMM/TE_contribution.txt		
#From chromHMM_class_all.txt
while read line; do grep $line all_chromHMM_TE_merge.txt | awk '{SUM+=$5}END{print SUM}' -; done < chromHMM_states.txt
while read line; do grep $line chromHMM_blocks.txt | awk '{SUM+=$1}END{print SUM}' -; done < chromHMM_states.txt

# (In Excel)	 
#TE_landscape/chromHMM/TEother_contribution.txt	
cat chromHMM_other/E*_15_coreMarks_mnemonics.bed_TEother_merge | awk -v OFS='\t' '{a[$7]+=$8;}END{for(i in a) {print i, a[i];}}' -
# Added Unconfident (total size), Unconfident-RC	 
awk -v OFS='\t' '{if(($5=="LTR?")||($5=="DNA?")||($5=="LINE?")||($5=="SINE?")||($5=="Unknown?")||($5=="Unknown")) print $0}' rmsk_other.txt | bedtools merge -i - | awk '{sum+=$3-$2}END{print sum}'

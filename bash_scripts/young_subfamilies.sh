# Young subfamilies enriched
# 8/2/2017

# Young subfamily members in state when enriched 	 
#TE_landscape/young_indv.txt	
awk -v OFS='\t' '{if(($4 == "AluYf5") && ($8 == "6_EnhG") && ($10 == "E030")) print $0}' chromHMM/all_chromHMM_TE_sorted.txt >> young_indv.txt
awk -v OFS='\t' '{if(($4 == "AluYf4") && ($8 == "1_TssA") && ($10 == "E091")) print $0}' chromHMM/all_chromHMM_TE_sorted.txt >> young_indv.txt
awk -v OFS='\t' '{if(($4 == "AluYf5") && ($8 == "1_TssA") && ($10 == "E091")) print $0}' chromHMM/all_chromHMM_TE_sorted.txt >> young_indv.txt
awk -v OFS='\t' '{if(($4 == "AluYd8") && ($8 == "1_TssA") && ($10 == "E043")) print $0}' chromHMM/all_chromHMM_TE_sorted.txt >> young_indv.txt

#TE_landscape/young_Gencode_genes.txt	
# Young subfamily members overlapping GENCODE genes	 
bedtools intersect -wo -a young_indv.txt -b ~/genic_features/Gencode_v19_genes.txt > young_Gencode_genes.txt

# Young subfamily members overlapping GENCODE promoters	 
#TE_landscape/young_Gencode_promoter.txt	
bedtools intersect -wo -a young_indv.txt -b ~/genic_features/Gencodev19_promoter.txt > young_Gencode_promoter.txt

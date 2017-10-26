# TEs in state when subfamily is enriched
# 9/30/2016, 5/4/2017, 5/24/2017, 5/29/2017, 5/30/2017, 5/31/2017

# Coordinates of subfamilies enriched in blood
#TE_landscape/enrichment/blood.txt	
# Coordinates of subfamilies enriched in mesenchymal cluster
#TE_landscape/enrichment/mesen_sub.txt	
# MER57E3 elements in the 1TssA state in at least 20 samples
#TE_landscape/enrichment/MER57E3_1TssA_20.bed	

# Individual TEs in the state from candidate subfamilies	 
#TE_landscape/enrichment/candidate_[state]_indv.txt [5 files]	
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "1_TssA")) print $0}' all_chromHMM_TE_sorted.txt >> candidate_1TssA_indv.txt; done < candidate_1TssA.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "2_TssAFlnk")) print $0}' all_chromHMM_TE_sorted.txt >> candidate_2TssAFlnk_indv.txt; done < candidate_2TssAFlnk.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "2_TssAFlnk")) print $0}' all_chromHMM_other_sorted.txt >> candidate_2TssAFlnk_indv.txt; done < candidate_2TssAFlnk.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "6_EnhG")) print $0}' all_chromHMM_TE_sorted.txt >> candidate_6EnhG_indv.txt; done < candidate_6EnhG.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "7_Enh")) print $0}' all_chromHMM_TE_sorted.txt >> candidate_7Enh_indv.txt; done < candidate_7Enh.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if(($4 == subfam) && ($8 == "7_Enh")) print $0}' all_chromHMM_other_sorted.txt >> candidate_7Enh_indv.txt; done < candidate_7Enh.txt

# Individual TEs from candidate subfamilies, filtered to enriched samples	 
#TE_landscape/enrichment/candidate_[state]_indv_enrich.txt [5 files]	
while read a b; do awk -v OFS='\t' -v subfam=$a -v sample=$b '{if(($4 == subfam) && ($10 == sample)) print $0}' candidate_1TssA_indv.txt; done < candidate_1TssA_enrich.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' > candidate_1TssA_indv_enrich.txt
while read a b; do awk -v OFS='\t' -v subfam=$a -v sample=$b '{if(($4 == subfam) && ($10 == sample)) print $0}' candidate_2TssAFlnk_indv.txt; done < candidate_2TssAFlnk_enrich.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' > candidate_2TssAFlnk_indv_enrich.txt
while read a b; do awk -v OFS='\t' -v subfam=$a -v sample=$b '{if(($4 == subfam) && ($10 == sample)) print $0}' candidate_3TxFlnk_indv.txt; done < candidate_3TxFlnk_enrich.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' > candidate_3TxFlnk_indv_enrich.txt
while read a b; do awk -v OFS='\t' -v subfam=$a -v sample=$b '{if(($4 == subfam) && ($10 == sample)) print $0}' candidate_6EnhG_indv.txt; done < candidate_6EnhG_enrich.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' > candidate_6EnhG_indv_enrich.txt 
while read a b; do awk -v OFS='\t' -v subfam=$a -v sample=$b '{if(($4 == subfam) && ($10 == sample)) print $0}' candidate_7Enh_indv.txt; done < candidate_7Enh_enrich.txt | awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' > candidate_7Enh_indv_enrich.txt

# TEs in state when enriched, by subfamily	 
#TE_landscape/enrichment/*_enrich1TssA.bed [40 files]	
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $1, $2, $3, $4, $8, $7}' candidate_1TssA_indv_enrich.txt > enrichment/$line\_enrich1TssA.bed; done < enrichment/candidate_1TssA.txt &			5/30/2017
#TE_landscape/enrichment/*_enrich2TssAFlnk.bed [93 files]		 
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $1, $2, $3, $4, $8, $7}' candidate_2TssAFlnk_indv_enrich.txt > enrichment/$line\_enrich2TssAFlnk.bed; done < candidate_2TssAFlnk.txt			5/31/2017
#TE_landscape/enrichment/*_enrich3TxFlnk.bed [14 files]		 
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $1, $2, $3, $4, $8, $7}' candidate_3TxFlnk_indv_enrich.txt > enrichment/$line\_enrich3TxFlnk.bed; done < candidate_3TxFlnk.txt &			5/30/2017
#TE_landscape/enrichment/*_enrich6EnhG.bed [38 files]		 
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $1, $2, $3, $4, $8, $7}' candidate_6EnhG_indv_enrich.txt > enrichment/$line\_enrich6EnhG.bed; done < candidate_6EnhG.txt			5/31/2017
#TE_landscape/enrichment/*_enrich7Enh.bed [166 files]		 
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $1, $2, $3, $4, $8, $7}' candidate_7Enh_indv_enrich.txt > enrichment/$line\_enrich7Enh.bed; done < candidate_7Enh.txt &			5/30/2017

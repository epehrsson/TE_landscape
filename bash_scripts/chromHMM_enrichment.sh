# chromHMM enrichment - bp and members in state
# 4/27/2016, 8/28/2016, 2/2/2017

# Number of bp in each subfamily x sample x state	 
#TE_landscape/chromHMM/subfamily/all_chromHMM_TE_subclass.txt	
q all_chromHMM_TE.txt

# Number of TE subfamily members in each state by sample (≥1bp)	 
#TE_landscape/chromHMM/subfamily/subfamily_state_sample_old.txt	
awk -v OFS='\t' '{if($9 > 0) a[$4, $5, $6, $8, $10]+=1;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], a[i];}}' all_chromHMM_TE.txt > subfamily_state_sample.txt

# Number of other TE subfamily members in each state by sample (≥1bp)	 
#TE_landscape/chromHMM/subfamily/other_subfamily_state_sample.txt	
awk -v OFS='\t' '{if($9 > 0) a[$4, $5, $6, $8, $10]+=1;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], a[i];}}' all_chromHMM_other.txt > other_subfamily_state_sample.txt

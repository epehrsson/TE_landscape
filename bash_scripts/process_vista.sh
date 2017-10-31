# Process VISTA enhancers
# 10/28/2016

# Human, positively validated VISTA enhancer headers
#TE_landscape/features/vista_enhancers/vista_enhancers_human.txt	

# Human, positively validated VISTA enhancers in bed format
#TE_landscape/features/vista_enhancers/vista_enhancers_human.bed	

cat features/intersect_features/vista_enhancers_rmsk_*.txt | awk -v OFS='\t' '{a[$4,$5,$6,$7,$8,$9,$10]+=1}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' -  > features/intersect_features/vista_enhancers_TE.txt
